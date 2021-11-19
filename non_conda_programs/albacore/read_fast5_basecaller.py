#!/home/user/miniconda3/envs/nanopore/bin/python3.6
"""
read_fast5_basecaller.py
Created on: February 12, 2016
Proprietary and confidential information of Oxford Nanopore Technologies, Limited
All rights reserved; (c)2016: Oxford Nanopore Technologies, Limited

Main entry point script for running read fast5 files through albacore.

Expected Input:

    *  Folder containing read fast5 files with raw data.
    *  Flow-cell and kit combination, or directly specify the configuration file
       to use. This will contain analysis parameters.
    *  Path to save log and intermediate files.

What It Does:

    * The save path is created if it does not exist.
    * A workspace folder is created under the save path.
    * Read files are copied into the workspace folder (only for fast5 output
      option).
    * If the recursive option is used, subfolders of the input path are searched
      as well.
    * Read files are processed, and basecall results are written back into the
      read fast5 files (only for fast5 output option).
    * Files that are completely finished are copied into batch subfolders
      of the workspace folder, maintaining an internal counter and starting
      a new folder when the counter reaches files_per_batch_folder (only for
      fast5 output option).
    * Fastq files are created containing the basecalled reads (only for fastq
      output option).
    * Log and summary files are created in the save folder as the program runs.

Notes:
    The --config parameter can be the full path to the config file, or just
    the filename of one of the standard config files provided with the
    software. If no path is provided, it will first look for the 'data' folder
    within the same path as this script. This is appropriate for checked out
    copies of the code. If this path does not exist, it assumed you are running
    an installed version of the code, and looks in the default installation
    path of the project, which is platform specific.

    The configuration file will specify model files by filename only. If you
    want to use a model that is not in the standard data folder, you will need
    to specify the path with the --data_path option. Note that ALL models used
    must be in the specified folder.

    This script spawns multiple threads, configured with the --worker_threads
    option. However, all IO and some of the work for aggregating basecall
    results is done in the master thread. So for processing large amounts of
    data on a system with lots of cores (ie more than 8), it is recommended
    that you distribute the data into multiple folders, and run multiple
    instances of this script, each with its own input folder.
"""
import os
import sys
import time
import shutil
import logging
from glob import glob
from argparse import ArgumentParser, SUPPRESS
from configparser import ConfigParser

from albacore.config_utils import (get_barcoding_options,
                                   parse_config_override, resume_options_1d,
                                   parse_telemetry_options)
from albacore.output_utils import output_summary, write_telemetry
from albacore.input_utils import list_input_files

from albacore import MIN_CFG_VERSION_I, VERSION_I, VERSION
from albacore.pipeline import Pipeline
from albacore.fast5_fastq_data_handler import Fast5FastQDataHandler
from albacore.summary import SummaryCSV
from albacore.config_selector import print_workflows, WORKFLOW_LINEAR
from albacore.log_utils import initialise_logger, log_submission
from albacore.path_utils import initialise_paths, get_default_path
from albacore.config_selector import choose_config
from albacore.telemetry import TelemetryAggregator
from albacore import telemetry_utils
from albacore.aligner_utils import initialise_aligner
from albacore.utilities.pipeline_utils import get_progress_bar


def main(opts, config_override):
    """ Prepares config values and logging before starting the main workflow
    processing.

    :param opts: Configuration options from the command line.
    :param config_override: dict whose key-value pair override the
        corresponding values from the configuration file.
    :returns: 0 on success.
    :rtype: int

    """

    # Intialise paths and loggers
    opts.read_path = initialise_paths(opts.save_path)
    logger = initialise_logger(opts.save_path, opts.debug, logfilename=opts.logfile, resume_mode=opts.resume)
    if opts.data_path is None:
        opts.data_path = opts.default_path

    # If config isn't set then we'll use flowcell and kit to find one
    if opts.config is None:
        opts.config, barcoding = choose_config(opts.data_path, opts.flowcell, opts.kit)
        if barcoding: # If the flowcell+kit includes barcoding we always run barcoding
            opts.barcoding = True
            logger.info('Kit {} includes barcodes "--barcoding" option set to True'.format(opts.kit))
    else: # Otherwise we look for it first as an abs path, then in the data path
        if not os.path.exists(opts.config):
            opts.config = os.path.join(opts.data_path, opts.config)
            if not os.path.exists(opts.config):
                raise Exception('Config file "{}" does not exist.'.format(opts.config))

    config = ConfigParser(interpolation=None)
    config.read(opts.config)

    # Clear out any compatibility information in the config file
    config.remove_section('compatibility')

    # The python loader expects a model path to be set, but we don't need it elsewhere.
    if config.has_section('basecaller') and not config.has_option('basecaller', 'model_path'):
        config.set('basecaller', 'model_path', opts.data_path)
    elif config.has_section('basecaller_hmm') and not config.has_option('basecaller_hmm', 'model_path'):
        config.set('basecaller_hmm', 'model_path', opts.data_path)
    elif config.has_section('basecaller_ocl') and not config.has_option('basecaller_ocl', 'model_path'):
        config.set('basecaller_ocl', 'model_path', opts.data_path)

    # in case a specific version is specified in the config file, check if this
    # version is in the range [min_config_version, current_version]
    config_version = config.get('pipeline', 'albacore_version', fallback=None)
    if config_version:
        config_version_info = tuple([int(num) for num in
                                     config_version.split('.')])
        if (config_version_info < MIN_CFG_VERSION_I
            or config_version_info > VERSION_I):
            msg = ("Current version of albacore ({}) is incompatible with the "
                   "config file's version ({}).".format(VERSION,
                                                        config_version))
            logger.error(msg)
            raise Exception(msg)
        config.remove_option('pipeline', 'albacore_version')

    # if --barcode_detector.config_file=... is given on the command line, we need
    # to write this to config *before* extracting option fields from it.
    if ('barcode_detector', 'config_file') in config_override:
        config.set('barcode_detector', 'config_file',
                   config_override[('barcode_detector', 'config_file')])
    get_barcoding_options(config, opts.barcoding, opts.data_path, logger,
                          opts.resume)
    # Now the config file is complete, all "usual" command-line arguments have
    # been parsed, and all sub-config files have been added. Now look for command-line overrides:
    for key, value in config_override.items():
        (section, parameter) = key
        if not config.has_section(section):
            config.add_section(section)
        config.set(section, parameter, value)
        logger.info('config parameter {}.{} overwritten by command-line value '
                    '{}'.format(section, parameter, value))

    # find external aligners and references; delete respective config section
    # if alignment or calibration strand detection is disabled
    references_loaded = \
        initialise_aligner(config, opts.align_ref, opts.data_path, 'aligner') + \
        initialise_aligner(config, None, opts.data_path, 'calib_detector')

    basecall_type = config.get('pipeline', 'basecall_type')
    raw_basecalling = not config.has_section('event_detector')
    ocl_basecall = config.has_section('basecaller_ocl')
    basecaller_section = 'basecaller_ocl' if ocl_basecall else 'basecaller'
    if basecall_type != 'linear' and raw_basecalling:
        raise Exception('Basecalling without event detection (a.k.a. raw basecalling) is currently '
                        'only available for 1D.')
    desc_base = {'linear': '1d',
                 'full_2d': '2d'}
    desc_file = 'layout_' + ('raw_' if raw_basecalling else '') + 'basecall_' \
        + desc_base[basecall_type] + ('_ocl' if ocl_basecall else '') + '.jsn'
    logger.info("Using config file: {}".format(desc_file))
    config.set('pipeline', 'desc_file', os.path.join(opts.data_path, desc_file))

    if config.has_option(basecaller_section, 'opencl_kernel_path'):
        config.set(basecaller_section, 'opencl_kernel_path',
                   os.path.join(opts.data_path, 
                                config.get(basecaller_section, 'opencl_kernel_path')))

    files, additional_summary_lines = list_input_files(opts) # will set opts.resume to False
                                                             # if no summary file is present

    # Initialize the summary file handler.
    include_calib = bool(config.get('calib_detector', 'method'))
    include_alignment = bool(config.get('aligner', 'method'))
    include_barcoding = config.has_option('barcode_detector',
                                          'arrangements_files')
    logger.info('include calibration strand detection: '
                '{}'.format(include_calib))
    logger.info('include alignment: {}'.format(include_alignment))
    logger.info('include barcoding: {}'.format(include_barcoding))

    telemetry_enabled = not opts.disable_pings
    (telemetry_urls,
     telemetry_segment_duration,
     telemetry_analysis_id) = parse_telemetry_options(config)
    telemetry_output_file = os.path.join(opts.save_path,
                                         'sequencing_telemetry.js')
    telemetry_aggregator = TelemetryAggregator(
        telemetry_segment_duration,
        albacore_opts=opts,
        include_1d_basecalling=True,
        include_1dsq_basecalling=False,
        include_calibration=include_calib,
        include_alignment=include_alignment,
        include_barcoding=include_barcoding,
        analysis_id=telemetry_analysis_id)
    logger.info('telemetry enabled: {}'.format(telemetry_enabled))
    summary_fields = output_summary(basecall_type, include_calib,
                                    include_barcoding, include_alignment)
    try:
        with SummaryCSV(os.path.join(opts.save_path, opts.summfile),
                        summary_fields,
                        start_new_file=not opts.resume) as summary:
            for filename in additional_summary_lines:
                summary.write({'filename': os.path.basename(filename)})
            process_pipeline(opts, config, files, summary, telemetry_aggregator)
        logger.info('Done.')
    except Exception as exc:
        logger.exception('Fatal error in pipeline')
        print('Fatal error in pipeline:', exc, file=sys.stderr)
        raise
    logger.info('Writing telemetry')
    write_telemetry(telemetry_aggregator, telemetry_output_file)
    if telemetry_enabled:
        for url in telemetry_urls:
            logger.info('Pinging telemetry to {}'.format(url))
            telemetry_utils.ping(telemetry_aggregator, url)
    if references_loaded:
        unload_commands = [config.get('aligner', 'unload_command').
                           format(executable=config.get('aligner',
                                                        'executable'),
                                  reference=ref)
                           for ref in references_loaded]
        print('\n    CAVE: Reference is still loaded. To unload, call:')
        for cmd in set(unload_commands):
            print('    ' + cmd)
        print('    ')
    return 0


def process_pipeline(opts, config, files, summary, telemetry_aggregator=None):
    logger = logging.getLogger('albacore')

    # Zero threads actually means turning off multithreading. But for purposes
    # of determining whether the pipeline is finished working, we need to know
    # the actual number of workers, which is still one.
    num_workers = opts.worker_threads if opts.worker_threads > 0 else 1
    # channels = list(range(1, num_workers + 1))
    reverse_direction = config.getboolean('basecaller', 'reverse_direction',
                                          fallback=False)
    u_substitution = config.getboolean('basecaller', 'u_substitution',
                                       fallback=False)
    new_config = os.path.join(opts.save_path, opts.cfgfile)
    with open(new_config, 'w') as cfg_out:
        config.write(cfg_out)

    # Initialize the basecall pipeline.
    data_handler = Fast5FastQDataHandler(new_config, opts)
    # Note that we pass opts.worker_threads (which may be zero), not num_workers.
    pipeline = Pipeline(opts.worker_threads, data_handler, new_config,
                        reverse_direction=reverse_direction,
                        u_substitution=u_substitution,
                        telemetry_aggregator=telemetry_aggregator)

    num_files = len(files)
    file_index = 0
    num_complete = 0

    # only copy files if Fast5 output is enabled
    if 'fast5' in data_handler.output_format:
        copied_files = []
        num_copied = 0
    else:
        copied_files = files
        num_copied = num_files
    no_more_files = False

    pbar = get_progress_bar(opts.debug, num_files)

    PASSING_DATA = 0
    FINISHING_JOBS = 1
    ALL_JOBS_COMPLETE = 2

    pipeline_state = PASSING_DATA
    while pipeline_state != ALL_JOBS_COMPLETE:
        if pipeline_state == FINISHING_JOBS:
            num_cached_reads = pipeline.finish_all_jobs()
            pipeline_state = ALL_JOBS_COMPLETE
        else:
            num_cached_reads = pipeline.cache_completed_calls()
            free_workers = pipeline.workers_ready()
            if free_workers == num_workers and no_more_files:
                # All workers are idle, and there is no more data to process. We're done.
                pipeline_state = FINISHING_JOBS
            will_sleep = True
        if free_workers == 0 and num_cached_reads == 0:
            if num_copied < num_files:
                try:
                    fname = files[num_copied]
                    destination = os.path.join(opts.read_path,
                                               os.path.basename(fname))
                    subpath = os.path.dirname(destination)
                    if not os.path.exists(subpath):
                        os.makedirs(subpath)
                    shutil.copyfile(fname, destination)
                    copied_files.append(destination)
                    num_copied += 1
                    will_sleep = False
                except Exception:
                    logger.error('Could not access file {}. '
                                 'File will be skipped.'.format(fname))
                    files.remove(fname)
                    num_files = len(files)

        if free_workers > 0 and file_index < num_files:
            # Submit new jobs before processing complete ones, so that our worker
            # threads are not left idle.
            submit = min(free_workers, num_files - file_index)
            for i in range(submit):
                if num_copied <= file_index:
                    try:
                        fname = files[file_index]
                        destination = os.path.join(opts.read_path,
                                                   os.path.basename(fname))
                        subpath = os.path.dirname(destination)
                        if not os.path.exists(subpath):
                            os.makedirs(subpath)
                        shutil.copyfile(fname, destination)
                        copied_files.append(destination)
                        num_copied += 1
                    except Exception:
                        logger.error('Could not access file {}. '
                                     'File will be skipped.'.format(fname))
                        files.remove(fname)
                        num_files = len(files)

                if num_copied > file_index:
                    short_name = os.path.basename(copied_files[file_index])
                    logger.info(log_submission(short_name))
                    try:
                        status = pipeline.submit_read(copied_files[file_index])
                        if not status:
                            logger.warning('Could not extract read data from '
                                           'file "{}".'.format(short_name))
                    except Exception as e:
                        logger.exception('Could not submit file "{}, '
                                         'error: {}".'.format(short_name, e))
                    file_index += 1
                    will_sleep = False

        if file_index == num_files:
            no_more_files = True

        if num_cached_reads > 0:
            # Process the data from the completed workers.
            summary_data = pipeline.process_cached_calls()
            for line in summary_data:
                line['filename'] = '{}.fast5'.format(line['label'])
                num_complete += 1
                pbar.update(num_complete)
                if 'error_message' in line:
                    logger.warning('Error processing '
                                   'file "{}".'.format(line['filename']))
                    logger.warning(line['error_message'])
                else:
                    summary.write(line)
                    logger.info('Finished processing '
                                'file "{}".'.format(line['filename']))
            will_sleep = False

        if will_sleep:
            # No workers are ready, and we haven't done any work. Take a nap.
            time.sleep(0.01)

    pbar.finish()


def this_folder():
    """ We need this in every script as the scripts may not live anywhere near
    the albacore package itself.

    """
    if getattr(sys, 'frozen', False):
        # The application is frozen
        return os.path.dirname(sys.executable)
    else:
        # The application is not frozen
        return os.path.dirname(__file__)


if __name__ == '__main__':
    # Check the data path.
    default_path = get_default_path(this_folder())

    # Check input parameters.
    list_parser = ArgumentParser(description='ONT Albacore Sequencing Pipeline '
                                             'Software', add_help=False)
    list_parser.add_argument('-l', '--list_workflows', action='store_true',
                             default=False,
                             help='List standard flowcell / kit combinations.')
    list_parser.add_argument('-v', '--version', action='version',
                             help='Print the software version.',
                             version='ONT Albacore Sequencing Pipeline '
                                     'Software (version {})'.format(VERSION), )
    opts, _ = list_parser.parse_known_args()

    if opts.list_workflows:
        print_workflows(default_path, WORKFLOW_LINEAR)

    parser = ArgumentParser(description='ONT Albacore Sequencing Pipeline '
                                        'Software', parents=[list_parser])
    parser.add_argument('-i', '--input', required=False, default=None,
                        help='Folder containing read fast5 files '
                             '(if not present, will expect file names on stdin).')
    parser.add_argument('-t', '--worker_threads', required=True, type=int,
                        help='Number of worker threads to use.')
    parser.add_argument('-s', '--save_path', required=True,
                        help='Path to save output.')
    parser.add_argument('--resume', metavar='X', nargs='?', const='failed', default='',
                        help='Resume previous run for the given save path. '
                             'Optional parameter X is for debugging purposes only.')
    parser.add_argument('-f', '--flowcell', required=False, default=None,
                        help='Flowcell used during the sequencing run.')
    parser.add_argument('-k', '--kit', required=False, default=None,
                        help='Kit used during the sequencing run.')
    parser.add_argument('--barcoding', required=False, default=False,
                        action='store_true', help='Search for barcodes to '
                                                  'demultiplex sequencing data.')
    parser.add_argument('-a', '--align_ref', required=False, default=None, help=SUPPRESS)
    # 'name of fasta genome to align against (will call bwa index if required).'
    parser.add_argument('-c', '--config', required=False, default=None,
                        help='Optional configuration file to use.')
    parser.add_argument('-d', '--data_path', required=False, default=None,
                        help='Optional path to model files.')
    parser.add_argument('-b', '--debug', action='store_true', default=False,
                        help='Output additional debug information to the log.')
    parser.add_argument('-r', '--recursive', action='store_true', default=False,
                        help='Recurse through subfolders for input data files.')
    parser.add_argument('-n', '--files_per_batch_folder', required=False,
                        type=int, default=4000,
                        help='Maximum number of files in each batch subfolder. '
                             'Set to 0 to disable batch subfolders.')
    parser.add_argument('-o', '--output_format', required=False, default='fastq',
                        help='desired output format, can be fastq,fast5 or only one of these.')
    parser.add_argument('-q', '--reads_per_fastq_batch', required=False,
                        type=int, default=4000,
                        help='number of reads per FastQ batch file. Set to 1 to'
                             ' receive one reads per file and file names which'
                             ' include the read ID. Set to 0 to have all reads'
                             ' per run ID written to one file.')
    parser.add_argument('--disable_filtering', action='store_true', default=False,
                        help='Disable filtering into pass/fail folders')
    parser.add_argument('--disable_pings', action='store_true', default=False,
                        help='Do not send summary information about the run')

    opts, unknown_args = parser.parse_known_args()
    opts.logfile = 'pipeline.log'
    opts.summfile = 'sequencing_summary.txt'
    opts.cfgfile = 'configuration.cfg'
    if opts.resume:
        resume_options_1d(opts, unknown_args) # if no log file or no summary file is found,
                                              # opts.resume will be reset to False
                                              # If config, log and summary file are found, unknown_args
                                              # will be cleared
    # We either need a config file or both a flowcell plus a kit
    if not opts.config and (not opts.flowcell or not opts.kit):
        parser.error('Either both flowcell and kit or a config file must be '
                     'specified.')
    opts.default_path = default_path
    config_override = parse_config_override(unknown_args)

    sys.exit(main(opts, config_override))
