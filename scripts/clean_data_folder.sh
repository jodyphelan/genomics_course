#! /usr/bin/env python
import sys
import subprocess
import os
import shutil

tokeep = [l.strip() for l in open("/usr/local/src/files.txt")]

for l in subprocess.Popen("find /home/user/data/",shell=True,stdout=subprocess.PIPE).stdout:
    f = l.decode().strip()
    if f not in tokeep:
        if os.path.isfile(f) or os.path.isdir(f):
            try:
                os.remove(f)
            except:
                shutil.rmtree(f)
