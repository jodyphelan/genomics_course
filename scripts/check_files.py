#! /usr/bin/env python
import sys
import subprocess
import os
import shutil

for l in open("/usr/local/src/files.txt"):
    f = l.strip()
    if not os.path.isfile(f) and not os.path.isdir(f):
        print(f)

