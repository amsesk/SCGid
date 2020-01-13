import subprocess
import os
import inspect
import sys
from scgid.module import Module
from scgid.modcomm import LoggingEntity, ErrorHandler, Head, pkgloc
from scgid.library import subprocessP, subprocessC

class SCGIDUpdate(Module, LoggingEntity, ErrorHandler, Head):
    def __init__(self):
        self.HOME, self.SCRIPTS = pkgloc()
        self.url = "https://www.github.com/amsesk/SCGid.git"
        self.local_branch = "dev"
        self.remote_branch = f"origin/{self.local_branch}"

    def check(self):
        pass

    def run(self):
        self.start_logging()

        self.logger.info(f"Entering SCGid package directory at `{self.HOME}`")
        os.chdir(self.HOME)

        p = subprocess.Popen(["git","remote","update"], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        local_tag, _ = subprocess.Popen(["git", "rev-parse", self.local_branch], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        remote_tag, _ = subprocess.Popen(["git", "rev-parse", self.remote_branch], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        local_tag = local_tag.strip()
        remote_tag = local_tag.strip()

        if local_tag == remote_tag:
            self.logger.info(f"SCGid is already up to date with {self.remote_branch}.")
        else:
            self.logger.info("SCGid can be updated.")

        ## A change! ##
        ## Another one ###
'''
## some variables
scgid_bin = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
scgid_home = os.path.dirname(scgid_bin)

## Confirm
while True:
    entry = raw_input("You are about to overwrite the scgid source code with updates from GitHub. Continue? [y/n]")
    entry = entry.lower()
    if entry not in ['y','n']:
        continue
    else:
        if entry == 'y':
            break
        else:
            sys.exit()

## navigate to local repo
os.chdir(scgid_home)
## save settings for hard reset from origin
os.rename(os.path.join(scgid_bin,'settings.py'), os.path.join(scgid_bin,'settings.py.sav'))

## Hard reset from origin
fetch = ['git','fetch','--all']
reset = ['git','reset','--hard','origin/master']

subprocess.call(fetch)
subprocess.call(reset)

## check if settings have changed. If not, restore old settings
with open(os.path.join(scgid_bin, 'settings.py.sav')) as old, open(os.path.join(scgid_bin, 'settings.py')) as new:
    new_names = [x.split("=")[0] for x  in new.readlines()]
    old_names = [x.split("=")[0] for x in old.readlines()]
    print (old_names)
    print (new_names)
os.rename(os.path.join(scgid_bin,'settings.py.sav'), os.path.join(scgid_bin,'settings.py'))
'''
