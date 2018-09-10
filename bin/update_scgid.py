import subprocess
import os
import inspect
import sys

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

## restore settings
os.rename(os.path.join(scgid_bin,'settings.py.sav'), os.path.join(scgid_bin,'settings.py'))
