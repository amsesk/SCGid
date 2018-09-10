import subprocess
import os
import inspect

## some variables
scgid_bin = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
scgid_home = os.path.dirname(scgid_bin)

## navigate to repo
os.chdir(scgid_home)
print os.getcwd()
## save settings for hard reset from origin
os.rename(os.path.join(scgid_bin,'settings.py'), os.path.join(scgid_bin,'settings.py.sav'))

## Hard reset from origin
fetch = ['git','fetch','--all']
reset = ['git','reset','--hard','origin/master']

subprocess.call(fetch)
subprocess.call(reset)

## restore settings
os.rename(os.path.join(scgid_bin,'settings.py.sav'), os.path.join(scgid_bin,'settings.py'))
