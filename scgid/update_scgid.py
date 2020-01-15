import subprocess
import os
import inspect
import sys
import signal
from scgid.module import Module
from scgid.modcomm import LoggingEntity, ErrorHandler, Head, pkgloc
from scgid.library import subprocessP, subprocessC, output_cols
from scgid.config import FileConfig

class SCGIDUpdate(Module, LoggingEntity, ErrorHandler, Head):
    def __init__(self, is_automated_update=False):
        self.HOME, self.SCRIPTS = pkgloc()
        self.url = "https://www.github.com/amsesk/SCGid.git"
        self.local_branch = "dev"
        self.remote_branch = f"origin/{self.local_branch}"
        self.is_automated_update = is_automated_update

    def is_updatable (self) -> bool:

        # Fetch origin to confirm updated or not
        fetch = ['git','fetch','--all']
        subprocess.call(fetch)

        local_tag, _ = subprocess.Popen(["git", "rev-parse", self.local_branch], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        remote_tag, _ = subprocess.Popen(["git", "rev-parse", self.remote_branch], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        local_tag = local_tag.strip()
        remote_tag = remote_tag.strip()

        if local_tag == remote_tag:
            return False
        else:
            return True

    def ask_user_to_update(self):

        # Wait 60 seconds for response below
        TIMEOUT = 120

        def timeout_handler():
            pass
            #self.logger.info(f"No response from user in {TIMEOUT} seconds. Continuing without update.")
    
        signal.signal(signal.SIGALRM, timeout_handler)

        response = None
        signal.alarm(TIMEOUT)
        try:
            response = input(f"{output_cols['GREEN']}Your version of SCGid is behind the official repository. Would you like to update SCGid now? {output_cols['RESET']}[y/n] ")

            while True:
                if response.lower() == 'y':

                    # Cancel alarm
                    signal.alarm(0)

                    # We're going to update
                    return True

                elif response.lower() == 'n':

                    # Cancel alarm
                    signal.alarm(0)

                    # We're not going to update
                    return False

                else:

                    # Wait for an accepted response
                    continue

        except:

            print("TIMEOUT")

            # Return False (i.e., Do Not Update) if we don't receive input within 120 seconds - the user is probably not running this interactively
            return False

    def update_scgid(self):
        
        # Load current config.yaml to preserve configuration settings after update
        config = FileConfig()
        config.load_yaml()

        os.chdir(os.getenv("HOME"))
        clone = ['git', 'clone', url]
        subprocess.call(clone)

        os.chdir("SCGid")
        install = ['python', 'setup.py', 'install']
        subprocess.call(install)

        '''
        ## Hard reset from origin
        fetch = ['git','fetch','--all']
        reset = ['git','reset','--hard',self.remote_branch]

        subprocess.call(fetch)
        subprocess.call(reset)
        '''

        config.write_yaml()

    def run(self):
        self.start_logging()

        here = os.getcwd()
        os.chdir(self.HOME)

        if not self.is_updatable():
            pass
            #self.logger.info(f"SCGid is already up to date with {self.remote_branch}.")
        else:
            if self.is_automated_update:
                conduct_update = self.ask_user_to_update()
                if conduct_update:
                    self.update_scgid()
            else:
                print("You asked for this so no need to ask to update.")
                self.update_scgid()

        os.chdir(here)

            


        ## A change! ##
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
