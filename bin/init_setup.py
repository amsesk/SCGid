#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 10:25:31 2018

@author: kevinamses
"""

### This script gets called after download and extraction of scgid-[version].tar.gz
### Automatically set-up settings.py and download databases

#%% imports
import os
import inspect
import sys
import re
import subprocess
import time
import ftplib
import urllib
import threading
import imp # to check for python-dependencies
from initlib import report_outcome, is_fasta, file_grep, ftp_retr_and_report, ftp_retr_progress, output_cols, CURSOR_UP_ONE, ERASE_LINE, spinner, ow_last_stdout

#%% some variables
scgid_bin = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
scgid_home = os.path.dirname(scgid_bin)
scgid_databases = os.path.abspath(os.path.join(scgid_home,"database"))
version = "0.1b"

#%% header

#print open("logo.txt").read()
print """
 ___   ___  ___  ____  ____  
/ __) / __)/ __)(_  _)(  _ \ 
\__ \( (__( (_-. _)(_  )(_) )
(___/ \___)\___/(____)(____/ 
"""
print "----- Welcome to scgid, version %s -----" % (version)
print "> scgid is provided under the GNU General Public License v3.0\n"
print "> This script is going to help you get scgid ready to run by defining some package-wide variables and downloading the necessary databases.\n"
if os.getcwd() != scgid_bin:
    os.chdir(scgid_bin)
    print "> Navigating to scgid bin directory at %s ...\n" % (scgid_bin)

#%% check python2.7 dependencies
print "> Scanning scgid python2.7 dependencies now..."
pydep = {
        'numpy': False,
        'pandas': False,
        'ete3': False,
        }
for pkg, check in pydep.iteritems():
    try:
        imp.find_module(pkg)
        pydep[pkg] = True
    except:
        while True:
            print "> Python module \'%s\' not installed..." % (pkg)
            entry = raw_input("Try installing %s with conda or pip? If in doubt, try conda first. [conda/pip/no] " %pkg)
            if entry not in ['conda','pip','no']:
                print "> " + output_cols['RED'] + "ERROR: Please enter only 'conda' or 'pip' or 'no'"+output_cols['RESET']
                continue
            else:
                if entry == 'conda':
                    try_install = subprocess.Popen(['conda','install',pkg], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                elif entry == 'pip':
                    subprocess.call(['pip', 'install', '--user', pkg])
                else:
                    print "> OK, please maually install %s." % (pkg)
                    break
            
                ## If you made it here, check to see if install worked and pkg is now available
                try:
                    imp.find_module(pkg)
                    pydep[pkg] = True
                    break
                except:
                    print "> " + outputcols['RED'] + "Something went wrong with automatic installation of %s. Please manually install it." % (pkg)+outputcols['RESET']
                    break
                
print "\npython2 Dependences\n"+"-"*24
for pkg, check in pydep.iteritems():
    if not check:
        print pkg+":\t"+output_cols["RED"]+"[NOT INSTALLED]"+output_cols["RESET"]
    else:
        print pkg+":\t"+output_cols["GREEN"]+"[INSTALLED]"+output_cols["RESET"]
print "-"*24
if not any (pydep.values()):
    print "\n> Finish installing python dependencies and then rerun to continue. Exitting...\n"
    sys.exit(1)
else:
    print "\n> All python dependencies installed and good to go. Continuing...\n"
                
#%% settings.py ###
settings = {} # variable: value
print "> Let's define the locations of your installed third-party dependencies"

#%% Locate third-party dependencies and add their paths to settings dict ##
files_to_find = {
        "esom_path": {
            'target': os.path.join("bin","esomstart"),
            'question': "Enter full path to ESOM install: "
            },
        "clams_path": {
            'target': "ClaMS-CLI.jar",
            'question': "Enter full path to folder containing ClaMS-CLI.jar: "
            },
        "path_to_Rscript": {
            'target': "Rscript",
            'question': "Enter full path to folder containing Rscript: "
            }
        }
for var, stuff in files_to_find.iteritems():
    while True:
        entry = raw_input(stuff['question'])
        if not os.path.isfile(os.path.join(entry, stuff['target'])):
            report_outcome (stuff['question'], "RED", "[NOT FOUND, TRY AGAIN]")
        else:
            report_outcome (stuff['question'], "GREEN", "[GOOD]")
            files_to_find[var] = entry
            break

settings.update(files_to_find)

#%% Get environmental variables, not really necessary for augustus config path? ##

settings['AUGUSTS_CONFIG_PATH'] = ""

#%% Figure out the swissprot protein database ##
print "\n> scgid is currently only compatible with swissprot-style protein databases. See README for more information on this."
while True:
    entry = raw_input("Download the current version of the swissprot database? [y/n] ")
    entry = entry.strip().lower()
    if entry not in ['y','n']:
        continue
    else:
        if entry == 'y': 
            while True:
                ## download database and update settings with version and path (both handled by update_swissprot.py
                entry = raw_input("Specify protein database download directory: [default = '{}'] ".format(scgid_databases))
                if len(entry) != 0:
                    if not os.path.isdir(entry):
                        print "> Directory does not exist. Try again."
                        continue
                    else:
                        scgid_databases = entry
                p = subprocess.call([sys.executable, 'update_swissprot.py', scgid_databases])
                ## get spdb_version and path_to_spdb from settings.py (find a less sloppy way to do this soon...)
                settings['path_to_spdb'] = file_grep ("path_to_spdb=", os.path.join(scgid_bin,"settings.py"),"first").strip().split("=")[1]
                settings['spdb_version'] = file_grep ("spdb_version=", os.path.join(scgid_bin,"settings.py"),"first").strip().split("=")[1]
                settings['path_to_spdb'] = settings['path_to_spdb'].replace('"','')
                settings['spdb_version'] = settings['spdb_version'].replace('"','')
                break

            break
        else: #entry.lower() is 'n'
            while True:
                question = "Path to your database (FASTA): "
                u_db = raw_input(question)
                
                ## Verify the existence and swissprot-like formatting of user-supplied database ##
                if not os.path.isfile(u_db):
                    report_outcome(question, "RED", "[NOT FOUND, TRY AGAIN]")
                else:
                    if not is_fasta(u_db):
                         report_outcome(question, "RED", "[SPECIFIED DB NOT IN FASTA FORMAT]")
                    else:
                        all_with_OS = True
                        for header in file_grep ("^>", u_db, 'multiple'):
                            if "OS=" not in header:
                                report_outcome(question, "RED", "[HEADERS NOT SWISSPROT-STYLE, see README]")
                                all_with_OS = False
                                break
                        if not all_with_OS:
                            continue
                        else:
                            report_outcome(question, "GREEN", "[GOOD]")
                            break
            ### add location of user-defined swissprot-style database to settings
            settings['path_to_spdb'] = u_db
            break

#%% Figure out the swissprot taxonomy database ##
#path_to_taxdb = "/home/aimzez/scgid/database/taxonomy-all-070318.tab" #local path
dwnld_link = "https://www.uniprot.org/taxonomy/?query=*&format=tab"
date = time.strftime("%m%d%y")
settings['taxonomy_all_tab'] = None

print "\n> scgid requires a copy of the swissprot taxonomy database"
entry = raw_input("Download the current version of the swissprot taxonomy database? [y/n] ")
entry = entry.strip().lower()
while True:    
    if entry not in ['y','n']:
        continue
    else:
        if entry == 'y':
            settings['taxonomy_all_tab'] = "uniprot-taxonomy-all-{}.tab".format(date)
            dwnld_to = os.path.join(scgid_databases, settings['taxonomy_all_tab'])
            if os.path.isfile(dwnld_to):
                os.remove(dwnld_to)
            print "\n> Downloading uniprot taxonomy database to: {}...\n".format(dwnld_to)
            spin = spinner()
            t = threading.Thread(target=urllib.urlretrieve, args=(dwnld_link, dwnld_to))
            t.start()
            while t.is_alive():
                print next(spin)
                t.join(0.2)
            print CURSOR_UP_ONE,ERASE_LINE  #erase last spinner char
            break 
        else: 
            print "> OK, you'll have to download and build the taxonomy database before using scgid. See README for more information."
            break

#%% Build taxonomy database
if settings['taxonomy_all_tab'] is not None:
    entry = raw_input("Build a taxonomy database for newly downloaded swissprot-style database? [y/n] ")
    while True:
        if entry.lower() not in ['y','n']:
            continue
        if entry.lower() == 'y':
            subprocess.call([sys.executable, os.path.join(scgid_bin,'build_taxdb.py'), '-db', settings['path_to_spdb'], '-t', '/home/aimzez/scgid/database/uniprot-taxonomy-all-090718.tab', '-o', '%s.taxdb' % (settings['path_to_spdb'])])
            settings['path_to_taxdb'] = settings['path_to_spdb']+'.taxdb'
            break
        else:
            settings['path_to_taxdb'] = None
            print "> You'll have to build the taxonomy datbase manually before using scgid. Try `scgid buildtaxdb [args...]"
            break
#%% write settings.py
print "\n> Writing your settings to %s" % (os.path.join(scgid_bin, "settings.py"))
print "\nYOUR SCGID SETTINGS\n" + "-"*24
with open('settings.py', 'w') as s:
    for var, value in settings.iteritems():
        print '%s = "%s"' % (var,value)
        s.write('%s="%s"\n' % (var,value))
print "-"*24+"\n"

#%% check $PATH
if scgid_bin not in os.environ['PATH']:
    print "> WARNING: scgid does not appear to be present in your PATH... \n\tAdd 'export PATH=$PATH:%s' to your .bashrc (or other file for your particular console) to add scgid to PATH" % (scgid_bin)
        
#%% DONE!     
print "\n" + output_cols["GREEN"] + "[ scgid init completed... happy squid-ing! ]"
