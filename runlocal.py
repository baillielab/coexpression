#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
''' run a pre-prepared coexpresion analysis '''

import sys, subprocess, threading, copy
from coexfunctions import *
#----------------------------------------
import argparse
# ----------------------------------------
import argparse
import copy
import subprocess
import sys
import threading

from coexfunctions import *

parser = argparse.ArgumentParser()

# essential
parser.add_argument('-d', '--working_files_dir',default='unset', help='input directory')
args = parser.parse_args()
module = sys.modules[__name__]
for name,value in vars(args).iteritems():
setattr(module, name, value)
#----------------------------------------
storage_dir_perm_results = os.path.join(working_files_dir, "perm_results_store")
jobfile = os.path.join(working_files_dir, "jobfile.txt")
collationcommandfile = os.path.join(working_files_dir, "collationcommandfile.txt")
#----------------------------------------
storedesettings=readsettings(working_files_dir)
numperms = storedesettings['numperms']
verbose = storedesettings['verbose']
#----------------------------------------
max_simultaneous = 20
num_started = 0
#----------------------------------------

def runcollation():
f=open(collationcommandfile)
cmd = string.strip(f.readlines()[0])
f.close()
print cmd
subprocess.call(cmd, shell=True) # execute a child program in a new process
print "done"

def getcompleted():
completedperms = sorted([int(string.split(x,'.')[-1]) for x in os.listdir(storage_dir_perm_results) if string.split(x,'.')[-2]=='indmeasure'])
return completedperms

def startjob(jobid):
global started
if verbose: print "starting job id:{}".format(jobid)
subprocess.Popen(commanddict[jobid], shell=True) # execute a child program in a new process
started.append(x)

def progressupdate():
print "numperms:{} len(completed):{} num_started:{} yet to start:{} running:{}".format(numperms, len(completed), len(started), len(yet_to_do), len(running))
print "running:", running
print "yet to start:", yet_to_do


def checkforcompletion():
global started, yet_to_do, running
completed   = getcompleted()
yet_to_do   = sorted([x for x in commanddict.keys() if x not in completed and x not in started])
running = list(set(started) - (set(completed) | set(already_done)))
if verbose: progressupdate()
if len(completed)>=numperms+1: # the job is finished. collate results
print "trigger fired: collation to begin now"
try:
t.cancel()
except:
pass #in case t hasn't started yet
runcollation()
else:
if len(running) < max_simultaneous and len(yet_to_do)>0:
startjob(yet_to_do[0])
t= threading.Timer(30, checkforcompletion).start()


f=open(jobfile)
localcommands = [string.strip(x) for x in f.readlines()]
f.close()
commanddict = {int(string.split(readoptions(x)['-mf'],'.')[-2]):x for x in localcommands}
#-----------------------
completed = getcompleted()
already_done = copy.copy(completed)
yet_to_do = sorted([x for x in commanddict.keys() if x not in completed])
started = []
running =[]
#-----------------------
print "already done:", completed
print "yet to do:", yet_to_do
for i,x in enumerate(yet_to_do[:max_simultaneous]):
startjob(x)
checkforcompletion()













