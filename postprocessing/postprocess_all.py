''' RUN ALL POST PROCESSING SCRIPTS OVER A SINGLE RESULTS WORKING DIRECTORY '''

import os
import subprocess
import ConfigParser
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-wd', '--workingdirectory',    help='full coexpression results directory')
parser.add_argument('-cc', '--config_file', default=os.path.join(scriptpath, '../app.cfg'), help='use a custom config file instead of app.cfg')
args = parser.parse_args()
#-----------------------------
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read(args.config_file)
pathtopython = config.get('directorypaths', 'pathtopython')
#-----------------------------
# run all of the following scripts on this directory
subprocess.call("{} {} -wd {}".format(pathtopython, os.path.join(scriptpath, "make_json.py"), args.workingdirectory), shell=True)