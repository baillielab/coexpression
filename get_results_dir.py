#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

'''Quick hack to get (generated) name of results directory
   - call with same arguments as when running 0-prepare-coex-py'''

import ConfigParser
import inspect
import os
import string
import sys

# -------------------------------------------------------------------------------
pathtocoexapp = inspect.stack()[0][1]
coexappdir = os.path.split(pathtocoexapp)[0]
config = ConfigParser.RawConfigParser()
config.read(os.path.join(coexappdir, 'app.cfg'))
sourcefilesdir = config.get('directorypaths', 'sourcefilesdir')
resdir = config.get('directorypaths', 'resdir')
pathtobedtools = config.get('directorypaths', 'pathtobedtools')
f5resource = config.get('directorypaths', 'f5resource')
# -------------------------------------------------------------------------------
if not os.path.isabs(sourcefilesdir):  #  then the user has probably made a relative path to the script.
    sourcefilesdir = os.path.join(coexappdir, sourcefilesdir)
if not os.path.isabs(resdir):  #  then the user has probably made a relative path to the script.
    resdir = os.path.join(coexappdir, resdir)
if not os.path.isabs(pathtobedtools):  #  then the user has probably made a relative path to the script.
    pathtobedtools = os.path.join(coexappdir, pathtobedtools)
# -------------------------------------------------------------------------------
import argparse

parser = argparse.ArgumentParser()

# essential
parser.add_argument('-f', '--snp_details_file', default='unset', help='input filepath')

# speed up for large datsets
parser.add_argument('-n', '--numperms', type=int, default=100, help='window size')
parser.add_argument('-p', '--permutationmode', type=str, choices=['circular', 'post', 'backcirc', 'background'],
                    default='circular', help='permutation mode')
parser.add_argument('-s', '--nsp', type=float, default=1.0,
                    help='node-specific pvalue in use; specify precision (1=perfect, 0=globaldistribution is used.)')
parser.add_argument('-t', '--numthreads', type=int, default=1, help='numthreads used')
parser.add_argument('-w', '--window', type=int, default=0, help='window size')

# useful
parser.add_argument('-b', '--backgroundfile', default='unset', help='backgroundfile')
parser.add_argument('-q', '--specialbedfile', default="not_specified", help='specialbedfile')
parser.add_argument('-x', '--specialexpressionfile', default="not_specified", help='specialexpressionfile')

# boring
parser.add_argument('-e', '--useremail', default="no_email", help='user email')
parser.add_argument('-y', '--version_requested', default="null", help='version of the software the user is expecting')

# occasional
parser.add_argument('-po', '--prepareonly', action="store_true", default=False, help='Stop after setup script')
parser.add_argument('-l', '--seedfile', default="none", help='seed file')
parser.add_argument('-v', '--verbose', action="store_true", default=False, help='increases verbosity')
parser.add_argument('-j', '--pval_to_join', type=float, default=0.1, help='pval_to_join')
parser.add_argument('-u', '--useaverage', action="store_true", default=False,
                    help='use the average coexpression score (eliminates position enrichment signal)')
parser.add_argument('-z', '--doitagain', action="store_true", default=False, help='create a whole new permutations set')
parser.add_argument('-cm', '--correlationmeasure', type=str, choices=['Spearman', 'Pearson'], default='Spearman',
                    help='filepath')

# obselete
parser.add_argument('-a', '--anticorrelation', action="store_true", default=False, help='anticorrelation')
parser.add_argument('-i', '--noiter', action="store_true", default=False, help='noiter')
parser.add_argument('-r', '--recordedges', action="store_true", default=False, help='recordedges')
parser.add_argument('-m', '--selectionthreshold_is_j', action="store_true", default=False,
                    help='match selectionthreshold_is_j')
parser.add_argument('-st', '--selectionthreshold', type=float, default=0, help='selectionthreshold')

args = parser.parse_args()
module = sys.modules[__name__]
for name, value in vars(args).iteritems():
    setattr(module, name, value)
# -------------------------------------------------------------------------------
nsp = float(nsp)
# -------------------------------------------------------------------------------
if snp_details_file == 'unset':
    snp_details_file = os.path.join(sourcefilesdir, 'test.bed')
common_label = os.path.split(snp_details_file)[1].replace(".bed", "")
# -------------------------------------------------------------------------------
# -###  SETTINGS
if selectionthreshold_is_j:
    selectionthreshold = -math.log(pval_to_join, 10)
data_start_col = 1  # nb zero  ex
data_start_row = 1  # 1084
chr_col = 1  # -address of each locus is in expression file in order to apply min_separation
start_col = 2
end_col = 3
gpl_len = 1000000  # number of correlations to include in the global correlation lists. #ONLY RELEVANT IF GENERATING NEW LIST!
soft_separation = 100000  # correlating promoters in this range will be joined. #IF pval_to_join >1 then this is a hard separation
progress_interval = 50
ramlimit = 200000000
# -------------------------------------------------------------------------------
# supplementary files
chrom_lengths_file = os.path.join(sourcefilesdir, "chromlengths.txt")
promoter_details = os.path.join(sourcefilesdir, "METADATA_U22_ENHANCERS")
feature_coordinates = os.path.join(sourcefilesdir, "f5ep300_100.bed")
expression_file = os.path.join(sourcefilesdir, "f5ep_ptt_condense.expression")
windowbed = os.path.join(pathtobedtools, "windowBed")
# schedulefiles 
schedulefilesdir = os.path.join(coexappdir, "schedulefiles")
recordfile = os.path.join(schedulefilesdir, 'already_run.sch')
alivefile = os.path.join(schedulefilesdir, 'alive.sch')
restartfile = os.path.join(schedulefilesdir, 'start.sch')
# correlation files
expfilename = string.split(expression_file, "/")[-1]
if correlationmeasure == "Pearson":
    correlationfile = "%s%s.correlationlist" % (sourcefilesdir, expfilename)
elif correlationmeasure == "Spearman":
    correlationfile = "%s%s.spearmanlist" % (sourcefilesdir, expfilename)
# -------------------------------------------------------------------------------
output_directory_label = common_label
if nsp == 1:
    output_directory_label += "_complete"
else:
    output_directory_label += "_s%s" % nsp
if permutationmode == 'post':
    output_directory_label += "_POST"
elif permutationmode == 'backcirc':
    output_directory_label += "_BACKCIRC"
elif permutationmode == 'circular':
    output_directory_label += "_CIRCULAR"
elif permutationmode == 'background':
    output_directory_label += "_BACKGROUND"
if selectionthreshold_is_j: output_directory_label += "_SIJ"
if anticorrelation: output_directory_label += "_ANTIin"
if seedfile != "none": output_directory_label += "_seed"
if noiter: output_directory_label += "_NOiter"
if useaverage: output_directory_label += "_useav"
if window > 0: output_directory_label += "_w%s" % window
output_directory_label += "_pj%s" % pval_to_join
expstring = string.split(specialexpressionfile, "/")[-1].replace('.expression', '')
if specialexpressionfile != "not_specified": output_directory_label += "_%s" % expstring
# -------------------------------------------------------------------------------
#  working files
working_files_dir = "%s%s" % (resdir, output_directory_label)

print working_files_dir
