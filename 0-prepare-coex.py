#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import copy
# -# get options
import datetime
import random
import resource
import subprocess
import sys
import time
import timeit
from ctypes import *
from subprocess import PIPE, STDOUT

from coexfunctions import *

version = getversion()
# -------------------------------------------------------------------------------
import inspect
import ConfigParser

pathtocoexapp = inspect.stack()[0][1]
coexappdir = os.path.split(pathtocoexapp)[0]
config = ConfigParser.RawConfigParser()
config.read(os.path.join(coexappdir, 'app.cfg'))
sourcefilesdir = config.get('directorypaths', 'sourcefilesdir')
resdir = config.get('directorypaths', 'resdir')
pathtobedtools = config.get('directorypaths', 'pathtobedtools')
f5resource = config.get('directorypaths', 'f5resource')
# -------------------------------------------------------------------------------
if not os.path.isabs(sourcefilesdir):  # then the user has probably made a relative path to the script.
    sourcefilesdir = os.path.join(coexappdir, sourcefilesdir)
if not os.path.isabs(resdir):  # then the user has probably made a relative path to the script.
    resdir = os.path.join(coexappdir, resdir)
if not os.path.isabs(pathtobedtools):  # then the user has probably made a relative path to the script.
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
    print "==reverting to test file because source file not set=="
    snp_details_file = os.path.join(sourcefilesdir, 'test.bed')
common_label = os.path.split(snp_details_file)[1].replace(".bed", "")
if verbose:
    print ("coexpression running", common_label, "permutationmode:", permutationmode, "numperms=", numperms, "email=",
           useremail, "nsp:", nsp, "seedfile", seedfile, "recordedges", recordedges, "doitagain", doitagain,
           "specialexpressionfile", specialexpressionfile)
# -------------------------------------------------------------------------------
tic = timeit.default_timer()
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
supplementary_label = 'coex'
# -------------------------------------------------------------------------------
os.environ['OMP_NUM_THREADS'] = str(numthreads)
# -------------------------------------------------------------------------------
# snp_details file format - BED FORMAT
snp_col = 3
chrom_col = 0  # chrom is an integer only
address_col = 1
snp_details_p_col = -1
# -------------------------------------------------------------------------------
# working files
working_files_dir = "%s%s" % (resdir, output_directory_label)
jobfile = os.path.join(working_files_dir, "jobfile.txt")
collationcommandfile = os.path.join(working_files_dir, "collationcommandfile.txt")
storage_dir_permutations = os.path.join(working_files_dir, "permutation_store")
thisqueue = os.path.join("%s/%s.queue" % (working_files_dir, supplementary_label))
settingsfile = os.path.join(working_files_dir, "settings.txt")
sofile = os.path.join(coexappdir, './coexpression_v2.so') # './' needed for cdll.LoadLibrary(sofile) below
if not os.path.exists(sofile):
    cmd = "gcc -shared  -O3 -fPIC -fopenmp %scoexpression_v2.c -o %scoexpression_v2.so" % (
        coexappdir, coexappdir)
    print "=== coexpression is not installed ==="
    try:
        os.system(cmd)
    except:
        print "please run this command:"
        print cmd
permstore = os.path.join(working_files_dir, "permstore_%s.txt" % (supplementary_label))
zerocountfile = os.path.join(working_files_dir, "zerocount_%s.txt" % (supplementary_label))
edgestore = os.path.join(working_files_dir, "edgestore_%s.txt" % (supplementary_label))
if seedfile != "none":
    seedfile = os.path.join(sourcefilesdir, seedfile)
# -------------------------------------------------------------------------------
for thisdir in [working_files_dir, storage_dir_permutations, resdir]:
    check_dir(thisdir)
t = open(thisqueue, "w")
t.write("Job started %s" % (str(datetime.datetime.now())))
t.close()
# -------------------------------------------------------------------------------
runcommand = "python \
0-run-coex.py -y %s -f %s -n %s -p %s -s %s -e %s -t %s -x %s -q %s -w %s" % (
version, snp_details_file, numperms, permutationmode, nsp, useremail, numthreads, specialexpressionfile, specialbedfile,
window)

runcommand += " -j %s" % (pval_to_join)
if verbose: runcommand += " -v"
if anticorrelation: runcommand += " -a"
if recordedges: runcommand += " -r"
if selectionthreshold_is_j: runcommand += " -m"
if useaverage: runcommand += " -u"
if verbose:
    print runcommand
# -------------------------------------------------------------------------------
if specialexpressionfile != "not_specified":
    expression_file = specialexpressionfile
    print "***"
    print "Using user-specified expressionfile:", expression_file
    print "***"
if specialbedfile != "not_specified":
    feature_coordinates = specialbedfile
    print "***"
    print "Using user-specified bedfile:", feature_coordinates
    print "***"
# -------------------------------------------------------------------------------
coexpression = cdll.LoadLibrary(sofile)
feature_label = string.split(string.split(feature_coordinates, "/")[-1], ".")[
    0]  # ie the "fantom5" in the coordinates file name
##################################################
# snps_mapped_filename=feature_label+"_"+common_label+".snps_mapped.bed"
# snps_mapped=os.path.join(storage_dir_permutations,snps_mapped_filename)
snps_mapped = os.path.join(storage_dir_permutations, feature_label + ".%s.bed" % (0))
statsfile = feature_label + "_" + supplementary_label + "_stats_snps.txt"
# -------------------------------------------------------------------------------
# -################################################################################
# - global objects used in functions.
globalcorrelationlist = []
nodespecificcorrelationlists = {}
shiftcount = 0
promindex_lookup = {}
alldata_c = []
result_all_c = []
dictlen = 0
entrylen = 0
exp_dict = {}

# -################################################################################

if verbose: print('verbose mode')
#seed random number generator to allow for testing
#random.seed(176202)

starttime = time.time()
# -######################
# read chromosome lengths
chrom_lengths = {}
genome_length = 0
f = open(chrom_lengths_file)
lines = [string.split(string.strip(x), ": ") for x in f.readlines()]
f.close()
for line in lines:
    chrom_lengths[line[0]] = int(line[1])
    genome_length += int(line[1])
print "====="
print "genome read, length=", genome_length
print "====="
# read snps
snp_bed_data = {}
snp_details = {}
snp_p_values = {}
all_snps = []
badlines = []
if verbose: print ("reading snp details (bed) file", snp_details_file)
f = open(snp_details_file)
lines = f.readlines()
f.close()
snp_count = 0
for line in lines:
    if line.startswith("#"):
        continue  # skipcomments
    if "\t" in line:
        line = string.split(string.strip(line), "\t")
    elif " " in line:
        line = string.split(string.strip(line), " ")
    else:
        line = string.split(string.strip(line))
    try:
        address = int(line[address_col])
    except:
        if verbose: print ("skipping line of snp_details file", line)
        try:
            badlines.append(line[0])
        except:
            pass
        continue
    if len(line) > 1:
        snp = line[snp_col]
        chrom = "chr%s" % line[chrom_col].replace("chr", "")
        if chrom not in chrom_lengths:
            if verbose: print ("skipping", chrom, "as length not known")
            continue  # skip chromosomes for which we don't know the length
        if snp not in all_snps:
            all_snps.append(snp)  # one big list of all snps included in this run
        try:
            snp_bed_data[chrom][address] = snp
        except:
            snp_bed_data[chrom] = {}
            snp_bed_data[chrom][address] = snp
        try:
            thisp = float(line[snp_details_p_col])  # 4th col of snp_details file
        except:
            thisp = 1
        try:
            snp_p_values[line[snp_col]]
            snp_p_values[line[snp_col]] = min(snp_p_values[line[snp_col]], thisp)
            if verbose: print ("snp p value already in dict. set to min: %s" % snp_p_values[line[snp_col]])
        except:
            snp_p_values[line[snp_col]] = thisp
        snp_details[snp_col] = [line[chr_col], line[address_col], line[snp_details_p_col]]
    else:
        if verbose: print ("line too short(", len(line), "):", line)
snp_count = len(all_snps)
if len(snp_bed_data) == 0:
    if verbose: print ("No snps in .snp_details file")
else:
    if verbose: print ("snps read successfully")
# -------------------------------------------------------------------------------
settingsdict = {
    "version": version,
    "numperms": numperms,
    "useremail": useremail,
    "seedfile": seedfile,
    "label": supplementary_label,
    "permutationmode": permutationmode,
    "verbose": verbose,
    "anticorrelation": anticorrelation,
    "nsp": nsp,
    "window": window,
    "gpl_len": gpl_len,
    "noiter": noiter,
    "selectionthreshold": selectionthreshold,
    "selectionthreshold_is_j": selectionthreshold_is_j,
    "pval_to_join": pval_to_join,
    "soft_separation": soft_separation,
    "progress_interval": progress_interval,
    "correlationmeasure": correlationmeasure,
    "data_start_col": data_start_col,
    "data_start_row": data_start_row,
    "chr_col": chr_col,
    "start_col": start_col,
    "end_col": end_col,
    "expression_file": expression_file,
    "feature_coordinates": feature_coordinates,
    "correlationfile": correlationfile,
    "supplementary_label": supplementary_label,
    "useaverage": useaverage,
    "genome_length": genome_length,
    "snp_count": snp_count,
    "runcommand": runcommand,
}
o = open(settingsfile, "w")
for key, value in settingsdict.iteritems():
    o.write("%s\t%s\n" % (key, value))
o.close()
# -------------------------------------------------------------------------------


# -#########################################################################################
# -#########################################################################################
# -######### ##########
# -#########START OF MAPPING ##########
# -#########  ##########
# -#########################################################################################
# -#########################################################################################
# -#########################################################################################

# --------------------------------------------------------
# DEFINE GENOME-WIDE ADDRESSES FOR TRANSLATION FUNCTIONS
chromrange = {}
chromnums = range(1, 23) + ['X', "Y", "M"]
chromlist = []
start = 0
for chromnum in chromnums:
    chrom = "chr%s" % chromnum
    chromlist.append(chrom)
    thislen = chrom_lengths[chrom]
    end = start + thislen
    chromrange[chrom] = [start, end]
    start += thislen
# --------------------------------------------------------
# feature_coordinates is already bed file from v28 onwards...
feature_bed = feature_coordinates
# -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# CHECK WHETHER WORTH CONTINUING...
temp_snp_bed = working_files_dir + "/temporary.bed"
o = open(temp_snp_bed, "w")
for chrom in snp_bed_data:
    for address in snp_bed_data[chrom]:
        snpname = snp_bed_data[chrom][address]
        try:
            o.write("%s\t%s\t%s\t%s\n" % (chrom.replace("chr", ""), int(address) - 1, address, snpname))
        except:
            if verbose: print ("error writing bed file:")
o.close()
command = windowbed + " -w %s" % (window) + " -a " + temp_snp_bed + " -b " + feature_bed + " > " + snps_mapped
if verbose: print (command)
p = subprocess.Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()
if verbose: print (output)  # so that errors are detected in progress file
if "Permission denied" in output:
    chmod_cmd = 'chmod 777 {}'.format(windowbed)
    print chmod_cmd
    p = subprocess.Popen(chmod_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    print p.stdout.read()
    print command
    p = subprocess.Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    print p.stdout.read()

f = open(snps_mapped)
lines = f.readlines()
f.close()
numhits = len(lines)
if numhits < 2 and len(badlines) < 2:
    give_up("not enough snps mapped - could there be a chr in the bedfile?")

# -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
# RANDOMISATION PRE-MAPPING
shifts = []
mapfiles = {0: snps_mapped}
if numhits > 1:  # CAN ONLY PROCEED IF THERE ARE ENOUGH SNPS TO BE WORTH MAPPING
    if permutationmode == 'circular':
        if verbose: print ("CIRCULAR randomisation before mapping")
        for i in range(numperms):
            shifts.append(int(random.random() * 10000000000))
            permfile = os.path.join(storage_dir_permutations, feature_label + ".%s.bed" % (i + 1))
            mapfiles[shifts[i]] = permfile
    elif permutationmode == 'background':
        backgroundname = os.path.split(backgroundfile)[-1]
        if verbose: print ("BACKGROUND randomisation before mapping")
        for i in range(numperms):
            shifts.append(i + 1)  # so that 0 is reserved for real data
            permfile = os.path.join(storage_dir_permutations, feature_label + ".%s.bed" % (i + 1))
            mapfiles[shifts[i]] = permfile
    elif permutationmode == 'backcirc':
        backgroundname = os.path.split(backgroundfile)[-1]
        if verbose: print ("BACKGROUND CIRCULAR randomisation before mapping")
        for i in range(numperms):
            shifts.append(int(random.random() * 10000000000))
            permfile = os.path.join(storage_dir_permutations, feature_label + ".%s.bed" % (i + 1))
            mapfiles[shifts[i]] = permfile
    shifts.append(0)  # the last one is the real data.

    # prepare for RANDOMISATION AGAINST BACKGROUND
    if permutationmode == 'background' or permutationmode == 'backcirc':
        f = open(backgroundfile)
        backlines = f.readlines()
        f.close()
    if permutationmode == 'backcirc':
        back_bed_data = {}  # create snp_bed_data equivalent object for background file
        backdict = {}
        print "creating background snp dict for BACKCIRC permutations..."
        for line in backlines:
            line = string.split(string.strip(line), '\t')
            chrom = 'chr' + line[0].replace('chr', '')
            address = int(line[1])
            snp = line[3]
            try:
                back_bed_data[chrom][address] = snp
            except:
                back_bed_data[chrom] = {}
                back_bed_data[chrom][address] = snp
            try:
                backdict["%s_%s" % (chrom, address)]
                print "duplicate address in background file", chrom, address, backdict["%s_%s" % (chrom, address)], snp
            except:
                try:
                    backdict["%s_%s" % (chrom, address)] = get_gwad(chrom, address, chromrange)
                except:
                    continue
                    # 
        # NOW ADD ALL THE REAL HITS TOO, JUST IN CASE THE BACKGROUND FILE ISN'T PERFECT
        # addition 26 Feb 2016 to allow immunochip background
        for chrom in snp_bed_data:
            for address in snp_bed_data[chrom]:
                snpaddress = "%s_%s" % (chrom, address)
                try:
                    backdict[snpaddress]
                except:
                    backdict["%s_%s" % (chrom, address)] = get_gwad(chrom, address, chromrange)
                    # 

    for shiftindex in range(len(shifts)):
        # MAKE PERMUTATIONS
        shift = shifts[shiftindex]
        # ===============
        if os.path.exists(mapfiles[shift]) and shift != 0:
            if doitagain:
                if verbose: print (
                "permutation file", mapfiles[shift], "exists but will not be used as doitagain is True")
            else:
                if verbose: print ("using existing permutation file", mapfiles[shift])
                continue  # don't replace existing permutation files. Saves time. Risk of error if the same common_label is used.
                # ===============
                # ------------------#  this is the randomisation code
        current_snp_bed_data = {}
        if permutationmode == 'circular' or shift == 0:
            # read through snp_bed_data, assigning new address to each snp
            if shift == 0:
                print 'writing bed file for real data'
            elif verbose:
                print ("CIRCULAR randomising shift = ", shift, permutationmode)
            for chrom in snp_bed_data:
                for address in snp_bed_data[chrom]:
                    snp = snp_bed_data[chrom][address]
                    if shift != 0:
                        snp = "rand_" + snp
                    this_gwad = get_gwad(chrom, address, chromrange)
                    this_gwad = this_gwad + shift
                    newchrom, newaddress = get_address(this_gwad, chromrange, chromlist)
                    if newchrom not in current_snp_bed_data:
                        current_snp_bed_data[newchrom] = {}
                    current_snp_bed_data[newchrom][newaddress] = snp
        elif permutationmode == 'background':
            if verbose: print ("BACKGROUND randomising shift = ", shift)
            current_snp_bed_data = {}
            for linenum in random.sample(range(len(backlines)), snp_count):
                line = string.split(string.strip(backlines[linenum]), '\t')
                chrom = line[0].replace('chr', '')
                address = line[2]
                try:
                    snp = line[4]
                except:
                    snp = "snpfromline_" + linenum
                if chrom not in current_snp_bed_data:
                    current_snp_bed_data[chrom] = {}
                current_snp_bed_data[chrom][address] = snp
        elif permutationmode == 'backcirc' and shift != '0':
            # read through snp_bed_data, assigning new address to each snp
            if verbose: print ("CIRCULAR randomising THE BACKGROUND FILE. shift = ", shift)
            sortedbackground = [key for key, value in sorted(backdict.iteritems(), key=lambda (k, v): (v, k),
                                                             reverse=False)]  # sorted by genome wide address
            for chrom in snp_bed_data:
                for address in snp_bed_data[chrom]:
                    snpaddress = "%s_%s" % (chrom, address)
                    newchoice = listitem(sortedbackground, sortedbackground.index(snpaddress) + shifts[shiftindex])
                    new = string.split(newchoice, '_')
                    newchrom = new[0]
                    newaddress = new[1]
                    if newchrom not in current_snp_bed_data:
                        current_snp_bed_data[newchrom] = {}
                    current_snp_bed_data[newchrom][newaddress] = snp
                    # thisperm = [listitem(sortedbackground, sortedbackground.index(snpaddress)+shifts[i]) for snpaddress in realdata]
                    # sortedrealdata = [key for key, value in sorted(realdata.iteritems(), key=lambda (k,v): (v,k), reverse=False)] #sorted by genome wide address

                    # ------------------#
                    # snp_bed_data to bed
        temp_snp_bed = working_files_dir + "/temporary.bed"
        o = open(temp_snp_bed, "w")
        for chrom in current_snp_bed_data:
            for address in current_snp_bed_data[chrom]:
                snpname = current_snp_bed_data[chrom][address]
                try:
                    o.write("%s\t%s\t%s\t%s\n" % (chrom.replace("chr", ""), address, int(address) + 1, snpname))
                except:
                    if verbose: print ("error writing bed file:")
        o.close()
        cmd = windowbed + " -w %s" % (window) + " -a " + temp_snp_bed + " -b " + feature_bed + " > " + mapfiles[shift]
        print cmd
        p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        print output  # so that errors are detected in progress file

# -#########################################################################################
# -#########################################################################################
# -#########################################################################################
# -########    ###########
# -######## START OF COEXPRESSION ANALYSIS        ###########
# -########    ###########
# -#########################################################################################
# -#########################################################################################
# -#########################################################################################

# ----------------------------------------------------#
# create dict of promoter addresses
if verbose: print ("reading feature coordinates...")
prom_addresses, feature_dict = readfeaturecoordinates(feature_coordinates)
# -------------------------------------#
# READ WHOLE EXPRESSION FILE INTO LARGE DICT - faster overall and ESSENTIAL FOR ESTIMATION OF P-VALUE.
if verbose: print ("reading expression file", expression_file)
if verbose: print("Now using %s bytes" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
exp_dict, header = read_expression_file(expression_file)
if verbose: print ("read")
# *************************************#
# READ METADATA IF AVAILABLE
try:
    f = open(promoter_details)
    lines = f.readlines()
    f.close()
except:
    lines = []
    pass
metadata = {}
knownnames = {}
for line in lines:
    line = string.split(string.strip(line), "\t")
    if len(line) < 2:
        continue
    metadata[line[0]] = line[1]
    try:
        name = (string.split(string.strip(line[1]), "@")[1])
    except:
        name = ""
    if name in knownnames:
        knownnames[name].append(line[0])
    elif name != "":
        knownnames[name] = [line[0]]
lines = []

# READ BONUS GENE NAMES/PROMOTERS IF AVAILABLE - these are additional gene or promoter names in the input set.
# ****DEPENDENT ON METADATA***
bonuspromoters = []
associatedgenes = {}
for name in badlines:  # lines that didn't have an integer in the right place in input file
    if name in knownnames.keys():  # gene name
        for prom in knownnames[name]:
            bonuspromoters.append(prom)
            associatedgenes[prom] = name
    elif name in exp_dict.keys():  # promoter name
        bonuspromoters.append(name)
        associatedgenes[name] = name
if len(bonuspromoters) > 0:
    if verbose: print("******* %s ADDITIONAL PROMOTERS ADDED*********" % len(bonuspromoters))
    if permutationmode != 'post':
        permutationmode = 'post'
        if verbose: print(
        "Permutations reset to post-mapping - additional promoters added manually so pre-mapping permutations are not applicable")

# *************************************#
# Make lists of correlations between random pairs
# This speeds calculation of empirical p-values considerably
if verbose: print ("reading correlation lists:", correlationmeasure)
if verbose: print (correlationfile)
if os.path.exists(correlationfile):
    f = open(correlationfile)
    lines = f.readlines()
    f.close()
    if len(lines) > gpl_len:
        lines = random.sample(lines, gpl_len)  # just use first million to save time...
    try:
        existingcorrelations = [float(x) for x in lines if string.strip(x) != ""]
    except:
        existingcorrelations = []
else:
    existingcorrelations = []
new_to_add = int(min(gpl_len, gpl_len - len(existingcorrelations)) * 1.2)  # make too many just to be sure.
if new_to_add <= 1:
    new_to_add = 1  # we need to initiate the AG system anyway here.
if verbose: print ("adding", new_to_add, "correlations to existing", len(existingcorrelations))

ptcf = exp_dict.keys()  # promoters to choose from

# NEWLY DEVELOPMENTS BY AG
dictlen = len(exp_dict)
entrylen = len(exp_dict[ptcf[0]])

# ctype structure for all promoter pair data
alldata_c = (c_double * (dictlen * entrylen))()

if correlationmeasure == "Spearman":  # creates an equivalent but ranked dictionary to 'exp_dict'
    exp_dict_rank = rankdict(exp_dict)

promindex_lookup = {}
ki = 0
for ikey in xrange(len(ptcf)):
    promindex_lookup[ptcf[ikey]] = ikey
    for entry in xrange(entrylen):
        if correlationmeasure == "Spearman":
            alldata_c[ki] = float(exp_dict_rank[ptcf[ikey]][entry])
        else:
            alldata_c[ki] = float(exp_dict[ptcf[ikey]][entry])
        ki = ki + 1

# ctype data structures to hold randomly generated indices
idxa_all_c = (c_int * new_to_add)()
idxb_all_c = (c_int * new_to_add)()

iPromPairs = 0
for x in range(new_to_add):
    idxa_candidate = random.randint(0, len(ptcf) - 1)
    idxb_candidate = random.randint(0, len(ptcf) - 1)

    # only consider correlations on different chromosomes in global distribution
    if ptcf[idxa_candidate][:5] != ptcf[idxb_candidate][:5]:
        idxa_all_c[iPromPairs] = idxa_candidate
        idxb_all_c[iPromPairs] = idxb_candidate
        iPromPairs = iPromPairs + 1

nPromPairs = iPromPairs
# ctype data structure to hold results
result_all_c = (c_double * nPromPairs)()

# get all the coefficients
coexpression.getPearson(byref(alldata_c), \
                        byref(result_all_c), \
                        byref(idxa_all_c), \
                        byref(idxb_all_c), \
                        c_int(nPromPairs), \
                        c_int(entrylen))

# set the NANs to -99, and store all in python list
for x in range(nPromPairs):
    if math.isnan(result_all_c[x]):
        correlation = -99
    else:
        correlation = result_all_c[x]
    existingcorrelations.append(correlation)
print len(existingcorrelations)

globalcorrelationlist = existingcorrelations
globalcorrelationlist.sort()  # essential if bisect is to be used
o = open(correlationfile, "w")
for r in globalcorrelationlist:
    o.write("%s\n" % r)
o.close()

for x in [0.05, 0.01, 0.001, 0.0005, 0.0001]:
    if verbose: print ("Correlation Centile at ", x, "=", centile(globalcorrelationlist, 1 - x))

# *************************************#
scatter_pairs = []
scatter_correlation = []
scatter_globalp = []
scatter_nsp = []
# *************************************#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RANDOMISATION POST MAPPING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if permutationmode == 'post':
    shifts = []
    mapfiles = {0: snps_mapped}
    # then randomise the existing real map file by rotation
    if verbose: print ("randomising by mapping")
    primary_dict = {}
    # READ SNPS_MAPPED
    f = open(snps_mapped)
    lines = f.readlines()
    f.close()
    for line in lines:
        line = string.split(string.strip(line), "\t")
        promoter = line[-1]
        primary_dict[promoter] = line[:-1]
        # -------------------------------------#
    for promoter in bonuspromoters:
        try:
            primary_dict[promoter] = [prom_addresses[promoter][0].replace('chr', ''), prom_addresses[promoter][1],
                                      prom_addresses[promoter][2], associatedgenes[promoter], "XXX", "XXX",
                                      "XXX"]  # put gene name in place of SNP in bed files
        except:
            primary_dict[promoter] = ["XXX", "XXX", "XXX", associatedgenes[promoter], "XXX", "XXX",
                                      "XXX"]  # put gene name in place of SNP in bed files
            # -------------------------------------#
            # rewrite snps_mapped with bonus promoters included
    o = open(snps_mapped, "w")
    for promoter in primary_dict:
        for entry in primary_dict[promoter]:  # same data, different label
            o.write("%s\t" % entry)
        o.write("%s\n" % promoter)  # a random promoter
    o.close()
    # -------------------------------------#

    # -------------------------------------#
    real_to_join = join_nearby(primary_dict.keys(), prom_addresses, globalcorrelationlist, exp_dict, pval_to_join,
                               correlationmeasure, soft_separation, anticorrelation)
    realdist = sorted([len(real_to_join[x]) for x in real_to_join])

    if verbose: print ("***", len(primary_dict.keys()), "promoters hit, in", len(real_to_join),
                       "distinct regions in the input dataset ***")
    if verbose: print ('realdist', realdist)

    # all_labels=prom_addresses.keys()
    all_labels = exp_dict.keys()
    all_labels.sort()
    shifts = [int(random.random() * 4 * len(all_labels)) for x in range(numperms)]

    permtimer = timeit.default_timer()

    for i in range(numperms):  # add some space for wastage from failed perms
        if verbose: print ("creating %s of %s post-mapping permutations" % (i + 1, numperms))
        # permfile=os.path.join(storage_dir_permutations,snps_mapped_filename+".%s"%i)
        permfile = os.path.join(storage_dir_permutations, feature_label + ".%s.bed" % (i + 1))
        shift = shifts[i]
        mapfiles[shift] = permfile
        # ===============
        if os.path.exists(mapfiles[shift]) and shift != 0:
            if doitagain:
                if verbose: print (
                "permutation file", mapfiles[shift], "exists but will not be used as doitagain is True")
            else:
                if verbose: print ("using existing permutation file", mapfiles[shift])
                continue  # don't replace existing permutation files. Saves time. Risk of error if the same common_label is used.
                # ===============

        if permutationmode == 'postcirc':
            # -----------circular post-mapping permutation
            thisperm = [listitem(all_labels, all_labels.index(promoter) + shift) for promoter in
                        primary_dict.keys()]  # A RANDOM SET OF PROMOTERS


            # -----------match real distribution
        else:
            real = sorted([len(real_to_join[a]) for a in real_to_join])

            bigcount = 0
            notready = True
            biglimit = 1
            nestlimit = 8

            while notready:
                bigcount += 1
                if bigcount > biglimit:
                    if verbose: print ("breaking big post permutations loop")
                    break  # just go ahead and use the last permutation anyway
                thisperm = []
                runcount = 0
                needed = copy.copy(realdist)  # need all
                while len(needed) > 0:
                    runcount += 1
                    if runcount > nestlimit:
                        if verbose: print ("breaking nested post permutations loop...")
                        break
                    if runcount == 0:  # first run
                        newlist = [listitem(all_labels, all_labels.index(promoter) + shift) for promoter in
                                   primary_dict.keys()]
                    elif runcount < 4:  # next few runs - try different shifts, select the proms that are needed.
                        newlist = [listitem(all_labels, all_labels.index(promoter) + int(
                            random.random() * len(primary_dict.keys() * 3))) for promoter in primary_dict.keys()]
                    else:  # shifts have failed, now create the right distribution by seeding from a random promoter
                        newlist = []
                        for n in needed:
                            r = int(n * 3)
                            randindex = random.choice(range(len(all_labels) - r))
                            t = [all_labels[randindex + x] for x in range(r)]  # transient
                            tg = join_nearby(t, prom_addresses, globalcorrelationlist, exp_dict, pval_to_join,
                                             correlationmeasure, soft_separation, anticorrelation)
                            b = [x for x in tg if len(tg[x]) >= n]
                            if len(b) > 0:
                                newlist += tg[random.choice(b)][:n]  # add in n adjacent promoters, chosen at random
                    newgroups = join_nearby(newlist, prom_addresses, globalcorrelationlist, exp_dict, pval_to_join,
                                            correlationmeasure, soft_separation, anticorrelation)
                    for new in newgroups:
                        if len(newgroups[new]) in needed:
                            thisperm += newgroups[new]
                            needed.remove(len(newgroups[new]))
                    perm_to_join = join_nearby(thisperm, prom_addresses, globalcorrelationlist, exp_dict, pval_to_join,
                                               correlationmeasure, soft_separation, anticorrelation)
                    # count how many are still outstanding
                    needed = []
                    perm = sorted([len(perm_to_join[b]) for b in perm_to_join])
                    for x in list(set(real)):
                        if real.count(x) > perm.count(x):
                            needed += [x for y in range(real.count(x) - perm.count(x))]
                if (len(perm_to_join) == len(real_to_join)) and (
                    len(thisperm) == len(primary_dict.keys())):  # we have the right distribution
                    notready = False
                else:
                    if verbose: print ("failed perm", perm)
            if notready:  # the perfect distribution was not found despite many searches. Now choose the right distribution from the permuted excess:
                nearestperm = []
                remaining = perm_to_join.keys()
                for desiredgrouplen in list(set(real)):
                    if verbose: print("desiredgrouplen", desiredgrouplen)
                    finished = False
                    acceptable_difference = 0
                    maxruns = desiredgrouplen * 2
                    runcount = 0
                    while not finished:
                        runcount += 1
                        if runcount > maxruns:
                            print "breaking post perm loop 2. "
                            break
                        searchpool = [x for x in remaining if
                                      abs(len(perm_to_join[x]) - desiredgrouplen) < acceptable_difference]
                        if verbose: print ("searchpool lengths = ", [len(perm_to_join[x]) for x in remaining if abs(
                            len(perm_to_join[x]) - desiredgrouplen) <= acceptable_difference])
                        if len(searchpool) < real.count(desiredgrouplen):
                            acceptable_difference += 1
                            if verbose: print ("broadening acceptable_difference:", acceptable_difference)
                        else:
                            for thisgrouplabel in random.sample(searchpool, real.count(desiredgrouplen)):
                                nearestperm += perm_to_join[thisgrouplabel]
                                remaining.remove(thisgrouplabel)
                            finished = True
                thisperm = nearestperm
                perm_to_join = join_nearby(thisperm, prom_addresses, globalcorrelationlist, exp_dict, pval_to_join,
                                           correlationmeasure, soft_separation, anticorrelation)

                # ---------------------------------------------
            if verbose: print ("perm_to_join:", sorted([len(perm_to_join[x]) for x in perm_to_join]))
            if verbose: print (
            len(primary_dict.keys()), "promoters hit, in", len(real_to_join), "distinct regions in input dataset")
            if verbose: print (
            len(thisperm), "promoters hit, in", len(perm_to_join), "distinct regions for this permutation")
            # ---------------------------------------------
        o = open(permfile, "w")
        for i in range(len(thisperm)):
            try:
                realpromoter = primary_dict.keys()[i]
            except:
                continue  # occasional error if more perm promoters get through than real promoters
            for entry in primary_dict[realpromoter]:  # same data, different label
                o.write("%s\t" % entry)
            o.write("%s\n" % thisperm[i])  # a random promoter
        o.close()
        # -----------
    shifts.append(0)  # the last one is the real data.

    if verbose: print (" perms created in: ", timeit.default_timer() - permtimer)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -------------------------------------------------------------------------------
shiftsdone = 0
if verbose: print ("searching for permstore:")
if os.path.exists(permstore) and os.path.exists(zerocountfile):
    if verbose: print ("permstore exists")
    if doitagain:
        if verbose: print ("but we're not using it because -z (doitagain) was set...")
    else:
        f = open(permstore)
        lines = f.readlines()
        f.close()
        shiftsdone = len(lines)
        if verbose: print ("Stored values available from %s previous permutations" % (shiftsdone))
# -------------------------------------------------------------------------------
# create directory for results to be stored in
storage_dir_perm_results = os.path.join(working_files_dir, "perm_results_store")
check_dir(storage_dir_perm_results)
# -------------------------------------------------------------------------------
jobs = []
for permnum, shift in enumerate(shifts[shiftsdone:]):
    # pcmd = "python 1-make-network.py -ef {} -fc {} -cf {} -mf {} -wd {} -sl {} -t {}".format(expression_file, feature_coordinates, correlationfile, mapfiles[shift], working_files_dir, supplementary_label, numthreads)
    pcmd = "python 1-make-network.py -mf {} -wd {}".format(mapfiles[shift], working_files_dir)
    jobs.append(pcmd)
o = open(jobfile, 'w')
o.write("{}".format('\n'.join(jobs[::-1])))  # reverse list so that biggest job (real data) is started first
o.close()

o = open(collationcommandfile, 'w')
# ccmd =  "python 2-collate-results.py -ef {} -fc {} -cf {} -mf {} -wd {} -sl {} -sc {} -gl {} -e {}".format(expression_file, feature_coordinates, correlationfile, mapfiles[shift], working_files_dir, supplementary_label, snp_count, genome_length, useremail)
ccmd = "python 2-collate-results.py -wd {}".format(working_files_dir)
o.write("{}".format(ccmd))
o.close()
# -------------------------------------------------------------------------------

if not prepareonly:
    print os.getcwd()
    cmd = "python runlocal.py -d {}".format(working_files_dir)
    print cmd
    subprocess.call(cmd, shell=True)
# -------------------------------------------------------------------------------
