#!/opt/local/bin/python

#-# get options 
import getopt, sys, os, time, threading, resource, gc, datetime
from subprocess import Popen, PIPE, STDOUT
from ctypes import *

import string, copy, stat, random, math, inspect, itertools, operator
from bisect import *
from scipy import stats
from scipy.stats import spearmanr
from scipy.stats import pearsonr

import pandas as pd
import networkx as nx
import timeit

import json
import numpy as np

###############################################################################
version = 'v1.0'
###############################################################################

try:
	opts, args = getopt.getopt(sys.argv[1:], "c: f: n: s: e: d: t: x: j: b: k: q: p v a i r z g h u ")
except:
	err = getopt.GetoptError
	print (str(err))
	sys.exit(2)
numperms = 10
useremail = "no_email"
seedfile = "none"
specialexpressionfile = "not_specified"
specialbedfile = "not_specified"
numthreads = 4
verbose = False
uberverbose = True
nsp = 1 #node-specific pvalue in use; specify precision (1=perfect, 0=globaldistribution is used.)
permutationmode = 'circular'
anticorrelation = False
noiter = False
recordedges=False
maxiperms=False
doitagain=False #create a whole new permutations set
common_label="null"
common_label2="null"
pval_to_join = 0.1 # *** 0.1 OPTIMAL *** >1 ==> join all
snp_details_file="supportingfiles/test.bed" #'unset'
backcirc = False
squarecoex = False
selectionthreshold_is_j=False
useaverage = False
for o, a in opts:
	if o == "-f":
		snp_details_file = a
	elif o == "-n":
		numperms = int(a)
	elif o == "-e":
		useremail = a
	elif o == "-t":
		numthreads = a
	elif o == "-x":
		specialexpressionfile = a
	elif o == "-j":
		pval_to_join = float(a)
	elif o == "-b":
		permutationmode = 'background'
		backgroundfile = a
	elif o == "-g":
		backcirc = True
	elif o == "-p":
		permutationmode = 'post'
	elif o == "-v":
		verbose = True
	elif o == "-a":
		anticorrelation = True
	elif o == "-i":
		noiter = True
	elif o == "-z":
		doitagain = True
if backcirc:
	permutationmode = 'backcirc' #reset this here so that the order of options doesn't matter.
if common_label=="null":
	if common_label2=="null":
		common_label2 = os.path.split(snp_details_file)[1].replace(".bed","")
	common_label = common_label2
	if common_label=="null": print "common_label not set"
elif common_label2 != "null": print "%s filename being ignored as common_label %s is set"%(common_label2, common_label)
print ("coexpression running", common_label, "permutationmode:", permutationmode, "numperms=", numperms, "email=", useremail, "nsp:", nsp, "seedfile", seedfile, "recordedges", recordedges, "doitagain", doitagain, "specialexpressionfile", specialexpressionfile)
nsp = float(nsp)
#_________________________________________________________________________________
tic = timeit.default_timer()
#_________________________________________________________________________________
#-###  SETTINGS
if selectionthreshold_is_j:
	selectionthreshold = -math.log(pval_to_join,10)
else:
	selectionthreshold=0
data_start_col=1 #nb zero  ex
data_start_row=1 #1084
chr_col=1  # -address of each locus is in expression file in order to apply min_separation
start_col=2
end_col=3
gpl_len = 1000000 # number of pearsons to include in the global pearson lists. #ONLY RELEVANT IF GENERATING NEW LIST!
write_layout_file = "no" #"yes" or anything
soft_separation=100000 #correlating promoters in this range will be joined. #IF pval_to_join >1 then this is a hard separation
progress_interval=50
ramlimit = 200000000
#_________________________________________________________________________________
#_________________________________________________________________________________
supplementary_label = common_label
if nsp==1:
	supplementary_label += "_complete"
else:
	supplementary_label += "_s%s"%nsp
if permutationmode == 'post': supplementary_label += "_POST"
elif permutationmode == 'backcirc': supplementary_label += "_BACKCIRC"
elif permutationmode == 'circular': supplementary_label += "_CIRCULAR"
elif permutationmode == 'background': supplementary_label += "_BACKGROUND"
if selectionthreshold_is_j: supplementary_label += "_SIJ"
if squarecoex: supplementary_label += "_SQUARECOEX"
if anticorrelation: supplementary_label += "_ANTIin"
if seedfile != "none": supplementary_label += "_seed"
if noiter: supplementary_label += "_NOiter"
if useaverage: supplementary_label += "_useav"
supplementary_label += "_pj%s"%pval_to_join
expstring = string.split(specialexpressionfile,"/")[-1].replace('.expression','')
if specialexpressionfile!="not_specified": supplementary_label += "_%s"%expstring
#_________________________________________________________________________________
#_________________________________________________________________________________
os.environ['OMP_NUM_THREADS'] = str(numthreads)
correlationmeasure = "Spearman"  #"Pearson" or "Spearman"
if correlationmeasure == "Pearson":
	doSpearman=0
elif correlationmeasure == "Spearman":
	doSpearman=1
#_________________________________________________________________________________
#snp_details file format - BED FORMAT
snp_col=3
chrom_col=0 #chrom is an integer only
address_col=1
snp_details_p_col=-1
#log_p_threshold=4 # only include snps with a negative_log_p ABOVE this threshold.
#_________________________________________________________________________________
############ settings to change when using a different machine ##################
#====
sourcefilesdir= "supportingfiles/"
chrom_lengths_file = sourcefilesdir+"chromlengths.txt" 
promoter_details   = sourcefilesdir+"METADATA_U22_ENHANCERS"
feature_coordinates= sourcefilesdir+"f5ep300_100.bed"
expression_file    = sourcefilesdir+"f5ep_ptt_condense.expression"
intersectBedlocation = sourcefilesdir+"intersectBed"
#====
resdir = "outputfiles/"
recordfile = 'already_run.sch'
alivefile = 'alive.sch'
restartfile = 'start.sch'
if snp_details_file == "unset":
	snp_details_file=os.path.join("/groups2/fantom5/coexpression_upload",common_label+".bed")
working_files_dir="%s%s"%(resdir, supplementary_label)
#====
storage_dir_premapping_randomisation = os.path.join(working_files_dir, "premapping_permutation_store")
storage_dir_post_mapping = os.path.join(working_files_dir, "postmapping_permutation_store")
thisqueue = os.path.join("%s/%s.queue"%(working_files_dir, supplementary_label))
ramfile = os.path.join(working_files_dir,"ramusage.txt")
settingsfile = os.path.join(working_files_dir,"settings.txt")
#====
sofile = sourcefilesdir+'coexpression_v2.so'
if not os.path.exists(sofile): 
	os.system("gcc -shared  -O3 -fPIC -fopenmp %scoexpression_v2.c -o %scoexpression_v2.so"%(sourcefilesdir, sourcefilesdir))
#====
permstore={}
zerocountfile={}
edgestore={}
#_*_*_*_*_*_*_*_*_
permstore = os.path.join(working_files_dir, "permstore_%s.txt"%(supplementary_label))
zerocountfile = os.path.join(working_files_dir, "zerocount_%s.txt"%(supplementary_label))
edgestore = os.path.join(working_files_dir, "edgestore_%s.txt"%(supplementary_label))
if seedfile != "none":
	seedfile=os.path.join(sourcefilesdir,seedfile)
print working_files_dir
print permstore
#_________________________________________________________________________________
runcommand = "python \
run_coexpression_\
%s.py -f %s -n %s -s %s -e %s -t %s -x %s -q %s"%(version, snp_details_file, numperms, nsp, useremail, numthreads, specialexpressionfile, specialbedfile)
runcommand += " -j %s"%(pval_to_join)
if permutationmode == 'post': runcommand += " -p"
elif permutationmode == 'background': runcommand += " -b %s"%(backgroundfile)
elif permutationmode == 'backcirc': runcommand += " -b %s -g"%(backgroundfile)
if verbose: runcommand += " -v"
if anticorrelation: runcommand += " -a"
if recordedges: runcommand += " -r"
if selectionthreshold_is_j: runcommand += " -h"
if useaverage: runcommand += " -u"
print runcommand
#_________________________________________________________________________________
if specialexpressionfile != "not_specified":
	expression_file = specialexpressionfile
	print "***"
	print "Using user-specified expressionfile:", expression_file
	print "***"
if specialbedfile!="not_specified":
	feature_coordinates = specialbedfile
	print "***"
	print "Using user-specified bedfile:", feature_coordinates
	print "***"
#_________________________________________________________________________________
gps = os.path.join(resdir,"globalperm.%s_%s"%(common_label, expstring))
if os.path.exists(gps):
	storage_dir_post_mapping = os.path.join(resdir,"globalperm.%s"%common_label,"postmapping_permutation_store")
	storage_dir_premapping_randomisation = os.path.join(resdir,"globalperm.%s"%common_label,"premapping_permutation_store")
	if not(doitagain):
		print "=======   using %s global postmapping permutation map files  ======="%(len(os.listdir(storage_dir_post_mapping)))
		print "=======   using %s global premapping permutation map files   ======="%(len(os.listdir(storage_dir_premapping_randomisation)))
	else:
		print "=======  Perms will be written to global perm store  ======="
#_________________________________________________________________________________
coexpression = cdll.LoadLibrary(sofile)
feature_label=string.split(string.split(feature_coordinates,"/")[-1],".")[0] #ie the "fantom5" in the coordinates file name
##################################################
snps_mapped_filename=feature_label+"_"+common_label+".snps_mapped.bed"
snps_mapped=os.path.join(working_files_dir,snps_mapped_filename)
statsfile = feature_label+"_"+supplementary_label+"_stats_snps.txt"
f5resource = "http://fantom.gsc.riken.jp/zenbu/gLyphs/#config=ne92nJ20PhPv5ziW90qnND;loc=hg19::"
#_________________________________________________________________________________
#-################################################################################
#- global objects used in functions.

#NB ALSO MUST PUT THESE AT THE START OF RUN ALL

globalpearsonlist=[]
nodespecificpearsonlists={}
shiftcount = 0
promindex_lookup={}
alldata_c = []
result_all_c = []
dictlen=0
entrylen=0
exp_dict={}

#-#################################################################################
#-######################		SCRIPT 2 FUNCTIONS	   ##########################
#-#################################################################################
def printv (text1, text2="", text3="", text4="", text5="", text6="", text7="", text8=""):
	if verbose:
		print text1, text2, text3, text4, text5, text6, text7, text8

def print_uberv (text1, text2="", text3="", text4="", text5="", text6="", text7="", text8=""):
	try:
		global uberverbose
	except:
		uberverbose = False
	if uberverbose:
		print text1, text2, text3, text4, text5, text6, text7, text8

def read_expression_file(filename, datastartcol='autodetect', datastartrow='autodetect', labelcol=0, labelrow=0):
	'''returns dict of {name:[list_of_floats], ...}, header row, startcol, startrow'''
	f=open(filename)
	lines = [string.split(string.strip(x),'\t') for x in f.readlines()]
	f.close()
	print_uberv("exp file opened %s: %s bytes"%('', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))	
	checklines = lines[:10] # these will be used to determine dsc and dsr.
	firstfloatfailcols=[]
	for r, row in enumerate(checklines):
		for c in range(len(row)-1, -1, -1):
			try:
				float(row[c])
			except:
				if row[c] != "NA":
					#print "unfloatable", "row",r,"col", c, checklines[r][c]
					firstfloatfailcols.append(c+1)
					break
	print_uberv("rows checked %s: %s bytes"%('', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))	
	if len(firstfloatfailcols)>0:
		dsc = firstfloatfailcols[-1]
		for r in range(len(firstfloatfailcols)-1,-1,-1):
			if firstfloatfailcols[r] != dsc:
				#print "first floatfailcol inconsistency", r, firstfloatfailcols[r]
				dsr = r + 1
				break
	else:
		dsr=0 # necessary in case firstfloatfailcols is empty
		dsc=0
	print_uberv("columns checked %s: %s bytes"%('', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))	
	if datastartcol=='autodetect':
		datastartcol = dsc
	else:
		datastartcol = int(datastartcol)
		if datastartcol != dsc:
			print "datastartcol (%s) discrepancy with autodetect (%s)"%(datastartcol, dsc)
	if datastartrow=='autodetect':
		datastartrow = dsr
	else:
		datastartrow = int(datastartrow)
		if datastartrow != dsr:
			print "datastartrow (%s) discrepancy with autodetect (%s)"%(datastartrow, dsr)
	print "reading expression file from column:%s and row:%s"%(dsc, dsr)
	exdic = {}
	progdict = {x:1 for x in range(0,len(lines),10000)}
	for i, line in enumerate(lines[datastartrow:]):
		try:
			progdict[i]
			print_uberv("%s lines read: %s bytes"%(i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
		except:
			pass
		try:
			exdic[line[labelcol]]
		except:
			exdic[line[labelcol]]=[]
		for x in line[datastartcol:]:
			try:
				exdic[line[labelcol]].append(float(x))
			except:
				exdic[line[labelcol]].append('NA')
	print "expression file read"
	return exdic, lines[labelrow][datastartcol:]

#_________________________________________________________________________________
#--- function to ensure that if script fails, this can be easily detected by another script.
end = False
def confirm_alive():
	global end, alivefile, shiftcount
	if not threading.enumerate()[0].is_alive():
		print "ending because MainThread has died"
		end = True
	wp = "w" #write mode for alivefile
	if os.path.exists(alivefile):
		wp = "a"
	a=open(alivefile,wp)
	a.write("%s\t%s\n"%(runcommand, shiftcount))
	a.close()
	if end:
		print "stop sending alive signal now", end
	else:
		# call f() again in 30 seconds
		threading.Timer(30, confirm_alive).start()
# start calling f now and every 30 sec thereafter
confirm_alive()
#_________________________________________________________________________________

def dir_list(pathname):
	target=os.path.abspath(pathname)
	file_list=[]
	initial_list=[]
	for root,dirs, files in os.walk(target):
		if root==target:
			initial_list=files
			for filename in initial_list:
				if string.split(filename,".")[-1]=="fa":
					file_list.append(str(filename))
	return file_list

def fix_permissions(this_path):
	os.system("/bin/chmod 755 %s"%(this_path))

def check_dir(this_dir):
	if not os.path.isdir(this_dir):
		os.mkdir(this_dir)
	fix_permissions(this_dir)

def get_gwad(chrom, address, masterchromrange):
	'''get a single number to identify this position on a concatentated genome'''
	try:
		masterchromrange[chrom]
	except:
		print "chrom not in masterchromrange", chrom
		return "chrom not in masterchromrange"
	if (masterchromrange[chrom][0] + address) < masterchromrange[chrom][1] or address < 0:
		return masterchromrange[chrom][0]+address
	else:
		printv ("error finding address - address is not on chromosome", chrom, address)

def get_address(gwad, masterchromrange, masterchromlist):
	'''go from gwad to chrom and position'''
	genome_max=masterchromrange[masterchromlist[-1]][1] # the highest value of all
	# first bring the value into range
	while gwad > genome_max:
		gwad = gwad - genome_max
	for chrom in masterchromlist:
		start=masterchromrange[chrom][0]
		end=masterchromrange[chrom][1]
		if gwad > start and gwad <= end:
			address = gwad - start
			return chrom, address
	printv ("***** ERROR NEEDS ATTENTION: address not found for gwad:", gwad)
	return "chr1", 0

def remove_overlap(thisdict):
	'''thisdict has format thisdict[chrom][start]=end'''
	newdict={}
	for chrom in thisdict:
		newdict[chrom]={}
		current_chrom_starts = sorted(thisdict[chrom].keys())
		last_start = current_chrom_starts[0]
		last_end = thisdict[chrom][last_start]
		breakfound = "no"
		for this_start in current_chrom_starts[1:]:
			this_end=thisdict[chrom][this_start]
			if this_start <= last_end: #overlap exists with last window
				last_start = min(this_start,last_start)
				last_end = max(this_end,last_end)
				breakfound = "no"
			else:  # there is no overlap. add concatentated range to dict
				newdict[chrom][last_start]=last_end
				last_start = this_start
				last_end = this_end
				breakfound = "yes"
		#add in the last one
		#if breakfound == "no":
		# peroxide j birthday
		newdict[chrom][last_start]=last_end
	return newdict

#-#################################################################################
#-######################		SCRIPT 5 FUNCTIONS	     ##########################
#-#################################################################################


def listitem(thislist, index): #FOR RANDOMISATION AT MAPPING
	'''correct list index to be within list range'''
	a=float(index)/len(thislist)
	b=a-int(round(a,0))
	c=b*len(thislist)
	d=int(round(c,0)) #to correct python int output
	return thislist[d]

def _ecdf(x):
	nobs = len(x)
	return np.arange(1, nobs + 1) / float(nobs)

def fdr_correction(pvals, alpha=0.05, method='indep'):
	pvals = np.asarray(pvals)
	shape_init = pvals.shape
	pvals = pvals.ravel()
	pvals_sortind = np.argsort(pvals)
	pvals_sorted = pvals[pvals_sortind]
	sortrevind = pvals_sortind.argsort()
	if method in ['i', 'indep', 'p', 'poscorr']:
		ecdffactor = _ecdf(pvals_sorted)
	elif method in ['n', 'negcorr']:
		cm = sum(1. / np.arange(1, len(pvals_sorted) + 1))
		ecdffactor = _ecdf(pvals_sorted) / cm
	else:
		raise ValueError("Method should be 'indep' and 'negcorr'")
	reject = pvals_sorted < (ecdffactor * alpha)
	if reject.any():
		rejectmax = max(np.nonzero(reject)[0])
	else:
		rejectmax = 0
	reject[:rejectmax] = True
	pvals_corrected_raw = pvals_sorted / ecdffactor
	pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
	pvals_corrected[pvals_corrected > 1.0] = 1.0
	pvals_corrected = pvals_corrected[sortrevind].reshape(shape_init)
	reject = reject[sortrevind].reshape(shape_init)
	return reject, pvals_corrected

def centile(list, percent, tot="len"):
	if len(list)<1:
		return "empty_list"
	if tot=="len":
		tot=len(list)
	if tot < len(list): #impossible if list is subset of tot!
		tot=len(list)
	index=int(percent * float(tot)) # the index in an imaginary complete list
	index=index-(tot-len(list)) #subtract the non-existent values
	if index<0:
		return "<%s"%list[0]
	try:
		list.sort()
		return list[index]
	except:
		return "centile_error"

def empiricalp(thisvalue, empiricalcollection):
	if len(empiricalcollection) == 0:
		return 1
	empiricalcollection.sort()
	return 1-float(bisect_left(empiricalcollection,thisvalue))/len(empiricalcollection)

def get_pearson(pa, pb, ed):
	if len(ed[pa])==0 or len(ed[pb])==0:
		print "empty list in get_pearson", pa, pb
	if sum(ed[pa])>0 and sum(ed[pb])>0:
		result = spearmanr(ed[pa],ed[pb])[0]
		if correlationmeasure == "Pearson":
			try:
				result = pearsonr(ed[pa],ed[pb])[0]
			except:
				print_uberv("pearsonr fail for ", pa, pb)
				result = -99
		elif correlationmeasure == "Spearman":
			try:
				result = spearmanr(ed[pa],ed[pb])[0]
			except:
				print_uberv("spearmanr fail for ", pa, pb)
				result = -99
	else:
		print_uberv("sum zero for one of these:", pa, pb)
		result = -99
	if math.isnan(result):
		print_uberv("isnan for ",correlationmeasure, pa, pb)
		return -99
	return result

def merge_ranges(ranges):
	''' Merges overlapping / adjacent ranges. 'ranges' must be sorted BEFORE
	running this function. Call with 'np.fromiter(merge_ranges(ranges),
	dtype='i8,i8')'. '''
	ranges = ranges + [0,+soft_separation] # add soft_separation
	ranges = iter(ranges) # generate iterator
	current_start, current_stop = next(ranges)
	for start, stop in ranges:
		if start > current_stop: # Gap between ranges
			yield current_start, current_stop # Output current range
			current_start, current_stop = start, stop # Start a new one.
		else: # Range overlapping (or adjacent) therefore merge.
			current_stop = max(current_stop, stop) # 'stop' not necessarily larger than 'current_stop'
	yield current_start, current_stop

def join_nearby(thesepromoters, theseaddresses):
	''' Creates a dictionary of the form (sorted lists):
	{prom10: [prom10, prom11 with link to prom10 etc... ],
	prom20: [prom20, prom21 with link to prom20 etc... ]} '''
	ptj={}
	promGrpGraph = nx.Graph() # creates graph promGrpGraph
	# creates DataFrame: index = promoter name, columns = chromosome, start, end
	a = pd.DataFrame({k: theseaddresses[k] for k in thesepromoters},index=['chromosome','start','end']).T
	a.sort(['chromosome','start','end'],inplace=True) # sorts DataFrame inplace - VERY IMPORTANT
	for chrom in pd.unique(a.chromosome):
		grps = np.fromiter(merge_ranges(a[a.chromosome==chrom][['start','end']].values), dtype='i8,i8')
		for grp_start, grp_end in grps:
			promGrpGraph.clear() # clears all the nodes and edges from promGrpGraph
			promGrp = list(a[(a.chromosome==chrom) & (a.start >= grp_start) & (a.start < grp_end)].index) # list of promoters in the group
			promGrpGraph.add_nodes_from(promGrp) # add the nodes to the graph promGrpGraph
			# only need to calculate triangular array. (Pearson / Spearman-r are symmetric corr(X,Y)==corr(Y,X))
			for i in range(len(promGrp)):
				for j in range(i+1,len(promGrp)): # +1 as do not need to calculate the diagaonal itself
					# ----
					r = get_pearson(promGrp[i],promGrp[j], exp_dict)
					# WHOLE-DISTRIBUTION P VALUE
					p1 = 1-float(bisect_left(globalpearsonlist,r))/len(globalpearsonlist)
					if anticorrelation:
						p2 = float(bisect_left(globalpearsonlist,r))/len(globalpearsonlist)
						p = min(p1,p2)
					else: p = p1
					# ----
					if p < pval_to_join:
						promGrpGraph.add_edge(promGrp[i], promGrp[j]) # add edges to the graph
			for i in list(nx.connected_components(promGrpGraph)):
				new = sorted(list(i)) # sort list. Is this necessary?
				ptj[new[0]]=new
	return ptj

#-#################################################################################
import matplotlib.pyplot as plt

def qq_plot(real, rand, filename):
	if len(rand) == 0:
		printv ("empty rand list for qq")
		return "empty"
	o = open(filename,"w")
	o.write("rand\t%s\trand_again\n"%filename)
	limit = min(len(rand),len(real))
	real.sort()
	rand.sort()
	expected=[]
	observed=[]
	for i in range(limit):
		indexa = int(len(rand)*float(i)/limit)
		indexb = int(len(real)*float(i)/limit)
		a = rand[indexa]
		b = real[indexb]
		o.write("%s\t%s\t%s\n"%(a,b,a))
		expected.append(a)
		observed.append(b)
	o.close()

	#----------------------------
	fig=plt.figure(figsize=(5, 5), dpi=300)
	ax=plt.axes() #[0,0,1,1], frameon=False) #fig.add_subplot(1,1,1)
	plt.box("off")
	#ax.set_title("Q:Q plot %s"%(supplementary_label))
	#ax.set_yticklabels(ticklabels, size=12, color="red") #use the list defined above
	ax.xaxis.set_ticks_position('bottom') #remove ticks from top and right
	ax.yaxis.set_ticks_position('left') #remove ticks from top and right
	ax.set_xlabel("Expected coexpression scores", size=12)
	ax.set_ylabel("Observed coexpression scores", size=12)
	#ax.set_xscale("log") #use logarithmic x axis
	#ax.set_yscale("log") #use logarithmic x axis
	#steps=2
	#ax.set_xticks(range(0, max(expected), steps))
	#ax.set_yticks(range(0,max(observed), steps))
	#ax.set_xticklabels([0,50,100], size=12)
	for i in range(len(expected)):
		ax.scatter([expected[i]],[expected[i]],color=["red"],alpha=0.7,marker="o", s=12)
		ax.scatter([expected[i]],[observed[i]],color=["blue"],alpha=0.7,marker="o", s=12)
	#axis limits
	if len(expected)>0 and len(observed) >0:
		plt.xlim(0,max(expected)+1)
		plt.ylim(0,max(observed)+1)
		fig.savefig("%s.svg"%(filename))
	#----------------------------

#-#################################################################################
#-#################################################################################
def give_up(message="giving up..."):
	global end
	end = True
	print message
	try:
		q = open(thisqueue,"w")
		q.write("Not enough hits to be worth pursuing")
		q.close()
		fix_permissions(thisqueue)
		f=open(recordfile,"a")
		f.write("%s\n"%common_label)
		f.close()
	except:
		pass
	sys.exit()

def fail(reason = "not sure"):
	global end
	end = True #instruct the alive function to cease.
	r = open(recordfile,"a")
	r.write("%s\n"%(common_label))
	r.close()
	print "script 6 failed:", common_label, reason, sys.exc_info()
	q = open(thisqueue,"w")
	q.write("We are sorry to say that your job id %s has failed. Please <a href='contact.php'>contact the administrator by email</a> to let us know, and we'll do our best to work out why and get it working.<br><br><br>\n\n\nError info: <br>\n %s<br> \n %s"%(common_label, sys.exc_info(), reason))
	q.close()
	fix_permissions(thisqueue)
	sys.exit()

def runall():
	#seed random number generator to allow for testing
	#random.seed(176202)
	global shiftcount
	global permutationmode
	global promindex_lookup
	global alldata_c
	global result_all_c
	global entrylen
	global globalpearsonlist
	global exp_dict
	starttime = time.time()
	#-######################
	print "resdir", resdir
	print "working_files_dir", working_files_dir

	#-######################


	for thisdir in [resdir, working_files_dir, storage_dir_premapping_randomisation, storage_dir_post_mapping]:
		check_dir(thisdir)

	#_________________________________________________________________________________
	settingsdict={"version":version, "numperms":numperms, "useremail":useremail, "seedfile":seedfile, "label":supplementary_label, "permutationmode":permutationmode, "verbose":verbose, "anticorrelation":anticorrelation, "nsp":nsp, "gpl_len":gpl_len, "write_layout_file":write_layout_file, "noiter":noiter, "selectionthreshold":selectionthreshold, "pval_to_join":pval_to_join, "soft_separation":soft_separation, "progress_interval":progress_interval, "correlationmeasure":correlationmeasure, "data_start_col":data_start_col, "data_start_row":data_start_row, "chr_col":chr_col, "start_col":start_col, "end_col":end_col, "runcommand":runcommand}
	o=open(settingsfile,"w")
	for key, value in settingsdict.iteritems():
		o.write("%s\t%s\n"%(key,value))
	o.close()
	#-#########################################################################################
	#-#########################################################################################
	#-#########################################################################################
	#-#########																		 ##########
	#-#########				START OF MAPPING										 ##########
	#-#########																	 	 ##########
	#-#########################################################################################
	#-#########################################################################################
	#-#########################################################################################

	#read chromosome lengths
	chrom_lengths={}
	genome_length=0
	f=open(chrom_lengths_file)
	lines=[string.split(string.strip(x),": ") for x in f.readlines()]
	f.close()
	for line in lines:
		chrom_lengths[line[0]]=int(line[1])
		genome_length+=int(line[1])
	print "====="
	print "genome length=", genome_length
	print "====="
	#read snps
	snp_bed_data={}
	snp_details={}
	snp_p_values={}
	all_snps=[]
	badlines = []
	printv ("reading snp details (bed) file", snp_details_file)
	f=open(snp_details_file)
	lines=f.readlines()
	f.close()
	snp_count=0
	for line in lines:
		if "\t" in line:
			line=string.split(string.strip(line),"\t")
		elif " " in line:
			line=string.split(string.strip(line)," ")
		else:
			line=string.split(string.strip(line))
		try:
			address=int(line[address_col])
		except:
			printv ("skipping line of snp_details file", line)
			try:
				badlines.append(line[0])
			except:
				pass
			continue
		if len(line)>1:
			snp=line[snp_col]
			chrom="chr%s"%line[chrom_col].replace("chr","")
			if chrom not in chrom_lengths:
				printv ("skipping", chrom, "as length not known")
				continue # skip chromosomes for which we don't know the length
			if snp not in all_snps:
				all_snps.append(snp) # one big list of all snps included in this run
			try:
				snp_bed_data[chrom][address]=snp
			except:
				snp_bed_data[chrom]={}
				snp_bed_data[chrom][address]=snp
			try:
				thisp=float(line[snp_details_p_col]) #4th col of snp_details file
			except:
				thisp=1
			try:
				snp_p_values[line[snp_col]]
				snp_p_values[line[snp_col]] = min(snp_p_values[line[snp_col]], thisp)
			except:
				snp_p_values[line[snp_col]] = thisp
			snp_details[snp_col]=[line[chr_col], line[address_col], line[snp_details_p_col]]
		else:
			printv ("line too short(",len(line),"):", line)
	snp_count=len(all_snps)
	if len(snp_bed_data)==0:
		printv ("No snps in .snp_details file")
	else:
		printv ("snps read successfully")

	#--------------------------------------------------------
	#DEFINE GENOME-WIDE ADDRESSES FOR TRANSLATION FUNCTIONS
	chromrange={}
	chromnums = range(1,23)+['X',"Y","M"]
	chromlist=[]
	start=0
	for chromnum in chromnums:
		chrom="chr%s"%chromnum
		chromlist.append(chrom)
		thislen=chrom_lengths[chrom]
		end=start+thislen
		chromrange[chrom]=[start,end]
		start+=thislen
	#--------------------------------------------------------
	#feature_coordinates is already bed file from v28 onwards...
	feature_bed = feature_coordinates
	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	#CHECK WHETHER WORTH CONTINUING...
	temp_snp_bed = working_files_dir+"/temporary.bed"
	o=open(temp_snp_bed,"w")
	for chrom in snp_bed_data:
		for address in snp_bed_data[chrom]:
			snpname = snp_bed_data[chrom][address]
			try:
				o.write("%s\t%s\t%s\t%s\n"%(chrom.replace("chr",""),int(address)-1,address,snpname))
			except:
				printv ("error writing bed file:")
	o.close()
	command = intersectBedlocation+" -a "+temp_snp_bed+" -b "+feature_bed+" -wa -wb > "+snps_mapped
	printv (command)
	p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
	output = p.stdout.read()
	printv (output) # so that errors are detected in progress file
	f = open(snps_mapped)
	lines=f.readlines()
	f.close()
	numhits = len(lines)
	if numhits < 2 and len(badlines)<2:
		give_up("not enough snps mapped - could there be a chr in the bedfile?")
	#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	#RANDOMISATION PRE-MAPPING
	shifts=[]
	mapfiles={0:snps_mapped}
	if numhits >1: #CAN ONLY PROCEED IF THERE ARE ENOUGH SNPS TO BE WORTH MAPPING
		if permutationmode == 'circular':
			printv ("CIRCULAR randomisation before mapping")
			for i in range(numperms):
				shifts.append(int(random.random()*10000000000))
				permfile=os.path.join(storage_dir_premapping_randomisation,feature_label+"_"+common_label+".snps_mapped.%s.bed"%(i+1))
				mapfiles[shifts[i]]=permfile
		elif permutationmode == 'background':
			backgroundname = os.path.split(backgroundfile)[-1]
			printv ("BACKGROUND randomisation before mapping")
			for i in range(numperms):
				shifts.append(i+1) # so that 0 is reserved for real data
				permfile=os.path.join(storage_dir_premapping_randomisation,feature_label+"_"+backgroundname+"_"+common_label+".snps_mapped.%s.bed"%(i+1))
				mapfiles[shifts[i]]=permfile
		elif permutationmode == 'backcirc':
			backgroundname = os.path.split(backgroundfile)[-1]
			printv ("BACKGROUND CIRCULAR randomisation before mapping")
			for i in range(numperms):
				shifts.append(int(random.random()*10000000000))
				permfile=os.path.join(storage_dir_premapping_randomisation,feature_label+"_"+backgroundname+"_backcirc2_"+common_label+".snps_mapped.%s.bed"%(i+1))
				mapfiles[shifts[i]]=permfile
		shifts.append(0) # the last one is the real data.

		#prepare for RANDOMISATION AGAINST BACKGROUND
		if permutationmode == 'background' or permutationmode == 'backcirc':
			f=open(backgroundfile)
			backlines=f.readlines()
			f.close()
		if permutationmode == 'backcirc':
			back_bed_data={} # create snp_bed_data equivalent object for background file
			backdict = {}
			print "creating background snp dict for BACKCIRC permutations..."
			for line in backlines:
				line=string.split(string.strip(line), '\t')
				chrom = 'chr'+line[0].replace('chr','')
				address = int(line[1])
				snp = line[3]
				try:
					back_bed_data[chrom][address]=snp
				except:
					back_bed_data[chrom]={}
					back_bed_data[chrom][address]=snp
				try:
					backdict["%s_%s"%(chrom,address)]
					print "duplicate address in background file", chrom, address, backdict["%s_%s"%(chrom,address)], snp
				except:
					try:
						backdict["%s_%s"%(chrom,address)]=get_gwad(chrom, address, chromrange)
					except:
						continue
			#•••••••••••
			# NOW ADD ALL THE REAL HITS TOO, JUST IN CASE THE BACKGROUND FILE ISN'T PERFECT
			# addition 26 Feb 2016 to allow immunochip background
			for chrom in snp_bed_data:
				for address in snp_bed_data[chrom]:
					snpaddress = "%s_%s"%(chrom,address)
					try:
						backdict[snpaddress]
					except:
						backdict["%s_%s"%(chrom,address)]=get_gwad(chrom, address, chromrange)
			#•••••••••••



		for i in range(len(shifts)):
			shift = shifts[i]
			#===============
			if os.path.exists(mapfiles[shift]) and shift != 0:
				if doitagain:
					printv ("permutation file", mapfiles[shift], "exists but will not be used as doitagain is True")
				else:
					printv ("using existing permutation file", mapfiles[shift])
					continue #don't replace existing permutation files. Saves time. Risk of error if the same common_label is used.
			#===============
			#------------------#  this is the randomisation code
			current_snp_bed_data={}
			if permutationmode == 'circular' or shift == 0: # this is where the real data gets written to bed file
				#read through snp_bed_data, assigning new address to each snp
				if shift == 0: print 'writing bed file for real data'
				else: printv ("CIRCULAR randomising shift = ", shift, permutationmode)
				for chrom in snp_bed_data:
					for address in snp_bed_data[chrom]:
						snp=snp_bed_data[chrom][address]
						if shift != 0:
							snp = "rand_" + snp
						try:
							this_gwad = get_gwad(chrom, address, chromrange)
						except:
							continue
						this_gwad = this_gwad + shift
						newchrom, newaddress = get_address(this_gwad, chromrange, chromlist)
						if newchrom not in current_snp_bed_data:
							current_snp_bed_data[newchrom]={}
						current_snp_bed_data[newchrom][newaddress]=snp
			elif permutationmode=='background':
				printv ("BACKGROUND randomising shift = ", shift)
				current_snp_bed_data={}
				for linenum in random.sample(range(len(backlines)),snp_count):
					line=string.split(string.strip(backlines[linenum]),'\t')
					chrom = line[0].replace('chr','')
					address = line[2]
					try:
						snp = line[4]
					except:
						snp = "snpfromline_"+linenum
					if chrom not in current_snp_bed_data:
						current_snp_bed_data[chrom]={}
					current_snp_bed_data[chrom][address]=snp
			elif permutationmode=='backcirc' and shift != '0':
				#read through snp_bed_data, assigning new address to each snp
				printv ("CIRCULAR randomising THE BACKGROUND FILE. shift = ", shift)
				sortedbackground = [key for key, value in sorted(backdict.iteritems(), key=lambda (k,v): (v,k), reverse=False)] #sorted by genome wide address
				for chrom in snp_bed_data:
					for address in snp_bed_data[chrom]:
						snpaddress = "%s_%s"%(chrom,address)
						newchoice = listitem(sortedbackground, sortedbackground.index(snpaddress)+shifts[i])
						new = string.split(newchoice,'_')
						newchrom = new[0]
						newaddress = new[1]
						if newchrom not in current_snp_bed_data:
							current_snp_bed_data[newchrom]={}
						current_snp_bed_data[newchrom][newaddress]=snp

			#------------------#
			#snp_bed_data to bed
			temp_snp_bed = working_files_dir+"/temporary.bed"
			o=open(temp_snp_bed,"w")
			for chrom in current_snp_bed_data:
				for address in current_snp_bed_data[chrom]:
					snpname = current_snp_bed_data[chrom][address]
					try:
						o.write("%s\t%s\t%s\t%s\n"%(chrom.replace("chr",""),address,int(address)+1,snpname))
					except:
						printv ("error writing bed file:")
			o.close()

			command = intersectBedlocation+" -a "+temp_snp_bed+" -b "+feature_bed+" -wa -wb > "+mapfiles[shift]
			p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
			output = p.stdout.read()
			printv(output) # so that errors are detected in progress file
			#os.system(intersectBedlocation+" -a "+temp_snp_bed+" -b "+feature_bed+" -wa -wb > "+mapfiles[shift] )

	#-#########################################################################################
	#-#########################################################################################
	#-#########################################################################################
	#-########																	    ###########
	#-########				 START OF COEXPRESSION ANALYSIS					        ###########
	#-########																	    ###########
	#-#########################################################################################
	#-#########################################################################################
	#-#########################################################################################
	#-#########################################################################################

	#----------------------------------------------------#
	#create dict of promoter addresses
	#read feature - same procedure as 2_map_snps_to_feature
	printv ("reading feature coordinates...")
	f=open(feature_coordinates)
	lines=f.readlines() #range_features file
	f.close()
	go="no"
	prom_addresses={}
	feature_dict={}
	#prog=range(0,len(lines),max(int(float(len(lines))/10),1))
	for i in range(len(lines)): #each line is a feature (ie a TSS region)
		#if i in prog: print i, "of", len(lines)
		line=lines[i]
		if go=="no":
			#skip past header information
			try:
				test=string.split(string.strip(line),"\t")
				int(test[2])
				go="yes"
			except:
				continue
		line=string.split(string.strip(line),"\t")
		label=line[-1]
		chrom=line[0]
		start=int(line[1])
		end=int(line[2])
		prom_addresses[label]=[chrom,start,end]
		st = min(start, end)
		en = max(start, end)
		try:
			feature_dict[chrom]
		except:
			feature_dict[chrom]={}
		try:
			feature_dict[chrom][st]
			if feature_dict[chrom][st] < en : #choose the greater end
				feature_dict[chrom][st] = en
		except:
			feature_dict[chrom][st]=en

	#-------------------------------------#
	#READ WHOLE EXPRESSION FILE INTO LARGE DICT - faster overall and ESSENTIAL FOR ESTIMATION OF P-VALUE.
	printv ("reading expression file", expression_file)
	print_uberv("%s: %s bytes"%(i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
	exp_dict, header = read_expression_file(expression_file)
	printv ("read")
	#*************************************#
	#READ METADATA IF AVAILABLE
	try:
		f=open(promoter_details)
		lines=f.readlines()
		f.close()
	except:
		lines=[]
		pass
	metadata={}
	knownnames={}
	for line in lines:
		line=string.split(string.strip(line),"\t")
		if len(line)<2:
			continue
		metadata[line[0]]=line[1]
		try:
			name = (string.split(string.strip(line[1]),"@")[1])
		except:
			name = ""
		if name in knownnames:
			knownnames[name].append(line[0])
		elif name != "":
			knownnames[name] = [line[0]]
	lines=[]

	#READ BONUS GENE NAMES/PROMOTERS IF AVAILABLE - these are additional gene or promoter names in the input set.
	#****DEPENDENT ON METADATA***
	bonuspromoters=[]
	associatedgenes={}
	for name in badlines: #lines that didn't have an integer in the right place in input file
		if name in knownnames.keys(): #gene name
			for prom in knownnames[name]:
				bonuspromoters.append(prom)
				associatedgenes[prom]=name
		elif name in exp_dict.keys(): #promoter name
			bonuspromoters.append(name)
			associatedgenes[name]=name
	if len(bonuspromoters) > 0:
		printv("******* %s ADDITIONAL PROMOTERS ADDED*********"%len(bonuspromoters))
		if permutationmode != 'post':
			permutationmode = 'post'
			printv("Permutaions reset to post-mapping - additional promoters added manually so pre-mapping permutations are not applicable")

	####################  READ SEEDFILE   #########################
	# this will seed every network with the contents of seedfile, then remove them all before calculating
	lines=[]
	if seedfile != "none":
		f = open(seedfile)
		lines=f.readlines()
		f.close()
	seedpromoters = []
	for line in lines:
		seedpromoters.append(string.split(string.strip(line),"\t")[-1]) # last entry in line is promoter (eg bedfile)

	#*************************************#
	#Make lists of correlations between random pairs
	#This speeds calculation of empirical p-values considerably
	printv ("reading correlation lists:", correlationmeasure)
	expfilename = string.split(expression_file,"/")[-1]
	if correlationmeasure == "Pearson":
		correlationfile = "%s%s.pearsonlist"%(sourcefilesdir, expfilename)
	elif correlationmeasure == "Spearman":
		correlationfile = "%s%s.spearmanlist"%(sourcefilesdir, expfilename)
	printv (correlationfile)
	if os.path.exists(correlationfile):
		f=open(correlationfile)
		lines=f.readlines()
		f.close()
		if len(lines)>gpl_len:
			lines = random.sample(lines,gpl_len) #just use first million to save time...
		try:
			existingcorrelations = [float(x) for x in lines if string.strip(x) != ""]
		except:
			existingcorrelations = []
	else:
		existingcorrelations = []
	new_to_add = int(min(gpl_len,gpl_len-len(existingcorrelations))*1.2) #make too many just to be sure.
	if new_to_add <=1:
		new_to_add = 1 #we need to initiate the AG system anyway here.
	printv ("adding" , new_to_add, "correlations to existing", len(existingcorrelations))

	ptcf = exp_dict.keys()  #promoters to choose from
	# NEWLY DEVELOPMENTS BY AG
	dictlen=len(exp_dict)
	entrylen=len(exp_dict[ptcf[0]])

	# ctype structure for all promoter pair data
	alldata_c= (c_double * (dictlen*entrylen))()

	if doSpearman==1: # creates an equivalent but ranked dictionary to 'exp_dict'
		exp_dict_rank={}
		for i in exp_dict.keys():
			exp_dict_rank[i]=stats.rankdata(exp_dict[i])

	promindex_lookup={}
	ki=0
	for ikey in xrange(len(ptcf)):
		promindex_lookup[ptcf[ikey]]=ikey
		for entry in xrange(entrylen):
			if doSpearman == 1:
				alldata_c[ki]=float(exp_dict_rank[ptcf[ikey]][entry])
			else:
				alldata_c[ki]=float(exp_dict[ptcf[ikey]][entry])
			ki=ki+1

	#ctype data structures to hold randomly generated indices
	idxa_all_c= (c_int * new_to_add)()
	idxb_all_c= (c_int * new_to_add)()

	iPromPairs=0
	for x in range(new_to_add):
		idxa_candidate=random.randint(0,len(ptcf)-1)
		idxb_candidate=random.randint(0,len(ptcf)-1)

		#only consider correlations on different chromosomes in global distribution
		if ptcf[idxa_candidate][:5] != ptcf[idxb_candidate][:5]:
			idxa_all_c[iPromPairs]=idxa_candidate
			idxb_all_c[iPromPairs]=idxb_candidate
			iPromPairs=iPromPairs+1

	nPromPairs=iPromPairs
	#ctype data structure to hold results
	result_all_c= (c_double * nPromPairs)()

	# get all the coefficients
	coexpression.getPearson(byref(alldata_c), \
							byref(result_all_c), \
							byref(idxa_all_c), \
							byref(idxb_all_c), \
							c_int(nPromPairs),  \
							c_int(entrylen))


	# set the NANs to -99, and store all in python list
	for x in range(nPromPairs):
		if math.isnan(result_all_c[x]):
			correlation=-99
		else:
			correlation=result_all_c[x]
		existingcorrelations.append(correlation)
	print len(existingcorrelations)

	globalpearsonlist = existingcorrelations
	globalpearsonlist.sort() # essential if bisect is to be used
	o=open(correlationfile,"w")
	for r in globalpearsonlist:
		o.write("%s\n"%r)
	o.close()

	for x in [0.05,0.01,0.001,0.0005,0.0001]:
		printv ("Correlation Centile at ",x,"=",centile(globalpearsonlist, 1-x))

	#*************************************#
	scatter_pairs=[]
	scatter_pearson=[]
	scatter_globalp=[]
	scatter_nsp=[]
	#*************************************#
	#RANDOMISATION POST MAPPING
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if permutationmode == 'post':
		shifts = []
		mapfiles={0:snps_mapped}
		#then randomise the existing real map file by rotation
		printv ("randomising by mapping")
		primary_dict={}
		#READ SNPS_MAPPED
		f=open(snps_mapped)
		lines=f.readlines()
		f.close()
		for line in lines:
			line=string.split(string.strip(line),"\t")
			promoter = line[-1]
			primary_dict[promoter]=line[:-1]
		#-------------------------------------#
		for promoter in bonuspromoters:
			try:
				primary_dict[promoter]=[prom_addresses[promoter][0].replace('chr',''),prom_addresses[promoter][1],prom_addresses[promoter][2],associatedgenes[promoter],"XXX","XXX","XXX"] #put gene name in place of SNP in bed files
			except:
				primary_dict[promoter]=["XXX","XXX","XXX",associatedgenes[promoter],"XXX","XXX","XXX"] #put gene name in place of SNP in bed files
		#-------------------------------------#
		#rewrite snps_mapped with bonus promoters included
		o=open(snps_mapped,"w")
		for promoter in primary_dict:
			for entry in primary_dict[promoter]: #same data, different label
				o.write("%s\t"%entry)
			o.write("%s\n"%promoter) # a random promoter
		o.close()
		#-------------------------------------#

		#-------------------------------------#
		real_to_join = join_nearby(primary_dict.keys(), prom_addresses)
		realdist = sorted([len(real_to_join[x]) for x in real_to_join])

		printv ("***", len(primary_dict.keys()), "promoters hit, in",len(real_to_join),"distinct regions in the input dataset ***")
		printv ('realdist', realdist)

		#all_labels=prom_addresses.keys()
		all_labels=exp_dict.keys()
		all_labels.sort()
		shifts=[int(random.random()*4*len(all_labels)) for x in range(numperms)]

		permtimer = timeit.default_timer()

		for i in range(numperms): #add some space for wastage from failed perms
			printv ("creating %s of %s post-mapping permutations"%(i+1,numperms))
			permfile=os.path.join(storage_dir_post_mapping,snps_mapped_filename+".%s"%i)
			shift = shifts[i]
			mapfiles[shift]=permfile
			#===============
			if os.path.exists(mapfiles[shift]) and shift != 0:
				if doitagain:
					printv ("permutation file", mapfiles[shift], "exists but will not be used as doitagain is True")
				else:
					printv ("using existing permutation file", mapfiles[shift])
					continue #don't replace existing permutation files. Saves time. Risk of error if the same common_label is used.
			#===============

			if permutationmode == 'postcirc':
				#-----------circular post-mapping permutation
				thisperm = [listitem(all_labels,all_labels.index(promoter)+shift) for promoter in primary_dict.keys()] # A RANDOM SET OF PROMOTERS


				#-----------match real distribution
			else:
				real = sorted([len(real_to_join[a]) for a in real_to_join])

				bigcount = 0
				notready = True
				biglimit = 1
				nestlimit = 8

				while notready:
					bigcount+=1
					if bigcount > biglimit:
						printv ("breaking big post permutations loop")
						break #just go ahead and use the last permutation anyway
					thisperm = []
					runcount=0
					needed = copy.copy(realdist) # need all
					while len(needed) > 0:
						runcount+=1
						if runcount > nestlimit:
							printv ("breaking nested post permutations loop...")
							break
						if runcount == 0: # first run
							newlist = [listitem(all_labels,all_labels.index(promoter)+shift) for promoter in primary_dict.keys()]
						elif runcount < 4: # next few runs - try different shifts, select the proms that are needed.
							newlist = [listitem(all_labels,all_labels.index(promoter)+int(random.random()*len(primary_dict.keys()*3))) for promoter in primary_dict.keys()]
						else: # shifts have failed, now create the right distribution by seeding from a random promoter
							newlist=[]
							for n in needed:
								r = int(n*3)
								randindex = random.choice(range(len(all_labels)-r))
								t=[all_labels[randindex + x] for x in range(r)] #transient
	 							tg = join_nearby(t, prom_addresses) #transient groups
	 							b=[x for x in tg if len(tg[x]) >= n]
	 							if len(b)>0:
									newlist+=tg[random.choice(b)][:n] #add in n adjacent promoters, chosen at random
						newgroups = join_nearby(newlist, prom_addresses)
						for new in newgroups:
							if len(newgroups[new]) in needed:
								thisperm += newgroups[new]
								needed.remove(len(newgroups[new]))
						perm_to_join = join_nearby(thisperm, prom_addresses)
						# count how many are still outstanding
						needed = []
						perm = sorted([len(perm_to_join[b]) for b in perm_to_join])
						for x in list(set(real)):
							if real.count(x) > perm.count(x):
								needed += [x for y in range(real.count(x) - perm.count(x))]
					if (len(perm_to_join) == len(real_to_join)) and (len(thisperm) == len(primary_dict.keys())): # we have the right distribution
						notready = False
					else:
						printv ("failed perm", perm)
				if notready: # the perfect distribution was not found despite many searches. Now choose the right distribution from the permuted excess:
					nearestperm = []
					remaining = perm_to_join.keys()
					for desiredgrouplen in list(set(real)):
						printv("desiredgrouplen", desiredgrouplen)
						finished = False
						acceptable_difference = 0
						maxruns = desiredgrouplen*2
						runcount = 0
						while not finished:
							runcount += 1
							if runcount > maxruns:
								print "breaking post perm loop 2. "
								break
							searchpool = [x for x in remaining if abs(len(perm_to_join[x])-desiredgrouplen) < acceptable_difference]
							printv ("searchpool lengths = ", [len(perm_to_join[x]) for x in remaining if abs(len(perm_to_join[x])-desiredgrouplen) <= acceptable_difference])
							if len(searchpool) < real.count(desiredgrouplen):
								acceptable_difference += 1
								printv ("broadening acceptable_difference:", acceptable_difference)
							else:
								for thisgrouplabel in random.sample(searchpool,real.count(desiredgrouplen)):
									nearestperm += perm_to_join[thisgrouplabel]
									remaining.remove(thisgrouplabel)
								finished = True
					thisperm = nearestperm
					perm_to_join = join_nearby(thisperm, prom_addresses)

				#---------------------------------------------
				printv ("perm_to_join:", sorted([len(perm_to_join[x]) for x in perm_to_join]))
				printv (len(primary_dict.keys()), "promoters hit, in",len(real_to_join),"distinct regions in input dataset")
				printv (len(thisperm), "promoters hit, in",len(perm_to_join),"distinct regions for this permutation")
				#---------------------------------------------
			o=open(permfile,"w")
			for i in range(len(thisperm)):
				try:
					realpromoter = primary_dict.keys()[i]
				except:
					continue #occasional error if more perm promoters get through than real promoters
				for entry in primary_dict[realpromoter]: #same data, different label
					o.write("%s\t"%entry)
				o.write("%s\n"%thisperm[i]) # a random promoter
			o.close()
			#-----------
		shifts.append(0) # the last one is the real data.

		printv (" perms created in: ", timeit.default_timer() - permtimer)

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#-------------------------------------------------------------------------------
	#PREPARE
	randompearson=[]
	already_calculated={}
	snp_map={}
	network_total={}
	num_nodes=0
	num_perm_nodes=0
	perm_edges=[]
	perm_indm=[]
	zerocount=0
	perm_pear=[]
	randompvals=[]
	randompearsons=[]
	#-------------
	shiftsdone = 0
	printv ("searching for permstore:")
	if os.path.exists(permstore) and os.path.exists(zerocountfile):
		printv ("permstore exists")
		if doitagain:
			printv ("but we're not using it because -z (doitagain) was set...")
		else:
			f=open(permstore)
			lines=f.readlines()
			f.close()
			shiftsdone = len(lines)
			printv ("Stored values available from %s previous permutations"%(shiftsdone))
	#-------------
	printv ("shifts:",shifts)
	if shiftsdone >= len(shifts):
		shiftsdone = len(shifts)-1
	shiftcount = shiftsdone
	permtimer = timeit.default_timer()

	for shift in shifts[shiftsdone:]:
		#-----
		print "last perm time = ", timeit.default_timer() - permtimer
		permtimer=timeit.default_timer()
		#-----

		shiftcount += 1
		if shift==0: #the last one
			printv ("real data this time...")
		#READ SNPS_MAPPED
		printv (mapfiles[shift])
		f=open(mapfiles[shift])
		lines=f.readlines()
		f.close()
		printv ("reading snp map(BED TOOLS): %s lines"%len(lines))

		snp_map={}
		for i in range(len(lines)):
			line=string.split(string.strip(lines[i]),"\t")
			if len(line)<4:
				continue
			this_promoter=line[-1]
			try:
				exp_dict[this_promoter]
			except:
				print "skipping", this_promoter, "as not in expression file"
				continue
			snp=line[3]
			#store p-value for each snp mapped
			try:
				snp_map[this_promoter]
				snp_map[this_promoter]["snps"] += "|" + snp
			except:
				snp_map[this_promoter]={"p":-9999, "snps":snp}
			try:
				logp=snp_p_values[snp] #always log p in snp_details file
					#if this is a script2 randomisation, "rand_" is appended to the snp name
			except:
				try:
					logp=snp_p_values[snp.replace("rand_","")]
				except:
					#ie if there is no p-value for this "snp" - might be a new bonus promoter
					logp=1000
			if logp > snp_map[this_promoter]["p"]:
				snp_map[this_promoter]["p"] = logp
		inputpromoterlist = snp_map.keys()
		for new in seedpromoters:
			if new not in snp_map.keys():
				snp_map[new]={"p":1000.125} #specific p-value to label these additions
		printv ('snp map read.')
		#-------------------------------------#
		edges_sought=0

		all_proms=[]
		individualscores={}
		promoters=snp_map.keys()
		num_nodes=len(promoters)
		num_perm_nodes+=len(promoters)
		if len(promoters)==0:
			continue
		prog=range(0,len(promoters),int(len(promoters)/(100/progress_interval))+1)

		#AG developments
		promaidx_all=[]
		prombidx_all=[]

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		use_specific_p = True
		if nsp <= 1.0/len(ptcf):
			use_specific_p = False
			print "nsp", nsp, "is too small because there are only",len(ptcf),"promoters. It must be greater than", 1.0/len(ptcf)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if use_specific_p:
			if nsp < 1:
				ptcf_sample = random.sample(ptcf, int(float(len(ptcf))*nsp))
			else:
				ptcf_sample = ptcf
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#CALCULATE ALL CORRELATIONS
			for i in range(len(promoters)):
				jrange = range(len(ptcf_sample)) # asymmetrical matrix for node specific pval
				for j in jrange:
					if i!=j:
						promaidx_all.append(promindex_lookup[promoters[i]])
						prombidx_all.append(promindex_lookup[ptcf_sample[j]])
			printv ("added", len(prombidx_all), "correlations to correlation index")

			nPromPairs=len(promaidx_all)
			idxa_all_c= (c_int *  nPromPairs)(*promaidx_all)
			idxb_all_c= (c_int *  nPromPairs)(*prombidx_all)

			#ctype data structure to hold results
			result_all_c= (c_double * nPromPairs)()

			tta=timeit.default_timer()

			# get all the coefficients
			coexpression.getPearson(byref(alldata_c), \
							byref(result_all_c), \
							byref(idxa_all_c), \
							byref(idxb_all_c), \
							c_int(nPromPairs),  \
							c_int(entrylen))

			

			ttb=timeit.default_timer()

			printv ("correlations complete")


			# set the NANs to -99, and store all in python list
			AGcorrelations=[]
			for x in range(nPromPairs):
				if math.isnan(result_all_c[x]):
					correlation=-99
				else:
					correlation=result_all_c[x]
				AGcorrelations.append(correlation)

			#++++++++++++++++++++++++++++++++++++
			ip=0
			for i in range(len(promoters)):
				nodespecificpearsonlists[promoters[i]] = []
				for j in jrange:
					if i!=j:
						nodespecificpearsonlists[promoters[i]].append(AGcorrelations[ip])
						ip=ip+1
				nodespecificpearsonlists[promoters[i]].sort() #this is essential so that p-vals can be calculated.
			#++++++++++++++++++++++++++++++++++++
			if recordedges and shift==0:
				scatter_pairs.append("realdatanow")
				scatter_pearson.append("realdatanow")
				scatter_globalp.append("realdatanow")
				scatter_nsp.append("realdatanow")
			#++++++++++++++++++++++++++++++++++++
		ip=0
		for i in range(len(promoters)):
			if i in prog: printv ("%s of %s promoters"%(i, len(promoters)))
			for j in range(len(promoters)): # NB must do this from both sides for nsp (ie not for j in range(i,len...))
				if i!=j:
					prom1 = promoters[i]
					prom2 = promoters[j]

					#pearson=AGcorrelations[ip]
					#pearson = AGcordict[prom1][prom2]
					pearson = get_pearson(prom1,prom2,exp_dict)

					#======================================
					if use_specific_p:
						if pearson == -99:
							p1=1
						else:
							p1 = 1-float(bisect_left(nodespecificpearsonlists[prom1],pearson))/len(nodespecificpearsonlists[prom1])
							if anticorrelation:
								p2 = float(bisect_left(nodespecificpearsonlists[prom1],pearson))/len(nodespecificpearsonlists[prom1])
								p1 = min(p1,p2)
					else: #legacy code. Enables calculation of whole-dist p
						# WHOLE-DISTRIBUTION P VALUE
						if pearson == -99:
							p1=1
						else:
							p1 = 1-float(bisect_left(globalpearsonlist,pearson))/len(globalpearsonlist)
							if anticorrelation:
								p2 = float(bisect_left(globalpearsonlist,pearson))/len(globalpearsonlist)
								p1 = min(p1,p2)
					#======================================
					if recordedges:
						scatter_pairs.append("%s\t%s"%(prom1,prom2))
						scatter_pearson.append(pearson)
						scatter_globalp.append(1-float(bisect_left(globalpearsonlist,pearson))/len(globalpearsonlist))
						scatter_nsp.append(1-float(bisect_left(nodespecificpearsonlists[prom1],pearson))/len(nodespecificpearsonlists[prom1]))
					#======================================

					try:
						log_edge_p = -math.log(p1,10)
					except:
						log_edge_p = 0

					#STORE P VALUES
					try:
						individualscores[prom1][prom2]=log_edge_p
					except:
						individualscores[prom1]={}
						individualscores[prom1][prom2]=log_edge_p
					if not(nsp):
						try:
							individualscores[prom2][prom1]=log_edge_p
						except:
							individualscores[prom2]={}
							individualscores[prom2][prom1]=log_edge_p
					ip=ip+1
		#-------------------------------------#
		for prom in individualscores:
			#just to make sure no promoters are left out
			all_proms.append(prom)
			for other_prom in individualscores[prom]:
				all_proms.append(other_prom)
		all_proms=list(set(all_proms))
		#•••••••••••••••••••••••••••••••••••••
		# check if this made a difference
		commontoboth = set(all_proms) & set(promoters)
		print "promoters:%s common:%s all_proms:%s "%(len(promoters), len(commontoboth), len(all_proms))
		#•••••••••••••••••••••••••••••••••••••

		#-------------------------------------#
		#work out which promoters should be joined as one.
		prom_to_join = join_nearby(all_proms, prom_addresses)
		joining_instructions={}
		for group_label in prom_to_join:
			for prom in prom_to_join[group_label]:
				joining_instructions[prom]=group_label
		printv (len(all_proms), "promoters hit, in",len(prom_to_join),"distinct regions")
		#-------------------------------------#
		printv ("merging ind scores")
		# NOW MERGE INDIVIDUAL SCORES.
		indjoined={}
		indmeasure={}
		whittledindmeasure={}
		weights={}
		winners={}

		allpromoterscores = {}

		for group_label in prom_to_join:
			scores={}
			#choose the contender from each group with the most unlikely coexpression anywhere in the network of __all__ promoters 
			for prom in prom_to_join[group_label]:
				if joining_instructions[prom] == group_label:
					scores[prom]=sum([individualscores[prom][otherprom] for otherprom in\
					individualscores[prom] if individualscores[prom][otherprom] >= selectionthreshold\
					and joining_instructions[otherprom] != group_label])
			if len(scores)==0:
				continue
			winners[group_label]=max(scores.iteritems(), key=operator.itemgetter(1))[0]
		#____________________________________________________________________________________
		# new in v0.63: repeat choice of winners using the new consensus to reduce the influence of large groups of promoters/enhancers
		for group_label in prom_to_join:
			scores={}
			checkscorelens=[]
			#choose the contender from each group with the most unlikely coexpression with another group
			for prom in prom_to_join[group_label]:
				if joining_instructions[prom] == group_label:
					thisscore = []
					for otherprom in winners.values():
						if prom != otherprom and joining_instructions[otherprom] != group_label:
							try:
								thisscore.append(individualscores[prom][otherprom])
							except:
								#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
								#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
								#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
								print "••••••••• prom lookup failure. output dump follows •••••••••"
								try:
									individualscores[prom]
								except:
									print "**** %s **** not in individual_scores"%prom
								print prom, otherprom
								print joining_instructions[prom]
								print joining_instructions[otherprom]
								print "====individual scores==="
								for x in individualscores:
									for y in individualscores[x]:
										print x,y,individualscores[x][y] 
								print "==== winners ===="
								for gl in winners:
									print gl, winners[gl]
								print "==== all_proms ===="
								print sorted(all_proms)
								print "==== promoters ===="
								print sorted(promoters)
								print "==== all_proms not promoters ===="
								print set(all_proms) - set(promoters)
								print "==== promoters not all_proms ===="
								print set(promoters) - set(all_proms)
								print "===== groups ======"
								for group_label2 in prom_to_join:
									print group_label2, prom_to_join[group_label2]
								print "================"
								#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
								#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
								#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
					scores[prom] = sum(thisscore)
					allpromoterscores[prom] = scores[prom] # store for output later
					checkscorelens.append(len(thisscore))
			#---------
			if checkscorelens.count(checkscorelens[0]) != len(checkscorelens): 
				printv("checkscorelens not equal!", set(checkscorelens) ) # each group should be scored against the same number of other groups. 
			if len(scores)==0:
				continue
			new_winner = max(scores.iteritems(), key=operator.itemgetter(1))[0]
			if new_winner != winners[group_label]:
				printv("==> new winner:%s (overtook %s)"%(new_winner, winners[group_label]))
			winners[group_label]=new_winner
			#---------

		#____________________________________________________________________________________
		for prom in individualscores:
			group_label = joining_instructions[prom]
			chosen = winners[group_label] #the chosen promoter for this group
			if group_label not in indjoined:
				indjoined[group_label]={}
			for other_prom in individualscores[prom]:
				other_group_label = joining_instructions[other_prom]
				other_chosen = winners[other_group_label]
				if chosen != other_chosen:
					try:
						indjoined[group_label][other_group_label] = individualscores[chosen][other_chosen]
					except:
						if nsp:
							print "error - no pval found", chosen, other_chosen
						else:
							indjoined[group_label][other_group_label] = individualscores[other_chosen][chosen]
		#_______________________________________________________________________________________

		# now calculate the coexpression scores
		for prom in indjoined:
			if useaverage:
				indmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in\
					indjoined[prom] if indjoined[prom][otherprom] >= selectionthreshold]) / len(indjoined)
			else:
				indmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in\
						indjoined[prom] if indjoined[prom][otherprom] >= selectionthreshold])

		cutlist = [] # list of promoters to remove from analysis
		if noiter == False: #cycle down the list removing top hit and recalculating indmeasure
			lastscore = -1
			for prom, score in sorted(indmeasure.iteritems(), key=lambda (k,v): (v,k), reverse=True):
				if score == lastscore:
					whittledindmeasure[prom] = lastresult
				else:
					if useaverage:
						whittledindmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in\
							indjoined[prom] if (indjoined[prom][otherprom] >= selectionthreshold and\
							otherprom not in cutlist)]) / len(indjoined)
					else:
						whittledindmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in\
							indjoined[prom] if (indjoined[prom][otherprom] >= selectionthreshold and\
							otherprom not in cutlist)])


				cutlist.append(prom)
				lastscore = score
				lastresult = whittledindmeasure[prom]
			indmeasure = whittledindmeasure

		#-------------------------------------#
		seedstoprune = [x for x in seedpromoters if x not in inputpromoterlist]
		for prom in seedstoprune:
			printv ('PRUNING SEED:', prom)
			puregroup = [x for x in prom_to_join[joining_instructions[prom]] if x not in seedstoprune]
			if len(puregroup) > 0:
				printv ("puregroup", puregroup)
			if len(puregroup) == 0:
				#then there are no input promoters in this group
				[joining_instructions[prom]]
				try:
					del indmeasure[joining_instructions[prom]]
				except:
					pass
				try:
					del indjoined[joining_instructions[prom]]
				except:
					pass
				for innocent in indjoined:
					try:
						del indjoined[innocent][joining_instructions[prom]]
					except:
						pass
		#-------------------------------------#

		#_______________________________________________________________________________________
		# correct coexpression scores for number of distinct observations made
		# • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • 
		# • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • 
		#indmeasure = {x:indmeasure[x]/len(indmeasure) for x in indmeasure} # eliminates the spatial enrichment signal
		#• • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • 
		#• • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • • 
		#_______________________________________________________________________________________

		#count up network totals
		network_total=0
		for group in indmeasure:
			network_total += indmeasure[group] #sum of all ind scores for this network.
		#-------------------------------------#
		printv ("**********s**********")
		printv ("SHIFT:",shift,[round(x,1) for x in sorted(indmeasure.values())])
		printv ("**********s**********")
		#-------------------------------------#
		#now collate the permuted results [all other code for permutations has been identical to the real code until here]
		if shift != 0: #then this is a permutation. Store permutations in case this script runs out of memory and is killed
			if os.path.exists(permstore):o=open(permstore,"a")
			else:o=open(permstore,"w")
			for prom in indmeasure:
				if indmeasure[prom] == 0:
					zerocount += 1
				else:
					perm_indm.append(indmeasure[prom])
					o.write("%s\t"%indmeasure[prom])
			o.write("\n")
			o.close()
			z=open(zerocountfile,"w")
			z.write("%s"%zerocount)
			z.close()
			if recordedges:
				if os.path.exists(edgestore):
					o=open(edgestore,"a")
				else:
					o=open(edgestore,"w")
				for prom1 in indjoined:
					for prom2 in indjoined[prom1]:
						edge = indjoined[prom1][prom2]
						o.write("%s\t"%edge)
				o.write("\n")
				o.close()

		#----------------
		#output ram usage
		ram = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		duration = time.time() - starttime
		try: o=open(ramfile,"a")
		except:
			o.open(ramfile,"w")
			o.write("PID: %s\n"%os.getpid())
		o.write("%s bytes, %s seconds\n"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, duration))
		printv ("%s bytes, %s seconds\n"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, duration))
		if ram > ramlimit:
			o.write("Process will die now as over ram limit.\n")
		o.close()
		if ram > ramlimit:
			try:
				o=open(restartfile,"a")
			except:
				o=open(restartfile,"w")
			o.write("\n%s\n"%runcommand)
			o.close()
			fail("Too much memory used: restart imminent")
		#----------------
		printv ("shift", shifts.index(shift), "of", numperms, "done.")
		printv ('--------------------------------')

	exp_dict={} #clear RAM - doesn't appear to work...
	gc.collect() # collect garbage
	printv ("exp_dicts cleared")
	if shift != 0:
		for i in range(5):
			printv ("************REAL DATA NOT INCLUDED - RANDOM DATA ONLY*************")
	printv (edges_sought, "total edges searched for shift=", shift)


	#-####################################################################################
	#-####################################################################################

	# write output files
	prefix=""
	#---- to make sure permutation source is clear in output files
	suffix = "%s_%s"%(correlationmeasure, selectionthreshold)
	if permutationmode == 'circular':
		suffix += "_%sCIRCULARperms"%(numperms)
	elif permutationmode == 'post':
		suffix += "_%sPOSTperms"%(numperms)
	elif permutationmode == 'background':
		suffix += "_%sBACKGROUNDperms"%(numperms)

	#--------
	# cover the contingency in which there are no promoters
	try:
		indjoined
	except:
		indjoined={}
		indmeasure={}
		individualscores={}
		prom_to_join={}
	#--------
	#lr = open(os.path.join(working_files_dir,prefix+common_label+"_%s_"%s+".qq_plot_specifics"),"w")
	realpvals=[]
	realmeasures=[]
	already=[]
	for prom1 in indjoined:
		realmeasures.append(indmeasure[prom1])
		for other_prom in indjoined[prom1]:
			combo = sorted([prom1,joining_instructions[other_prom]])
			combo = "%s_%s"%(combo[0],combo[1])
			#lr.write("%s\t%s\t%s\n"%(prom1,combo,indjoined[prom1][other_prom]))
			if combo not in already:
				already.append(combo)
				realpvals.append(indjoined[prom1][other_prom])
	#------------put zeros back into perm_indm-----------------
	if os.path.exists(permstore) and os.path.exists(zerocountfile) and not doitagain:
		f=open(permstore)
		lines=f.readlines()
		f.close()
		perm_indm = [] ###### delete existing contents of perm_indm
		for line in lines:
			perm_indm += [float(x) for x in string.split(string.strip(line),"\t") if x != ""]
		printv ("%s nonzero values in total from %s previous permutations"%(len(perm_indm), len(lines)))
		f=open(zerocountfile)
		line = f.readline() # single line
		f.close()
		zerocount = int(string.strip(line))
		fullrange_perm_indm = perm_indm + [0 for x in range(zerocount)]
	else: # no previous permutations have been run...
		fullrange_perm_indm = perm_indm
	#------------save qq_plot files-----------------
	#qq_plot(realpvals, randompvals, os.path.join(working_files_dir,prefix+common_label+"_%s_.%s.qq_plot_pvals"%(s,suffix)))
	qqfile = prefix+supplementary_label+"_%s.qq_plot"%(suffix)
	try:
		qq_plot(realmeasures,fullrange_perm_indm, os.path.join(working_files_dir,qqfile))
		qq_plot(realmeasures,fullrange_perm_indm, os.path.join(working_files_dir,qqfile)) #duplicate in offline directory
	except:
		fail("qq plot not working")
	#------------calculate FDR----------------------
	permuted_indm_p = sorted([empiricalp(x,fullrange_perm_indm) for x in realmeasures], reverse=False)
	fdr_permuted_indm = list(fdr_correction(permuted_indm_p,0.05,"negcorr")[1])
	fdr_BH_permuted_indm = list(fdr_correction(permuted_indm_p)[1])
	#lr.close()
	fdrstore={}

	printv (network_total)
	#output individual scores
	pseudonyms={}
	htmlresults = prefix+supplementary_label+"_%s.html"%(suffix)
	h = open(os.path.join(working_files_dir,htmlresults),"w")
	h.write("<table id=\"results\">\n<tr>\n<th>Top promoter (click for details in new tab)</th>\n<th>SNPs in top_promoter</th>\n<th>Coexpression score</th>\n<th>FDR</th>\n</tr>")
	fullresults = prefix+supplementary_label+"_%s.individual_scores"%(suffix)
	outputfilecontents=''
	outputfilecontents+="Top_promoter[SNPs_in_top_promoter]\tpromoter_snps_in_region\tchosen_promoter_id\tpromoters\toriginal_p\tcoexpression_score\traw_p(perm)\tBenjamini-Hochberg\n"
	#for key in indjoined:
	i=0
	highlight=[]
	for key, value in sorted(indmeasure.iteritems(), key=lambda (k,v): (v,k), reverse=True):
		snps=[]
		group_label = joining_instructions[key]
		for prom in prom_to_join[group_label]:
			try:
				snp_list=string.split(snp_map[prom]["snps"],"|")
			except:
				snp_list=["no_snps.seed_promoter"]
			for snp in snp_list:
				if snp not in snps:
					snps.append(snp)
		snp_string=""
		for snp in snps:
			snp_string+="|%s"%snp
		if snp_string[0]=="|":
			snp_string = snp_string[1:]
		prom_string=""
		m = []
		for promoter in prom_to_join[joining_instructions[key]]:
			prom_string+="|%s"%promoter
			if promoter in metadata:
				m.append("%s"%metadata[promoter])
			else:
				m.append(promoter)
		m = "|".join(m)
		if prom_string[0]=="|":
			prom_string = prom_string[1:]
		wname=winners[group_label]
		if winners[group_label] in metadata:
			wname = metadata[winners[group_label]]
			wname.replace(",+","[+]")
			wname.replace(",-","[-]")
			wname.replace(",",", ")
			wname.replace("[+]",",+ ")
			wname.replace("[-]",",- ")
		pseudonyms[winners[group_label]]=wname
		theseps=[0]
		try:
			winning_snps = snp_map[winners[group_label]]["snps"]
		except:
			winning_snps = ""
		wslist = string.split(winning_snps,"|")
		wslist = [w for w in wslist if w!=""]
		for snp in wslist:
			try:
				theseps.append(snp_p_values[snp])
			except:
				printv ("failed to find original p value for ", snp)
		winning_snps = winning_snps.replace("||","|")
		try:
			if winning_snps[-1]=="|":
				winning_snps = winning_snps[:-1]
		except:
			pass
		realmeasures.sort()
		permuted_indm_p.sort(reverse=False)
		fdr_permuted_indm.sort(reverse=True)
		outputfilecontents+= "%s [%s]\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
		(wname, winning_snps, snp_string, winners[group_label], m, min(theseps), indmeasure[key],\
		permuted_indm_p[i],fdr_BH_permuted_indm[i])
		fdrstore[key] = fdr_BH_permuted_indm[i]

		if fdr_BH_permuted_indm[i] < 0.05:
			#highlight this hit in the output html file
			h.write("<tr class=\"sig\">\n")
			#and highlight all promoters in this group in the biolayout network
			highlight.append(winners[group_label])
			for prom in prom_to_join[group_label]:
				highlight.append(prom)
		else:
			h.write("<tr>\n")
		h.write("<td><a href=\"%s%s\" target=\"_blank\">%s</a></td><td>%s</td><td>%s</td><td>%s</td></tr>\n"\
		%(f5resource, winners[group_label], wname, winning_snps.replace("|",", "),\
		round(indmeasure[key],1), round(fdr_BH_permuted_indm[i],2)))
		i+=1
	h.write("</table>")
	h.close()
	o=open(os.path.join(working_files_dir,fullresults),"w")
	o.write(outputfilecontents)
	o.close()
	o2=open(os.path.join(working_files_dir,fullresults),"w") #duplicate in offline directory
	o2.write(outputfilecontents)
	o2.close()

	#-------------------------------------#
	#write a summary of the coex evidence at every location considered
	allscoresfile = os.path.join(working_files_dir, "allscores.txt")
	allscoresoutlist=[]
	for prom in joining_instructions:
		allscoresoutlist.append([prom_addresses[prom][0], prom_addresses[prom][1],prom_addresses[prom][2], prom, joining_instructions[prom], allpromoterscores[prom], indmeasure[joining_instructions[prom]], fdrstore[joining_instructions[prom]] ])
	allscoresoutlist.sort(key = lambda x: int(x[1]) ) #sort by start position	
	allscoresoutlist.sort(key = lambda x: int(x[0]) ) #sort by chrom
	o=open(allscoresfile,'w')
	o.write("chrom\tstart\tend\tname\tgroup_name\traw_coex\tgroup_coex\tgroup_FDR\n")
	for thisline in allscoresoutlist:
		o.write("%s\n"%('\t'.join([str(x) for x in thisline])))
	o.close()
    
	
	#-------------------------------------#
	# write big network file (for circos)
	networkname = "%s_%s.network.php"%(supplementary_label, suffix)
	networkfile = os.path.join(working_files_dir, networkname)
	if write_layout_file == "yes":
		nftext=''
		already={}
		#write a line for each pair of promoters with edge between them
		for prom1 in individualscores:
			for prom2 in individualscores[prom1]:
				if prom1 != prom2:
					prom1sub = joining_instructions[prom1] #substitute group label
					prom2sub = joining_instructions[prom2] #substitute group label
					if prom1sub != prom2sub:
						c="%s"%("\t".join(sorted([prom1,prom2])))
						try:
							edge = indjoined[prom1sub][prom2sub] #-math.log(joined_ps[prom1][prom2][7],10)
							try:
								already[c]
							except:
								already[c]=[]
							already[c].append(edge)
						except:
							try:
								indjoined[prom1sub]
								print prom1sub, "in indjoined"
								try:
									indjoined[prom1sub][prom2sub]
									print prom2sub, "in indjoined[%s]"%prom1sub
								except:
									print prom2sub, "not in indjoined[%s]"%prom1sub
							except:
								print prom1sub, "not in indjoined"
		o=open(networkfile,'w')
		o.write('<? header("Access-Control-Allow-Origin: http://coexpression.net")?>\nsource\ttarget\tvalue\n')
		for c in already:
			avedge = float(sum(already[c]))/len(already[c])
			o.write("%s\t%s\n"%(c,avedge))
		o.close()

	#-------------------------------------#
	#write   layout   file # also used by d3
	layoutname = "%s_%s.layout"%(supplementary_label, suffix)
	layoutfile = os.path.join(working_files_dir, layoutname)
	if write_layout_file == "yes":
		lftext = ''
		lftext+="//BioLayout Express 3D Version 2.1  Layout File\n"
		already=[]
		included=[]
		#write a line for each pair of promoters with a strong edge between them
		for prom1 in indjoined:
			for prom2 in indjoined[prom1]:
				prom1sub = joining_instructions[prom1] #substitute group label
				prom2sub = joining_instructions[prom2] #substitute group label
				prom1winner = winners[prom1sub]
				prom2winner = winners[prom2sub]
				edgelist = []
				try:
					edgelist.append(indjoined[prom1][prom2]) #-math.log(joined_ps[prom1][prom2][7],10)
				except:
					pass
				try:
					edgelist.append(indjoined[prom2][prom1]) #-math.log(joined_ps[prom1][prom2][7],10)
				except:
					pass
				edge = float(sum(edgelist))/len(edgelist) #average of 2 edges (these are different for nsp)
				#printv (prom1, prom2sub, edge, selectionthreshold, edge >= selectionthreshold)
				selthreshforbiolayout = -1 #selectionthreshold
				if edge >= selthreshforbiolayout:
					combo=[prom1,prom2sub]
					combo.sort()
					c="%s_%s"%(combo[0],combo[1])
					if c not in already:
						lftext+=('\"%s\"\t\"%s\"\t%s\n'%(pseudonyms[prom1winner],pseudonyms[prom2winner],edge))
						included.append(prom1winner)
						included.append(prom2winner)
		included=set(included)
		#WRITE NODESIZES
		for entry in included:
			theseps2=[]
			try:
				snps = snp_map[entry]["snps"]
				wslist = string.split(snps,"|")
				wslist = [w for w in wslist if w!=""]
				for snp in wslist:
					theseps2.append(snp_p_values[snp])
			except:
				# if any snp
				theseps2.append(1)
			nodesize = math.sqrt(max(theseps2))*10
			lftext+=('//NODESIZE\t\"%s\"\t%s\n'%(pseudonyms[entry],nodesize))
			#WRITE NODECOLOURS
			status = "FALSE"
			if entry in highlight:
				status = "TRUE"
			lftext+=('//NODECLASS\t\"%s\"\t%s\t\"significantly_coexpressed\"\n'%(pseudonyms[entry],status))

		lftext+=("//EDGECOLOR\t\"#FF0033\"\n")
		lftext+=("//EDGESIZE\t2\n")
		lftext+=("//NODECLASSCOLOR\t\"TRUE\"\t\"significantly_coexpressed\"\t\"#0000FF\"\n")
		lftext+=("//NODECLASSCOLOR\t\"FALSE\"\t\"significantly_coexpressed\"\t\"#BAAAD1\"\n")
		lftext+=("//CURRENTCLASSSET\t\"significantly_coexpressed\"\n")
		lftext+=("//DEFAULTSEARCH\t\"%s\"\n"%(f5resource))
		o=open(layoutfile,"w")
		o.write(lftext)
		o.close()
		o=open(os.path.join(working_files_dir, layoutname),"w")
		o.write(lftext)
		o.close()
		os.system("zip %s.zip %s"%(layoutfile,layoutfile))


	#------------------------
	#calculate stats for shift = 0
	f=open(mapfiles[0])
	lines=f.readlines()
	f.close()
	tss_snps = []
	proms = []
	for line in lines:
		line=string.split(string.strip(line),"\t")
		proms.append(line[-1])
		tss_snps.append(line[3])
	proms = list(set(proms))
	tss_snps = list(set(tss_snps))
	hpm_genome = (snp_count*1000000.0/genome_length*1.0)
	min_feature_dict = remove_overlap(feature_dict)
	range_length = 0
	for chrom in min_feature_dict:
		for start in min_feature_dict[chrom]:
			range_length += min_feature_dict[chrom][start] - start
	try:
		hpm=(float(len(proms))*1000000) / (float(range_length))
	except:
		hpm="div0"

	htmlstatsfile = supplementary_label+"_stats.html"
	o=open(os.path.join(working_files_dir,htmlstatsfile),"w")
	o.write("snps_searched\t%s\n"%(snp_count))
	o.write("tss_snps\t%s\n"%(len(tss_snps)))
	o.write("hpm_genome\t%s\n"%(hpm_genome))
	o.write("hpm_tss\t%s\n"%(hpm))
	o.write("permutationmode\t%s\n"%(permutationmode))
	o.write("numperms\t%s\n"%(numperms))
	o.write("nsp\t%s\n"%(nsp))
	# per range results
	o.write("allproms\t%s\n"%(len(all_proms)))
	o.write("distinct\t%s\n"%(len(prom_to_join)))
	o.close()
	#------------------------
	o=open(os.path.join(working_files_dir,"%s.scatter"%(supplementary_label)),"w")
	o.write("Spearman_r\tGlobal_p\tSpecific_p\n")
	for i in range(len(scatter_pearson)):
		o.write("%s\t%s\t%s\n"%(scatter_pearson[i],scatter_globalp[i],scatter_nsp[i]))
	o.close()
	#------------------------
	filelist = "%s.filelist"%(supplementary_label)
	o=open(os.path.join(working_files_dir,"%s.filelist"%(supplementary_label)),"w")
	o.write("qqfile\t%s.svg\n"%(qqfile))
	o.write("layoutfile\t%s\n"%(layoutname))
	o.write("statsfile\t%s\n"%(htmlstatsfile))
	o.write("fullresults\t%s\n"%(fullresults))
	o.write("htmlresults\t%s\n"%(htmlresults))
	o.write("email\t%s\n"%(useremail))
	o.close()

	for filename in [htmlstatsfile, htmlresults, filelist]:
		#TEXT FILES READABLE BY PHP GO IN WORKING FILES DIR
		filename = os.path.join(working_files_dir, filename)
		fix_permissions(filename)

	for filename in ["%s.svg"%(qqfile), "%s.zip"%(layoutname), fullresults]:
		#IMAGE FILES AND ZIP FILES GO IN IMAGE DIR
		filename = os.path.join(working_files_dir, filename)
		fix_permissions(filename)

	r = open(recordfile,"a")
	r.write("%s\n"%(common_label))
	r.close()

	import smtplib
	from email.mime.text import MIMEText
	msg = MIMEText("Your script (%s) at coexpression.net has completed. Follow this link to see the results:\n http://coexpression.net/view_results.php?id=%s"%(common_label, supplementary_label))
	sender = "autosend@coexpression.net"
	recipient = useremail
	msg['Subject'] = 'Script completed'
	msg['From'] = sender
	msg['To'] = recipient

	try:
		s = smtplib.SMTP('localhost')
		s.sendmail(sender, [recipient], msg.as_string())
		s.sendmail(sender, [sender], msg.as_string())
		s.quit()
	except:
		pass

	#clean up randomisation files
	for shift in mapfiles:
		permfile = mapfiles[shift]
		#os.system("rm %s"%(permfile))

	#remove queue file
	try:
		os.system("rm %s"%(thisqueue))
	except:
		pass
	toc = timeit.default_timer()
	totaltime = toc - tic
	o=open(os.path.join(working_files_dir,'time.txt'),'w')
	o.write("%s\t%s\t%s\n"%(runcommand, numthreads,totaltime))
	o.close()


#------------
check_dir(resdir)
check_dir(working_files_dir)
t=open(thisqueue, "w")
t.write("Job started %s"%(str(datetime.datetime.now())))
t.close()
#------------
if verbose:
	runall()
else:
	try:
		runall()
		f=open(recordfile,"a")
		f.write("%s\n"%common_label)
		f.close()
	except:
		fail("runall failed")

end = True #instruct the alive function to cease.
