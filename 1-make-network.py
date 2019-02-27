#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from coexfunctions import *
version = getversion()

import sys, timeit, string, os, operator
from ctypes import *
import json
import resource

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import argparse

parser = argparse.ArgumentParser()
# always needed
parser.add_argument('-mf', '--mapfile', type=str, default='null', help='filepath')
parser.add_argument('-wd', '--working_files_dir', type=str, default='null', help='filepath')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# usually default
parser.add_argument('-t', '--numthreads', type=int, default=1, help='numthreads used')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# optional / unnecessary
parser.add_argument('-r', '--recordedges', action="store_true", default=False, help='recordedges')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
args = parser.parse_args()
module = sys.modules[__name__]
for name, value in vars(args).iteritems():
    print name, value
    setattr(module, name, value)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GET SETTINGS
storedesettings = readsettings(working_files_dir)
# numperms = storedesettings['numperms']
verbose = storedesettings['verbose']
feature_coordinates = storedesettings['feature_coordinates']
correlationfile = storedesettings['correlationfile']
expression_file = storedesettings['expression_file']
supplementary_label = storedesettings['supplementary_label']
correlationmeasure = storedesettings['correlationmeasure']
soft_separation = storedesettings['soft_separation']
nsp = storedesettings['nsp']
seedfile = storedesettings['seedfile']
pval_to_join = storedesettings['pval_to_join']
useaverage = storedesettings['useaverage']
anticorrelation = storedesettings['anticorrelation']
noiter = storedesettings['noiter']
selectionthreshold = storedesettings['selectionthreshold']
selectionthreshold_is_j = storedesettings['selectionthreshold_is_j']
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import inspect
import ConfigParser
pathtocoexapp = inspect.stack()[0][1]
coexappdir = os.path.split(os.path.abspath(pathtocoexapp))[0]
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read(storedesettings['config_file'])
sourcefilesdir = config.get('directorypaths', 'sourcefilesdir')
resdir = config.get('directorypaths', 'resdir')
pathtobedtools = config.get('directorypaths', 'pathtobedtools')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NAME OUTPUT DIR
storage_dir_perm_results = os.path.join(working_files_dir, "perm_results_store")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# READ INPUT FILES
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sofile = os.path.join(coexappdir, 'coexpression_v2.so')
coexpression = cdll.LoadLibrary(sofile)
if verbose: print "sofile read"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exp_dict, header = read_expression_file(expression_file)
if verbose: print "expression file read"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
prom_addresses, feature_dict, sl = readfeaturecoordinates(feature_coordinates)
if verbose: print "feature_coordinates read"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
f = open(correlationfile)
lines = f.readlines()
f.close()
globalcorrelationlist = [float(x) for x in lines if string.strip(x) != ""]
if verbose: print "global correlations read"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# read seedfile
# this will seed every network with the contents of seedfile, then remove them all before calculating
lines = []
if seedfile != "none":
    if os.path.exists(seedfile):
        f = open(seedfile)
        lines = f.readlines()
        f.close()
        if verbose: print "seedfile read"
seedpromoters = []
for line in lines:
    seedpromoters.append(string.split(string.strip(line), "\t")[-1])  # last entry in line is promoter (eg bedfile)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if selectionthreshold_is_j:
    selectionthreshold = -math.log(pval_to_join, 10)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
promindex_lookup = {}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ptcf = exp_dict.keys()  # promoters to choose from
# NEWLY DEVELOPMENTS BY AG
dictlen = len(exp_dict)
entrylen = len(exp_dict[ptcf[0]])

# ctype structure for all promoter pair data
alldata_c = (c_double * (dictlen * entrylen))()

if correlationmeasure == "Spearman":  # creates an equivalent but ranked dictionary to 'exp_dict'
    exp_dict_rank = {}
    for i in exp_dict.keys():
        exp_dict_rank[i] = stats.rankdata(exp_dict[i])

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

if verbose: print "promindex made"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
nodespecificcorrelationlists = {}
# DO COEXPRESSION
permtimer = timeit.default_timer()

# READ SNPS_MAPPED
f = open(mapfile)
lines = f.readlines()
f.close()
if verbose: print ("reading snp map(BED TOOLS): %s lines" % len(lines))

snp_map = {}
for i in range(len(lines)):
    line = string.split(string.strip(lines[i]), "\t")
    if len(line) < 4:
        continue
    this_promoter = line[-1]
    try:
        exp_dict[this_promoter]
    except:
        print "skipping", this_promoter, "as not in expression file"
        continue
    snp = line[3]
    # store p-value for each snp mapped
    try:
        snp_map[this_promoter]
        snp_map[this_promoter]["snps"] += "|" + snp
    except:
        snp_map[this_promoter] = {"p": -9999, "snps": snp}
    try:
        thisp = snp_p_values[snp]
        # if this is a script2 randomisation, "rand_" is appended to the snp name
    except:
        try:
            thisp = snp_p_values[snp.replace("rand_", "")]
        except:
            # ie if there is no p-value for this "snp" - might be a new bonus promoter
            thisp = 1000
    if thisp < snp_map[this_promoter]["p"]:
        snp_map[this_promoter]["p"] = thisp

if verbose: print ('snp map read.')
# -------------------------------------#
inputpromoterlist = snp_map.keys()
for new in seedpromoters:
    if new not in snp_map.keys():
        snp_map[new] = {"p": 1000.125}  # specific p-value to label these additions
# -------------------------------------#
edges_sought = 0

all_proms = []
individualscores = {}
promoters = snp_map.keys()

# if len(promoters)==0:
# print "not enough promoters"
# sys.exit()

# AG developments
promaidx_all = []
prombidx_all = []

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use_specific_p = True
if nsp <= 1.0 / len(ptcf):
    use_specific_p = False
    print "nsp", nsp, "is too small because there are only", len(ptcf), "promoters. It must be greater than", 1.0 / len(
        ptcf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tta = timeit.default_timer()
if use_specific_p:
    tt_p1 = timeit.default_timer()
    if nsp < 1:
        ptcf_sample = random.sample(ptcf, int(float(len(ptcf)) * nsp))
    else:
        ptcf_sample = ptcf
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CALCULATE ALL CORRELATIONS
    if verbose: print "starting phase 1"
    tt_c1 = timeit.default_timer()
    for i in range(len(promoters)):
        jrange = range(len(ptcf_sample))  # asymmetrical matrix for node specific pval
        for j in jrange:
            if i != j:
                promaidx_all.append(promindex_lookup[promoters[i]])
                prombidx_all.append(promindex_lookup[ptcf_sample[j]])
    if verbose: print ("added", len(prombidx_all), "correlations to correlation index")
    tt_c2 = timeit.default_timer()
    if verbose: print "phase 1 done", tt_c2 - tt_c1

    if verbose: print "starting phase 2"
    nPromPairs = len(promaidx_all)
    idxa_all_c = (c_int * nPromPairs)(*promaidx_all)
    idxb_all_c = (c_int * nPromPairs)(*prombidx_all)
    tt_c3 = timeit.default_timer()
    if verbose: print "phase 2 done", tt_c3 - tt_c2

    # ctype data structure to hold results
    result_all_c = (c_double * nPromPairs)()

    # get all the coefficients
    coexpression.getPearson(byref(alldata_c), \
                            byref(result_all_c), \
                            byref(idxa_all_c), \
                            byref(idxb_all_c), \
                            c_int(nPromPairs), \
                            c_int(entrylen))
    tt_c4 = timeit.default_timer()
    if verbose: print "phase 3 done", tt_c4 - tt_c3

    if verbose: print "starting phase 4"
    # set the NANs to -99, and store all in python list
    AGcorrelations = []
    for x in range(nPromPairs):
        if math.isnan(result_all_c[x]):
            correlation = -99
        else:
            correlation = result_all_c[x]
        AGcorrelations.append(correlation)
    tt_c5 = timeit.default_timer()
    if verbose: print "phase 4 done", tt_c5 - tt_c4

    # ++++++++++++++++++++++++++++++++++++
    ip = 0
    for i in range(len(promoters)):
        nodespecificcorrelationlists[promoters[i]] = []
        for j in jrange:
            if i != j:
                nodespecificcorrelationlists[promoters[i]].append(AGcorrelations[ip])
                ip = ip + 1
        nodespecificcorrelationlists[promoters[i]].sort()  # this is essential so that p-vals can be calculated.
    tt_c6 = timeit.default_timer()
    if verbose: print "phase 5 done", tt_c6 - tt_c5
    # ++++++++++++++++++++++++++++++++++++
    if recordedges:
        scatter_pairs.append("realdatanow")
        scatter_pearson.append("realdatanow")
        scatter_globalp.append("realdatanow")
        scatter_nsp.append("realdatanow")
        # ++++++++++++++++++++++++++++++++++++
    tt_p2 = timeit.default_timer()
    if verbose: print "pval prep done in:", tt_p2 - tt_p1
ttb = timeit.default_timer()
ip = 0

for i in range(len(promoters)):
    tt_l1 = timeit.default_timer()
    for j in range(len(promoters)):  # NB must do this from both sides for nsp (ie not for j in range(i,len...))
        if i != j:
            prom1 = promoters[i]
            prom2 = promoters[j]
            # correlation=AGcorrelations[ip]
            # correlation = AGcordict[prom1][prom2]
            correlation = get_correlation(prom1, prom2, exp_dict, verbosefunction=True)
            # ======================================
            if use_specific_p:
                if correlation == -99:
                    p1 = 1
                else:
                    p1 = 1 - float(bisect_left(nodespecificcorrelationlists[prom1], correlation)) / len(
                        nodespecificcorrelationlists[prom1])
                    if anticorrelation:
                        p2 = float(bisect_left(nodespecificcorrelationlists[prom1], correlation)) / len(
                            nodespecificcorrelationlists[prom1])
                        p1 = min(p1, p2)
            else:  # legacy code. Enables calculation of whole-dist p
                # WHOLE-DISTRIBUTION P VALUE
                if correlation == -99:
                    p1 = 1
                else:
                    p1 = 1 - float(bisect_left(globalcorrelationlist, correlation)) / len(globalcorrelationlist)
                    if anticorrelation:
                        p2 = float(bisect_left(globalcorrelationlist, correlation)) / len(globalcorrelationlist)
                        p1 = min(p1, p2)
                        # ======================================
            if recordedges:
                scatter_pairs.append("%s\t%s" % (prom1, prom2))
                scatter_correlation.append(correlation)
                scatter_globalp.append(
                    1 - float(bisect_left(globalcorrelationlist, correlation)) / len(globalcorrelationlist))
                scatter_nsp.append(1 - float(bisect_left(nodespecificcorrelationlists[prom1], correlation)) / len(
                    nodespecificcorrelationlists[prom1]))
                # ======================================
            try:
                log_edge_p = -math.log(p1, 10)
            except:
                log_edge_p = 0
                # STORE P VALUES
            try:
                individualscores[prom1][prom2] = log_edge_p
            except:
                individualscores[prom1] = {}
                individualscores[prom1][prom2] = log_edge_p
            if not (nsp):
                try:
                    individualscores[prom2][prom1] = log_edge_p
                except:
                    individualscores[prom2] = {}
                    individualscores[prom2][prom1] = log_edge_p
            ip = ip + 1
    tt_l2 = timeit.default_timer()
    if verbose: print 'loop done. [{} of {}] [{} of {}] time={}'.format(i, len(promoters), j, len(promoters), tt_l2 - tt_l1)

ttc = timeit.default_timer()
# -------------------------------------#
for prom in individualscores:
    # just to make sure no promoters are left out
    all_proms.append(prom)
    for other_prom in individualscores[prom]:
        all_proms.append(other_prom)
all_proms = list(set(all_proms))
# 
# check if this made a difference
commontoboth = set(all_proms) & set(promoters)
if verbose: print "promoters:%s common:%s all_proms:%s " % (len(promoters), len(commontoboth), len(all_proms))
# 

# -------------------------------------#
# work out which promoters should be joined as one.
prom_to_join = join_nearby(all_proms, prom_addresses, globalcorrelationlist, exp_dict, pval_to_join, correlationmeasure,
                           soft_separation, anticorrelation)
joining_instructions = {}
for group_label in prom_to_join:
    for prom in prom_to_join[group_label]:
        joining_instructions[prom] = group_label
if verbose: print (len(all_proms), "promoters hit, in", len(prom_to_join), "distinct regions")
# -------------------------------------#
if verbose: print ("merging ind scores")
# NOW MERGE INDIVIDUAL SCORES.
indjoined = {}
indmeasure = {}
whittledindmeasure = {}
weights = {}
winners = {}

allpromoterscores = {}

for group_label in prom_to_join:
    scores = {}
    # choose the contender from each group with the most unlikely coexpression anywhere in the network of __all__ promoters
    for prom in prom_to_join[group_label]:
        if joining_instructions[prom] == group_label:
            scores[prom] = sum([individualscores[prom][otherprom] for otherprom in \
                                individualscores[prom] if individualscores[prom][otherprom] >= selectionthreshold \
                                and joining_instructions[otherprom] != group_label])
    if len(scores) == 0:
        continue
    winners[group_label] = max(scores.iteritems(), key=operator.itemgetter(1))[0]
# ____________________________________________________________________________________
# new in v0.63: repeat choice of winners using the new consensus to reduce the influence of large groups of promoters/enhancers
for group_label in prom_to_join:
    scores = {}
    checkscorelens = []
    # choose the contender from each group with the most unlikely coexpression with another group
    for prom in prom_to_join[group_label]:
        if joining_instructions[prom] == group_label:
            thisscore = []
            for otherprom in winners.values():
                if prom != otherprom and joining_instructions[otherprom] != group_label:
                    try:
                        thisscore.append(individualscores[prom][otherprom])
                    except:
                        # 
                        # 
                        # 
                        print "prom lookup failure. output dump follows "
                        try:
                            individualscores[prom]
                        except:
                            print "**** %s **** not in individual_scores" % prom
                        print prom, otherprom
                        print joining_instructions[prom]
                        print joining_instructions[otherprom]
                        print "====individual scores==="
                        for x in individualscores:
                            for y in individualscores[x]:
                                print x, y, individualscores[x][y]
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
                        # 
                        # 
                        # 
            scores[prom] = sum(thisscore)
            allpromoterscores[prom] = scores[prom]  # store for output later
            checkscorelens.append(len(thisscore))
    if checkscorelens.count(checkscorelens[0]) != len(checkscorelens):
        if verbose: print("checkscorelens not equal!",
                          set(checkscorelens))  # each group should be scored against the same number of other groups.
    if len(scores) == 0:
        continue
    new_winner = max(scores.iteritems(), key=operator.itemgetter(1))[0]
    if new_winner != winners[group_label]:
        if verbose: print("==> new winner:%s (overtook %s)" % (new_winner, winners[group_label]))
    winners[group_label] = new_winner

# ____________________________________________________________________________________
for prom in individualscores:
    group_label = joining_instructions[prom]
    chosen = winners[group_label]  # the chosen promoter for this group
    if group_label not in indjoined:
        indjoined[group_label] = {}
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
# _______________________________________________________________________________________

# now calculate the coexpression scores
for prom in indjoined:
    if useaverage:
        indmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in \
                                indjoined[prom] if indjoined[prom][otherprom] >= selectionthreshold]) / len(indjoined)
    else:
        indmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in \
                                indjoined[prom] if indjoined[prom][otherprom] >= selectionthreshold])

cutlist = []  # list of promoters to remove from analysis
if noiter == False:  # cycle down the list removing top hit and recalculating indmeasure
    lastscore = -1
    for prom, score in sorted(indmeasure.iteritems(), key=lambda (k, v): (v, k), reverse=True):
        if score == lastscore:
            whittledindmeasure[prom] = lastresult
        else:
            if useaverage:
                whittledindmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in \
                                                indjoined[prom] if (indjoined[prom][otherprom] >= selectionthreshold and \
                                                                    otherprom not in cutlist)]) / len(indjoined)
            else:
                whittledindmeasure[prom] = sum([indjoined[prom][otherprom] for otherprom in \
                                                indjoined[prom] if (indjoined[prom][otherprom] >= selectionthreshold and \
                                                                    otherprom not in cutlist)])
        cutlist.append(prom)
        lastscore = score
        lastresult = whittledindmeasure[prom]
    indmeasure = whittledindmeasure
# _______________________________________________________________________________________
seedstoprune = [x for x in seedpromoters if x not in inputpromoterlist]
for prom in seedstoprune:
    if verbose: print ('PRUNING SEED:', prom)
    puregroup = [x for x in prom_to_join[joining_instructions[prom]] if x not in seedstoprune]
    if len(puregroup) > 0:
        if verbose: print ("puregroup", puregroup)
    if len(puregroup) == 0:
        # then there are no input promoters in this group
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
# _______________________________________________________________________________________

permtimer2 = timeit.default_timer()
if verbose: 
    print "prep", tta - permtimer
    print "loop1", ttb - tta
    print "loop2", ttc - ttb
    print "end", permtimer2 - ttc
    print "all", permtimer2 - permtimer

if verbose: print ("********** 1-make-network.py results **********")
if verbose: print ("indmeasures:", [round(x, 1) for x in sorted(indmeasure.values())])

# _______________________________________________________________________________________

# output the findings for this network run - to be collated by 2-collate-results.py

objectdump = {
    'indmeasure': indmeasure,
    'individualscores': individualscores,
    'prom_to_join': prom_to_join,
    'indjoined': indjoined,
    'joining_instructions': joining_instructions,
    'winners': winners,
    'allpromoterscores': allpromoterscores,
    'snp_map': snp_map,
}

for label, var in objectdump.iteritems():
    filename = string.split(os.path.split(mapfile)[-1].replace('.bed', ''), '.')
    filename = filename[:-1] + [label] + filename[-1:]
    interimresultsfile = os.path.join(storage_dir_perm_results, '.'.join(filename))
    o = open(interimresultsfile, 'w')
    json.dump(objectdump[label], o)
    o.close()
