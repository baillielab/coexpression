#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from coexfunctions import *

version = getversion()

import sys, string, os
import json

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import argparse

parser = argparse.ArgumentParser()
# always needed
parser.add_argument('-wd', '--working_files_dir', type=str, default='null', help='filepath')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# usually default
parser.add_argument('-t', '--numthreads', type=int, default=4, help='numthreads used')
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
useremail = storedesettings['useremail']
snp_count = int(storedesettings['snp_count'])
genome_length = int(storedesettings['genome_length'])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import inspect
import ConfigParser

pathtocoexapp = inspect.stack()[0][1]
coexappdir = os.path.split(os.path.abspath(pathtocoexapp))[0]
config = ConfigParser.RawConfigParser()
config.read(os.path.join(coexappdir, 'app.cfg'))
sourcefilesdir = config.get('directorypaths', 'sourcefilesdir')
resdir = config.get('directorypaths', 'resdir')
pathtobedtools = config.get('directorypaths', 'pathtobedtools')
f5resource = config.get('directorypaths', 'f5resource')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
feature_label = string.split(string.split(feature_coordinates, "/")[-1], ".")[
    0]  # ie the "fantom5" in the coordinates file name
storage_dir_perm_results = os.path.join(working_files_dir, "perm_results_store")
storage_dir_permutations = os.path.join(working_files_dir, "permutation_store")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
write_layout_file = "yes"  # "yes" or anything
snps_mapped = os.path.join(storage_dir_permutations, feature_label + ".%s.bed" % (0))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# read perm indmeasures
interimfiles = os.listdir(storage_dir_perm_results)
perm_indm = []
for permresultsfile in [x for x in interimfiles if ".indmeasure." in x and not x.endswith('.0')]:
    f = open(os.path.join(storage_dir_perm_results, permresultsfile))
    indmeasure = json.load(f)
    f.close()
    perm_indm += indmeasure.values()
numperms = len(perm_indm)  # the actual number of perms completed
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# read real data
realdataresults = [x for x in interimfiles if x.endswith('.0')]
filestems = [x for x in realdataresults if x.endswith('.indmeasure.0')]
if len(filestems) > 1:
    print "Error: more than one filestem in {}".format(storage_dir_perm_results)
    sys.exit()
filestem = filestems[0].replace('.indmeasure.0', '')

f = open(os.path.join(storage_dir_perm_results, filestem + '.indmeasure.0'))  # main results
indmeasure = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.individualscores.0'))  # for circos
individualscores = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.prom_to_join.0'))
prom_to_join = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.indjoined.0'))
indjoined = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.joining_instructions.0'))
joining_instructions = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.winners.0'))
winners = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.allpromoterscores.0'))
allpromoterscores = json.load(f)
f.close()
f = open(os.path.join(storage_dir_perm_results, filestem + '.snp_map.0'))
snp_map = json.load(f)
f.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

# write output files
prefix = ""
suffix = ""
outputfilelabel = "_".join([x for x in [prefix, supplementary_label, suffix] if x != ''])
# --------
realmeasures = [indmeasure[prom1] for prom1 in indjoined]
# ------------save qq_plot file-----------------
try:
    qqfile = outputfilelabel + ".qq_plot"
    o = open(os.path.join(working_files_dir, qqfile), "w")
    o.write("rand\t%s\n" % qqfile)
    limit = min(len(fullrange_perm_indm), len(realmeasures))
    realmeasures.sort()
    fullrange_perm_indm.sort()
    for i in range(limit):
        indexa = int(len(fullrange_perm_indm) * float(i) / limit)
        indexb = int(len(realmeasures) * float(i) / limit)
        a = fullrange_perm_indm[indexa]
        b = realmeasures[indexb]
        o.write("%s\t%s\n" % (a, b))
    o.close()
except:
    pass
# ------------calculate FDR----------------------
permuted_indm_p = sorted([empiricalp(x, perm_indm) for x in realmeasures], reverse=False)
fdr_permuted_indm = list(fdr_correction(permuted_indm_p, 0.05, "negcorr")[1])
fdr_BH_permuted_indm = list(fdr_correction(permuted_indm_p)[1])
# lr.close()
fdrstore = {}

# output individual scores
pseudonyms = {}
htmlresults = outputfilelabel + ".html"
h = open(os.path.join(working_files_dir, htmlresults), "w")
h.write(
    "<table id=\"results\">\n<tr>\n<th>Top promoter (click for details in new tab)</th>\n<th>SNPs in top_promoter</th>\n<th>Coexpression score</th>\n<th>FDR</th>\n</tr>")
fullresults = outputfilelabel + ".individual_scores"
outputfilecontents = ''
outputfilecontents += "Top_promoter[SNPs_in_top_promoter]\tpromoter_snps_in_region\tchosen_promoter_id\tpromoters\toriginal_p\tcoexpression_score\traw_p(perm)\tBenjamini-Hochberg\n"
# for key in indjoined:
i = 0
highlight = []

for key, value in sorted(indmeasure.iteritems(), key=lambda (k, v): (v, k), reverse=True):
    snps = []
    group_label = joining_instructions[key]
    for prom in prom_to_join[group_label]:
        try:
            snp_list = string.split(snp_map[prom]["snps"], "|")
        except:
            snp_list = ["no_snps.seed_promoter"]
        for snp in snp_list:
            if snp not in snps:
                snps.append(snp)
    snp_string = ""
    for snp in snps:
        snp_string += "|%s" % snp
    if snp_string[0] == "|":
        snp_string = snp_string[1:]
    prom_string = ""
    m = []
    for promoter in prom_to_join[joining_instructions[key]]:
        prom_string += "|%s" % promoter
        if promoter in metadata:
            m.append("%s" % metadata[promoter])
        else:
            m.append(promoter)
    m = "|".join(m)
    if prom_string[0] == "|":
        prom_string = prom_string[1:]
    wname = winners[group_label]
    if winners[group_label] in metadata:
        wname = metadata[winners[group_label]]
        wname.replace(",+", "[+]")
        wname.replace(",-", "[-]")
        wname.replace(",", ", ")
        wname.replace("[+]", ",+ ")
        wname.replace("[-]", ",- ")
    pseudonyms[winners[group_label]] = wname

    try:
        winning_snps = snp_map[winners[group_label]]["snps"]
    except:
        winning_snps = ""
    wslist = string.split(winning_snps, "|")
    wslist = [w for w in wslist if w != ""]

    minp = 1

    winning_snps = winning_snps.replace("||", "|")
    try:
        if winning_snps[-1] == "|":
            winning_snps = winning_snps[:-1]
    except:
        pass
    realmeasures.sort()
    permuted_indm_p.sort(reverse=False)
    fdr_permuted_indm.sort(reverse=True)
    outputfilecontents += "%s [%s]\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                          (wname, winning_snps, snp_string, winners[group_label], m, minp, indmeasure[key], \
                           permuted_indm_p[i], fdr_BH_permuted_indm[i])
    fdrstore[key] = fdr_BH_permuted_indm[i]

    if fdr_BH_permuted_indm[i] < 0.05:
        # highlight this hit in the output html file
        h.write("<tr class=\"sig\">\n")
        # and highlight all promoters in this group in the biolayout network
        highlight.append(winners[group_label])
        for prom in prom_to_join[group_label]:
            highlight.append(prom)
    else:
        h.write("<tr>\n")
    h.write("<td><a href=\"%s%s\" target=\"_blank\">%s</a></td><td>%s</td><td>%s</td><td>%s</td></tr>\n" \
            % (f5resource, winners[group_label], wname, winning_snps.replace("|", ", "), \
               round(indmeasure[key], 1), round(fdr_BH_permuted_indm[i], 2)))
    i += 1
h.write("</table>")
h.close()
o = open(os.path.join(working_files_dir, fullresults), "w")
o.write(outputfilecontents)
o.close()
o2 = open(os.path.join(working_files_dir, fullresults), "w")  # duplicate in offline directory
o2.write(outputfilecontents)
o2.close()

# -------------------------------------#
if True:
    # write a summary of the coex evidence at every location considered
    prom_addresses, feature_dict = readfeaturecoordinates(feature_coordinates)
    allscoresfile = os.path.join(working_files_dir, "allscores.txt")
    allscoresoutlist = []
    for prom in joining_instructions:
        allscoresoutlist.append([prom_addresses[prom][0], prom_addresses[prom][1], prom_addresses[prom][2], prom,
                                 joining_instructions[prom], allpromoterscores[prom],
                                 indmeasure[joining_instructions[prom]], fdrstore[joining_instructions[prom]]])
    allscoresoutlist.sort(key=lambda x: int(x[1]))  # sort by start position
    allscoresoutlist.sort(key=lambda x: int(x[0]))  # sort by chrom
    o = open(allscoresfile, 'w')
    o.write("chrom\tstart\tend\tname\tgroup_name\traw_coex\tgroup_coex\tgroup_FDR\n")
    for thisline in allscoresoutlist:
        o.write("%s\n" % ('\t'.join([str(x) for x in thisline])))
    o.close()

# -------------------------------------#
# write big network file (for circos)
networkname = outputfilelabel + ".network.php"
networkfile = os.path.join(working_files_dir, networkname)
if write_layout_file == "yes":
    nftext = ''
    already = {}
    # write a line for each pair of promoters with edge between them
    for prom1 in individualscores:
        for prom2 in individualscores[prom1]:
            if prom1 != prom2:
                prom1sub = joining_instructions[prom1]  # substitute group label
                prom2sub = joining_instructions[prom2]  # substitute group label
                if prom1sub != prom2sub:
                    c = "%s" % ("\t".join(sorted([prom1, prom2])))
                    try:
                        edge = indjoined[prom1sub][prom2sub]  # -math.log(joined_ps[prom1][prom2][7],10)
                        try:
                            already[c]
                        except:
                            already[c] = []
                        already[c].append(edge)
                    except:
                        try:
                            indjoined[prom1sub]
                            print prom1sub, "in indjoined"
                            try:
                                indjoined[prom1sub][prom2sub]
                                print prom2sub, "in indjoined[%s]" % prom1sub
                            except:
                                print prom2sub, "not in indjoined[%s]" % prom1sub
                        except:
                            print prom1sub, "not in indjoined"

o = open(networkfile, 'w')
o.write('<? header("Access-Control-Allow-Origin: http://baillielab.net")?>\nsource\ttarget\tvalue\n')
for c in already:
    avedge = float(sum(already[c])) / len(already[c])
    o.write("%s\t%s\n" % (c, avedge))
o.close()

# -------------------------------------#
# write   layout   file # also used by d3
layoutname = outputfilelabel + ".layout"
layoutfile = os.path.join(working_files_dir, layoutname)
if write_layout_file == "yes":
    lftext = ''
    lftext += "//BioLayout Express 3D Version 2.1  Layout File\n"
    already = []
    included = []
    # write a line for each pair of promoters with a strong edge between them
    for prom1 in indjoined:
        for prom2 in indjoined[prom1]:
            prom1sub = joining_instructions[prom1]  # substitute group label
            prom2sub = joining_instructions[prom2]  # substitute group label
            prom1winner = winners[prom1sub]
            prom2winner = winners[prom2sub]
            edgelist = []
            try:
                edgelist.append(indjoined[prom1][prom2])  # -math.log(joined_ps[prom1][prom2][7],10)
            except:
                pass
            try:
                edgelist.append(indjoined[prom2][prom1])  # -math.log(joined_ps[prom1][prom2][7],10)
            except:
                pass
            edge = float(sum(edgelist)) / len(edgelist)  # average of 2 edges (these are different for nsp)
            # if verbose: print (prom1, prom2sub, edge, selectionthreshold, edge >= selectionthreshold)
            selthreshforbiolayout = -1  # selectionthreshold
            if edge >= selthreshforbiolayout:
                combo = [prom1, prom2sub]
                combo.sort()
                c = "%s_%s" % (combo[0], combo[1])
                if c not in already:
                    lftext += ('\"%s\"\t\"%s\"\t%s\n' % (pseudonyms[prom1winner], pseudonyms[prom2winner], edge))
                    included.append(prom1winner)
                    included.append(prom2winner)
    included = set(included)

    # WRITE NODECOLOURS
    for entry in included:
        status = "FALSE"
        if entry in highlight:
            status = "TRUE"
        lftext += ('//NODECLASS\t\"%s\"\t%s\t\"significantly_coexpressed\"\n' % (pseudonyms[entry], status))

    lftext += ("//EDGECOLOR\t\"#FF0033\"\n")
    lftext += ("//EDGESIZE\t2\n")
    lftext += ("//NODECLASSCOLOR\t\"TRUE\"\t\"significantly_coexpressed\"\t\"#0000FF\"\n")
    lftext += ("//NODECLASSCOLOR\t\"FALSE\"\t\"significantly_coexpressed\"\t\"#BAAAD1\"\n")
    lftext += ("//CURRENTCLASSSET\t\"significantly_coexpressed\"\n")
    lftext += ("//DEFAULTSEARCH\t\"%s\"\n" % (f5resource))
    o = open(layoutfile, "w")
    o.write(lftext)
    o.close()
    o = open(os.path.join(working_files_dir, layoutname), "w")
    o.write(lftext)
    o.close()
    os.system("zip %s.zip %s" % (layoutfile, layoutfile))

# ------------------------
# calculate stats for shift = 0
f = open(snps_mapped)
lines = f.readlines()
f.close()
tss_snps = []
proms = []
for line in lines:
    line = string.split(string.strip(line), "\t")
    proms.append(line[-1])
    tss_snps.append(line[3])
proms = list(set(proms))
tss_snps = list(set(tss_snps))
hpm_genome = (snp_count * 1000000.0 / genome_length * 1.0)
min_feature_dict = remove_overlap(feature_dict)
range_length = 0
for chrom in min_feature_dict:
    for start in min_feature_dict[chrom]:
        range_length += min_feature_dict[chrom][start] - start
try:
    hpm = (float(len(proms)) * 1000000) / (float(range_length))
except:
    hpm = "div0"

htmlstatsfile = supplementary_label + "_stats.html"
o = open(os.path.join(working_files_dir, htmlstatsfile), "w")
o.write("snps_searched\t%s\n" % (snp_count))
o.write("tss_snps\t%s\n" % (len(tss_snps)))
o.write("hpm_genome\t%s\n" % (hpm_genome))
o.write("hpm_tss\t%s\n" % (hpm))
o.write("numperms\t%s\n" % (numperms))
o.write("nsp\t%s\n" % (nsp))
# per range results
o.write("allproms\t%s\n" % (len(allpromoterscores)))
o.write("distinct\t%s\n" % (len(prom_to_join)))
o.close()
# ------------------------
filelist = "%s.filelist" % (supplementary_label)
o = open(os.path.join(working_files_dir, "%s.filelist" % (supplementary_label)), "w")
o.write("qqfile\t%s.svg\n" % (qqfile))
o.write("layoutfile\t%s\n" % (layoutname))
o.write("statsfile\t%s\n" % (htmlstatsfile))
o.write("fullresults\t%s\n" % (fullresults))
o.write("htmlresults\t%s\n" % (htmlresults))
o.write("email\t%s\n" % (useremail))
o.close()
