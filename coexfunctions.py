#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
version = 'v1.21'

# duplicated imports
import math
import os
import string
from bisect import *

import matplotlib.pyplot as plt
# unique imports
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats


def getversion():
    return version


def fix_permissions(this_path):
    os.system("/bin/chmod 755 %s" % (this_path))


def check_dir(this_dir):
    if not os.path.isdir(this_dir):
        os.mkdir(this_dir)
    fix_permissions(this_dir)


# ------------------------------------
def readsettings(wfd):
    f = open(os.path.join(wfd, 'settings.txt'))
    lines = [string.split(string.strip(x), '\t') for x in f.readlines()]
    f.close()
    thesesettings = {}
    for line in lines:
        if len(line) == 0:
            continue
        if len(line) == 1:
            thesesettings[line[0]] = ''
        if len(line) > 1:
            thesesettings[line[0]] = '\t'.join(line[1:])
    for s in thesesettings:
        try:
            thesesettings[s] = float(thesesettings[s])
        except:
            pass
    for s in thesesettings:
        try:
            thesesettings[s] = str2bool(thesesettings[s])
        except:
            pass
    return thesesettings


def str2bool(v):
    if v.lower() == "false":
        return False
    elif v.lower() == "true":
        return True
    else:
        return v


def readoptions(cmd):
    cmd = [x for x in string.split(cmd, ' ') if x != ""]
    if 'python' in cmd[0]:
        cmd = cmd[1:]
    if cmd[0].endswith('.py'):
        cmd = cmd[1:]
    opts = {}
    i = -1
    while i < (len(cmd) - 1):
        i += 1
        if cmd[i].startswith('-'):
            if i >= (len(cmd) - 1):
                opts[cmd[i]] = True
                continue
            elif cmd[i + 1].startswith('-'):
                opts[cmd[i]] = True
                continue
            else:
                opts[cmd[i]] = cmd[i + 1]
                i += 1
        else:
            opts[cmd[i]] = True
    return opts


# ------------------------------------

def get_gwad(chrom, address, masterchromrange):
    '''get a single number to identify this position on a concatentated genome'''
    try:
        masterchromrange[chrom]
    except:
        print "chrom not in masterchromrange", chrom
        return "chrom not in masterchromrange"
    if (masterchromrange[chrom][0] + address) < masterchromrange[chrom][1] or address < 0:
        return masterchromrange[chrom][0] + address
    else:
        if verbose: print ("error finding address - address is not on chromosome", chrom, address)


def get_address(gwad, masterchromrange, masterchromlist):
    '''go from gwad to chrom and position'''
    genome_max = masterchromrange[masterchromlist[-1]][1]  # the highest value of all
    # first bring the value into range
    while gwad > genome_max:
        gwad = gwad - genome_max
    for chrom in masterchromlist:
        start = masterchromrange[chrom][0]
        end = masterchromrange[chrom][1]
        if gwad > start and gwad <= end:
            address = gwad - start
            return chrom, address
    if verbose: print ("***** ERROR NEEDS ATTENTION: address not found for gwad:", gwad)
    return "chr1", 0


def read_expression_file(filename, datastartcol='autodetect', datastartrow='autodetect', labelcol=0, labelrow=0):
    '''returns dict of {name:[list_of_floats], ...}, header row, startcol, startrow'''
    f = open(filename)
    lines = [string.split(string.strip(x), '\t') for x in f.readlines()]
    f.close()
    checklines = lines[:10]  # these will be used to determine dsc and dsr.
    firstfloatfailcols = []
    for r, row in enumerate(checklines):
        for c in range(len(row) - 1, -1, -1):
            try:
                float(row[c])
            except:
                if row[c] != "NA":
                    # print "unfloatable", "row",r,"col", c, checklines[r][c]
                    firstfloatfailcols.append(c + 1)
                    break
    if len(firstfloatfailcols) > 0:
        dsc = firstfloatfailcols[-1]
        for r in range(len(firstfloatfailcols) - 1, -1, -1):
            if firstfloatfailcols[r] != dsc:
                # print "first floatfailcol inconsistency indicating data starts above row", r, firstfloatfailcols[r]
                dsr = r + 1
                break
        try:
            dsr
        except:
            dsr = 0  # if data starts in row 0 then this won't be detected above
    else:
        dsr = 0  # necessary in case firstfloatfailcols is empty
        dsc = 0
    if datastartcol == 'autodetect':
        datastartcol = dsc
    else:
        datastartcol = int(datastartcol)
        if datastartcol != dsc:
            print "datastartcol (%s) discrepancy with autodetect (%s)" % (datastartcol, dsc)
    if datastartrow == 'autodetect':
        datastartrow = dsr
    else:
        datastartrow = int(datastartrow)
        if datastartrow != dsr:
            print "datastartrow (%s) discrepancy with autodetect (%s)" % (datastartrow, dsr)
    print "reading expression file from column:%s and row:%s" % (dsc, dsr)
    exdic = {}
    linelens = []
    for i, line in enumerate(lines[datastartrow:]):
        linelens.append(len(line))
        try:
            exdic[line[labelcol]]
            print "duplicate label in expression file! ({}) taking first value".format(line[labelcol])
            continue
        except:
            exdic[line[labelcol]] = []
        for j, x in enumerate(line[datastartcol:]):
            try:
                exdic[line[labelcol]].append(float(x))
            except:
                exdic[line[labelcol]].append('NA')
    linelens = sorted(list(set(linelens)))
    if len(linelens) > 1:
        print "Be aware: {} different line lengths in expressin file ({})".format(len(linelens),
                                                                                  ','.join([str(x) for x in linelens]))
    return exdic, lines[labelrow][datastartcol:]


def readfeaturecoordinates(featfile):
    f = open(featfile)
    lines = f.readlines()  # range_features file
    f.close()
    go = "no"
    pa = {}
    fd = {}
    # prog=range(0,len(lines),max(int(float(len(lines))/10),1))
    for i in range(len(lines)):  # each line is a feature (ie a TSS region)
        # if i in prog: print i, "of", len(lines)
        line = lines[i]
        if go == "no":
            # skip past header information
            try:
                test = string.split(string.strip(line), "\t")
                int(test[2])
                go = "yes"
            except:
                continue
        line = string.split(string.strip(line), "\t")
        label = line[-1]
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        pa[label] = [chrom, start, end]
        st = min(start, end)
        en = max(start, end)
        try:
            fd[chrom]
        except:
            fd[chrom] = {}
        try:
            fd[chrom][st]
            if fd[chrom][st] < en:  # choose the greater end
                fd[chrom][st] = en
        except:
            fd[chrom][st] = en
    return pa, fd


def join_nearby(thesepromoters, theseaddresses, globcors, this_exp_dict, pvtj=0.1, cormeasure="Spearman",
                softsep=100000, anticor=False):
    ''' Creates a dictionary of the form (sorted lists):
    {prom10: [prom10, prom11 with link to prom10 etc... ],
    prom20: [prom20, prom21 with link to prom20 etc... ]} '''
    ptj = {}
    promGrpGraph = nx.Graph()  # creates graph promGrpGraph
    # creates DataFrame: index = promoter name, columns = chromosome, start, end
    a = pd.DataFrame({k: theseaddresses[k] for k in thesepromoters}, index=['chromosome', 'start', 'end']).T
    a.sort(['chromosome', 'start', 'end'], inplace=True)  # sorts DataFrame inplace - VERY IMPORTANT
    for chrom in pd.unique(a.chromosome):
        grps = np.fromiter(merge_ranges(a[a.chromosome == chrom][['start', 'end']].values, softsep), dtype='i8,i8')
        for grp_start, grp_end in grps:
            promGrpGraph.clear()  # clears all the nodes and edges from promGrpGraph
            promGrp = list(a[(a.chromosome == chrom) & (a.start >= grp_start) & (
            a.start < grp_end)].index)  # list of promoters in the group
            promGrpGraph.add_nodes_from(promGrp)  # add the nodes to the graph promGrpGraph
            # only need to calculate triangular array. (Pearson / Spearman-r are symmetric corr(X,Y)==corr(Y,X))
            for i in range(len(promGrp)):
                for j in range(i + 1, len(promGrp)):  # +1 as do not need to calculate the diagaonal itself
                    # ----
                    r = get_correlation(promGrp[i], promGrp[j], this_exp_dict, cormeasure)
                    # WHOLE-DISTRIBUTION P VALUE
                    p1 = 1 - float(bisect_left(globcors, r)) / len(globcors)
                    if anticor:
                        p2 = float(bisect_left(globcors, r)) / len(globcors)
                        p = min(p1, p2)
                    else:
                        p = p1
                    # ----
                    if p < pvtj:
                        promGrpGraph.add_edge(promGrp[i], promGrp[j])  # add edges to the graph
            for i in list(nx.connected_components(promGrpGraph)):
                new = sorted(list(i))  # sort list. Is this necessary?
                ptj[new[0]] = new
    return ptj


def listitem(thislist, index):  # FOR RANDOMISATION AT MAPPING
    '''correct list index to be within list range'''
    a = float(index) / len(thislist)
    b = a - int(round(a, 0))
    c = b * len(thislist)
    d = int(round(c, 0))  # to correct python int output
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
    if len(list) < 1:
        return "empty_list"
    if tot == "len":
        tot = len(list)
    if tot < len(list):  # impossible if list is subset of tot!
        tot = len(list)
    index = int(percent * float(tot))  # the index in an imaginary complete list
    index = index - (tot - len(list))  # subtract the non-existent values
    if index < 0:
        return "<%s" % list[0]
    try:
        list.sort()
        return list[index]
    except:
        return "centile_error"


def empiricalp(thisvalue, empiricalcollection):
    if len(empiricalcollection) == 0:
        return 1
    empiricalcollection.sort()
    return 1 - float(bisect_left(empiricalcollection, thisvalue)) / len(empiricalcollection)


def get_correlation(pa, pb, ed, cm='Spearman', verbosefunction=False):
    if len(ed[pa]) == 0 or len(ed[pb]) == 0:
        print "empty list in get_correlation", pa, pb
    if sum(ed[pa]) > 0 and sum(ed[pb]) > 0:
        result = stats.spearmanr(ed[pa], ed[pb])[0]
        if cm == "Pearson":
            try:
                result = stats.pearsonr(ed[pa], ed[pb])[0]
            except:
                if verbosefunction: print("stats.pearsonr fail for ", pa, pb)
                result = -99
        elif cm == "Spearman":
            try:
                result = stats.spearmanr(ed[pa], ed[pb])[0]
            except:
                if verbosefunction: print("spearmanr fail for ", pa, pb)
                result = -99
    else:
        if verbosefunction: print("EXPRESSION DATA is zero in all samples for one of these:", pa, pb)
        result = -99
    if math.isnan(result):
        if verbosefunction: print("isnan for ", cm, pa, pb)
        return -99
    return result


def merge_ranges(ranges, ss=100000):
    ''' Merges overlapping / adjacent ranges. 'ranges' must be sorted BEFORE
    running this function. Call with 'np.fromiter(merge_ranges(ranges),
    dtype='i8,i8')'. '''
    ranges = ranges + [0, +ss]  # add soft_separation
    ranges = iter(ranges)  # generate iterator
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > current_stop:  # Gap between ranges
            yield current_start, current_stop  # Output current range
            current_start, current_stop = start, stop  # Start a new one.
        else:  # Range overlapping (or adjacent) therefore merge.
            current_stop = max(current_stop, stop)  # 'stop' not necessarily larger than 'current_stop'
    yield current_start, current_stop


def rankdict(thisdict):
    thisdict_rank = {}
    for i in thisdict.keys():
        thisdict_rank[i] = stats.rankdata(thisdict[i])
    return thisdict_rank


# -#################################################################################

def remove_overlap(thisdict):
    '''thisdict has format thisdict[chrom][start]=end'''
    newdict = {}
    for chrom in thisdict:
        newdict[chrom] = {}
        current_chrom_starts = sorted(thisdict[chrom].keys())
        last_start = current_chrom_starts[0]
        last_end = thisdict[chrom][last_start]
        breakfound = "no"
        for this_start in current_chrom_starts[1:]:
            this_end = thisdict[chrom][this_start]
            if this_start <= last_end:  # overlap exists with last window
                last_start = min(this_start, last_start)
                last_end = max(this_end, last_end)
                breakfound = "no"
            else:  # there is no overlap. add concatentated range to dict
                newdict[chrom][last_start] = last_end
                last_start = this_start
                last_end = this_end
                breakfound = "yes"
                # add in the last one
                # if breakfound == "no":
        # peroxide j birthday
        newdict[chrom][last_start] = last_end
    return newdict


def qq_plot(real, rand, filename):
    if len(rand) == 0:
        if verbosefunction: print ("empty rand list for qq")
        return "empty"
    o = open(filename, "w")
    o.write("rand\t%s\trand_again\n" % filename)
    limit = min(len(rand), len(real))
    real.sort()
    rand.sort()
    expected = []
    observed = []
    for i in range(limit):
        indexa = int(len(rand) * float(i) / limit)
        indexb = int(len(real) * float(i) / limit)
        a = rand[indexa]
        b = real[indexb]
        o.write("%s\t%s\t%s\n" % (a, b, a))
        expected.append(a)
        observed.append(b)
    o.close()

    # ----------------------------
    fig = plt.figure(figsize=(5, 5), dpi=300)
    ax = plt.axes()  # [0,0,1,1], frameon=False) #fig.add_subplot(1,1,1)
    plt.box("off")
    # ax.set_title("Q:Q plot %s"%(supplementary_label))
    # ax.set_yticklabels(ticklabels, size=12, color="red") #use the list defined above
    ax.xaxis.set_ticks_position('bottom')  # remove ticks from top and right
    ax.yaxis.set_ticks_position('left')  # remove ticks from top and right
    ax.set_xlabel("Expected coexpression scores", size=12)
    ax.set_ylabel("Observed coexpression scores", size=12)
    # ax.set_xscale("log") #use logarithmic x axis
    # ax.set_yscale("log") #use logarithmic x axis
    # steps=2
    # ax.set_xticks(range(0, max(expected), steps))
    # ax.set_yticks(range(0,max(observed), steps))
    # ax.set_xticklabels([0,50,100], size=12)
    for i in range(len(expected)):
        ax.scatter([expected[i]], [expected[i]], color=["red"], alpha=0.7, marker="o", s=12)
        ax.scatter([expected[i]], [observed[i]], color=["blue"], alpha=0.7, marker="o", s=12)
        # axis limits
    if len(expected) > 0 and len(observed) > 0:
        plt.xlim(0, max(expected) + 1)
        plt.ylim(0, max(observed) + 1)
        fig.savefig("%s.svg" % (filename))
        # ----------------------------
