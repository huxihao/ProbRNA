## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

## Learning domain patten from data
## Author: Xihao Hu <huxihao@gmail.com>
## Version 1.0 - October 26, 2011

from tools import *
from ParsScores import *
from gene_mapping import map2genome, get3length

from multiprocessing import Pool
from random import seed, choice, sample, random
import math, string

MODEL = 'libsvm'

def ss_to_lr(ss):
    """ Transform secondary structure to marking stem loop regions """
    lr = list(ss)
    last = -1
    for i in range(len(lr)):
        if lr[i] == '(':
            last = i
        elif lr[i] == ')':
            if last != -1 and i-last-1 >= 5:
                for j in range(last+1, i):
                    lr[j] = 'o'
                last = -1
    return ''.join(lr)

def read_table(filename, select):
    import csv
    data_table = csv.reader(open(filename, "rb"), delimiter=',')
    data = {}
    head = data_table.next()
    gcol, scol = 0, 1
    for col in xrange(len(head)):
        if head[col] == "gene":
            gcol = col
        if head[col] == select:
            scol = col
    head[gcol] = "[" + head[gcol] + "]"
    head[scol] = "<" + head[scol] + ">"
    print "-".join(head),
    index = 0
    gene = ""
    values = []
    for row in data_table:
        if row[gcol] == "" or row[gcol] == gene:
            values.append(row[scol])
            continue
        ## update a new gene
        index += 1
        if gene != "":
            data[gene] = values
        gene = row[gcol]
        values = [row[scol]]
        ## Test
        #if int(row[0]) != index:
        #    warning("Wrong indexes! %s != %s\n"%(row[0], index))
        #    exit()
    data[gene] = values
    print "with", len(data), "genes, eg:", values[0]
    return data

def sigmoid(x):
    if x > 100: return 1
    if x < -100: return 0
    return 1 / (1 + math.exp(-x))

def format_data_near(outfile, zipcodes, encode, fitFile="win40.csv", window=40, near=100):
    """ Positive set from positions within the region
        Negtive set from positions near the region """
    span = window/2 ## define half the window size
    pars_data = import_pars()
    genome = map2genome()
    output = open(outfile, 'w')

    probv, probs = 0, 0
    if "ProbVS" in encode:
        probv = read_table(fitFile, "pbv")
        probs = read_table(fitFile, "pbs")

    ## headers
    output.write("label,gene,pos")
    if "SeqIndex" in encode:
        for j in xrange(-span, span):
            output.write(",char%s"%j)
    if "SeqBinary" in encode:
        for j in xrange(-span, span):
            output.write(",A%s,U%s,C%s,G%s"%(j,j,j,j))
    if "SeqGC" in encode:
        output.write(",GC,AU")
    if "SeqDiNu" in encode:
        for nu1 in ['A','U','C','G']:
            for nu2 in ['A','U','C','G']:
                output.write(",%s%s"%(nu1, nu2))
    if "SeqRatio" in encode:
        output.write(",A,U,C,G")
    if "PredSS3" in encode:
        for j in xrange(-span, span):
            output.write(",SpL%s,SpR%s,SpU%s"%(j,j,j))
    if "PredSS2" in encode:
        for j in xrange(-span, span):
            output.write(",pP%s,pU%s"%(j,j))
    if "PARS" in encode:
        for j in xrange(-span, span):
            output.write(",pars%s"%j)
    if "PARS2" in encode:
        for j in xrange(-span, span):
            output.write(",pars2%s"%j)
    if "LogVS" in encode:
        for j in xrange(-span, span):
            output.write(",lV%s,lS%s"%(j,j))
    if "ProbVS" in encode:
        for j in xrange(-span, span):
            output.write(",pV%s,pS%s"%(j,j))
    output.write("\n")

    data_size = 0
    for gene, zipcode, region in zipcodes:
        lens = get3length(genome[gene])
        pars_gene = pars_data[gene]
        seq = pars_gene["FOLD_SEQ"]
        ss = pars_gene["FOLD_SS"]
        lr = ss_to_lr(ss)
        prob_l = pars_gene["FOLD_PROB_L"]
        prob_r = pars_gene["FOLD_PROB_R"]
        score = pars_gene["PARS"]
        v1 = pars_gene["V1"]
        s1 = pars_gene["S1"]

        pv, ps = 0, 0
        if "ProbVS" in encode:
            pv = [float(val) for val in probv[gene]]
            ps = [float(val) for val in probs[gene]]

        ## the index of region is begin with 1 and close on both end
        region_begin, region_end = [int(val)+lens[0] for val in region.split('~')]
        print gene, zipcode, region_begin, region_end, len(seq)==sum(lens)
        for i in xrange(region_begin - near, region_end + near):
            if i < span or i >= len(seq) - span:
                continue
            ## region [RL, i, RR); span [WL, i, WR)
            RL = region_begin - 1; RR = region_end
            WL = i - span; WR = i + span
            #label = max(0, min(RR,WR)-max(RL,WL))/float(min(RR-RL,WR-WL))
            #label = 2*(label - 0.5)
            if RL <= i and i <= RR:
                label = 1
            else:
                label = -1
            ele_list = [label, "%s-%s"%(gene,zipcode), i+1]

            if "SeqIndex" in encode:
                for j in xrange(WL, WR):
                    ele_list.append("ACGU".find(seq[j]) + 1) ## return index
            if "SeqBinary" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(seq[j]=='A'), int(seq[j]=='U'),
                        int(seq[j]=='C'), int(seq[j]=='G')])
            if "SeqGC" in encode:
                ele_list.append((seq.count('G',WL,WR)+seq.count('C',WL,WR))/float(window))
                ele_list.append((seq.count('A',WL,WR)+seq.count('U',WL,WR))/float(window))
            if "SeqDiNu" in encode:
                for nu1 in ['A','U','C','G']:
                    for nu2 in ['A','U','C','G']:
                        ele_list.append(sum([int(seq[i]==nu1 and seq[i+1]==nu2)
                                        for i in xrange(WL,WR-1)])/float(window-1))
            if "SeqRatio" in encode:
                for nu in ['A','U','C','G']:
                    ele_list.append(seq.count(nu,WL,WR)/float(window))
            if "PredSS3" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([prob_l[j], prob_r[j], (1-prob_l[j]-prob_r[j])])
            if "PredSS2" in encode:
                for j in xrange(WL, WR):
                    #ele_list.extend([int(ss[j]!='.'), int(ss[j]=='.')])
                    ele_list.extend([prob_l[j]+prob_r[j], 1-prob_l[j]-prob_r[j]])
            if "PARS" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j]+7)/14.0) ## normalize
            if "PARS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j])**2/49.0) ## normalize
            if "LogVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([math.log(v1[j]+1,2), math.log(s1[j]+1,2)])
            if "ProbVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([pv[j], ps[j]])
            output.write(",".join([str(ele) for ele in ele_list])+"\n")
            data_size += 1
    output.close()
    return data_size

def format_data_perm(outfile, zipcodes, encode, fitFile="win40.csv", window=40, fold=1, rep=10):
    """ Positive set from positions within the region
        Negative set from a base permutation on the sample in the positive set """
    seed(2012)
    span = window/2 ## define half the window size
    pars_data = import_pars()
    genome = map2genome()
    output = open(outfile, 'w')

    probv, probs = 0, 0
    if "ProbVS" in encode:
        probv = read_table(fitFile, "pbv")
        probs = read_table(fitFile, "pbs")

    ## headers
    output.write("label,gene,pos")
    if "SeqIndex" in encode:
        for j in xrange(-span, span):
            output.write(",char%s"%j)
    if "SeqBinary" in encode:
        for j in xrange(-span, span):
            output.write(",A%s,U%s,C%s,G%s"%(j,j,j,j))
    if "SeqGC" in encode:
        output.write(",GC,AU")
    if "SeqDiNu" in encode:
        for nu1 in ['A','U','C','G']:
            for nu2 in ['A','U','C','G']:
                output.write(",%s%s"%(nu1, nu2))
    if "SeqRatio" in encode:
        output.write(",A,U,C,G")
    if "PredSS3" in encode:
        for j in xrange(-span, span):
            output.write(",SpL%s,SpR%s,SpU%s"%(j,j,j))
    if "PredSS2" in encode:
        for j in xrange(-span, span):
            output.write(",pP%s,pU%s"%(j,j))
    if "PARS" in encode:
        for j in xrange(-span, span):
            output.write(",pars%s"%j)
    if "PARS2" in encode:
        for j in xrange(-span, span):
            output.write(",pars2%s"%j)
    if "LogVS" in encode:
        for j in xrange(-span, span):
            output.write(",lV%s,lS%s"%(j,j))
    if "ProbVS" in encode:
        for j in xrange(-span, span):
            output.write(",pV%s,pS%s"%(j,j))
    output.write("\n")

    data_size = 0
    for gene, zipcode, region in zipcodes:
        lens = get3length(genome[gene])
        pars_gene = pars_data[gene]
        seq = pars_gene["FOLD_SEQ"]
        ss = pars_gene["FOLD_SS"]
        lr = ss_to_lr(ss)
        prob_l = pars_gene["FOLD_PROB_L"]
        prob_r = pars_gene["FOLD_PROB_R"]
        score = pars_gene["PARS"]
        v1 = pars_gene["V1"]
        s1 = pars_gene["S1"]

        split_name = gene
        if fold > 1:
            split_name = "fold_%s"%(choice(range(fold))+1)
        else:
            split_name = "%s-%s"%(gene,zipcode)

        pv, ps = 0, 0
        if "ProbVS" in encode:
            pv = [float(val) for val in probv[gene]]
            ps = [float(val) for val in probs[gene]]

        ## the index of region is begin with 1 and close on both end
        region_begin, region_end = [int(val)+lens[0] for val in region.split('~')]
        print gene, zipcode, region_begin, region_end, len(seq)==sum(lens)
        for i in xrange(region_begin-1, region_end):
            if i < span or i >= len(seq) - span:
                continue
            ## region [RL, i, RR); span [WL, i, WR)
            RL = region_begin - 1; RR = region_end
            WL = i - span; WR = i + span
            if RL <= i and i <= RR:
                label = 1
            else:
                label = -1
            ele_list = [label, split_name, i+1]

            ## permuate `rep' times to generate negative set
            neg_list = [[-1, split_name, i+1] for k in xrange(rep)]
            neg_idx = [sample(range(WL,WR),WR-WL) for k in xrange(rep)]

            if "SeqIndex" in encode:
                for j in xrange(WL, WR):
                    ele_list.append("ACGU".find(seq[j]) + 1) ## return index
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].append("ACGU".find(seq[j]) + 1) ## return index
            if "SeqBinary" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(seq[j]=='A'), int(seq[j]=='U'), int(seq[j]=='C'), int(seq[j]=='G')])
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].extend([int(seq[j]=='A'), int(seq[j]=='U'), int(seq[j]=='C'), int(seq[j]=='G')])
            if "SeqGC" in encode:
                ele_list.append((seq.count('G',WL,WR)+seq.count('C',WL,WR))/float(window))
                ele_list.append((seq.count('A',WL,WR)+seq.count('U',WL,WR))/float(window))
            if "SeqDiNu" in encode:
                for nu1 in ['A','U','C','G']:
                    for nu2 in ['A','U','C','G']:
                        ele_list.append(sum([int(seq[i]==nu1 and seq[i+1]==nu2)
                                        for i in xrange(WL,WR-1)])/float(window-1))
            if "SeqRatio" in encode:
                for nu in ['A','U','C','G']:
                    ele_list.append(seq.count(nu,WL,WR)/float(window))
                for k in xrange(rep):
                    neg_list[k].extend(ele_list[-4:])
            if "PredSS3" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([prob_l[j], prob_r[j], (1-prob_l[j]-prob_r[j])])
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].extend([prob_l[j], prob_r[j], (1-prob_l[j]-prob_r[j])])
            if "PredSS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(ss[j]!='.'), int(ss[j]=='.')])
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].extend([int(ss[j]!='.'), int(ss[j]=='.')])
            if "PARS" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j]+7)/14.0) ## normalize
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].append((score[j]+7)/14.0) ## normalize
            if "PARS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j])**2/49.0) ## normalize
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].append((score[j])**2/49.0) ## normalize
            if "LogVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([math.log(v1[j]+1,2), math.log(s1[j]+1,2)])
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].extend([math.log(v1[j]+1,2), math.log(s1[j]+1,2)])
            if "ProbVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([pv[j], ps[j]])
                for k in xrange(rep):
                    for j in neg_idx[k]:
                        neg_list[k].extend([pv[j], ps[j]])
            output.write(",".join([str(ele) for ele in ele_list])+"\n")
            for k in xrange(rep):
                output.write(",".join([str(ele) for ele in neg_list[k]])+"\n")
            data_size += 2
    output.close()
    return data_size

def format_data_gene(outfile, zipcodes, encode, fitFile="win40.csv", window=40, isGroup=True, place="All", fold=1):
    """ Position set from positions within all the regions in a gene
        Negative set from positions remaining in that gene """
    seed(2012)
    span = window/2 ## define half of the window size
    pars_data = import_pars()
    genome = map2genome()
    output = open(outfile, 'w')
    data_size = 0

    probv, probs = 0, 0
    if "ProbVS" in encode:
        probv = read_table(fitFile, "pbv")
        probs = read_table(fitFile, "pbs")

    ## headers
    output.write("label,gene,pos")
    if "SeqIndex" in encode:
        for j in xrange(-span, span):
            output.write(",char%s"%j)
    if "SeqBinary" in encode:
        for j in xrange(-span, span):
            output.write(",A%s,U%s,C%s,G%s"%(j,j,j,j))
    if "SeqGC" in encode:
        output.write(",GC,AU")
    if "SeqDiNu" in encode:
        for nu1 in ['A','U','C','G']:
            for nu2 in ['A','U','C','G']:
                output.write(",%s%s"%(nu1, nu2))
    if "SeqRatio" in encode:
        output.write(",A,U,C,G")
    if "PredSS3" in encode:
        for j in xrange(-span, span):
            output.write(",SpL%s,SpR%s,SpU%s"%(j,j,j))
    if "PredSS2" in encode:
        for j in xrange(-span, span):
            output.write(",pP%s,pU%s"%(j,j))
    if "PARS" in encode:
        for j in xrange(-span, span):
            output.write(",pars%s"%j)
    if "PARS2" in encode:
        for j in xrange(-span, span):
            output.write(",pars2%s"%j)
    if "LogVS" in encode:
        for j in xrange(-span, span):
            output.write(",lV%s,lS%s"%(j,j))
    if "ProbVS" in encode:
        for j in xrange(-span, span):
            output.write(",pV%s,pS%s"%(j,j))
    output.write("\n")

    for gene, zipcode, region in zipcodes:
        lens = get3length(genome[gene])
        pars_gene = pars_data[gene]
        seq = pars_gene["FOLD_SEQ"]
        ss = pars_gene["FOLD_SS"]
        lr = ss_to_lr(ss)
        prob_l = pars_gene["FOLD_PROB_L"]
        prob_r = pars_gene["FOLD_PROB_R"]
        score = pars_gene["PARS"]
        v1 = pars_gene["V1"]
        s1 = pars_gene["S1"]

        pv, ps = 0, 0
        if "ProbVS" in encode:
            pv = [float(val) for val in probv[gene]]
            ps = [float(val) for val in probs[gene]]

        ## the index of region is begin with 1 and close on both end
        region_begin, region_end = [int(val)+lens[0] for val in region.split('~')]
        print gene, zipcode, region_begin, region_end, len(seq)==sum(lens)

        ## generate postive and negative sample index list
        region_list = range(max(span, region_begin-1), min(region_end, len(seq)-span))

        if len(region_list) == 0:
            print zipcode, "is removed, please use smaller window size."
            continue

        neg_list = range(span, region_begin-1) + range(region_end, len(seq)-span)

        if isGroup: ## group by gene
            counter = 0
            for gene1, zipcode1, region1 in zipcodes:
                if gene1 != gene: ## not the same RNA, so go on
                    continue
                counter += 1 ## count the number of zipcodes on this RNA
                if zipcode1 == zipcode and counter != 1: ## not the first one
                    neg_list = [] ## clear negtive list
                    break
                r_b1, r_e1 = [int(val)+lens[0] for val in region1.split('~')]
                for i in range(r_b1-1, r_e1):
                    if i in neg_list:
                        neg_list.remove(i)

        #region_list.extend(random.sample(neg_list, len(region_list)))
        region_list.extend(neg_list)

        split_name = gene
        if fold > 1:
            split_name = "fold_%s"%(choice(range(fold))+1)
        elif not isGroup:
            split_name = "%s-%s"%(gene,zipcode)

        for i in region_list:
            if place == "All":
                pass
            elif place == "5UTR":
                if i >= lens[0]:
                    continue
            elif place == "CDS":
                if i < lens[0] or i >= lens[0]+lens[1]:
                    continue
            elif place == "3UTR":
                if i < lens[0] + lens[1]:
                    continue

            ## region [RL, i, RR); span [WL, i, WR)
            RL = region_begin - 1; RR = region_end
            WL = i - span; WR = i + span

            if RL <= i and i <= RR:
                label = 1
            else:
                label = -1

            ## begin package
            ele_list = [label, split_name, i+1]

            if "SeqIndex" in encode:
                for j in xrange(WL, WR):
                    ele_list.append("ACGU".find(seq[j]) + 1) ## return index
            if "SeqBinary" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(seq[j]=='A'), int(seq[j]=='U'),
                        int(seq[j]=='C'), int(seq[j]=='G')])
            if "SeqGC" in encode:
                ele_list.append((seq.count('G',WL,WR)+seq.count('C',WL,WR))/float(window))
                ele_list.append((seq.count('A',WL,WR)+seq.count('U',WL,WR))/float(window))
            if "SeqDiNu" in encode:
                for nu1 in ['A','U','C','G']:
                    for nu2 in ['A','U','C','G']:
                        ele_list.append(sum([int(seq[i]==nu1 and seq[i+1]==nu2)
                                        for i in xrange(WL,WR-1)])/float(window-1))
            if "SeqRatio" in encode:
                for nu in ['A','U','C','G']:
                    ele_list.append(seq.count(nu,WL,WR)/float(window))
            if "PredSS3" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([prob_l[j], prob_r[j], (1-prob_l[j]-prob_r[j])])
            if "PredSS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(ss[j]!='.'), int(ss[j]=='.')])
                    #ele_list.extend([prob_l[j]+prob_r[j], 1-prob_l[j]-prob_r[j]])
            if "PARS" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j]+7)/14.0) ## normalize
            if "PARS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j])**2/49.0) ## normalize
            if "LogVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([math.log(v1[j]+1,2), math.log(s1[j]+1,2)])
            if "ProbVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([pv[j], ps[j]])
            output.write(",".join([str(ele) for ele in ele_list])+"\n")
            data_size += 1
    output.close()
    return data_size


def format_data_gene_ratio(outfile, zipcodes, encode, fitFile="win40.csv", window=40, isGroup=True, place="All", fold=1, cutoff=0.05):
    """ Position set from positions within all the regions in a gene
        Negative set from positions remaining in that gene 
        Choose negative set with the same nucleotide ratios as positive set """
    seed(2012)
    span = window/2 ## define half of the window size
    pars_data = import_pars()
    genome = map2genome()
    output = open(outfile, 'w')
    data_size = 0

    probv, probs = 0, 0
    if "ProbVS" in encode:
        probv = read_table(fitFile, "pbv")
        probs = read_table(fitFile, "pbs")

    ## headers
    output.write("label,gene,pos")
    if "SeqIndex" in encode:
        for j in xrange(-span, span):
            output.write(",char%s"%j)
    if "SeqBinary" in encode:
        for j in xrange(-span, span):
            output.write(",A%s,U%s,C%s,G%s"%(j,j,j,j))
    if "SeqGC" in encode:
        output.write(",GC,AU")
    if "SeqDiNu" in encode:
        for nu1 in ['A','U','C','G']:
            for nu2 in ['A','U','C','G']:
                output.write(",%s%s"%(nu1, nu2))
    if "SeqRatio" in encode:
        output.write(",A,U,C,G")
    if "PredSS3" in encode:
        for j in xrange(-span, span):
            output.write(",SpL%s,SpR%s,SpU%s"%(j,j,j))
    if "PredSS2" in encode:
        for j in xrange(-span, span):
            output.write(",pP%s,pU%s"%(j,j))
    if "PARS" in encode:
        for j in xrange(-span, span):
            output.write(",pars%s"%j)
    if "PARS2" in encode:
        for j in xrange(-span, span):
            output.write(",pars2%s"%j)
    if "LogVS" in encode:
        for j in xrange(-span, span):
            output.write(",lV%s,lS%s"%(j,j))
    if "ProbVS" in encode:
        for j in xrange(-span, span):
            output.write(",pV%s,pS%s"%(j,j))
    output.write("\n")

    pos_cc = neg_cc = 0
    for gene, zipcode, region in zipcodes:
        lens = get3length(genome[gene])
        pars_gene = pars_data[gene]
        seq = pars_gene["FOLD_SEQ"]
        ss = pars_gene["FOLD_SS"]
        lr = ss_to_lr(ss)
        prob_l = pars_gene["FOLD_PROB_L"]
        prob_r = pars_gene["FOLD_PROB_R"]
        score = pars_gene["PARS"]
        v1 = pars_gene["V1"]
        s1 = pars_gene["S1"]

        pv, ps = 0, 0
        if "ProbVS" in encode:
            pv = [float(val) for val in probv[gene]]
            ps = [float(val) for val in probs[gene]]

        ## the index of region is begin with 1 and close on both end
        region_begin, region_end = [int(val)+lens[0] for val in region.split('~')]
        print gene, zipcode, region_begin, region_end, len(seq)==sum(lens),

        ## generate postive and negative sample index list
        region_list = range(max(span, region_begin-1), min(region_end, len(seq)-span))

        if len(region_list) == 0:
            print zipcode, "is removed, please use smaller window size."
            continue

        neg_list = range(span, region_begin-1) + range(region_end, len(seq)-span)

        if isGroup: ## group by gene
            counter = 0
            for gene1, zipcode1, region1 in zipcodes:
                if gene1 != gene: ## not the same RNA, so go on
                    continue
                counter += 1 ## count the number of zipcodes on this RNA
                if zipcode1 == zipcode and counter != 1: ## not the first one
                    neg_list = [] ## clear negtive list
                    break
                r_b1, r_e1 = [int(val)+lens[0] for val in region1.split('~')]
                for i in range(r_b1-1, r_e1):
                    if i in neg_list:
                        neg_list.remove(i)

        ## get nucleotides ratios in the postive set
        positive_seq = seq[(region_begin-1) : region_end]
        pos_seq_len = float(len(positive_seq))
        ratios = {'A':positive_seq.count('A')/pos_seq_len,
                  'U':positive_seq.count('U')/pos_seq_len,
                  'C':positive_seq.count('C')/pos_seq_len,
                  'G':positive_seq.count('G')/pos_seq_len}

        #region_list.extend(sample(neg_list, len(region_list)))
        new_pos_list = []
        for i in region_list:
        #    WL = i - span; WR = i + span
        #    dis = 0
        #    for char in list("AUCG"):
        #        dis += (ratios[char] - seq[WL:WR].count(char)/float(window))**2
        #    if dis**0.5 <= cutoff:
            ## include all postive regions
            new_pos_list.append(i)
        new_neg_list = []
        for i in neg_list:
            if place == "All":
                pass
            elif place == "5UTR":
                if i >= lens[0]:
                    continue
            elif place == "CDS":
                if i < lens[0] or i >= lens[0]+lens[1]:
                    continue
            elif place == "3UTR":
                if i < lens[0] + lens[1]:
                    continue
            #negative_seq = seq[(i-span) : (i+span)]/float(window)
            ## composition on negative measured by similar length
            negative_seq = seq[max(0, i-len(positive_seq)/2) : min(i+len(positive_seq)/2, len(seq))]
            dis = 0
            for char in list("AUCG"):
                dis += (ratios[char] - negative_seq.count(char)/float(len(negative_seq)))**2
            if dis**0.5 <= cutoff:
                new_neg_list.append(i)
        print len(new_pos_list), len(new_neg_list)
        if new_pos_list == [] or new_neg_list == []:
            continue ## no positive or negative sample

        pos_cc += len(new_pos_list)
        neg_cc += len(new_neg_list)

        region_list = new_pos_list + new_neg_list

        split_name = gene
        if fold > 1:
            split_name = "fold_%s"%(choice(range(fold))+1)
        elif not isGroup:
            split_name = "%s-%s"%(gene,zipcode)

        for i in region_list:
            ## region [RL, i, RR); span [WL, i, WR)
            RL = region_begin - 1; RR = region_end
            WL = i - span; WR = i + span

            if RL <= i and i <= RR:
                label = 1
            else:
                label = -1

            ## begin package
            ele_list = [label, split_name, i+1]

            if "SeqIndex" in encode:
                for j in xrange(WL, WR):
                    ele_list.append("ACGU".find(seq[j]) + 1) ## return index
            if "SeqBinary" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(seq[j]=='A'), int(seq[j]=='U'),
                        int(seq[j]=='C'), int(seq[j]=='G')])
            if "SeqGC" in encode:
                ele_list.append((seq.count('G',WL,WR)+seq.count('C',WL,WR))/float(window))
                ele_list.append((seq.count('A',WL,WR)+seq.count('U',WL,WR))/float(window))
            if "SeqDiNu" in encode:
                for nu1 in ['A','U','C','G']:
                    for nu2 in ['A','U','C','G']:
                        ele_list.append(sum([int(seq[i]==nu1 and seq[i+1]==nu2)
                                        for i in xrange(WL,WR-1)])/float(window-1))
            if "SeqRatio" in encode:
                for nu in ['A','U','C','G']:
                    ele_list.append(seq.count(nu,WL,WR)/float(window))
            if "PredSS3" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([prob_l[j], prob_r[j], (1-prob_l[j]-prob_r[j])])
            if "PredSS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(ss[j]!='.'), int(ss[j]=='.')])
                    #ele_list.extend([prob_l[j]+prob_r[j], 1-prob_l[j]-prob_r[j]])
            if "PARS" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j]+7)/14.0) ## normalize
            if "PARS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j])**2/49.0) ## normalize
            if "LogVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([math.log(v1[j]+1,2), math.log(s1[j]+1,2)])
            if "ProbVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([pv[j], ps[j]])
            output.write(",".join([str(ele) for ele in ele_list])+"\n")
            data_size += 1
    output.close()
    show("#Pos %d, #Neg %d\n"%(pos_cc, neg_cc))
    return data_size

def format_data_all(outfile, candidate, encode, window, pars_data, genome, probv, probs, place="All"):
    """ Don't give labels, but just encode all positions in the gene for prediction.  """
    span = window/2 ## define half of window size

    output = open(outfile, 'w')
    ## headers
    output.write("label,gene,pos")
    if "SeqIndex" in encode:
        for j in xrange(-span, span):
            output.write(",char%s"%j)
    if "SeqBinary" in encode:
        for j in xrange(-span, span):
            output.write(",A%s,U%s,C%s,G%s"%(j,j,j,j))
    if "SeqGC" in encode:
        output.write(",GC,AU")
    if "SeqDiNu" in encode:
        for nu1 in ['A','U','C','G']:
            for nu2 in ['A','U','C','G']:
                output.write(",%s%s"%(nu1, nu2))
    if "SeqRatio" in encode:
        output.write(",A,U,C,G")
    if "PredSS3" in encode:
        for j in xrange(-span, span):
            output.write(",SpL%s,SpR%s,SpU%s"%(j,j,j))
    if "PredSS2" in encode:
        for j in xrange(-span, span):
            output.write(",pP%s,pU%s"%(j,j))
    if "PARS" in encode:
        for j in xrange(-span, span):
            output.write(",pars%s"%j)
    if "PARS2" in encode:
        for j in xrange(-span, span):
            output.write(",pars2%s"%j)
    if "LogVS" in encode:
        for j in xrange(-span, span):
            output.write(",lV%s,lS%s"%(j,j))
    if "ProbVS" in encode:
        for j in xrange(-span, span):
            output.write(",pV%s,pS%s"%(j,j))
    output.write("\n")

    for gene in pars_data:
        if candidate != [] and gene not in candidate:
            continue
        lens = get3length(genome[gene])
        pars_gene = pars_data[gene]
        seq = pars_gene["FOLD_SEQ"]
        ss = pars_gene["FOLD_SS"]
        lr = ss_to_lr(ss)
        prob_l = pars_gene["FOLD_PROB_L"]
        prob_r = pars_gene["FOLD_PROB_R"]
        score = pars_gene["PARS"]
        v1 = pars_gene["V1"]
        s1 = pars_gene["S1"]

        pv, ps = 0, 0
        if "ProbVS" in encode:
            pv = [float(val) for val in probv[gene]]
            ps = [float(val) for val in probs[gene]]

        ## generate postive and negative sample index list
        for i in xrange(span, (len(seq)-span)):
            if place == "All":
                pass
            elif place == "5UTR":
                if i >= lens[0]:
                    continue
            elif place == "CDS":
                if i < lens[0] or i >= lens[0]+lens[1]:
                    continue
            elif place == "3UTR":
                if i < lens[0] + lens[1]:
                    continue

            WL = i - span; WR = i + span
            ele_list = [0, gene, i+1] ## label, name, pos
            if "SeqIndex" in encode:
                for j in xrange(WL, WR):
                    ele_list.append("ACGU".find(seq[j]) + 1) ## return index
            if "SeqBinary" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(seq[j]=='A'), int(seq[j]=='U'),
                        int(seq[j]=='C'), int(seq[j]=='G')])
            if "SeqGC" in encode:
                ele_list.append((seq.count('G',WL,WR)+seq.count('C',WL,WR))/float(window))
                ele_list.append((seq.count('A',WL,WR)+seq.count('U',WL,WR))/float(window))
            if "SeqDiNu" in encode:
                for nu1 in ['A','U','C','G']:
                    for nu2 in ['A','U','C','G']:
                        ele_list.append(sum([int(seq[i]==nu1 and seq[i+1]==nu2)
                                        for i in xrange(WL,WR-1)])/float(window-1))
            if "SeqRatio" in encode:
                for nu in ['A','U','C','G']:
                    ele_list.append(seq.count(nu,WL,WR)/float(window))
            if "PredSS3" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([prob_l[j], prob_r[j], (1-prob_l[j]-prob_r[j])])
            if "PredSS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([int(ss[j]!='.'), int(ss[j]=='.')])
                    #ele_list.extend([prob_l[j]+prob_r[j], 1-prob_l[j]-prob_r[j]])
            if "PARS" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j]+7)/14.0) ## normalize
            if "PARS2" in encode:
                for j in xrange(WL, WR):
                    ele_list.append((score[j])**2/49.0) ## normalize
            if "LogVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([math.log(v1[j]+1,2), math.log(s1[j]+1,2)])
            if "ProbVS" in encode:
                for j in xrange(WL, WR):
                    ele_list.extend([pv[j], ps[j]])
            output.write(",".join([str(ele) for ele in ele_list])+"\n")
    output.close()


def motif_build(align_file, start, end):
    from Bio import Motif
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    motif = Motif.Motif(alphabet=IUPAC.unambiguous_rna)
    with open(align_file) as align:
        for line in align:
            sequence = line.strip().upper().replace('T','U')
            if len(sequence) == 0 or sequence[0] == '#':
                continue
            motif.add_instance(Seq(sequence[start:end], motif.alphabet))
    show("Consensus Motif: %s\n"%motif.consensus())
    for row in motif.log_odds():
        for char in "AUCG":
            show(row[char])
        show("\n")
    with open("motif_fasta.txt", 'w') as tmpfile:
        tmpfile.write(motif._to_fasta())
    #m.weblogo("motif_logo.gif", "GIF")
    return motif
    
def motif_score(seq, motif):
    from Bio import Motif
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    max_score = 0
    ps = motif.search_pwm(Seq(seq, motif.alphabet), normalized=True, threshold=max_score, both=False)
    for pos, score in ps:
        if score > max_score:
            max_score = score
    return max_score

def learn_domain_process((zipcode, encode, scanWin, model, fitWin, fold, format)):
    fitted = "PredAll_%s_%s_PARS_V1_S1.csv"%(model, fitWin)
    if 'ProbVS' in encode:
        features = 'Domain%s_%s%s_'%(len(zipcode),format,fold)+'_'.join(encode)+'_ScanWin%s_%s_FitWin%s.csv'%(scanWin, model, fitWin)
    else:
        features = 'Domain%s_%s%s_'%(len(zipcode),format,fold)+'_'.join(encode) + "_ScanWin%s.csv"%(scanWin)
    results = {"ScanWin":scanWin, "Encode":'-'.join(encode), 'FitModel':model, 'FitWin':fitWin}
    if not os.path.exists(features + ".log") or sum([l.startswith('Total') for l in open(features + ".log")])==0:
        if format == 'whole':
            data_size = format_data_gene(features, zipcode, encode, fitted, scanWin, True, "All", fold)
        elif format == 'permute':
            data_size = format_data_perm(features, zipcode, encode, fitted, scanWin, fold, rep=5)
        elif format == 'ratio':
            data_size = format_data_gene_ratio(features, zipcode, encode, fitted, scanWin, True, 'All', fold)
        else:
            raise ValueError('Unknown format to generate features: %s'%format)
        os.system("Rscript ../src/learn_domain_"+MODEL+".R CV %s "%features)
        #os.remove(features)
    else:
        print 'Read existing log file:', features
    outfile = open(features + ".log", 'r')
    for line in outfile:
        ele = line.strip().split('\t')
        if len(ele) < 2:
            break
        item, value = ele[:2]
        results[item] = value
    outfile.close()
    print "Process End:", features
    return results

def main(action, parameter):
    if action == "help":
        print """
        Welcome to use learn_domain.py script

The command format are:
  WorkPath=[a_valid_path]
  Action=[figure3 |figure4 |figure6| predict] | [_kim] | test
  Paras=[ThrdNum_w_MaxK |ThrdNum_K_MaxW |Place_ThrdNum_K_MaxW |->]
         -> w_SeqBinary-ProbVS-LogVS-PARS-PARS2-PredSS2-PredSS3

Examples:
  python src/learn_domain.py WorkPath=work Action=figure3 Paras=9_40_50
  python src/learn_domain.py WorkPath=work Action=figure3 Paras=9_100_50
  python src/learn_domain.py WorkPath=work Action=figure4 Paras=9_2_300
  python src/learn_domain.py WorkPath=work Action=figure6 Paras=9_2_300
  python src/learn_domain.py WorkPath=work Action=predict Paras=All_100_PARS-ProbVS
  python src/learn_domain.py WorkPath=work Action=predict_kim Paras=All_40_SeqBinary
        """
        return

    import_pars() ## initilization
    zipcodeList = [ #List of zipcodes
        ("YKL185W","E1min","635~683"), #Jambhekar
        ("YKL185W","E2Bmin","1279~1314"), #Jambhekar
        ("YKL185W","Umin","1766~1819"), #Jambhekar
        #("YKL185W","Other","1684~1719"), #Jambhekar
        ("YNL283C","WSC2N","418~471"), #Jambhekar
        ("YNL283C","WSC2C","1354~1384"), #Jambhekar
        ("YMR202W","ERG2N","180~250"), #Jambhekar
        #("YLL001W","DNM1N","605~805"), #Jambhekar
        #("YLL001W","DNM1C","1656~1752"), #Jambhekar
        ("YOR247W","SRL1C","419~596"), #Jambhekar
        ("YLL028W","TPO1N","2~178"), #Jambhekar
        #("YJL172W","CPS1CR","1305~1456"), #Jambhekar
        ("YKL185W","E2A","1109~1185"), #Olivier
        #("YBR086C","IST2","2716~2777"), #Olivier
        ("YMR171C","EAR1","1572~1621")] #Olivier

    puf3sites = [ # doi: 10.1261/rna.7270204 Figure 1
        ("YLL009C","SiteC","215~232"),
        ("YLL009C","SiteA","233~252"),
        ("YLL009C","SiteB","284~300")]

    if action == "test":
        learn_domain_process((zipcodeList, ["SeqBinary"], 10, "test.csv", 2, 1, 'whole'))
        return
        
    fold = 1 ## 1 is leave-one-gene-out validation
    format = 'whole' ## way to select negative sets
    if action.endswith("kim"):
        with open("candidate_list.txt", 'r') as region_file:
            zipcodeList = []
            for line in region_file:
                gene, name, region = line.strip().split("\t")[0:3]
                zipcodeList.append((gene, name, region))
        fold = 5 ## 5-fold cross validation
        if 'ratio' in action:
            format = 'ratio' ## avoid the sequence bias
        elif 'permute' in action:
            format = 'permute'

    if action.startswith("figure"):
        #####################################################################
        ## Cross Validation
        NumThrd, WinSize, MaxWin = parameter.split('_')
        pool = Pool(int(NumThrd)) ## Number of workers
        pool_para = []
        if action.startswith("figure3"):
            encode = ["ProbVS"]
            scanWin = int(WinSize) ## default is 40
            for model in ["MixPoiSep", "MixPoiLin", "MPLComSam", "MPLComOpp"]:
                for fitWin in xrange(2, int(MaxWin)+1, 2): 
                    pool_para.append((zipcodeList, encode, scanWin, model, fitWin, fold, format))
        elif action.startswith("figure4"): ## individual features
            model = "MixPoiLin"
            fitWin = int(WinSize) ## default is 2
            encode_list = [["SeqBinary"],["SeqRatio"],["SeqGC"],["SeqDiNu"],["ProbVS"],["LogVS"],["PARS"],["PARS2"],["PredSS2"],["PredSS3"]]
            if format == 'permute':
                encode_list = [["SeqBinary"],["ProbVS"],["LogVS"],["PARS"],["PARS2"],["SeqRatio"]]
            for encode in encode_list:
                for scanWin in range(10, int(MaxWin)+1, 10): 
                    pool_para.append((zipcodeList, encode, scanWin, model, fitWin, fold, format))
        elif action.startswith("figure5"): ## selected features
            model = "MixPoiLin"
            fitWin = int(WinSize) ## default is 2
            encode_list = [["SeqBinary"],["ProbVS"],["PARS"],["SeqRatio"],["SeqBinary","ProbVS"],["SeqBinary","PARS"]]
            for encode in encode_list:
                for scanWin in range(10, int(MaxWin)+1, 10): 
                    pool_para.append((zipcodeList, encode, scanWin, model, fitWin, fold, format))
        elif action.startswith("figure6"): ## features with sequences
            model = "MixPoiLin"
            fitWin = int(WinSize) ## default is 2
            encode_list = [[],["ProbVS"],["LogVS"],["PARS"],["PARS2"],["PredSS2"],["PredSS3"]]
            if format == 'permute':
                encode_list = [[],["SeqBinary"],["ProbVS"],["LogVS"],["PARS"],["PARS2"],["SeqRatio"]]
            for encode in encode_list:
                encode = ["SeqBinary"] + encode
                for scanWin in range(10, int(MaxWin)+1, 10): 
                    pool_para.append((zipcodeList, encode, scanWin, model, fitWin, fold, format))
        results = pool.map(learn_domain_process, pool_para)
        items = results[0].keys()
        items.sort()
        for item in items:
            show(item)
        show('\n')
        for result in results:
            for item in items:
                show(result[item])
            show('\n')
        return

    if action.startswith("predict"):
        ###########################################################################
        ## Screening using trained model
        place = "All" ## default
        para = "1" ## default
        Place, WinSize, EncodeString = parameter.split('_') ## read in
        place = Place
        scanWin = int(WinSize)
        encode = EncodeString.split('-')
        fitted_pars = "PredAll_MixPoiLin_2_PARS_V1_S1.csv"
        if action.endswith("kim"):
            format_data_perm("domain_encode_temp.csv", zipcodeList, encode, fitted_pars, scanWin, 5, rep=5)
        else:
            format_data_gene("domain_encode_temp.csv", zipcodeList, encode, fitted_pars, scanWin, True, place, 5)

        os.system("Rscript ../src/learn_domain_"+MODEL+".R Train domain_encode_temp.csv " + para)

        ## read huge files
        pars_data = import_pars()
        genome = map2genome()
        full_list, groups = pars_gene_groups(pars_data)
        ## predict subset if given
        if os.path.exists("gene_list.txt"):
            full_list = []
            with open("gene_list.txt",'r') as infile:
                for line in infile:
                    full_list.append(line.strip())

        probv, probs = 0, 0
        if "ProbVS" in encode:
            probv = read_table(fitted_pars, "pbv")
            probs = read_table(fitted_pars, "pbs")

        predicted_label = open("predicted_domain_labels.csv", 'w')
        predicted_label.write("gene,pos,pred\n")
        start = end = 0
        while end < len(full_list): ## not cover the full list
            start = end
            end += 30
            if end > len(full_list):
                end = len(full_list)

            gene_list = full_list[start:end]
            format_data_all("domain_encode_temp.csv", gene_list, encode, scanWin, pars_data, genome, probv, probs, place)
            os.system("Rscript ../src/learn_domain_"+MODEL+".R Predict domain_encode_temp.csv")
            gene_label = read_table("domain_encode_temp.csv", "label")
            gene_pos = read_table("domain_encode_temp.csv", "pos")
            for gene in gene_label:
                label = gene_label[gene]
                pos = gene_pos[gene]
                for i in xrange(len(label)):
                    if i == 0:
                        predicted_label.write("%s,%s,%s\n"%(gene, pos[i], label[i]))
                    else:
                        predicted_label.write(",%s,%s\n"%(pos[i], label[i]))
        predicted_label.close()
        print len(read_table("predicted_domain_labels.csv", "pred"))
        return
    ## End of main

if __name__ == "__main__":
    ## Default settings
    WORK_PATH = "."
    ACTION = "help"
    PARAMETER = "default"

    ## Read parameters
    for arg in sys.argv:
        ## example: python AA.py WorkPath=C:/aa/
        if arg.startswith("WorkPath"):
            WORK_PATH = arg[(string.find(arg,"=")+1):]
        if arg.startswith("Action"):
            ACTION = arg[(string.find(arg,"=")+1):]
        if arg.startswith("Paras"):
            PARAMETER = arg[(string.find(arg,"=")+1):]
    if ACTION == "help":
        main(ACTION, PARAMETER)
        sys.exit(0)

    ## Set work path
    os.chdir(WORK_PATH)
    ## Delete old log file
    if os.path.exists(LOG_FILE):
        os.remove(LOG_FILE)
    show("# Current work path is %s\n"%os.path.abspath("."))
    show("# Action is %s with parameter %s\n"%(ACTION, PARAMETER))
    ## Set timer
    START_TIME = time.time()
    main(ACTION, PARAMETER)
    END_TIME = time.time()
    show(time.strftime("# Start  at: %c\n", time.gmtime(START_TIME)))
    show(time.strftime("# Finish at: %c\n", time.gmtime(END_TIME)))
    show(time.strftime("# Total Time: %j day %H:%M:%S\n",\
            time.gmtime(END_TIME-START_TIME))) 

