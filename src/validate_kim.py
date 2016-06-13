## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

## Validate the PAR-CLIP results from Kim lab
## Author: Xihao Hu <huxihao@gmail.com>
## Version 1.0 - March 7, 2012

from tools import *
from ParsScores import import_pars
from gene_mapping import map2genome, get3length, get_splice, global2local
import math

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

def get_region_list(Type, kim_data='../data/KimLab/120214FREEZE.txt'):
    ## read pars data set
    pars = import_pars()

    ## read genome map
    genome = map2genome()

    ## read Kim's data
    region_list = []
    data_file = open(kim_data, 'r')

    ## read in all data from Kim's Lab
    su3utr, su4mid, su5utr = 0, 0, 0
    cc3utr, cc4mid, cc5utr = 0, 0, 0
    cc_in, cc_use = 0, 0
    cc_hit, cc_miss = 0, 0
    gene_list = set()
    print data_file.readline() ## header
    len_dist = [] ## binding region length distribution
    rpm_dist = [] ## RPM distribution
    for line in data_file:
        eid, chro, strand, start, stop, gene, pos, rpm = line.strip().split('\t')
        cc_in += 1
        if gene not in genome:
            continue
        if gene not in pars:
            continue
        if pars[gene]["TYPE"] != "mRNA":
            continue
        pos = pos.replace("'",'') ## remove the prim
        if Type != "All" and Type != "Splice" and pos != Type: ## select by region position
            continue
        genome_mark = genome[gene]
        start = global2local(genome_mark, int(start))
        stop = global2local(genome_mark, int(stop))
        lengths = get3length(genome_mark)
        splice = get_splice(genome_mark) 

        if stop < start:
            start, stop = stop, start
        len_dist.append((stop-start))
        rpm_dist.append(float(rpm))

        ## filtering
        if stop <0 or start <0 or stop-start<10 or stop-start>40:
            continue

        if Type == "Splice":
            if splice == []:
                continue
            isNear = False
            for site in splice:
                if abs(start+stop-2*site) < 100:
                    cc_hit += 1
                    isNear = True
                else:
                    cc_miss += 1
            if not isNear:
                continue

        su3utr += lengths[0] 
        su4mid += lengths[1] 
        su5utr += lengths[2] 
        if start+stop < 2*lengths[0]:
            cc3utr += 1
        elif start+stop < 2*(lengths[0] + lengths[1]):
            cc4mid += 1
        else:
            cc5utr += 1
        region_list.append((gene, start, stop, pos, rpm))
        gene_list.add(gene)
        cc_use += 1
    data_file.close()
    with open("gene_list.txt",'w') as outfile:
        for gene in gene_list:
            outfile.write(gene+"\n")
    show('There are %s binding regions from %s genes.\n'%(len(region_list), len(gene_list)))
#    show("%d are imported while using %d with ratio %.3f\n"%(cc_in, cc_use, cc_use/float(cc_in)))
#    show("%d hit splicing sites, while %s are missed\n"%(cc_hit, cc_miss))
#    show("Total  -> 5UTR: %d, CDS: %d, 3UTR: %d\n"%(su5utr, su4mid, su3utr))
#    show("Region -> 5UTR: %d, CDS: %d, 3UTR: %d\n"%(cc5utr, cc4mid, cc3utr))
#    show("Ratio  -> 5UTR: %.5f, CDS: %.5f, 3UTR: %.5f\n"%(cc5utr/float(su5utr), cc4mid/float(su4mid), cc3utr/float(su3utr)))
#    import numpy as np
#    show('Bin\tRegion Length\tRPM*0.1\n')
#    count1, bins = np.histogram(np.array(len_dist), bins=range(100))
#    count2, bins = np.histogram(np.array(rpm_dist)*0.1, bins=range(100))
#    for bin in zip(bins, count1, count2):
#        for e in bin:
#            show(e)
#        show('\n')
    return region_list

def choose_candidate(region_list, minRPM):
    ## read pars data set
    pars = import_pars()

    ## read genome map
    genome = map2genome()

    ## create candidate list for training and compare score distributions
    ValRpm, ValMean, ValFull = [], [], []
    candidate_file = open("candidate_list.txt", 'w')
    for gene in pars:
        ValFull.extend([float(val) for val in pars[gene]["PARS"]])
    cc_can = 0

    caA = caC = caG = caU = 0
    crA = crC = crG = crU = 0
    for gene, start, stop, pos, rpm in region_list:
        if float(rpm) < minRPM:
            continue
        seq = pars[gene]["FOLD_SEQ"]
        reg = seq[start:(stop+1)]
        caA += seq.count('A')
        caC += seq.count('C')
        caG += seq.count('G')
        caU += seq.count('U')
        crA += reg.count('A')
        crC += reg.count('C')
        crG += reg.count('G')
        crU += reg.count('U')

        scores = pars[gene]["PARS"]
        region_scores = [float(val) for val in scores[start:(stop+1)]]

        ValRpm.append(float(rpm))
        ValMean.append(mean_std(region_scores)[0])

        lengths = get3length(genome[gene])
        candidate_file.write("%s\tRPM%s\t%s~%s\t%s-%s-%s\n"%(gene, rpm, 
            start+1-lengths[0], stop+1-lengths[0], lengths[0], lengths[1], lengths[2]))
        cc_can += 1
    candidate_file.close()
    if len(ValMean) > 0:
        show("Nucle\tA\tC\tG\tU\n")
        show("Whole")
        show(round(caA/float(caA+caC+caG+caU),4))
        show(round(caC/float(caA+caC+caG+caU),4))
        show(round(caG/float(caA+caC+caG+caU),4))
        show(round(caU/float(caA+caC+caG+caU),4))
        show("\n")
        show("Region")
        show(round(crA/float(crA+crC+crG+crU),4))
        show(round(crC/float(crA+crC+crG+crU),4))
        show(round(crG/float(crA+crC+crG+crU),4))
        show(round(crU/float(crA+crC+crG+crU),4))
        show("\n")
        show("Coef: %.3f, Pearson: %.3f\n"%correlation(ValRpm, ValMean))
        show("Whole RNA is %.3f + %.3f\n"%mean_std(ValFull))
        show("Binding Region is %.3f + %.3f\n"%mean_std(ValMean))
    return cc_can

def validate_prediction(region_list, maxRPM):
    ## read predicted scores
    scoreIdx = read_table("predicted_domain_labels.csv", "pos")
    scoreVal = read_table("predicted_domain_labels.csv", "pred")

    ## check prediction results gene by gene
    validate_file = open("validate_table.txt", 'w')
    validate_file.write("index,pos,label,score\n")

    cc_gene = 1
    for gene in scoreIdx:
        idx = [int(val)-1 for val in scoreIdx[gene]]
        val = [float(val) for val in scoreVal[gene]]
        lab = [0 for i in idx]
        for r_gene, r_start, r_stop, r_pos, r_rpm, in region_list:
            if gene != r_gene:
                continue
            if maxRPM >0 and float(r_rpm) >= maxRPM:
                continue
            for i in xrange(r_start, r_stop+1):
                if i in idx:
                    lab[idx.index(i)] = 1
        if sum(lab) == 0 or sum(lab) == len(lab):
            continue
        for i in xrange(len(idx)):
            validate_file.write("%s,%s,%s,%s\n"%(gene, idx[i], lab[i], val[i]))
        cc_gene += 1
    validate_file.close()
    os.system("Rscript ../src/validate_kim.R")

def main(action, parameter):
    if action == "help":
        print """
Welcome to use validate_kim.py script

The command format are:
  WorkPath=[a_valid_path]
  Action=[test| choose_list| validate_list| choose_validate]
  Paras=[minRPM_Place | maxRPM_Place| RPM_Place_paras]
  Place=[All| 5UTR| CDS| 3UTR]

Examples:
  python src/validate_kim.py WorkPath=work Action=choose_list Paras=Puf3p_0_All
  python src/validate_kim.py WorkPath=work Action=validate_list Paras=gPAR_1000_All

Or using one command:
  python src/validate_kim.py WorkPath=work Action=choose_validate 
                             Paras=10000_All_40_SeqBinary
            """
        return
    kim_data = '../data/KimLab/120214FREEZE.txt' ## old data set

    if action == "test":
        kim_data = '../data/KimLab/Puf3p_PAR-CLIP_crosslinking_sit.txt'
        choose_candidate(get_region_list("All", kim_data), 100000)
        return

    ## read parameters
    para_list = parameter.split('_')
    Data, RPM, Place = para_list[:3]
    if Data == 'Puf3p':
        kim_data = '../data/KimLab/Puf3p_PAR-CLIP_crosslinking_sit.txt'
    elif Data == 'gPAR':
        kim_data = '../data/KimLab/gPAR-CLIP_crosslinking_sites_top.txt'
    print 'Use data set', kim_data

    ## run all
    if action == "all_cases":
        import learn_domain
        for model in ["40_SeqBinary", "40_SeqRatio", "40_ProbVS", "40_PARS", "40_PredSS3"]:
            for RPM, Place in [("1000", "5UTR"), ("1000", "CDS",),\
                               ("10000", "3UTR"), ("10000", "All")]:
                choose_candidate(get_region_list(Place, kim_data), float(RPM))
                learn_domain.main("predict_kim", Place + "_" + model)
                validate_prediction(get_region_list(Place, kim_data), float(RPM))
        return
    if action == "data_sizes":
        for Place in ['All', '5UTR', 'CDS', '3UTR']:
            show(Place)
            data = get_region_list(Place, kim_data)
            for RPM in [0, 10, 100, 1000, 2000]:
                show(choose_candidate(data, float(RPM)))
            show('\n')
        return

    if action == "choose_list" or action == "choose_validate":
        choose_candidate(get_region_list(Place, kim_data), float(RPM))
    if action == "choose_validate":
        import learn_domain
        learn_domain.main("predict_kim", "_".join(para_list[2:]))
    if action == "validate_list" or action == "choose_validate":
        validate_prediction(get_region_list(Place, kim_data), float(RPM))

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

