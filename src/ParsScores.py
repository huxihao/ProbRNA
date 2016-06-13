## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

## Deal with PARS score data files
## Author: Xihao Hu <huxihao@gmail.com>
## Version 1.0 - January 22, 2011

import os,sys,string,time,math
from tools import *

def read_pars(DataPath="../data/PARS/"):
    ## Function: reading data set from Kertesz etal. Nature 2010
    ## Defination of data structure:
    ## 1. saved in a dictionary, each element represents a gene
    ##    ID: ID from sce_genes.fasta
    ## 2. for each gene id, attributes save in a dictionary, we have:
    ##    "SEQ": DNA sequence from sce_genes_fasta (string)
    ##    "TYPE": type of transcript (mRNA, rRNA, snoRNA, tRNA, snRNA, ncRNA)
    ##    "V1": a list of V1 scores from sce_V1.tab (list of int)
    ##    "S1": a list of S1 scores from sce_S1.tab (list of int)
    ##    "PARS": a list of PARS scores from sce_S1.tab (list of float)
    ##    "FOLD_SEQ": RNA sequence from sce_genes_folded.tab
    ##    "FOLD_SS": predicted folding results from sce_genes_rnafold_ss.tab
    ##    "FOLD_PROB_L"/"FOLD_PROB_R": from sce_genes_rnafold_prob.tab
    ## In checking the constance, the reading order is:
    ##  Gene List File -> PARS Score File -> Folding File -> Type File
    ## Unmatched data will been given warnings.
    data = {}
    records = read_fasta(DataPath + "sce_genes.fasta")
    for eid, seq in records:
        if data.has_key(eid):
            warning("Warning: duplicated keys for gene id of %s"%eid)
        else:
            one_gene = {}
            one_gene["SEQ"] = seq
            data[eid] = one_gene
    ## read PARS scores
    infile = open(DataPath + "sce_Score.tab",'r')
    for line in infile:
        (Id, Len, Scores) = line.split('\t')
        if data.has_key(Id):
            data[Id]["PARS"] = [float(val) for val in Scores.strip().split(';')]
            ## check the length of sequence and length of scores
            if len(data[Id]["SEQ"]) != len(data[Id]["PARS"]):
                warning("unmatched lengthes for %s"%Id)
        else:
            warning("no gene record for %s"%Id)
            one_gene = {}
            one_gene["PARS"] = [float(val) for val in Scores.strip().split(';')]
            data[Id] = one_gene
    infile.close()
    infile = open(DataPath + "sce_V1.tab",'r')
    for line in infile:
        (Id, Len, Scores) = line.split('\t')
        data[Id]["V1"] = [int(val) for val in Scores.strip().split(';')]
    infile.close()
    infile = open(DataPath + "sce_S1.tab",'r')
    for line in infile:
        (Id, Len, Scores) = line.split('\t')
        data[Id]["S1"] = [int(val) for val in Scores.strip().split(';')]
    infile.close()
    ## read predicted probability of pairing from RNAstructure
    infile = open(DataPath + "sce_genes_rnafold_prob.tab",'r')
    for line in infile:
        (Id, Len, ProbLeft, ProbRight) = line.split('\t')
        data[Id]["FOLD_PROB_L"] = [float(val) for val in ProbLeft.strip().split(';')]
        data[Id]["FOLD_PROB_R"] = [float(val) for val in ProbRight.strip().split(';')]
    infile.close()
    ## read predicted folding from PARS
    #infile = open(DataPath + "sce_genes_folded.tab",'r')
    infile = open(DataPath + "sce_genes_rnafold_ss.tab",'r')
    for line in infile:
        (Id, Seq, Fold) = line.strip().split('\t')
        if data.has_key(Id):
            one_gene = data[Id]
            one_gene["FOLD_SEQ"] = "%s"%Seq
            one_gene["FOLD_SS"] = Fold
            if not one_gene.has_key("SEQ"):
                warning("missing DNA sequence but contain RNA for %s"%Id)
            elif len(one_gene["SEQ"]) != len(Seq):
                warning("unmatched sequences for %s\n_%s\n_%s"%(Id,
                        one_gene["SEQ"], Seq))
        else:
            warning("missing record for %s"%Id)
    infile.close()
    ## read category of genes
    type_count = {}
    infile = open(DataPath + "type_of_genes.txt",'r')
    for line in infile:
        (Id, gene_type) = line.strip().split('\t')
        type_count[gene_type] = type_count.setdefault(gene_type, 0) + 1
        if data.has_key(Id):
            one_gene = data[Id]
            one_gene["TYPE"] = gene_type
        else:
            hasMatchedGeneName = False
            for gene_name in data:
                if gene_name.startswith(Id):
                    data[gene_name]["TYPE"] = gene_type
                    hasMatchedGeneName = True
            if not hasMatchedGeneName:
                warning("missing name of %s with type of %s"%(Id,gene_type))
    for type_name in type_count:
        print "Type of genes:", type_name, type_count[type_name]
    infile.close()
    ## fix control sample
    data["RPR1"] = data.pop("RPR1-di")
    data["RPR1"]["TYPE"] = 'ncRNA'
    data["YKL185W"] = data.pop("YKL185W-di")
    data["YKL185W"]["TYPE"] = 'mRNA'
    data["YNL229C"] = data.pop("YNL229C-di")
    data["YNL229C"]["TYPE"] = 'mRNA'
    ## return the data structure
    print "# Read in PARS with %s genes"%len(data)
    return data

def export_pars(pars, gene_list=[], max_length=-1, outfilename="PARS_DATA_SET"):
    outfile = open(outfilename, 'w')
    cc = 0
    if len(gene_list) == 0:
        gene_list = pars.keys()
    for gene in gene_list:
        if not pars[gene].has_key("FOLD_SEQ"):
            continue
        seq = pars[gene]["FOLD_SEQ"]
        if max_length > 0 and len(seq) > max_length:
            outfile.write("%s\t*\t.\n"%(gene))
            outfile.write("0\n0\n0\n")
            continue
        cc += 1
        pars_gene = pars[gene]
        typ = pars_gene["TYPE"]
        ss = pars_gene["FOLD_SS"]
        score = [str(val) for val in pars_gene["PARS"]]
        v1 = [str(val) for val in pars_gene["V1"]]
        s1 = [str(val) for val in pars_gene["S1"]]
        fold_l = [str(round(val,3)) for val in pars_gene["FOLD_PROB_L"]]
        fold_r = [str(round(val,3)) for val in pars_gene["FOLD_PROB_R"]]
        outfile.write("%s\t%s\t%s\t%s\n"%(gene, typ, seq, ss))
        outfile.write("%s\n%s\n%s\n"%(';'.join(score), ';'.join(v1),';'.join(s1)))
        outfile.write("%s\n%s\n"%(';'.join(fold_l),';'.join(fold_r)))
    outfile.close()
    print "Export", cc, "PARS data"

def import_pars(infilename="PARS_DATA_SET"):
    pars = {}
    if not os.path.exists(infilename):
        export_pars(read_pars(), outfilename=infilename)
    infile = open(infilename, 'r')
    cc = 0
    gene_name = ""
    one_gene = {}
    for line in infile:
        if cc%6 == 0:
            gene_name, gene_type, gene_seq, gene_ss = line.strip().split('\t')
            one_gene = {}
            one_gene["TYPE"] = gene_type
            one_gene["FOLD_SEQ"] = gene_seq
            one_gene["FOLD_SS"] = gene_ss
        elif cc%6 == 1:
            one_gene["PARS"] = [float(val) for val in line.strip().split(';')]
        elif cc%6 == 2:
            one_gene["V1"] = [int(val) for val in line.strip().split(';')]
        elif cc%6 == 3:
            one_gene["S1"] = [int(val) for val in line.strip().split(';')]
        elif cc%6 == 4:
            one_gene["FOLD_PROB_L"] = [float(val) for val in line.strip().split(';')]
        elif cc%6 == 5:
            one_gene["FOLD_PROB_R"] = [float(val) for val in line.strip().split(';')]
            pars[gene_name] = one_gene
        cc += 1
    print "Import", cc/6, "PARS data"
    return pars

def descriptions_from_GCD(filename="../data/PARS/SGD_features.tab"):
    infile = open(filename, 'r')
    data = {}
    for line in infile:
        ele = line[:-1].split('\t')
        gene_name = ele[3] # 4.   Feature name (optional)
        standard = ele[4] # 5.   Standard gene name (optional)
        description = ele[15] # 16.  Description (optional)
        if gene_name not in data:
            data[gene_name] = "%s\t%s"%(standard, description)
    infile.close()
    return data

def poisson_cutoff(scores, window):
    ## histogram
    import numpy
    hist, bin_edges = numpy.histogram(scores, range(0,max(scores)+window,window))
    start, end = 0, 0
    for i in range(len(hist)):
        if hist[i] == max(hist):
            start = bin_edges[i]
            end = bin_edges[i+1]
            break
    mean, cc = 0, 0
    for s in scores:
        if start <= s and s < end:
            mean += s
            cc += 1
    mean = mean/float(cc)
    from scipy.stats import poisson
    return poisson.isf(q, mean), mean

def poisson_pvalue(scores, window):
    ## histogram
    import numpy
    hist, bin_edges = numpy.histogram(scores, range(0,max(scores)+window,window))
    start, end = 0, 0
    for i in range(len(hist)):
        if hist[i] == max(hist):
            start = bin_edges[i]
            end = bin_edges[i+1]
            break
    mean, cc = 0, 0
    for s in scores:
        if start <= s and s < end:
            mean += s
            cc += 1
    mean = mean/float(cc)
    from scipy.stats import poisson
    pvalues = []
    for s in scores:
        pvalues.append(1.0-poisson.cdf(s, mean))
    return pvalues, mean

def poisson_chi2(obs):
    ## histogram
    window = 2
    import numpy
    hist, bin_edges = numpy.histogram(obs, \
            range(0,max(obs)+window,window))
    start, end = 0, 0
    for i in range(len(hist)):
        if hist[i] == max(hist):
            start = bin_edges[i]
            end = bin_edges[i+1]
            break
    mean, cc = 0, 0
    for ob in obs:
        if start <= ob and ob < end:
            mean += ob
            cc += 1
    mean = mean/float(cc) ## expected to be lambda
    from scipy.stats import poisson
    f_obs = []
    f_exp = []
    last_cdf = 0.0
    for v in range(max(obs)+1):
        f_obs.append(obs.count(v))
        now_cdf = poisson.cdf(v, mean)
        f_exp.append(now_cdf-last_cdf)
        last_cdf = now_cdf
        if last_cdf > 0.99: ## ignore outliers
            break
    obs_sum = sum(f_obs)
    for i in range(len(f_obs)):
        f_exp[i] *= obs_sum
        #print '%d    %.2f\t%.2f'%(i, f_obs[i], f_exp[i])
    from scipy.stats import chisquare
    import scipy
    ## degree: k-1-p, default p=0
    chi2, pvalue = chisquare(scipy.array(f_obs), scipy.array(f_exp), 0)
    return pvalue, mean

def pars_gene_groups(pars):
    gene_list = []
    for key in pars.keys():
        if pars[key].has_key("FOLD_SEQ"): ## has valid RNA sequence
            gene_list.append(key)
    gene_list.sort()

    ## search genes with same RNA sequences
    groups = []
    cc1 = 0
    for geneA in gene_list:
        parsA = pars[geneA]
        for geneB in gene_list:
            if geneA == geneB:
                break
            parsB = pars[geneB]
            if parsA["FOLD_SEQ"] == parsB["FOLD_SEQ"]:
                s1 = parsA["PARS"]
                s2 = parsB["PARS"]
                for k in range(len(s1)):
                    if s1[k] != s2[k]: ## different PARS score
                        cc1 += 1
                #print "%s*%s"%(geneA, geneB),
                isContain = False
                for group in groups:
                    if geneA in group or geneB in group:
                        isContain = True
                        if geneA not in group:
                            group.append(geneA)
                        if geneB not in group:
                            group.append(geneB)
                        break
                if not isContain:
                    groups.append([geneA, geneB])
    print "Number of Mismatch PARS scores", cc1

    ## APPLYING SOME ANALYSIS
    if False: ## save group information
        for group in groups:
            group.sort()
            for member in group:
                show(member)
            show('\n')

    if True: ## combine group and export pars
        ## only use the first member in each group
        for group in groups:
            group.sort()
            for member in group[1:]: ## remove others in gene_list
                gene_list.remove(member)
    return gene_list, groups

def get_pars_score(v1, s1):
    sum_v = 0.834
    sum_s = 1
    kv = float(sum_v + sum_s)/float(sum_v + sum_v)
    ks = float(sum_v + sum_s)/float(sum_s + sum_s)
    #kv = 2**1.07 - 1
    #ks = 2**0.93 - 1
    return math.log((kv*v1+1)/(ks*s1+1),2)

def main():
    pars = read_pars()
    gene_list = []
    for key in pars.keys():
        if pars[key].has_key("FOLD_SEQ"):
            gene_list.append(key)
    gene_list.sort()
    print "First 10 data points are:"
    for gene in gene_list[0:100:10]:
        attr_list = pars[gene]
        print "%7s"%gene,
        for attr in attr_list:
            print attr, #len(pars[gene][attr]),
        print ''
    #print pars["RDN18-1"]["PARS"]

    ## search genes with same RNA sequences
    gene_list, groups = pars_gene_groups(pars)
    export_pars(pars)
    pars = import_pars()

    if True: ## 3' term counts and pvalues for V1 and S1
        ## get global distribution
        v1cc = [0 for i in range(50)]
        s1cc = [0 for i in range(50)]
        sum_v1, sum_s1 = 0, 0
        cc_all_pos = 0
        cc_use_pos = 0
        for gene in gene_list:
            if "FOLD_SEQ" not in pars[gene]:
                continue
            gene_type = pars[gene].get("TYPE")
            seq = pars[gene]["FOLD_SEQ"]
            scores = pars[gene]["PARS"]
            v1 = pars[gene]["V1"]
            s1 = pars[gene]["S1"]

            ## check RNA structure and expression level
            if False:
                show(math.log(sum(v1)/len(seq) + 1))
                show(math.log(sum(s1)/len(seq) + 1))
                show(sum([abs(val) for val in scores])/len(seq))
                show("\n")

            ## calculate global counts
            for i in range(50):
                v1cc[i] += v1[len(seq)-i-1]
                s1cc[i] += s1[len(seq)-i-1]
            for i in range(len(seq)):
                cc_all_pos +=1
                #if v1[i] == 0 and s1[i] == 0:
                if scores[i] > 5 or scores[i] <-5:
                    cc_use_pos +=1
            sum_v1 += sum(v1)
            sum_s1 += sum(s1)

            if False: ## check the calculation of PARS score
                for j in range(len(scores)):
                    if v1[j] > 50 or s1[j] > 50:
                        print scores[j], "%.2f"%get_pars_score(v1[j],s1[j])
            elif False: ## p-values
                ## remove 0 counts near 3' term
                ## re-define the range of all elements
                seq = seq[0:-24]
                v1 = v1[0:-24]
                s1 = s1[0:-24]
                v1_chi2, v1_mean = poisson_chi2(v1)
                s1_chi2, s1_mean = poisson_chi2(s1)
                ## calculate p-values
                pos_list = []
                v1_pvalue, v1_mean = poisson_pvalue(v1, 2)
                s1_pvalue, s1_mean = poisson_pvalue(s1, 2)
                for i in range(len(seq)):
                    if v1_pvalue[i] < 0.01 and s1_pvalue[i] < 0.01\
                            and scores[i] < 3 and scores[i] > -3:
                        pos_list.append(i)
                if len(pos_list):
                    show(gene)
                    show(len(pos_list))
                    show('\n')
        ## save global counts
        #for i in range(50):
        #    show('%s\t%s\t%s\n'%(i+1, v1cc[i], s1cc[i]))
        print cc_all_pos, cc_use_pos, cc_use_pos/float(cc_all_pos)
      #  show(sum_v1)
      #  show(sum_s1)
      #  show(sum_v1/float(sum_s1))
      #  show('\n')

if __name__ == "__main__":
    ## Default settings
    WORK_PATH = "../work/"
    ## Read parameters
    for arg in sys.argv:
        ## example: python AA.py WorkPath=C:/aa/
        if arg.startswith("WorkPath"):
            WORK_PATH = arg[(string.find(arg,"=")+1):]
    ## Setting work path
    os.chdir(WORK_PATH)
    show("# Current work path is %s\n"%os.path.abspath("."))
    ## Delete old log file
    if os.path.exists(LOG_FILE):
        os.remove(LOG_FILE)
    START_TIME = time.time()
    main()
    END_TIME = time.time()
    show(time.strftime("# Start  at: %c\n", time.gmtime(START_TIME)))
    show(time.strftime("# Finish at: %c\n", time.gmtime(END_TIME)))
    show(time.strftime("# Total Time: %j day %H:%M:%S\n",\
            time.gmtime(END_TIME-START_TIME))) 

