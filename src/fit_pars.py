## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

from ParsScores import *

from multiprocessing import Pool

def ss_to_lr(ss):
    "Secondary Structure to Loop Region"
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

def extract_data(outfilename, enzyme, spanLen=30, minRPKM=10, shortList=[], saveZero=True):
    pars_data = import_pars()
    from gene_mapping import map2genome
    genome = map2genome()
    gene_list, groups = pars_gene_groups(pars_data)
    read_sum = 0
    for gene in pars_data:
        read_sum += sum(pars_data[gene][enzyme])
    print "Number of total reads in", enzyme, "is", read_sum

    ## choose a list of genes
    choose_gene = []
    for gene in gene_list:
        gene_cell = pars_data[gene]
        counts = gene_cell[enzyme]
        rpkm = 10**9 * sum(counts) / float(read_sum * (len(counts)))
        if gene_cell.get("TYPE") != "mRNA": ## don't use non-coding RNAs to train
            rpkm = 0
        #show(gene);show(rpkm);show("\n")
        choose_gene.append([gene, rpkm, counts])

    ## sort by RPKM
    #choose_gene.sort(lambda x,y: cmp(y[1],x[1]))

    ## extract pars data of selected genes
    outfile = open(outfilename, 'w')
    cc = 0
    outfile.write("index,tag,seq,count,gene\n")
    for gene, rpkm, counts in choose_gene:
        gene_all = pars_data[gene]
        seq = gene_all["FOLD_SEQ"].replace('U','T') ## Use DNA
        pred_ss = gene_all["FOLD_SS"]
        genome_mark = []
        if gene in genome:
            genome_mark = genome[gene]
        lens = [0,0,0] ## [5UTR, CDS, 3UTR]
        for g_mark, g_start, g_end in genome_mark[1:]:
            g_len = g_end - g_start
            if g_len >= 0: g_len = 1 + g_len
            else: g_len = 1 - g_len
            if g_mark == "5UTR": lens[0] += g_len
            elif g_mark == "Exon": lens[1] += g_len
            elif g_mark == "3UTR": lens[2] += g_len
        ## validate whether the lengths are valid
        #if lens[1] % 3 != 0:
        #    warning("Wrong CDS length of %s in %s"%(lens[1],gene))
        if sum(lens) != len(seq):
            warning("Unmatched lengths: %s != %s"%(sum(lens),len(seq)))

        ## mseq format
        cc += 1
        pred_lr = ss_to_lr(pred_ss)
        for i in range(len(seq)): ## move for each position
            tag = "1" ## CDS
            if i < lens[0]: ## 5UTR
                tag = "5"
            elif i >= lens[0]+lens[1]: ## 3UTR
                tag = "3"
            ## 24 is the length of blind tail
            if i < spanLen or i >= len(seq)-max(24, spanLen):
                tag = "-" + tag ## marked as negative
            if rpkm < minRPKM:
                tag = "0" ## gene with low expression
            if shortList != [] and gene not in shortList:
                tag = "0" ## not the choosen one
            ## Encoding Scheme
            seq_code = ""
            seq_code += seq[i].upper()
        #    if pred_ss[i] == '.': seq_code += "S" 
        #    else: seq_code += "P"
        #    if pred_lr[i] == 'o': seq_code += "O"
        #    else: seq_code += "N"

            out_gene = ""
            if i == 0:
                out_gene = gene
            if tag == "0" and saveZero == False:
                continue
            outfile.write(",".join([str(cc), tag, seq_code, str(counts[i]), out_gene])+"\n")
    outfile.close()

def read_table(datafilename, gcol=4, vcol=5):
    import csv
    data_table = csv.reader(open(datafilename, "rb"), delimiter=',')
    data = {}
    head = data_table.next()
    head[gcol] = "[" + head[gcol] + "]"
    head[vcol] = "<" + head[vcol] + ">"
    print "Read:", "-".join(head),
    index = 0
    gene = ""
    values = []
    for row in data_table:
        if row[gcol] != "":
            index += 1
            if gene != "":
                data[gene] = values
            gene = row[gcol]
            values = []
        if int(row[0]) != index:
            warning("Wrong indexes! %s != %s\n"%(row[0], index))
            exit()
        values.append(row[vcol])
    data[gene] = values
    print "with", len(data), "genes."
    return data

def validate_regions(v1file, s1file=""):
    import math
    if s1file == "": ## combined file
        TAG = read_table(v1file, 5,1)
        SEQ = read_table(v1file, 5,2)
        V1count = read_table(v1file, 5,3)
        S1count = read_table(v1file, 5,4)
        V1pred = read_table(v1file, 5,6)
        S1pred = read_table(v1file, 5,7)
        V1prob = read_table(v1file, 5,8)
        S1prob = read_table(v1file, 5,9)
    else: ## seperate file
        TAG = read_table(v1file, 4,1)
        SEQ = read_table(v1file, 4,2)
        V1count = read_table(v1file, 4,3)
        S1count = read_table(s1file, 4,3)
        V1pred = read_table(v1file, 4,6)
        S1pred = read_table(s1file, 4,6)
        V1prob = read_table(v1file, 4,7)
        S1prob = read_table(s1file, 4,7)

    outfile = open("roc_table.csv", 'w')
    outfile.write("domain,real,PARS,V1,S1,Prob.V1,Prob.S1,Test\n")
    for gene, base, start, end, real_ss, pred_ss, pars, v1, s1 in read_real_domains():
        ## has both V1 and S1
        if gene not in V1count or gene not in S1count:
            continue
        tag = TAG[gene]
        seq = SEQ[gene]
        for i in xrange(start-1, end):
            outfile.write("%s[%s-%s],"%(gene, start, end))
        #    outfile.write("%s,"%gene)
#            outfile.write("%s,"%seq[i])
            if real_ss[i-start+1] == '.':
                outfile.write("0,")
            else:
                outfile.write("1,")
            outfile.write("%s,"%(pars[i-start+1]))
            outfile.write("%s,"%(v1[i-start+1]))
            outfile.write("%s,"%(s1[i-start+1]))
            outfile.write("%s,"%(float(V1prob[gene][i])))
            outfile.write("%s,"%(float(S1prob[gene][i])))
            outfile.write("%s\n"%(float(V1prob[gene][i])-float(S1prob[gene][i])))
    for id, gene, seq, real_ss, pred_ss, pars, v1, s1 in read_blast_map():
        ## has both V1 and S1
        if gene not in V1count or gene not in S1count:
            continue
        tag = TAG[gene]
        seq = SEQ[gene]
        for i in xrange(len(seq)):
            outfile.write("%s_%s,"%(id, gene))
#            outfile.write("%s,"%seq[i])
            if real_ss[i] == '.':
                outfile.write("0,")
            else:
                outfile.write("1,")
            outfile.write("%s,"%(pars[i]))
            outfile.write("%s,"%(v1[i]))
            outfile.write("%s,"%(s1[i]))
            outfile.write("%s,"%(float(V1prob[gene][i])))
            outfile.write("%s,"%(float(S1prob[gene][i])))
            outfile.write("%s\n"%(float(V1prob[gene][i])-float(S1prob[gene][i])))
    outfile.close()
    ## DO PLOT
    os.system("Rscript ../src/plot_roc.R")
    plot_file = "validate_regions_plot.pdf"
    if os.path.exists(plot_file):
        os.remove(plot_file)
    os.rename("plot_roc.pdf", plot_file)


def find_regions(v1file, s1file):
    TAG = read_table(v1file, 1)
    SEQ = read_table(v1file, 2)
    V1count = read_table(v1file, 3)
    S1count = read_table(s1file, 3)
    V1prob = read_table(v1file, 6)
    S1prob = read_table(s1file, 6)

    for line in open("../data/ASH1/mRNA_transport.txt",'r'):
        gene, oldName, description = line.strip().split('\t')
        ## has both V1 and S1
        if gene not in V1count or gene not in S1count:
            continue

        fiveUTR = 0
        if True: ## get the length of 5UTR
            for tag in TAG[gene]:
                if int(tag) <= -2:
                    fiveUTR += 1
                else:
                    break

        window = 100
        scores = []
        if True:
            v1 = [float(score) for score in V1prob[gene]]
            s1 = [float(score) for score in S1prob[gene]]
            scores.append([sum(v1[0:window]), sum(s1[0:window])])
            for i in xrange(0, len(v1) - window):
                sum_v1 = scores[i][0] - v1[i] + v1[i+window]
                sum_s1 = scores[i][1] - s1[i] + s1[i+window]
                scores.append([sum_v1, sum_s1])
        scores_sort = sorted(scores, cmp=lambda x,y: cmp(combine(y), combine(x)))
        cutoff = combine(scores_sort[int(0.05 * len(scores_sort))])
    #    if cutoff < 30:
    #        continue
        show(gene)
        show(oldName)
        show(fiveUTR)
        show(len(v1))
        show(round(cutoff, 1))

        ## imerge regions
        left_left = -1
        left_right = -1
        for left in xrange(0, len(scores)):
            if combine(scores[left]) >= cutoff:
            #    print gene, left+1, left+window
                if left_left < 0:
                    left_left = left
                    left_right = left
                    continue
                if left <= left_right + window:
                    left_right = left
                else:
                    show("%s-%s"%(left_left+1 -fiveUTR, left_right+window -fiveUTR))
                    left_left = -1
        if left_left >= 0: ## last one
            show("%s-%s"%(left_left+1 -fiveUTR, left_right+window -fiveUTR))
        show("\n")
            
def read_real_domains(filename = "../data/ASH1/known_domains.txt"):
    infile = open(filename, 'r')
    data = []
    gene = ''
    start, end = 0, 0 
    seq, real_ss, pred_ss = '','',''
    pars, v1, s1 = [],[],[]
    for line in infile:
        ele = line.strip().split('\t')
        if len(ele) == 0:
            continue
        if line.startswith("Domain"):
            if gene != '':
                data.append([gene, seq, start, end, real_ss, pred_ss, pars, v1, s1])
            gene = ele[1]
            start = int(ele[2])
            end = int(ele[3])
            seq, real_ss, pred_ss = '','',''
            pars, v1, s1 = [],[],[]
        elif line.startswith("Sequence"):
            pass
        else:
            seq += ele[0]
            real_ss += ele[2]
            pred_ss += ele[3]
            pars.append(float(ele[4]))
            v1.append(int(ele[5]))
            s1.append(int(ele[6]))
    return data

def read_blast_map(filename="rna_blast_map.txt"):
    infile = open(filename, 'r')
    data = []
    infile.readline().replace('\t', '\n')
    cc = 0
    for line in infile:
        ele = line.strip().split('\t')
        id = ele[0]
        real_ss = ele[3]
        gene = ele[4].split(';')[0]
        seq = ele[5].replace('U','T')
        pred_ss = ele[6]
        pars = [float(val) for val in ele[7].split(';')]
        v1 = [int(val) for val in ele[8].split(';')]
        s1 = [int(val) for val in ele[9].split(';')]
        data.append((id, gene, seq, real_ss, pred_ss, pars, v1, s1))
    return data


####################################################################
## Fit V1 and S1 counts together, use union data
def extract_data_combine(V1, S1, V1_S1):
    v1file = open(V1, 'r')
    s1file = open(S1, 'r')
    outfile = open(V1_S1, 'w')
    outfile.write("index,tag,seq,v1,s1,gene\n")
    cc, ccv, ccs = 0, 0, 0
    ## read headers
    v1line = v1file.readline()
    s1line = s1file.readline()
    while True:
        v1line = v1file.readline()
        s1line = s1file.readline()
        if v1line == "" or s1line == "":
            break
        index, tagv, seq, v1, gene = v1line.strip().split(',')
        index, tags, seq, s1, gene = s1line.strip().split(',')
        if tagv == "0" or tags == "0":
            if gene != "" and tagv == "0":
                ccv += 1
            if gene != "" and tags == "0":
                ccs += 1
            continue
        if gene != "":
            cc += 1
        tag = str(min(int(tagv), int(tags)))
        outfile.write(",".join([str(cc), tag, seq, v1, s1, gene]) + "\n")
    v1file.close()
    s1file.close()
    outfile.close()
    print "Combine", cc, "data points, and remove", ccv, "from V1", ccs, "from S1"

def validate_zipcode(fit_file, window=10, zipcode_file="../data/ASH1/zipCodeList.txt"):
    ''' validate zipcode regions '''
    import math
    ## Read from fitted data file
    TAG = read_table(fit_file, 5,1)
    SEQ = read_table(fit_file, 5,2)
    V1count = read_table(fit_file, 5,3)
    S1count = read_table(fit_file, 5,4)
    V1pred = read_table(fit_file, 5,6)
    S1pred = read_table(fit_file, 5,7)
    V1prob = read_table(fit_file, 5,8)
    S1prob = read_table(fit_file, 5,9)
    zipcodes = [line.strip().split('\t') for line in open(zipcode_file,'r')]
    outfile = open("roc_table.csv", 'w')
    outfile.write("domain,real,V1,S1,PARS,PARS^2,V1*S1,Prob\n")
    for gene, zipcode, region in zipcodes:
        if gene not in TAG:
            continue
        ## get the length of 5UTR CDS 3UTR
        lens = [0, 0, 0]
        for tag in TAG[gene]:
            if abs(int(tag)) == 5:
                lens[0] += 1
            elif abs(int(tag)) == 1:
                lens[1] += 1
            elif abs(int(tag)) == 3:
                lens[2] += 1
        print gene, zipcode, sum(lens),
        regionL, regionR = region.split('-')
        regionL = int(regionL) + lens[0]
        regionR = int(regionR) + lens[0]
        print regionL, regionR
        ## output region list in ROC format
        _v1 = [int(val) for val in V1count[gene]]
        _s1 = [int(val) for val in S1count[gene]]
        _pars = [math.log((_v1[i]+1.0)/(_s1[i]+1.0)) for i in xrange(len(_v1))]
        _pars2 = [val*val for val in _pars]
        _v1s1 = [_v1[i]*_s1[i] for i in xrange(len(_v1))]
        _prob = [float(val) for val in V1prob[gene]]
        sum_v1 = sum_s1 = sum_pars = sum_pars2 = sum_prob = 0
        for i in xrange(len(_v1)-window+1):
            if i == 0: ## init
                sum_v1 = sum(_v1[0:window])
                sum_s1 = sum(_s1[0:window])
                sum_pars = sum(_pars[0:window])
                sum_pars2 = sum(_pars2[0:window])
                sum_v1s1 = sum(_v1s1[0:window])
                sum_prob = sum(_prob[0:window])
            else:
                sum_v1 = sum_v1 - _v1[i-1] + _v1[i+window-1]
                sum_s1 = sum_s1 - _s1[i-1] + _s1[i+window-1]
                sum_pars = sum_pars - _pars[i-1] + _pars[i+window-1]
                sum_pars2 = sum_pars2 - _pars2[i-1] + _pars2[i+window-1]
                sum_v1s1 = sum_v1s1 - _v1s1[i-1] + _v1s1[i+window-1]
                sum_prob = sum_prob - _prob[i-1] + _prob[i+window-1]
            match = 0 ## not a match
            if (i+1-regionL)*(i+window-regionR) <= 0:
                match = 1 ## match
            #outfile.write("All,")
            outfile.write("%s-%s,"%(gene, zipcode))
            outfile.write("%s,"%match)
            outfile.write("%s,"%sum_v1)
            outfile.write("%s,"%sum_s1)
            outfile.write("%s,"%sum_pars)
            outfile.write("%s,"%sum_pars2)
            outfile.write("%s,"%sum_v1s1)
            outfile.write("%s\n"%sum_prob)
    outfile.close()
    ## DO PLOT
    os.system("Rscript ../src/plot_roc.R")
    plot_file = "validate_zipcode_plot.pdf"
    if os.path.exists(plot_file):
        os.remove(plot_file)
    os.rename("plot_roc.pdf", plot_file)

def compare_3m(WinSize_CutOff = "2_8000"):
    show("# Run compare_3m with parameters %s\n"%WinSize_CutOff)
    ## get R^2 from PL, MP, MPL from synthetic data
    WinSize, CutOff = WinSize_CutOff.split('_')
    extract_data("PARS_V1.csv", "V1", int(WinSize)/2, float(CutOff), saveZero=False)
    extract_data("PARS_S1.csv", "S1", int(WinSize)/2, float(CutOff), saveZero=False)
    os.system("Rscript ../src/fit_pars_3m.R")

def combine_fit_process((action, model, window)):
    short_list = ["YKL185W", "YNL327W", "YLR190W", "YOR247W", "YLL028W",
         "YNL283C", "YPL084W", "YPR119W", "YJL172W", "YLL001W", 
         "YMR202W", "YBR086C", "YGR040W", "YMR296C", "YNL103W",
         "YLR332W", "YGR023W", "YLR434C", "YMR301C", "YJR121W",
         "YEL052W", "YBR003W", "YGL143C", "YER154W", "YPL172C", 
         "YIL022W", "YMR171C"]
    #short_list = ["RDN25-1", "RDN58-1", "RDN18-1", "RDN5-1"]
    term_span = window/2
    v1_table = "%s_%s_%s_PARS_V1.csv"%(action, model, window)
    s1_table = "%s_%s_%s_PARS_S1.csv"%(action, model, window)
    comb_table = "%s_%s_%s_PARS_V1_S1.csv"%(action, model, window)
    if not os.path.exists(comb_table + ".log"):
        extract_data(v1_table, "V1", term_span, 1000)
        extract_data(s1_table, "S1", term_span, 1000)
        extract_data_combine(v1_table, s1_table, comb_table)
        if action == "CV": ## cross-validation
            os.system("Rscript ../src/fit_pars_all.R CV %s %s %s"%(model, window, comb_table))
            os.remove(comb_table)
        else: ## train and predict
            os.system("Rscript ../src/fit_pars_all.R Train %s %s %s"%(model, window, comb_table))
            if action == "PredShort":
                extract_data(v1_table, "V1", term_span, 0, short_list)
                extract_data(s1_table, "S1", term_span, 0, short_list)
            elif action == "PredAll":
                extract_data(v1_table, "V1", term_span, 0, [])
                extract_data(s1_table, "S1", term_span, 0, [])
            extract_data_combine(v1_table, s1_table, comb_table)
            os.system("Rscript ../src/fit_pars_all.R Predict %s %s %s"%(model, window, comb_table))
            os.remove(comb_table+".model")
        os.remove(v1_table)
        os.remove(s1_table)
    else:
        print 'Read existing log file for', comb_table
    ## read log and clean temp files
    with open(comb_table+".log") as log_file:
        return log_file.readline() ## read the first line

def combine_fit(ThrdNum_Model_MaxWin = "1_MixPoiSep_2"):
    show("# Run combine_fit with parameters %s\n"%ThrdNum_Model_MaxWin)
    ThrdNum, Model, MaxWin = ThrdNum_Model_MaxWin.split('_')
    pool = Pool(int(ThrdNum)) ## Number of workers
    pool_para = []
    for action in ["CV", "PredAll"]:
        for winSize in xrange(2, int(MaxWin)+1, 2):
            if Model == "AllModel":
                for model in ["MixPoiSep", "MixPoiLin", "MPLComSam", "MPLComOpp"]:
                    pool_para.append((action, model, winSize))
            else:
                pool_para.append((action, Model, winSize))
        if ThrdNum_Model_MaxWin == "1_MixPoiSep_2":
            break ## for testing the modular and so skip PredAll
    logs = pool.map(combine_fit_process, pool_para)
    for log in logs:
        if log != "":
            show(log)

def plot_data(Model_WinSize = "MixPoiLin_2"):
    show("# Run plot_data with parameters %s\n"%Model_WinSize)
    Model, WinSize = Model_WinSize.split('_')
    extract_data("PARS_V1.csv", "V1", int(WinSize)/2, 1000)
    extract_data("PARS_S1.csv", "S1", int(WinSize)/2, 1000)
    extract_data_combine("PARS_V1.csv", "PARS_S1.csv", "PARS_V1_S1.csv")
    os.system("Rscript ../src/fit_pars_all.R Plot %s %s PARS_V1_S1.csv"%(Model,WinSize))

def main(action, parameter):
    if action == "help":
        print """
Welcome to use fit_pars.py script

The command format are:
  WorkPath=[a_valid_path]
  Action=[test| compare_3m |combine_fit |plot_data]
  Paras=[default| WinSize_CutOff |ThrdNum_Model_MaxWin |Model_WinSize]
  Model=[AllModel| MixPoiSep |MixPoiLin |MPLComSam |MPLComOpp]:

Examples:
  python src/fit_pars.py WorkPath=work Action=compare_3m Paras=60_1000
  python src/fit_pars.py WorkPath=work Action=combine_fit Paras=9_MixPoiSep_50
  python src/fit_pars.py WorkPath=work Action=plot_data Paras=MixPoiLin_14

Or you can use the following command to test all actions with default settings
  python src/fit_pars.py WorkPath=work Action=test
            """
        return
    import_pars() ## initilization
    if action == "compare_3m":
        if parameter == "default": compare_3m()
        else: compare_3m(parameter)
    if action == "combine_fit" or action == "test":
        if parameter == "default": combine_fit()
        else: combine_fit(parameter)
    if action == "plot_data":
        if parameter == "default": plot_data()
        else: plot_data(parameter)

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

