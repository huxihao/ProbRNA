## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

## Detacting Mapping problems in PARS dataset
## Author: Xihao Hu <huxihao@gmail.com>
## Version 1.0 - March 24, 2011

from tools import *

def map2genome(codingfile = "../data/PARS/sce_transcriptome_global.tab"):
    infile = open(codingfile)
    genome = {}
    last_gene = ""
    last_ele = []
    for line in infile:
        chro, gene, start, end, mark = line.strip().split('\t')
        if mark not in ["5UTR", "Exon", "3UTR", "Transcript"]:
            continue
        if gene in genome:
            gene_map = genome[gene]
            if chro != gene_map[0]:
                warning("% is marked on chromasome %s and %s"%(gene,\
                    chro, gene_map[0]))
            if gene == last_gene: ## speed up for searching
                last_ele.append((mark, int(start), int(end)))
            else: ## because of overlapped genes
                genome[gene].append((mark, int(start), int(end)))
        else:
            last_ele = genome[gene] = [chro, (mark, int(start), int(end))]
            last_gene = gene
    #show("# Read in %s mapped genes\n"%len(genome))
    ## format of gemome:  save in a dictionary
    ## key -> [chromasome number, (mark, start, end), (mark, start, end), ...]
    return genome

def get3length(genome_mark):
    lens = [0,0,0] ## [5UTR, CDS, 3UTR]
    for g_mark, g_start, g_end in genome_mark[1:]:
        g_len = g_end - g_start
        if g_len >= 0: g_len = 1 + g_len
        else: g_len = 1 - g_len
        if g_mark == "5UTR": lens[0] += g_len
        elif g_mark == "Exon": lens[1] += g_len
        elif g_mark == "3UTR": lens[2] += g_len
    return lens

def get_splice(genome_mark):
    sites = []
    isFirst = True
    pos = 0
    for g_mark, g_start, g_end in genome_mark[1:]:
        g_len = g_end - g_start
        if g_len >= 0: g_len = 1 + g_len
        else: g_len = 1 - g_len
        if g_mark == "5UTR": 
            pos += g_len
        elif g_mark == "Exon": 
            if isFirst: isFirst = False
            else: sites.append(pos)
            pos += g_len
        elif g_mark == "3UTR": 
            pass
    return sites

def global2local(genome_mark, index, hasIntron = False):
    ## global index is count from 1
    ## local index is count from 0 to length-1
    left_count, right_count = 0, 0
    isIntron = True
    for g_mark, g_start, g_end in genome_mark[1:]:
        if (not hasIntron) and g_mark == "Transcript":
            continue
        elif (hasIntron) and g_mark != "Transcript":
            continue
        if g_start < g_end: ## positive strand
            ## >>>>g_start>>>index>>>g_end>>>[g_start,g_end]>>>
            if index < g_start:
                right_count += g_end - g_start + 1
            elif index > g_end:
                left_count += g_end - g_start + 1
            else:
                isIntron = False
                left_count += index - g_start
                right_count += g_end - index
        else: ## nagative strand
            ## <<<[g_end,g_start]<<<g_end<<<index<<<g_start<<<<
            if index > g_start:
                right_count += g_start - g_end + 1
            elif index < g_end:
                left_count += g_start - g_end + 1
            else:
                isIntron = False
                left_count += g_start - index
                right_count += index - g_end
        #print g_mark, g_start, g_end ,index,left_count ,right_count
    if (not hasIntron) and isIntron:
        return -1 ## error
    return left_count

def main():
    from ParsScores import import_pars, pars_gene_groups
    pars = import_pars()
    gene_list, groups = pars_gene_groups(pars)
    genome = map2genome()

    if True: ## 3' term counts and pvalues for V1 and S1
        ## get global distribution
        v1cc = [0 for i in range(50)]
        s1cc = [0 for i in range(50)]
        cc = 0
        for gene in gene_list:
            _seq = pars[gene]["FOLD_SEQ"]
            _v1 = pars[gene]["V1"]
            _s1 = pars[gene]["S1"]
            for i in range(50):
                v1cc[i] += int(_v1[len(_seq)-i-1])
                s1cc[i] += int(_s1[len(_seq)-i-1])
            splice_sites = get_splice(genome[gene])
            if splice_sites != []:
                show(gene)
                for site in splice_sites:
                    show(site)
                    cc += 1
                show("\n")
        show("Here are %d splicing sites\n"%cc)
        for i in range(50):
            show('%s\t%s\t%s\n'%(i+1, v1cc[i], s1cc[i]))

    if False: ## unique reads mapping
        ReadLen = 35
        uniqueSub = {}
        for gene in gene_list:
            gene_map = genome[gene]
            rna_type = pars[gene]["TYPE"]
            _seq = pars[gene]["FOLD_SEQ"]
            _v1 = pars[gene]["V1"]
            _s1 = pars[gene]["S1"]
            ## mapping
            coding_length = 0
            chro = gene_map[0]
            for mark, i_st, i_ed in gene_map[1:]:
                if i_st > i_ed:
                    coding_length += i_st - i_ed + 1
                else:
                    coding_length += i_ed - i_st + 1
            if True and len(_seq) != coding_length: ## check
                show(gene)
                show(rna_type)
                show(len(_seq))
                show(coding_length)
                show("\n")

            for i in range(1, len(_seq)-ReadLen+1): ## shift
                substring = _seq[i:i+ReadLen]
                uniqueSub[substring] = 1+uniqueSub.setdefault(substring, 0)
        overlapped = {}
        for gene in gene_list:
            gene_map = genome[gene]
            _seq = pars[gene]["FOLD_SEQ"]
            _v1 = pars[gene]["V1"]
            _s1 = pars[gene]["S1"]
            for i in range(1, len(_seq)-ReadLen+1): ## shift
                substring = _seq[i:i+ReadLen]
                if uniqueSub[substring] > 1:
                    chro = gene_map[0]
                    posi = i
                    for mark, i_st, i_ed in gene_map[1:]:
                        if i_st < i_ed:
                            regin_len = i_ed - i_st + 1
                            if posi < regin_len:
                                posi += i_st
                                break
                            else:
                                posi -= regin_len
                        else: ## reverse
                            posi = len(_seq) - ReadLen - posi
                            regin_len = i_st - i_ed + 1
                            if posi < regin_len:
                                posi += i_ed
                                break
                            else:
                                posi -= regin_len
                                posi = len(_seq) - ReadLen - posi
                    if False: ## swich
                        show(substring)
                        show(gene)
                        show(i+1)
                        show(_v1[i-1])
                        show(_s1[i-1])
                        show("Chro %s"%chro)
                        show(posi)
                        show('\n')
                    overlapped.setdefault(substring, [])
                    overlapped[substring].append((gene,i,chro,posi,\
                            int(_v1[i-1]),int(_s1[i-1]))) ## shift
        ### print overlapped subsequences
        if True: ## swich
            show("# Contain %s overlapped %s-mers\n"%(len(overlapped), ReadLen))
            for substring in overlapped:
                matched = overlapped[substring]
                #if list(matched[0])[4]==list(matched[1])[4] \
                #        and list(matched[0])[5]==list(matched[1])[5]:
                #    continue
                if list(matched[0])[2]==list(matched[1])[2] \
                        and list(matched[0])[3]==list(matched[1])[3]:
                    continue
                show(substring)
                for gene,i,chro,p,v,s in matched:
                    show("%s\t%s\t%s\t%s\t%s %s"%(gene,i,chro,p,v,s))
                show("\n")

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

