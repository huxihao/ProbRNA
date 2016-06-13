## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

import os,sys,time,string

## define log file name
LOG_FILE = "log.txt"

def warning(info):
    print "Warning:", info

def show(info):
    out = str(info)
    outfile = open(LOG_FILE,"a")
    ## output a table like file
    if str(out).endswith("\n"):
        outfile.write(str(out))
        print out[0:-1]
    else:
        outfile.write("%s\t"%out)
        print out,
    outfile.close()

def write_fasta(filename, data):
    outfile = open(filename,'w')
    for fid,seq in data:
        outfile.write(">%s\n%s\n"%(fid,seq))
    outfile.close()

def read_fasta(filename):
    data = []
    infile = open(filename,'r')
    fid = ""
    seq = ""
    for line in infile:
        if line[0] == '>':
            if len(seq) > 0:
                data.append((fid,seq))
                fid = ""
                seq = ""
            fid = line[1:-1] ## remove '>' and '\n'
        else:
            seq += line[0:-1]
    if len(seq) > 0:
        data.append((fid,seq))
    infile.close()
    return data

def norm_line(a, st=0, ed=1):
    a_max = float(max(a))
    a_min = float(min(a))
    r = []
    for e in a:
        if a_max > a_min:
            r.append(st+(ed-st)*(e-a_min)/(a_max-a_min))
        else:
            r.append((st+ed)*0.5)
    return r

def norm_rank(a, use_percent=True):
    pair = []        
    range_len = len(a)
    for i in range(range_len):
        pair.append([a[i], i])
    pair.sort(lambda x,y: cmp(x[0],y[0]))
    r = range(range_len) ## initialization
    for i in range(range_len):
        if use_percent:
            r[pair[i][1]] = (i+1)/float(range_len)
        else:
            r[pair[i][1]] = i
    return r

def mean_std(score):
    mean, std = 0,0
    range_len = len(score)
    for k in range(range_len):
        mean += score[k]
    mean /= float(range_len)
    for k in range(range_len):
        std += (score[k]-mean)**2
    if range_len == 1:
        std = 0
    else:
        std /= float(range_len - 1)
    std = std**0.5
    return mean, std

def roman2arabic(roman):
    d={"I":1,"V":5,"X":10,"L":50,"C":100,"D":500,"M":1000}
    l=map(lambda x:d[x],list(roman))
    for i in range(len(l)-1):
      if l[i] < l[i+1]: 
          l[i] = -l[i]
    return sum(l)

def arabic2roman(arabic):
    ones=["","I","II","III","IV","V","VI","VII","VIII","IX"]
    tens=["","X","XX","XXX","XL","L","LX","LXX","LXXX","XC"]
    hund=["","C","CC","CCC","CD","D","DC","DCC","DCCC","CM"]
    thou=["","M","MM","MMM"]
    N=str(format(arabic, '04n'))
    I=ones[int(N[-1])]
    X=tens[int(N[-2])]
    C=hund[int(N[-3])]
    M=thou[int(N[-4])]
    return M+C+X+I

def correlation(x,y):
    from scipy.stats import pearsonr, spearmanr
    return (pearsonr(x,y)[0], spearmanr(x,y)[0])

