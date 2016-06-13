''' 
Author: Xihao Hu <huxihao@gmail.com>
Example: create a new file and type

from tools import *

def main(para):
    print para

if __name__ == '__main__': main_fun(main)

'''

import os,sys,time,string

## Global Log File Name
DEFINE_LOG_FILE = "log.txt"

def show(info="\n", END=False):
    ''' Output information to be in a table-like form.
        Both print on screen and save into the file of DEFINE_LOG_FILE.
    '''
    out = str(info)
    outfile = open(DEFINE_LOG_FILE,"a")
    if isinstance(info, int):
        outfile.write("%s\t"%out)
        print string.ljust(out,5),
    elif isinstance(info, float):
        outfile.write("%s\t"%out)
        print "%.3f"%info,
    elif isinstance(info, list) or isinstance(info, set) or isinstance(info, tuple):
        for e in info:
            show(e)
    elif out.endswith("\n"):
        outfile.write(out)
        print out[:min(len(out)-1,78)]
    else:
        outfile.write(str(out)+"\t")
        print string.ljust(out[:min(len(out),25)],25),
    if END:
        print ''
        outfile.write("\n")
    outfile.close()
    return ''

def main_fun(main):
    ''' Running command: python template.py WorkPath=C:/work_path/
        Read parameters and set WorkPath to default '''
    para = {"ExeFile": os.path.basename(sys.argv[0]), 
            "LogFile": "log." + os.path.basename(sys.argv[0])+".txt", 
            "WorkPath":"./"}
    for arg in sys.argv[1:]:
        name, value = arg.split("=")
        para[name] = value
    ## Set paths
    os.chdir(para["WorkPath"])
    para["WorkPath"] = os.path.abspath(".")
    global DEFINE_LOG_FILE
    DEFINE_LOG_FILE = os.path.abspath(para["LogFile"])
    ## Delete old log file
    if os.path.exists(DEFINE_LOG_FILE):
        os.remove(DEFINE_LOG_FILE)
    START_TIME = time.time()
    main(para)
    END_TIME = time.time()
    show()
    for name in para:
        show("# Parameter: " + name + " =\t" + para[name] + "\n")
    show(time.strftime("# Begin  at: %Y/%m/%d %H:%M:%S\n", time.gmtime(START_TIME)))
    show(time.strftime("# Finish at: %Y/%m/%d %H:%M:%S\n", time.gmtime(END_TIME)))
    show(time.strftime("# Total Time: %d-th day %H:%M:%S\n", time.gmtime(END_TIME-START_TIME))) 


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
    for k in xrange(range_len):
        mean += score[k]
    mean /= float(range_len)
    for k in xrange(range_len):
        std += (score[k]-mean)**2
    if range_len == 1:
        std = 0
    else:
        std /= float(range_len - 1)
    std = std**0.5
    return [mean, std]

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

def randomword(length):
    import random
    return ''.join(random.choice(string.uppercase) for i in range(length))

def str2num(s):
    return ''.join(map(lambda x : '%.3d' % ord(x), s))

def num2str(s):
    return ''.join([chr(int(s[i:i+3])) for i in range(0, len(s), 3)])

def correlation(x,y,rank=False):
    from scipy.stats import pearsonr, spearmanr
    if rank: return spearmanr(x,y)[0]
    else: return pearsonr(x,y)[0]

def histogram(val, bins, prob=True):
    ''' bins = [a, b, c] --> [a,b), [b,c), others '''
    counts = [0 for i in bins]
    for v in val:
        find = False
        for j in xrange(len(bins)-1):
            if bins[j] <= v and v < bins[j+1]:
                counts[j] += 1
                find = True
                break
        if not find:
            counts[-1] += 1
    if prob:
        return [c/float(len(val)) for c in counts]
    return counts

def performance_slow(real, pred, x='FPR', y='TPR', eps=1e-50, trim=4):
    ''' Measure the performance by returning the area under the curve
        Ref: http://en.wikipedia.org/wiki/Receiver_operating_characteristic
        * 'eps' was used to avoid devided-by-zero error
        * 'trim' was used to merge close points on the curve
    '''
    from numpy import trapz
    curve = set()
    point = set()
    for cutoff in (list(set(pred))+[float('inf')]):
        TP, TN, FP, FN = 0, 0, 0, 0
        for r, p in zip(real, pred):
            if r:
                if p >= cutoff: TP += 1
                else: FN += 1
            else:
                if p >= cutoff: FP += 1
                else: TN += 1
        outcome = []
        for measure in [x, y]:
            if measure == 'TPR':
                outcome.append(round(TP/(TP+FN+eps),trim))
            elif measure == 'SPC':
                outcome.append(round(TN/(FP+TN+eps),trim))
            elif measure == 'PPV':
                outcome.append(round((TP+eps)/(TP+FP+eps),trim))
            elif measure == 'NPV':
                outcome.append(round(TN/(TN+FN+eps),trim))
            elif measure == 'FPR':
                outcome.append(round(FP/(FP+TN+eps),trim))
            elif measure == 'FDR':
                outcome.append(round(FP/(FP+TP+eps),trim))
            elif measure == 'FNR':
                outcome.append(round(FN/(FN+TP+eps),trim))
            else:
                raise ValueError('Unknown measurement %s'%measure)
        if tuple(outcome) in point: ## remove overlapped points
            continue
        point.add(tuple(outcome))
        outcome.append(cutoff) ## add cutoff
        curve.add(tuple(outcome))
    curve = list(curve)
    curve.sort()
    curve_x, curve_y, curve_c = zip(*curve)
    curve_area = trapz(curve_y, curve_x)
    return round(curve_area,trim), curve_x, curve_y, curve_c

def performance(real, pred, x='FPR', y='TPR', eps=1e-50, trim=4):
    ''' Measure the performance by returning the area under the curve
        Ref: http://en.wikipedia.org/wiki/Receiver_operating_characteristic
        * 'eps' was used to avoid devided-by-zero error
        * 'trim' was used to merge close points on the curve
    '''
    from numpy import trapz
    pairs = zip(pred, real)
    pairs.sort()
    curve = set()
    point = set()
    ## predict all to be positive
    TP, TN, FP, FN = real.count(True), 0, real.count(False), 0
    i = -1; j = -1;
    all_cutoff = sorted(list(set(pred)))
    for cutoff in [min(all_cutoff)-1] + all_cutoff: ## having less positive
        neg = 0; pos = 0
        while j >= 0 and  j < len(pairs): ## end normally
            p, r = pairs[j]
            if i == j:
                assert p == cutoff ## p[i] must be equal to the cutoff
            if r: pos += 1
            else: neg += 1
            if p != cutoff: ## roll back
                if r: pos -= 1
                else: neg -= 1
                j -= 1
                break
            j += 1
        ## points within [i,j] are now predicted to be negative
        TP -= pos; TN += neg; FP -= neg; FN += pos
        ## update i,j
        j += 1; i = j
        outcome = []
        for measure in [x, y]:
            if measure == 'TPR':
                outcome.append(round(TP/(TP+FN+eps),trim))
            elif measure == 'SPC':
                outcome.append(round(TN/(FP+TN+eps),trim))
            elif measure == 'PPV':
                outcome.append(round((TP+eps)/(TP+FP+eps),trim))
            elif measure == 'NPV':
                outcome.append(round(TN/(TN+FN+eps),trim))
            elif measure == 'FPR':
                outcome.append(round(FP/(FP+TN+eps),trim))
            elif measure == 'FDR':
                outcome.append(round(FP/(FP+TP+eps),trim))
            elif measure == 'FNR':
                outcome.append(round(FN/(FN+TP+eps),trim))
            else:
                raise ValueError('Unknown measurement %s'%measure)
        if tuple(outcome) in point: ## remove overlapped points
            continue
        point.add(tuple(outcome))
        outcome.append(cutoff) ## add cutoff
        curve.add(tuple(outcome))
    curve = list(curve)
    curve.sort()
    curve_x = []; curve_y = []; curve_c = []
    x1 = -1; y1 = -1
    for x,y,c in curve: ## remove points within the lines
        if x == x1 and y >= y1:
            if len(curve_x) >= 2 and curve_x[-2] == x:
                curve_y[-1] = y
                curve_c[-1] = c
                y1 = y
                continue
        elif y == y1 and x >= x1:
            if len(curve_y) >= 2 and curve_y[-2] == y:
                curve_x[-1] = x
                curve_c[-1] = c
                x1 = x
                continue
        curve_x.append(x)
        curve_y.append(y)
        curve_c.append(c)
        x1 = x; y1 = y
    curve_area = trapz(curve_y, curve_x)
    return round(curve_area,trim), curve_x, curve_y, curve_c

def main(para):
    d = []
    infile = open('gPAR-CLIP_crosslinking_sites.txt','r')
    infile.readline()
    for line in infile:
        ele = line.split('\t')
        d.append(abs(float(ele[3])-float(ele[4])))
    infile.close()
    show(min(d))
    show(max(d))
    show()
    bins = range(int(min(d)), int(max(d))+1)
    cc = histogram(d, bins)
    show(bins)
    show()
    show(cc)
    show()

if __name__ == "__main__": main_fun(main)
