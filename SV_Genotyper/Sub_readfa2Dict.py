import gzip
from collections import deque
from gzip import BadGzipFile
import re
def readfa2Dict(fa):
    bigFa = {}
    geneID = ''
    geneSeq = deque()
    with gzip.open(fa,'rb') as f_in:
        try:
            f_in.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzip.open(fa, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = line.strip().split(b">")[1].split(b' ')[0].decode()
                            geneseq = deque()
                        else:
                            geneID = line.strip().split(b'>')[1].split(b' ')[0].decode()
                    else:
                        geneSeq.append(line.strip().decode())
        else:
            with open(fa, 'r') as fin:
                for line in fin:
                    if ">" in line:
                        if geneID != '':
                            bigFa[geneID] = ''.join(geneSeq)
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                            geneSeq = deque()
                        else:
                            geneID = re.split('\s+', line.strip().split('>')[1])[0]
                    else:
                        geneSeq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    ####### the last line cant be iterated, so we should one more code to store it into dict ###########
    if geneID != '':
        bigFa[geneID] = ''.join(geneSeq)
    return bigFa
