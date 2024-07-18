import os
import argparse
from Bio import SeqIO
import time
import pandas as pd
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='gene_count')

parser.add_argument('--ref_genome', required=True, help='select reference genome(e.g. GRCh37)')
parser.add_argument('--input_dir', required=True, help='input data directory')
parser.add_argument('--output_dir', required=True, help='output data directory')
parser.add_argument('--threads', default = 10, type = int, help='threads to use in multiprocessing')

args = parser.parse_args()
















# Mutation Count
def changeValue(before, after):
    if (before == 'G' and after == 'C') or (before == 'A' and after == 'T'):
        return after, before
    if before == 'G' and after == 'A':
        return 'C', 'T'
    if before == 'G' and after == 'T':
        return 'C', 'A'
    if before == 'A' and after == 'G':
        return 'T', 'C'
    if before == 'A' and after == 'C':
        return 'T', 'G'
    return before, after

genes = './data/gene_'+str(args.ref_genome)+'.txt'
seqs = list(SeqIO.parse('./data/'+str(args.ref_genome)+'.fa', 'fasta'))
for i in seqs:
    i.id = i.id.replace('chr', '')
chrs = [s.id for s in seqs]

flist = []

input_folder = args.input_dir

for fname in os.listdir(input_folder): 
    if fname[-4:] == '.vcf':
        flist.append(fname)

g = pd.read_csv(genes, sep='\t', header=None)
t1 = pd.read_csv('./data/snvtype01.txt', header=None)
t2 = pd.read_csv('./data/snvtype02.txt', header=None)
types = t1.loc[:, 0] + t2.loc[:, 0].str.split('>').str.get(1)
columns = g.loc[:, 4]

geneList = []

with open(genes, 'r') as f1:
    for line1 in f1.readlines():
        tmp1 = line1.split()
        gene = {
            "startPos": int(tmp1[1]),
            "endPos": int(tmp1[2]),
            "chr": tmp1[3]
        }
        geneList.append(gene)

def calculate(i): # vcf 파일 반복
    start = time.time()
    df = pd.DataFrame(0, index=types, columns=columns)
    fname = flist[i]
    with open(input_folder + '/' + fname, 'r') as f2:
        for line2 in f2.readlines():
            if line2[0] == '#': continue
            tmp2 = line2.split()
            if tmp2[0] not in chrs:
                continue

            Pos = int(tmp2[1])
            loc = int(tmp2[1])
            subseq = str(seqs[chrs.index(tmp2[0])][loc - 2:loc + 1].seq)  # loc이 index인지 순서인지
            before, after = changeValue(subseq[1], tmp2[4])
            subseq2 = subseq[0] + before + subseq[2]
            if subseq2 + tmp2[4] not in df.index: continue

            j = 0
            for gene in geneList:
                colname = columns[j]
                j += 1
                startPos = gene["startPos"]
                endPos = gene["endPos"]
                Chr = gene["chr"]
                if tmp2[0] != Chr or len(tmp2[3]) > 1 or len(tmp2[4]) > 1 or Pos < startPos or Pos > endPos:
                    continue
                if subseq[1] == tmp2[3]:
                    df.loc[subseq2 + after, colname] += 1
                    break
                else:
                    print('Not match! Fname: {}, subseq: {}, vcf: {} to {}'.format(fname, subseq, tmp2[3], tmp2[4]))
    end = time.time()
    
    f = open("./log/gene_cnt_log.txt", "a")
    seconds = end - start
    minutes = int(seconds // 60)
    seconds = int(seconds % 60)
    hours = int(minutes // 60)
    minutes = minutes % 60
    print(fname, file = f)
    print(f"Time: {hours}:{minutes}:{seconds}", file = f)
    f.close()

    df.to_csv(f'{args.output_dir}/Gene_count/{fname[:-4]}_cnt.csv', sep=',', index=True, header=True)


















if __name__ == "__main__":
    pool = Pool(processes = args.threads)
    os.makedirs(f'{args.output_dir}/Gene_count', exist_ok = True)
    result = pool.map(calculate, range(len(flist)))
