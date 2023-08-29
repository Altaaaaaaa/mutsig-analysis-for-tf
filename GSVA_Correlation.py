import os
import pandas as pd
import numpy as np
import scipy.stats as stat
from numpy import dot
from numpy.linalg import norm
from sklearn.metrics.pairwise import cosine_similarity
from tqdm.notebook import tqdm
import argparse

parser = argparse.ArgumentParser(description="Performing GSVA and Calculate correlation between GSVA and mutation count",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-g", "--geneset", help="Path to geneset file")
parser.add_argument("-e", "--expression", help="Path to expression file")
parser.add_argument("-o1", "--gsva", help="Path to GSVA output file (.tsv)")
parser.add_argument("-m", "--mutation", help="Path to directory with mutation count files (.csv)")
parser.add_argument("-o2", "--correlation", help="Prefix of correlation output file (.csv)")
args = vars(parser.parse_args())

os.system(f"Rscript GSVA.R -g {args['geneset']} -e {args['expression']} -o {args['gsva']}")

##### Load GSVA data #####

gsva_d = pd.read_csv(args['gsva'], sep='\t', index_col=0) # gsva score

gsva_data = []
gsva_gene = []

df1 = pd.DataFrame(columns=gsva_d.columns)
df0 = pd.DataFrame(columns=gsva_d.columns)

for i in range(len(gsva_d)):
    if gsva_d.index[i][-2:] == '_1':
        df1.loc[gsva_d.index[i][:-2]] = gsva_d.loc[gsva_d.index[i]]
    if gsva_d.index[i][-2:] == '_0':
        df0.loc[gsva_d.index[i][:-2]] = gsva_d.loc[gsva_d.index[i]]
    gsva_gene.append(gsva_d.index[i][:-2])
        
gsva_data.append(df0)
gsva_data.append(df1)


##### Load mutation count data #####

pcawg_colon_folder = args['mutation'] + '/'
pcawg_colon = []
pcawg_colon_sample_name = []
pcawg_colon_C = []

for fName in tqdm(gsva_d.columns.tolist()):
    pcawg_colon.append(pd.read_csv(pcawg_colon_folder + fName + "_TF.csv", sep=',', index_col=0)) # 입력받는 파일
    pcawg_colon_C.append(pd.read_csv(pcawg_colon_folder + fName + "_C.csv", sep=',', index_col=0)) # 입력받는 파일
    pcawg_colon_sample_name.append(fName)

# Check the orders

colon_TF_symbol = list(pcawg_colon[0].columns)

for i in range(len(gsva_data)):
    gsva_data[i] = gsva_data[i][gsva_data[i].index.isin(colon_TF_symbol)]


##### Calculate correlation #####

column = ['gene', 'sig', 'c', 'p']

index = 0
sig = pcawg_colon_C[0].columns
types = pcawg_colon[0].index

cnt = 0

for gs in tqdm(range(len(gsva_data))):
    df = pd.DataFrame(columns=column)

    genes = gsva_data[gs].index[:]
    gsva_cal = gsva_data[gs]
    for gene in tqdm(genes):  
        log_ = []
        print(f'gene:{gene}_{gs}')
        for s in sig:
            df.loc[index, 'gene'] = gene
            df.loc[index, 'sig'] = str(s)
            v1 = []
            log_ = []

            for i in range(len(pcawg_colon)):
                sum = 0
                type_i = 0
                for type in types:
                    sum += pcawg_colon[i].loc[type, gene] * pcawg_colon_C[i].iloc[type_i, int(s)]
                    type_i += 1
                v1.append(sum)
            for i in pcawg_colon_sample_name:
                log_.append(gsva_cal[i][gene])

            if np.sum(v1) == 0:
                cnt+=1
            
                if cnt%10 == 0:
                    print(cnt)
            
            # print(f'v1:{v1}')
            # print(f'log:{log_}')

            r = np.nan
            p_value = np.nan
            if np.isfinite(log_[0]):
                r, p_value = stat.spearmanr(v1, log_) 

            df.loc[index, 'c'] = r
            df.loc[index, 'p'] = p_value
            index += 1

    df.to_csv(args['correlation']+str(gs)+'.csv')  # correlation 분석 결과 correlation -> gsva 결과 