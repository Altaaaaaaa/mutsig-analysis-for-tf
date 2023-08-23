import os
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stat
from numpy import dot
from numpy.linalg import norm
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns
from pdf2image import convert_from_path

# 인자값을 받을 수 있는 인스턴스 생성
parser = argparse.ArgumentParser(description='mutational signature analysis program')

# 입력받을 인자값 등록
parser.add_argument('--gsva_folder', required=True, help='gsva_result_folder')
parser.add_argument('--tf_file', required=True, help='tf_database_file')

# 입력받은 인자값을 args에 저장 (type: namespace)
args = parser.parse_args()

# Gene count 불러오기
yMatList = []

for fName in os.listdir("./Gene_count"):
    if fName[-4:] == '.csv':
        yMatList.append(pd.read_csv("./Gene_count/" + fName, index_col=0)) # [:11]
        

ext_P = pd.read_csv("./ext_data/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt", sep='\t', index_col=0) #Extractor 결과 중 process 사용
ext_E = pd.read_csv("./ext_data/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt", sep='\t', index_col=0) #Extractor 결과 중 exposure 사용

# Contribution 계산

def cal_contri(file_name,p_mat,e_mat):
    #E = E.T

    samples = e_mat.index

    sig = p_mat.columns
    types = p_mat.index

    # total -> 분자 합

    contri_list = []
    contri_name_list = []

    for s in range(len(samples)):
        df = pd.DataFrame()
        for j in range(len(types)):
            numerator = []
            sum = 0.0
            for i in range(len(sig)):
                p = p_mat.iloc[j, i]
                e = e_mat.iloc[s, i]
                result = p * e
                numerator.append(result)
                sum += result
            for ss in range(len(sig)):
                if sum == 0.0:
                    df.loc[j, ss] = 0
                else:
                    df.loc[j, ss] = numerator[ss] / sum
            #print(samples[s])
        df.to_csv(file_name+'C/' + samples[s] + '_C.csv') # contribution 결과 파일
        contri_list.append(df)
        contri_name_list.append(samples[s])

    return contri_list, contri_name_list

cList, contri_name = cal_contri('./output/', ext_P, ext_E)

contri_dict = dict(zip(contri_name,cList))
contri_dict = dict(sorted(contri_dict.items()))

# TF 추출 -> 저장

db_TF = pd.read_csv(args.tf_file, sep=',', index_col=0) # TF 데이터

db_TF_list = list(set(db_TF['name']))

with open('./output/tf/tf_name_list.txt', 'a') as file:
    file.write('\n'.join(db_TF_list))

# Gene count 불러오기

count_folder = "./Gene_count"
count = []
count_name = []

for fName in os.listdir(count_folder):
    if fName[-4:] == '.csv' and fName[-7:-4] == 'cnt':
        count.append(pd.read_csv(count_folder + "/" + fName, sep=',', index_col=0)) # 입력받는 파일
        count_name.append(fName[:-8])

count_dict = dict(zip(count_name,count))
count_dict = dict(sorted(count_dict.items()))

smp_name = list(count_dict.keys())

# Gene id가 TF인 것만 추출하여 저장

count_TF = []
for name in smp_name:
    sample_sym = count_dict[name][db_TF_list]
    count_TF.append(sample_sym)

for i in range(len(smp_name)):
    count_TF[i].to_csv('./TF_count/'+smp_name[i]+'_TF.csv')

count_TF_dict = dict(zip(smp_name,count_TF))
count_TF_dict = dict(sorted(count_TF_dict.items()))

# gsva 결과 파일 가져오기

gsva_d = pd.read_csv(args.gsva_folder, sep='\t', index_col=0) # gsva score

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

# count 파일의 gene id와 일치하는 것만 남김
gsva_TF = list(count_TF[0].columns)

for i in range(len(gsva_data)):
    gsva_data[i] = gsva_data[i][gsva_data[i].index.isin(gsva_TF)]


# correlation 분석

column = ['gene', 'sig', 'c', 'p']

index = 0
sig = contri_dict[smp_name[0]].columns
types = count_TF_dict[smp_name[0]].index

for gs in range(len(gsva_data)):
    df = pd.DataFrame(columns=column)

    genes = gsva_data[gs].index[:]
    gsva_cal = gsva_data[gs]
    for gene in genes:  
        log_ = []
        for s in sig:
            df.loc[index, 'gene'] = gene
            df.loc[index, 'sig'] = str(s)
            v1 = []
            log_ = []

            for i in range(len(smp_name)):
                sum = 0
                type_i = 0
                for type in types:
                    sum += count_TF_dict[smp_name].loc[type, gene] * contri_dict[smp_name].iloc[type_i, int(s)]
                    type_i += 1
                v1.append(sum)
            for i in smp_name:
                log_.append(gsva_cal[i][gene])

            r = np.nan
            p_value = np.nan
            if np.isfinite(log_[0]):
                r, p_value = stat.pearsonr(v1, log_) 

            df.loc[index, 'c'] = r
            df.loc[index, 'p'] = p_value
            index += 1

    df.to_csv('./output/Cor/Result_'+str(gs)+'.csv')  # correlation 분석 결과