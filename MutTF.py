import os
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stat

parser = argparse.ArgumentParser(description='mutational signature analysis program')

parser.add_argument('--ext_dir', required=True, help = 'Directory of signature extraction output')
parser.add_argument('--gsva_file', required=True, help='File name of GSVA output')
parser.add_argument('--count_dir', required=True, help='Directory of gene count output')
parser.add_argument('--tf_file', required=True, help='File name of TF-TG database')
parser.add_argument('--corr_dir', required=True, help = 'Output directory of correlation results')

args = parser.parse_args()
        

ext_P = pd.read_csv(f"./{args.ext_dir}/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt", sep='\t', index_col=0) #Extractor 결과 중 process 사용
ext_E = pd.read_csv(f"./{args.ext_dir}/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt", sep='\t', index_col=0) #Extractor 결과 중 exposure 사용


# Calculated contribution

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


# TF extraction in TF-TG database file

db_TF = pd.read_csv(str(args.tf_file), sep=',', index_col=0) 

db_TF_list = list(set(db_TF['name']))


# Load gene count

count_folder = args.count_dir
count = []
count_name = []

for fName in os.listdir(count_folder):
    if fName[-4:] == '.csv' and fName[-7:-4] == 'cnt':
        count.append(pd.read_csv(count_folder + "/" + fName, sep=',', index_col=0))
        count_name.append(fName[:-8])

count_dict = dict(zip(count_name,count))
count_dict = dict(sorted(count_dict.items()))

smp_name = list(count_dict.keys())


# Save only TF genes

# count_TF = []
# for name in smp_name:
#     sample_sym = count_dict[name][db_TF_list]
#     count_TF.append(sample_sym)

# for i in range(len(smp_name)):
#     count_TF[i].to_csv('./TF_count/'+smp_name[i]+'_TF.csv')

count_TF_dict = dict(zip(smp_name,count_TF))
count_TF_dict = dict(sorted(count_TF_dict.items()))


# Load GSVA results

gsva_d = pd.read_csv(str(args.gsva_file), sep='\t', index_col=0)

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


# Leave genes in GSVA file matching genes in gene count file
gsva_TF = list(count_TF[0].columns)

for i in range(len(gsva_data)):
    gsva_data[i] = gsva_data[i][gsva_data[i].index.isin(gsva_TF)]


# Correlation analysis

column = ['gene', 'sig', 'c', 'p']

index = 0
sig = contri_dict[smp_name[0]].columns
types = count_TF_dict[smp_name[0]].index

pn = ['neg','pos']

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

    df.to_csv(f'./{args.corr_dir}/Result_'+pn[gs]+'.csv') 


# Filtering
    df = df.dropna()

    PartRe_list = []
    sig_df = pd.DataFrame(columns=df.columns)
    cnt = 0
    for s in range(9):  
        cal_df = pd.DataFrame(columns=df.columns)  
        idx = 0
        for i in range(len(df)): 
            if df.iloc[i, 2] == '' or df.iloc[i, 3] == '': continue
            # Num = int(df.iloc[i,0])
            Sig = int(df.iloc[i, 1])
            Cor = float(df.iloc[i, 2])
            P_val = float(df.iloc[i, 3])
            if Sig == s and Cor >= 0.35 and P_val <= 0.05: 
                print(Sig, Cor, P_val) 
                cal_df.loc[idx, :] = df.iloc[i, :]
                sig_df.loc[cnt, :] = df.iloc[i, :]
                idx += 1
                cnt += 1
        cal_df.drop(cal_df.columns[1], axis=1)   

        PartRe_list.append(cal_df)
    sig_df.drop(cal_df.columns[1], axis=1)
    sig_df.to_csv(f'./{args.corr_dir}/Filt_Result_'+pn[gs]+'.csv', sep=',', index=True, header=True) 

