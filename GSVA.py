import os
import pandas as pd
import numpy as np
import scipy.stats as stat
from tqdm.notebook import tqdm
import operator
import argparse

parser = argparse.ArgumentParser(description="Performing GSVA based on TF-TG gene set",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-g", "--geneset", help="Path to original geneset file (.csv)")
parser.add_argument("-e", "--expression", help="Path to expression file (.tsv)")
parser.add_argument("-o", "--gsva", help="Path to GSVA output file (.tsv)")
args = vars(parser.parse_args())

##### Load input files #####
norm_count_sym = pd.read_csv(args['expression'], sep='\t', index_col = 0)
colon_TF_file = pd.read_csv(args['geneset'], sep=',', index_col = 0)

TF_name_list = list(set(colon_TF_file['name']))


##### Calculate correlation between TF and TG expression #####

print('Calcualte correlation between TF and TG expression')

colon_TF_TG_cor = pd.DataFrame(columns=['TF','target','tissue','r','p'])

for i in tqdm(range(len(TF_name_list))):
    tg = colon_TF_file[colon_TF_file['name']==TF_name_list[i]]

    tg_sym = list(tg['description'])
    
    tf_tg_match = norm_count_sym[norm_count_sym.index.isin(tg_sym)]
    tg_com = tf_tg_match
    tf_com = norm_count_sym[norm_count_sym.index==TF_name_list[i]]
    if len(tf_com) == 0 or len(tg_com) == 0:
        continue
    
    log_= list(tf_com.iloc[0])

    for j in range(len(tg_com)):
        r = np.nan
        p_value = np.nan

        v1 = list(tf_tg_match.iloc[j])
        
        if np.isfinite(log_[0]):
            r, p_value = stat.pearsonr(v1, log_)
        
        if r < 0:
            TF_name = TF_name_list[i]+'_0'
        else:
            TF_name = TF_name_list[i]+'_1'

        colon_TF_TG_cor.loc[len(colon_TF_TG_cor)] = [TF_name, tf_tg_match.index[j], 'colon',r,p_value]


##### Filter based on p.value #####

colon_TF_TG_corp = colon_TF_TG_cor[colon_TF_TG_cor['p']<=0.05]

colon_TF_TG_corp.to_csv(f"{args['geneset'][:-4]}_sep.txt", sep='\t', index=True, header=True)


##### Perform GSVA #####

print('Perform GSVA')

os.system(f"Rscript GSVA.R -g {args['geneset'][:-4]}_sep.txt -e {args['expression']} -o {args['gsva']}")
