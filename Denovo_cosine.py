import pandas as pd
import argparse
import matplotlib.pyplot as plt
from numpy import dot
from numpy.linalg import norm
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='de novo-cosmic cosine similarity')

parser.add_argument('--in_dir', required=True, help='Directory of vcf files')
parser.add_argument('--ext_dir', required=True, help='Directory of signature extraction results')
parser.add_argument('--ref_genome', required=True, help='select reference genome(e.g. GRCh37)')
parser.add_argument('--version', required=True, help='select cosmic signature version(e.g. 3.3)')

args = parser.parse_args()











def cos(P,COSMIC_P):
    def cos_sim(A, B):
        return dot(A, B)/(norm(A)*norm(B))

    cossim = pd.DataFrame(columns=COSMIC_P.columns, index=P.columns)

    for i in range(len(P.columns)):
        sim = []
        cosmax = 0
        p_vector = P[P.columns[i]].values # p_vector에 P행렬의 값 벡터로 변환
        for j in range(len(COSMIC_P.columns)):
            cosp_vector = COSMIC_P[COSMIC_P.columns[j]].values # cosp_vector에 COSMIC_P33행렬의 값 벡터로 변환
            simsim = cos_sim(p_vector, cosp_vector)
            sim.append(cos_sim(p_vector, cosp_vector)) # 코사인 유사도 계산한 값 sim 리스트에 추가
            cossim.iloc[i,j] = float(cos_sim(p_vector, cosp_vector))
            #cossim.iloc[i,j] = np.array((cos_sim(p_vector, cosp_vector)), dtype=float)
            if cos_sim(p_vector, cosp_vector) == max(sim): # 만약 계산한 코사인 유사도 값이 sim 리스트에 있는 가장 큰 값이면 j 값 저장
                cosmax = j
        if max(sim) >= 0: # 0.85이상으로 고치기     # sim 리스트의 가장 큰 값 출력
            print(P.columns[i],'is similar to', COSMIC_P.columns[cosmax], ' The similarity is', max(sim))

    # Draw heatmap
    cossim2 = cossim.astype(float) # datatype 변환

    plt.figure(figsize = (30, len(P.columns))) # (# of cosmic signature, # of de novo signature)
    sns.heatmap(cossim2, cmap='coolwarm', 
                    annot=True,
                    fmt=".3f",
                    annot_kws={'size':10},
                    cbar=True,
                    square=True)
    plt.title('Cosine Similarity', fontsize=20)

    cossim2 = cossim.astype(float)

    plt.show()
    plt.savefig('Heatmap.png',dpi=300) 









COSMIC_P = pd.read_csv(f"{args.in_dir}/COSMIC_v"+str(args.version)+'_SBS_'+str(args.ref_genome)+'.txt', sep='\t', index_col=0) # COSMIC의 Signature
P1 = pd.read_csv(f"{args.ext_dir}/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt", sep='\t', index_col=0).sort_values(by='MutationType') #Extractor 결과 중 process 사용





cos(P1,COSMIC_P)
