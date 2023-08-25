import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# 인자값을 받을 수 있는 인스턴스 생성
parser = argparse.ArgumentParser(description='de novo-cosmic cosine similarity')

# 입력받을 인자값 등록
parser.add_argument('--pos_neg', required=True, help='select pos/neg group(pos or neg)')
parser.add_argument('--tf_group_file', required=True, help='tf-tg data divided into two groups')
parser.add_argument('--sig_num', required=True, help='enter the number of signature')


# 입력받은 인자값을 args에 저장 (type: namespace)
args = parser.parse_args()











cor_output = pd.read_csv("./output/Cor/Filt_Result_"+str(args.pos_neg)+".csv", sep=',', index_col=0)
TF_list=list(set(cor_output.gene)) # TF의 symbol id 리스트로 만들어주기(중복제거)


##########바꿔주기
# TF데이터(GSVA 전, group별로 분리)
gsva_cor = pd.read_csv(str(args.tf_file), sep=',', index_col=0)











def make_edge(pn):
    tg_group = gsva_cor[gsva_cor['TF'].str.contains('_'+str(pn), na = False)]
    tg_group['TF'] = [x.split('_')[0] for x in list(tg_group['TF'])]

    tg_group = tg_group.loc[:,['TF','target','r']]
    tg_group.columns = ['protein1','protein2','combined_score']
    
    if pn == 0:
        tg_group['combined_score'] = tg_group['combined_score'].transform(lambda x : -x*1000)
    else:
        tg_group['combined_score'] = tg_group['combined_score'].transform(lambda x : x*1000)
    
    tflist = list(set(cor_output['gene']))
    tftg_net = tg_group[tg_group['protein1'].isin(tflist)] # protein1이 있던 부분이 TF임

    tglist = list(set(tftg_net['protein2'])-set(tflist))

    tftg_net.to_csv('./output/Node_classi/edge_'+str(args.pos_neg)+".csv", sep=',', index=True, header=True)

    return tflist, tglist, tftg_net











if str(args.pos_neg) == 'neg':
    TF_list, tglist, link_filt = make_edge(0)
else:
    TF_list, tglist, link_filt = make_edge(1)








# signature를 0,1로 표현한 데이터프레임 만들기
gene_sig_df = pd.DataFrame(index=TF_list,columns=['sig.'+str(i) for i in range(int(args.sig_num))])
gene_sig_df = gene_sig_df.fillna(0)

for i in range(len(cor_output)):
    gene_name = cor_output.gene[i]
    gene_sig_n = 'sig.'+str(cor_output.sig[i])
    gene_sig_df.loc[gene_name,gene_sig_n] = 1

for tgname in tglist:
    gene_sig_df.loc[tgname] = [0,0,0,0,0,0,0,0,0]

node_gene_sig = []

for i in range(len(gene_sig_df)):
    node_gene_sig.append(gene_sig_df.iloc[i].to_list())















def node_iter(n,node_info,mul_gene):

    #new_count = []

    node_gene_sig = node_info

    # 1번째로 TG들의 node label 설정

    #print(f"iter: {1}")

    for i in range(len(gene_sig_df)): # (len(gene_sig_df)): (100)   10)):#
        count = node_gene_sig[i]
        prot1 = gene_sig_df.index[i]
        
        if i >= 5: # tg들
            prot2 = link_filt[link_filt['protein2'] == prot1]
            #print(prot2) # 해당 gene과 연결된 edge 정보 출력

            for k in range(len(prot2)):
                node_gene = prot2.iloc[k][0]
                count+=prot2.iloc[k][2] * np.array(gene_sig_df.loc[node_gene])
                count = count.tolist()
            #print(f'count : {count}') # edge 정보(가중치) x label 결과

        idx = np.array(np.nonzero(count)).tolist() # count 배열 중 큰 label을 index로 설정하여 출력
        re_count = [0 for i in range(9)]
        for nz in idx[0]:
            re_count[nz] = 1
        node_gene_sig[i] = re_count


    # 두번째부터는 edge 정보를 이용해서
    for it in range(1,n):
        #print(f"iter: {it+1}")

        for i in range(len(gene_sig_df)): # (len(gene_sig_df)): (100)   10)):#

            if np.sum(node_gene_sig[i])>1 or (gene_sig_df.index[i] in mul_gene):
                count = [0 for i in range(9)]
                prot1 = gene_sig_df.index[i]
                if it == 1:
                    mul_gene.append(prot1)
                
                if i < 5: # TF들
                    prot2 = link_filt[link_filt['protein1'] == prot1]
                    #print(prot2) # 해당 gene과 연결된 edge 정보 출력

                    for k in range(len(prot2)):
                        node_gene = prot2.iloc[k][1]
                        count+=prot2.iloc[k][2] * np.array(gene_sig_df.loc[node_gene])
                        count = count.tolist()
                    #print(f'count : {count}') # edge 정보(가중치) x label 결과
                elif i >= 5: # tg들
                    prot2 = link_filt[link_filt['protein2'] == prot1]
                    for k in range(len(prot2)):
                        node_gene = prot2.iloc[k][0]
                        count+=prot2.iloc[k][2] * np.array(gene_sig_df.loc[node_gene])
                        count = count.tolist()
                    #print(f'count : {count}') # edge 정보(가중치) x label 결과
            
                idx = np.argmax(count) # count 배열 중 큰 label을 index로 설정하여 출력
                re_count = [0 for i in range(9)]
                re_count[idx] = 1
                node_gene_sig[i] = re_count


    return node_gene_sig, mul_gene









mul_gene = []
node_gene_sig, mul = node_iter(4,node_gene_sig,mul_gene)










df_index = gene_sig_df.index
df_columns = ['sig.'+str(i) for i in range(int(args.sig_num))]

gene_sig_df2 = pd.DataFrame(data=node_gene_sig, index=df_index, columns=df_columns, dtype=None, copy=False)

gene_sig_df2.to_csv('./output/Node_classi/node_result_'+str(args.pos_neg)+".csv", sep=',', index=True, header=True)









# 그래프를 그리기 위한 node 리스트 구하기

list_sig = []

tg_sig_df = gene_sig_df2[5:]

col = list(tg_sig_df.columns)

for col_name in col:
    sig = list(tg_sig_df[tg_sig_df[col_name] == 1].index)
    list_sig.append(sig)

tf_sig = []

tf_sig_df = gene_sig_df2[:5]

col = list(tf_sig_df.columns)

for col_name in col:
    sig = list(tf_sig_df[tf_sig_df[col_name] == 1].index)
    tf_sig.append(sig)





link_filt['combined_score'] = link_filt['combined_score'].transform(lambda x : 1000 - x)



# 그래프 생성
classi_G = nx.from_pandas_edgelist(link_filt, 'protein1', 'protein2', edge_attr='combined_score', create_using = nx.Graph())


# 컬러 리스트
c_list = ['skyblue', 'pink', 'Magenta', 'Yellow', 'Gray', 'Purple', 'Brown', '#cbbbf4', 'Teal']

c_list2 = ['#5CACEE', '#FF6B8B', 'Magenta', 'Yellow', 'Gray', 'Purple', 'Brown', '#6D66C2', 'Teal']

# 그래프 그리기
pos = nx.spring_layout(classi_G)

# 노드 리스트별로 그래프에 노드 추가
for i, node_list in enumerate(list_sig):
    nx.draw_networkx_nodes(classi_G, pos, nodelist=node_list, node_color=c_list[i], label=f'{col[i]}', node_size=50)
    
    # TF의 노드 그리기,레이블 표시
for i, tf_list in enumerate(tf_sig):
    nx.draw_networkx_nodes(classi_G, pos, nodelist=tf_list, node_color=c_list[i], edgecolors=c_list2[i], linewidths=2, node_size=2500) 

nx.draw_networkx_edges(classi_G, pos, style='solid', edge_color='black', alpha=0.1) # 

node_labels = {node: node if node in TF_list else '' for node in classi_G.nodes()}
nx.draw_networkx_labels(classi_G, pos, labels=node_labels, font_size=15, font_color='black')








plt.savefig('./output/Node_classi/node_figure_'+str(args.pos_neg)+'png')

plt.rcParams['figure.figsize'] = [20, 10]
plt.legend(loc='upper right', fontsize=15)
plt.show()
