# Mutational Signature Analysis tool

This framework can discover candidate transcription factors(TFs) that regulate target gene expression by mutational signatures.

## input data

분석이 필요한 데이터는 input_data에 넣기(vcf 형태)
우리가 제공하는 기본 데이터는 data 폴더에 있음(e.g. )


## <a name="How to execute code"></a> How to execute code

To install the current version of this Github repo, git clone this repo or download the zip file.
Unzip the contents of SigProfilerExtractor-master.zip or the zip file of a corresponding branch.

In the command line, please run the following:
```bash
$ python CODE.py --data_folder=[datafolder name] --deseq_folder=[deseq result foloder name]
```
[datafolder name] => yMat = mutation count(vcf file to csv), result of Sigprofiler  
[deseqfolder name] => result of DESeq(tsv)



| Category | Parameter | Variable Type | Parameter Description |
| --------- | --------------------- | -------- |-------- |
| **Input Data** |  |  | |
|  | **input_type** | String | The type of input:<br>"vcf": used for vcf format inputs. Input data is signature analysis data obtained as a result of sigprofiler. As input data of sigprofiler, a somatic mutation dataset was used.  |
|  | **output** | String | The name of the output folder. The output folder will be generated in the current working directory.  |
|  | **reference_genome** | String | The name of the reference genome. The default reference genome is "GRCh38". This parameter is applicable only if the input_type is "vcf". | 
|  | **DESeq data** | String | DESeq represents the result of DEG analysis by obtaining the LogFC value of the mutation dataset, and refers to gene expression data. | 
|  | **Used data** | String | DESeq represents the result of DEG analysis by obtaining the LogFC value of the mutation dataset, and refers to gene expression data. | 
|  | **data** | String | <ul><li>COSMIC signature: 30 signature data already analyzed</li><li>Gene-Gene link data: Gene-gen network data for signature classification</li></ul> |
| **Output Data** |  |  | |
|  | **result of sigprofiler** | String | Results of analyzing somatic mutation sample data using signature analysis tool called sigprofiler.(process and exposure data)  |
|  | **heatmap** | Image | Heatmap to express the similarity between the most likely signatures analyzed using sigprofiler and already studied cosmic signature data through cosine similarity. | 
|  | **gene data** | Table | The result of identifying the genes whose correlation and p-value between each contribution of the sample and the gene expression data satisfies a certain range. | 
|  | **Results of node classification** | Image and Table | Through node classification, genes identified by multiple signatures are classified as single signatures that are most similar. | 
