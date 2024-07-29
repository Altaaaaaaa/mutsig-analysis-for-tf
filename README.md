# Mutational Signature Analysis tool

**MutTF** is **a multi-omics analysis framework** that combines gene expression data with mutational signatures based on non-negative matrix decomposition and correlation analysis. MutTF can discover **candidate transcription factors(TFs) that regulate target gene expression by mutational signatures**.

<!--
나중에 여기에 논문 링크 넣기
글자에 색 입히기: ${\textsf{\color{green}Green}}$
-->

![Workflow of 'MutTF'](./readme_img/workflow_new.png)

## Input Data

1. `Variant files (.vcf)` <br>
  VCF files from whole genome sequencing are required in this analysis. You can obtain VCF files after steps of aligning the reads to a reference genome, marking duplicates, performing local realignment and base quality recalibration, and calling variants using variant calling software.
2. `Expression file (.tsv)` <br>
  Expression file from RNA sequencing are required in this analysis. You can obtain expression files after steps of aligning the reads to a reference genome or transcriptome, and quantifying the read counts per gene or transcript.
3. (Optional) `TF-TG geneset file (.txt)` <br>
  We provide TF-TG geneset file obtained from hTFTarget, but if you wish to use a manual TF-TG geneset, the format of the file should be like this:
   | name | description | ... |
   | --- | --- | --- |
   | TF_0 | TG_0 |
   | TF_0 | TG_1 |
   | ... | TG_2 |

    **"name"** column contains TF, and **"description"** column contains group of TGs regulated by the corresponding TF.

  Since VCF files used in this project require a dbGaP access request, only the RNA expression data are provided.


## Installation
1. Clone repository.
```
git clone https://github.com/BML-cbnu/MutTF
cd MutTF
```
2. Install the requirments.
```
pip install -r requirements.txt
```


## How to execute code

- [Step1) Mutational signature extraction](#Step1-Mutational-signature-extraction)  
- [Step2) Gene_count](#Step2-Gene_count)   
- [Step3) GSVA](#Step3-GSVA)   
- [Step4) MutTF](#Step4-MutTF)
- [Optional step](#Optional-Code)   

   
In the command line, please run the following:

### Step1. Mutational signature extraction

* `Input`: <br>
  VCF file per sample (.vcf)
* `Variable`: <br>
  * [reference genome] &rarr; Enter the reference genome you want to analyze (e.g. GRCh37).
  * [minimum] &rarr; Minimum number of signatures to extract
  * [maximum] &rarr; Maximum number of signatures to extract
  * [input directory] &rarr; Directory where vcf files are located (e.g. input_data).
  * [output directory] &rarr; Directory where the output data should be stored.
  * [threads] &rarr; Number of threads to use in signature extraction
* `Description`: <br>
  Used **SigprofilerMatrixGenerator** to convert vcf files into count matrix, and used **sigProfilerExtractor** to extract signatures based on the count matrix generated.
  The optimal number of signature will be selected and used for further analysis. (Refer to './[output directory]/SBS96/SBS96_selection_plot.pdf' for the best number of signature)
  In this project, we used SBS96-based signatures (96 types of mutations in Single Base Substitution) in further analysis.
  Refer to https://cancer.sanger.ac.uk/signatures/tools/.
* `Output`: <br>
  Directory including signature extraction results (./[output directory]) <br>
  The results are as shown in the tables below: <br>

  **Exposure Matrix** <br>
  (./[output directory]/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt)

   | Samples | SBS96A | SBS96B | ... |
   | --- | --- | --- | --- |
   | Sample 1 | 22 | 40 |
   | Sample 2 | 35 | 13 |
   | ... | 16 | 32 |

  **Process Matrix** <br>
  (./[output directory]/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt)

   | MutationType | SBS96A | SBS96B | ... |
   | --- | --- | --- | --- |
   | A[C>A]A | 0.024 | 0.014 |
   | A[C>A]C | 0.012 | 0.052 |
   | ... | 0.081 | 0.068 |

```bash
$ python Signature_extraction.py --ref_genome=[reference genome] --minimum=[minimum] --maximum=[maximum] --input_dir=[input data directory] --output_dir=[output data directory] --threads=[threads]
```

---

### Step2. Gene_count

* Input: <br>
  VCF file per sample (.vcf)
* Variable: <br>
  * [reference genome] &rarr; Enter the reference genome you want to analyze (e.g. GRCh37).
  * [input directory] &rarr; Directory where vcf files are located (e.g. input_data).
  * [output directory] &rarr; Directory where the output data should be stored.
  * [threads] &rarr; Number of threads to use in multiprocessing
* Description: <br>
  Before we calculate the contribution of signatures, we need **gene-specific mutation counts** calculated using the annotation file of reference genome.
* Output: <br>
  Gene count matrix per sample (./[output directory]/*_cnt.csv) <br>
  The results are as shown in the table below: <br>

   |  | Gene 1 | Gene 2 | ... |
   | --- | --- | --- | --- |
   | ACA>A | 2 | 0 |
   | ACC>A | 0 | 1 |
   | ... | 1 | 1 |

```bash
$ python Gene_count.py --ref_genome=[reference genome] --input_dir=[input directory] --output_dir=[output directory] --threads=[threads]
```

---

### Step3. GSVA

* Input: <br>
  TF-TG geneset file (.txt), Expression file (.tsv)
* Variable: <br>
  * [TF-TG geneset file] &rarr; TF-TG geneset file (e.g. ./hTFTarget/colon_TF-Target-information.txt)
  * [Expression file] &rarr; File name of gene expression file
  * [GSVA output file] &rarr; File name of GSVA output results
* Description: <br>
  Seperate TG into positively and negatively regulated groups based on correlation coefficient with corresponding TF expression value.
  Based on these groups, perform GSVA. 
* Output: <br>
  GSVA output file (./[GSVA output file].tsv) <br>
  The results are as shown in the table below: <br>

   | Genesets | Sample 1 | Sample 2 | ... |
   | --- | --- | --- | --- |
   | TF1_0 | 0.4 | 0.3 |
   | TF1_1 | -0.9 | -0.1 |
   | ... | 0.1 | -0.6 |      

```bash
$ python GSVA.py -g [TF-TG geneset file] -e [Expression file] -o [GSVA output file]
```

---

### Step4. MutTF

* Input: <br>
  Signature extraction results (dir), Gene count matrix per sample (.csv), TF-TG geneset file (.txt), GSVA results (.tsv)
* Variable: <br>
  * [Signature extraction directory] &rarr; Directory of signature extraction results (output from **Signature_extraction.py**)
  * [Gene count directory] &rarr; Directory with gene-wise mutation count files (output from **Gene_count.py**)
  * [TF-TG geneset file] &rarr; TF-TG geneset file used in **GSVA.py** (e.g. ./hTFTarget/colon_TF-Target-information.txt)
  * [GSVA output file] &rarr; File name of GSVA results (output from **GSVA.py**)
  * [Correlation output directory] &rarr; Directory of correlation results between signature-induced mutation count and GSVA
* Description: <br>
  Calculate the signature's contribution (by sample).
  Analyze the correlation between gene-specific counts by signature and the GSVA score.
* Output: <br>
  Correlation result matrix (All and Filtered) (./[Correlation output directory]/[All/Filt]\_Result\_[pos/neg].csv) <br>
  The results are as shown in the table below: <br>

   | No. | Gene | sig | r | p |
   | --- | --- | --- | --- | --- |
   | 0 | Gene id | signature id | correlation coefficient | p-value |

```bash
$ python MutTF.py --ext_dir=[Signature extraction directory] --count_dir=[Gene count directory] --tf_file=[TF-TG geneset file] --gsva_folder=[GSVA output file] --corr_dir=[Correlation output directory]
```

---

## Optional Code

**Node_classification**

* Input: <br>
  Correlation results, GSVA results
* Variable: <br>
  * [Correlation results] &rarr; File name of correlation results (output from **MutTF.py**)
  * [direction] &rarr; Enter the group for which you want to proceed node classification (pos or neg)
  * [GSVA results] &rarr; File name of GSVA results (output from **GSVA.py**)
  * [Number of signatures] &rarr; The optimal number of signatures used for analysis
  * [Output directory] &rarr; Directory of node classification results
* Output: <br>
  Files including result of node classification and visualization. <br>
  The visualized graph figure is saved as '[Output_directory]/node_figure_XXX.png'

```bash
$ python Node_classification.py --corr_dir=[Correlation results] --pos_neg=[direction] --gsva=[GSVA file] --sig_num=[Number of signatures] --out_dir=[Output directory]
```

---

**Denovo_cosine**

* Input: <br>
  Matrix P
* Variable: <br>
  * [Signature extraction directory] &rarr; Directory of signature extraction results
  * [reference genome] &rarr; Enter the reference genome you want to analyze (e.g. GRCh37).
  * [version] &rarr; Enter the version of cosmic signature you want to compare (e.g. 3.3.1)
* Output: <br>
  Image showing cosine similarity

* A heat map shows how the optimal signature extracted by De novo Signatures from **NMF.py** is similar to the cosmic signature.
* We referred from https://cancer.sanger.ac.uk/signatures/downloads/.
* Examples are as follows:
![Denovo_cosine](./readme_img/cosine.png) 

```bash
$ python Denovo_cosine.py --ext_dir=[Signature extraction directory] --ref_genome=[reference genome] --version=[version]
```
