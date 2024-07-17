# Mutational Signature Analysis tool

This framework can discover candidate transcription factors(TFs) that regulate target gene expression by mutational signatures.

<!--
나중에 여기에 논문 링크 넣기
-->

## Introduction of MutTF

>Background: Over the past decade, mutational signatures representing specific patterns of mutations have been studied to decipher the multiple causes of somatic mutations in cancer genome. However, previous studies have focused on identifying the novel mutational signatures or improving the accuracy of signature extraction. In this study, we propose a multi-omics analysis framework, MutTF, that discovers the candidate Mutational signatures that induce somatic mutations on Transcription Factors (TFs) and subsequently affect the regulation of target genes (TGs) of TFs.<br><br>
Results: In this study, 59 colon cancer samples from PCAWG whole-genome sequencing data were analyzed with the proposed framework. Based on the correlation analysis between signature-wise mutation counts of TF and Gene Set Variation Analysis scores of TGs, we identified 14 significant Signature-TF pairs. Most of the pairs were further validated with VIPER and statistical tests to show that expressions of TGs or TF itself were significantly different between samples with and without signature-induced mutations.<br><br>
Conclusions: The proposed framework enabled the selection of mutational signatures that influenced TF’s regulation on TGs in colorectal cancer, including three pairs that showed differential expression in both TF and TGs: SBS10a-ZBTB7B, SBS10b-GMEB2, and SBS94-RORA. Most of the signatures and TF genes detected by MutTF were relevant to colorectal cancer based on literature analysis, showing the reliability of the proposed framework.

<br>

![Workflow of 'MutTF'](./readme_img/workflow_new.png)

## Input Data

Input data should be put in **vcf files** representing **mutation counts**.


If everything is ready, please put the mutation data that needs analysis in directory named './input_data'

>The sample data is in directory named './sample'

## Installation
Clone repository.
```
git clone https://github.com/BML-cbnu/MutTF
cd MutTF
```

## How to execute code

- [Step1) Matrix Generator](#Step1-Matrix-Generator)   
- [Step2) NMF](#Step2-NMF)   
- [Step3) Gene_count](#Step3-Gene_count)   
- [Step4) gsva](#Step4-gsva)   
- [Step5) MutTF](#Step5-MutTF)
- [Optional step](#Optional-Code)   
   

In the command line, please run the following:

### Step1. Matrix Generator

* input: vcf file
* output: Matrix file (.all)
* variable:
  * [reference genome] => Enter the reference genome you want to analyze (e.g. GRCh37).
* Use **SigprofilerMatrixGenerator** to convert input files into count matrix(M).
* We referred from https://cancer.sanger.ac.uk/signatures/tools/.
* The results are as shown in the table below:

| MutationType | Sample 1 | Sample 2 | ... |
| --- | --- | --- | --- |
| A[C>A]A | 200 | 143 |
| A[C>A]C | 21 | 131 |
| ... | 214 | 654 |

```bash
$ python MatGen.py --ref_genome=[reference genome]
```

---
### Step2. NMF

* input: count matrix(M)
* variable:
  * [reference genome] => Enter the reference genome you want to analyze (e.g. GRCh37).
  * [minimum] => Minimum number of signatures to extract
  * [maximum] => Maximum number of signatures to extract
* The results made through **MatGen.py** were used as input data.
* We use SBS96.all(96 types of mutations in Single Base Substation).
* After execution, the optimal number of signature will be selected and used for analysis (Refer to './ext_data/SBS/SBS96_selection_plot.pdf' for the best number of signature).
* The results are as shown in the tables below: <br>

> Exposure Matrix

| Samples | SBS96A | SBS96B | ... |
| --- | --- | --- | --- |
| Sample 1 | 22 | 40 |
| Sample 2 | 35 | 13 |
| ... | 16 | 32 |

> Process Matrix

| MutationType | SBS96A | SBS96B | ... |
| --- | --- | --- | --- |
| A[C>A]A | 0.024 | 0.014 |
| A[C>A]C | 0.012 | 0.052 |
| ... | 0.081 | 0.068 |

```bash
$ python NMF.py --ref_genome=[reference genome] --min=[minimum] --max=[maximum]
```

---

### Step3. Gene_count

* input: vcf file, fasta, gtf file
* variable:
  * [reference genome] => Enter the reference genome you want to analyze (e.g. GRCh37).
* Before we calculate the contribution of signatures, we need **gene-specific counts** from a particular gene region.
* It is a code that calculates gene-specific counts using the annotation file of reference genome.
* The results are as shown in the table below:

|  | Gene 1 | Gene 2 | ... |
| --- | --- | --- | --- |
| ACA>A | 2 | 0 |
| ACC>A | 0 | 1 |
| ... | 1 | 1 |

```bash
$ python Gene_count.py --ref_genome=[reference genome]
```

---

### Step4. GSVA

* input: Original TF-TG geneset data, Expression file
* output : separated TF-TG geneset data, GSVA output file

```bash
$ python GSVA.py -g [Original TF-TG geneset data] -e [Expression file] -o1 [Separated TF-TG geneset data] -o2 [GSVA output file]
```

---

### Step5. MutTF

* Our main analysis model, **'MutTF'**
* input: Gene count Matrix, GSVA score, TF-TG  data
* variable:
  * [gsva_result_folder] => Enter the directory where the results file from **gsva.py** is located.
  * [tf_database_file] => Enter the TF-TG database file you want to analyze
* Calculate the signature's contribution (by sample).
* The correlation between gene-specific counts by signature and the GSVA score was analyzed.
* The result file is saved in the directory named './output/Cor/'.
* The results are as shown in the table below:

| No. | Gene | sig | r | p |
| --- | --- | --- | --- | --- |
| 0 | Gene id | signature id | correlation coefficient | p-value |

```bash
$ python MutTF.py --gsva_folder=[gsva_result_folder] --tf_file=[tf_database_file]
```

---

## Optional Code

**Node_classification**

* input: gene expression data
* variable:
  * [pos or neg] => Enter the group for which you want to proceed node classification (pos or neg)
  * [tg divided into two groups] => As a result of correlation analysis based on gene expression, it means tf-tg data divided into positive and negative.
  * [the number of signatures] => The optimal number of signatures used for analysis
* You can proceed with node classification for multi-signature gene based on the result file from **Step 5**.
* The result file is saved in the directory named './output/Node_classi/'.
* The visualized graph figure is saved as './output/Node_classi/node_figure_XXX.png'

```bash
$ python Node_classification.py --pos_neg=[pos or neg] --tf_group_file=[tg divided into two groups] --sig_num=[the number of signatures]
```

---

**Denovo_cosine**

* input: matrix P
* variable:
  * [reference genome] => Enter the reference genome you want to analyze (e.g. GRCh37).
  * [version] => Enter the version of cosmic signature you want to compare (e.g. 3.3.1)
* A heat map shows how the optimal signature extracted by De novo Signatures from **NMF.py** is similar to the cosmic signature.
* We referred from https://cancer.sanger.ac.uk/signatures/downloads/.
* Examples are as follows:
![Denovo_cosine](./readme_img/cosine.png) 

```bash
$ python Denovo_cosine.py --ref_genome=[reference genome] --version=[version]
```
