# bpsp 
BPSP is a method that can identify branch point based on only the intron sequence.  

### How BPSP works?

BPSP predicts the branch point sequence by integrating the degenerative motif of BPS and PPT characteristics. Specifically, BPSP uses a mixture model to infer the BPS motif and a set of weighted octanucleotides to estimate the contribution of the 65kDa subunit of U2AF (U2AF65). 

### Detailed information on BPSP and citation

A paper describing BPSP is under review.  


### Dependent libraries or software

- R(http://www.r-project.org/) : for generating the PPT scores and training data

### Quickstart

#### Step 1: 
R CMD BATCH --no-save --no-restore '--args dist_polyn=fin PPTscore=fout' get_ppt_score.r tmp.Rout
##### 
-  `--dist_polyn: FILE  input file with the counts of octanucleotide in background, branch point and PPT regions`
-  `--PPTscore:   FILE  output file with PPT scores`

#### Step 2: 
R CMD BATCH --no-save --no-restore '--args dist_polyn=fin train_data=fout' get_train_data.r tmp.Rout
#####
-  `--dist_polyn:  FILE input the counts of  heptanucleotide in background, branch point and PPT regions`
-  `--train_data:  FILE output file with the data set training the MM model`

#### Step 3: 
bpmotif.pl [options]

##### Required:
-  `--pmotif   FILE   prior motif (energy motif in our paper)`  
-  `--motif    FILE   inferred motif`  
-  `--bpolyn   FILE   training data from step 2`  

#### Step 3: 
bpsp.pl [options]

##### Required:
-  `--onlyM   INT   1 : only use motif; 0 : use both motif and PPT. [1]`
-  `--nBP     INT   reported mumber of BPs. [3]`
-  `--motif   FILE  motif file (inferred from last step)`
-  `--intron  FILE  intron file`
-  `--PPT     FILE  PPT score (inferred from step 1)`
-  `--out     FILE  output file`



##### Required:

##### Options:
-  `-m, --min          INT    minimal number of total reads covering candidate editing sites  [default: 5]`   
-  `-p, --paired-end   INT    1:paired-end RNA-Seq reads; 0:single-end [default: 1]`   
-  `-s, --strand       INT    0:non-strand specific RNA-Seq; 1: strand-specific RNA-Seq and read 1 (first read for the paired-end reads) is sense to RNA; 2: strand-specific RNA-Seq and read 1 is anti-sense to RNA [default: 0]`

###### NOTE:
Here, the BPS region is defined as the 21-34 nucleotides (nt) upstream of the 3SS where most branchpoints are located. We also define the PPT region as the 4-15 nt upstream of the 3SS; and the background region as the 187-200 nt upstream of 3SS because no general splice element is reported for this region. Noticeably, these regions are only relatively enriched or devoid of corresponding signals, and their contrast will provide a statistical clue on what the true signal looks like.

##### Required format of the dist_polyn file:

- column 1 : The name of the heptanucleotide/octanucleotide 
- column 2 : The numbers in the background region   
- column 3 : The numbers in the BPS region   
- column 4 : The numbers in the PPT region   

##### Required format of the prior motif file:

- column 1 : The name of the position 
- column 2 : The frequencies of A in all positions 
- column 3 : The frequencies of C in all positions 
- column 4 : The frequencies of G in all positions 
- column 5 : The frequencies of T in all positions 

##### Required format of the intron file:

- column 1 : The id of the intron 
- column 2 : The sequence of the intron 

##### Format of the output file:

- column 1 : The id of the intron 
- column 2 : The predicted branch point positions relative to the 3 pair end splicing site (multiple sites are separated by ",")
- column 3 : The predicted branch point sequences (multiple sites are separated by ",")
- column 4 : The scores of the predicted branch point sequences based on the BPS motif (multiple sites are separated by ",")
- column 5 : The scores of the PPTs fowllowing the predicted branch point sequences (multiple sites are separated by ",")
- column 6 : The scores of the BPS and PPTs (multiple sites are separated by ",")

#### Examples:

1. R CMD BATCH --no-save --no-restore '--args dist_polyn="example/dist_polyn_8" PPTscore="example/ppt.score.8"' get_ppt_score.r tmp.Rout
2. R CMD BATCH --no-save --no-restore '--args dist_polyn="example/dist_polyn_7" train_data="example/bp.polyn"' get_train_data.r tmp.Rout
3. bpmotif.pl --pmotif example/energy_motifBP --motif example/motifBP --bpolyn example/bp.polyn
4. bpsp.pl --motif example/motifBP --intron example/intron_example --PPT example/ppt.score.8 --out example/res
