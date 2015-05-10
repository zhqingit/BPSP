# bpsp 
BPSP is a method that can identify branch point based on only the intron sequence.  

### How BPSP works?

BPSP predicts the branch point sequence by integrating the degenerative motif of BPS and PPT characteristics. Specifically, BPSP uses a mixture model to infer the BPS motif and a set of weighted octanucleotides to estimate the contribution of the 65kDa subunit of U2AF (U2AF65). 

### Detailed information on GIREMI and citation

A paper describing GIREMI is published at Nature Methods (http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3314.html).  


### Dependent libraries or software

- R(http://www.r-project.org/) : for generating the PPT scores and training data

### Quickstart

#### Step 1: 
R CMD BATCH --no-save --no-restore '--args dist_polyn=fin PPTscore=fout' get_ppt_score.r tmp.Rout
##### 
-   fin:  input file with the counts of octanucleotide in background, branch point and PPT regions
-   fout: output file with PPT scores

#### Step 2: 
R CMD BATCH --no-save --no-restore '--args dist_polyn=fin train_data=fout' get_train_data.r tmp.Rout
#####
-   fin:  input the counts of  heptanucleotide in background, branch point and PPT regions
-   fout: output file with the data set training the MM model

#### Step 3: 
bpmotif.pl [options]

##### Required:
-  `--pmotif   FILE   prior motif (energy motif in our paper)`  
-  `--motif    FILE   inferred motif`  
-  `--bpolyn   FILE   training data from step 2`  

##### Options:
-  `-m, --min          INT    minimal number of total reads covering candidate editing sites  [default: 5]`   
-  `-p, --paired-end   INT    1:paired-end RNA-Seq reads; 0:single-end [default: 1]`   
-  `-s, --strand       INT    0:non-strand specific RNA-Seq; 1: strand-specific RNA-Seq and read 1 (first read for the paired-end reads) is sense to RNA; 2: strand-specific RNA-Seq and read 1 is anti-sense to RNA [default: 0]`

In this paper, the BPS region is defined as the 21-34 nucleotides (nt) upstream of the 3SS where most branchpoints are located[27]. We also define the PPT region as the 4-15 nt upstream of the 3SS; and the background region as the 187-200 nt upstream of 3SS because no general splice element is reported for this region. Noticeably, these regions are only relatively enriched or devoid of corresponding signals, and their contrast will provide a statistical clue on what the true signal looks like.

##### Required format of the file containing the list of SNVs (-l option):

- column 1 : The name of the chromosome or scaffold   
- column 2 : The starting position of the SNV in the chromosome or scaffold (0-based)   
- column 3 : The ending position of the SNV in the chromosome or scaffold (1-based)   
- column 4 : The name of the gene harboring this SNV; “Inte”: the SNV resides in the Intergenic region   
- column 5 : A flag, 1: the SNV belongs to dbSNP; 0: otherwise   
- column 6 : Strand (+ or -); “#” for “Inte” gene  


##### Format of the output file:

######NOTE: This output file includes a rich list of information about the SNVs. Not all sites in this file are predicted as RNA editing sites, see the ifRNAE field.

- chr                     : Name of the chromosome or scaffold     
- coordinate              : Position of the SNVs in the chromosome or scaffold (1-based)    
- strand                  : Strand information
- ifSNP                   : 1, If the SNV is included in dbSNP; 0: otherwise.
- gene                    : Name of the gene harboring this SNV
- reference_base          : The nucleotide of this SNV in the reference chromosome (+ strand)
- upstream_1base          : The upstream neighboring nucleotide of this SNV in the reference chromosome (+ strand)
- downstream_1base        : The downstream neighboring nucleotide of this SNV in the reference chromosome  (+ strand)
- major_base              : The major nucleotide of the SNV in the RNA-seq data     
- major_count             : Number of reads with the major nucleotide    
- tot_count               : Total number of reads covering this SNV in the RNA-Seq data   
- major_ratio             : The ratio of major nucleotide (major_count/tot_count)   
- MI                      : The mutual information of this SNV if a value exists   
- pvalue_mi               : P-value from the MI test if applicable    
- estimated_allelic_ratio : Estimated allelic ratio of the gene harboring this SNV
- ifNEG                   : 1: this SNV was a negative control in the training data  
- RNAE_t                  : Type of RNA editing or RNA-DNA mismatches (A-to-G, etc)
- A,C,G,T                 : Numbers of reads with specific nucleotides at this site
- ifRNAE                  : 1: the SNV is predicted as an RNA editing site based on MI analysis; 
						  2: the SNV is predicted as an RNA editing site based on GLM 
						  0: the SNV is not predicted as an RNA editing site
  
