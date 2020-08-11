# DRM

## Information about the repository

This repository contains the computer code to implement the Deviation Regression Model (DRM) on PLINK genotype data in R. The DRM tests for an association between a SNP's genotype and the variance of a quantitative phenotype. These SNPs can be referred to as vQTLs. 

If you have any questions, please contact:

Andrew Marderstein \
anm2868@med.cornell.edu

## Dependencies

The DRM.R script requires the following R packages:

	BEDMatrix
	data.table
	parallel
	
These packages load PLINK data into R, potentially large files into R, and parallelize the genome-wide analysis.

## Input arguments

### 1) genotype

This is the full path to the genotype file. The genetic data is in PLINK bed/bim/fam format. The path includes the prefix, but not the file suffix. For example:

```
/path/to/files/plinkFiles
```

where the PLINK genotype files are:

```
/path/to/files/plinkFiles.bim
/path/to/files/plinkFiles.fam
/path/to/files/plinkFiles.bed
```

Within the second column, the fam file contains the individual IDs for mapping to the phenotype file.

### 2) phenotype

This is the full path to the phenotype file. The first line in the file contains column names. The file has a column labeled IID, which contains the same information as present in the genotype file's second column (individual IDs). Genotype and phenotype data will be linked within the script, such that the rows do not to be aligned prior to running.

The remaining columns in the phenotype file contain phenotypes with the phenotype names shown in the first row.

Example:

```
(vQTL) [anm2868@cannes vqtl]$ head /athena/elementolab/scratch/anm2868/bmi_data.txt
FID     IID     bmi
1001 9888 -1.065152800531
1002 9889 2.76930734592663
1003 9890 0.293101015253496
1004 9891 0.424138249409053
1005 9892 1.27352707643111
...
...
...
```

### 3) phenoName

The phenotype to be used in analysis. This matches identically to the relevant column name in the phenotype file.

### 4) output

This is the full path and file name to where results should be saved.

### 5) num_cores

Number of cores to be used in parallelizing the analysis. Use 1 for no parallelization.

### 6) startInd

OPTIONAL: The first column index in the genotype file to run analysis on. Genotype indices start at 1. If not specified, then the analysis will run from the 1st SNP.

### 7) endInd
OPTIONAL: The final column index in the genotype file to run analysis on. The analysis will run from startInd to endInd. If not specified, then the analysis will run from startInd to the last column index.



## Running the analysis

The analysis can be run by entering the following in the command-line:

```
Rscript DRM.R $genotype $phenotype $phenoName $output $num_cores $startInd $endInd
```

## Example

The following analysis performs the DRM variance test on the first 100 SNPs found on chromosome 22. BMI is used as the phenotype and 4 cores are used for parallelization. 

```
CHR=22
startInd=1
endInd=100
genotype=/athena/elementolab/scratch/anm2868/ukbb.$CHR.impute
phenotype=/athena/elementolab/scratch/anm2868/bmi_data.txt
phenoName=bmi
output=/athena/elementolab/scratch/anm2868/bmi_results.vGWAS.txt
num_cores=4

Rscript DRM.R $genotype $phenotype $phenoName $output $num_cores $startInd $endInd
```

While we suggest parallelizing the analysis around a subset of SNPs (e.g., 5,000 SNPs at a time), all SNPs within a genotype file can be ran by omitting the final two arguments:

```
Rscript DRM.R $genotype $phenotype $phenoName $output $num_cores
```

