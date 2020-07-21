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

The fam file continues the individual IDs for mapping to the phenotype file in the second column.

### 2) phenotype

This is the full path to the phenotype file. The first line in the file contains column names. The file has a column labeled IID, which contains the same information as present in the genotype file's second column (individual IDs). Genotype and phenotype data will be linked within the script, such that the rows do not to be aligned prior to running.

The remaining columns in the phenotype file contain phenotypes with the phenotype names shown in the first row.

### 3) phenoName

The phenotype to be used in analysis. This matches identically to the relevant column name in the phenotype file.

### 4) output

This is the full path and file name to where results should be saved.

### 5) num_cores

Number of cores to be used in parallelizing the analysis. Use 1 for no parallelization.


## Running the analysis

The analysis can be run by entering the following in the command-line:

```
Rscript DRM.R $genotype $phenotype $phenoName $output $num_cores
```

