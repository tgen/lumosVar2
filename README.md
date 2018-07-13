# LumosVar2
Calls somatic SNVs, indels, and allelic copy number jointly across multiple samples from the same patient.  These can be standard tumor/normal pair, longitudinal samples, primary/met, etc.  Can also be used for tumor only calling, ideally with a high tumor content and a low tumor content sample.

### Citation
Joint analysis of matched tumor samples with varying tumor contents improves somatic variant calling in the absence of a germline sample
Rebecca F Halperin, Winnie S Liang, Sidharth Kulkarni, Erica E. Tassone, Jonathan Adkins, Daniel Enriquez, Nhan L Tran, Nicole C Hank, James Newell, Chinnappa Kodira, Ronald Korn, Michael E Berens, Seungchan Kim, Sara A Byron
bioRxiv 364943; doi: https://doi.org/10.1101/364943

## Prerequisites
### System Requirements
- Linux OS (tested on Centos 7 and Ubuntu 17)

### Dependencies
- bcftools and htslib (tested with 1.2-1.8)
http://www.htslib.org/download/
- Matlab Runtime (MCR 9.0)
http://www.mathworks.com/products/compiler/mcr/
- pyyaml and libYAML
https://pyyaml.org/wiki/PyYAMLDocumentation
- GSL - GNU Scientific Library
https://www.gnu.org/software/gsl/

### Bam preperation
- Bams should be created using bwa-mem
http://bio-bwa.sourceforge.net/
- Bams should be indexed by samtools

### VCFs
- SNP VCFs can be downloaded from 1000 genomes, Exac, or similar population genotyping projects
- SNP VCFs must have the population allele frequency as "AF" in the INFO field
- Cosmic VCFs can be downloaded from cosmic or other cancer mutation database
- Cosmic VCFs must have "CNT" in the info field indicating count of samples having mutation
- There must be one VCF per chromsome for both SNP and COSMIC VCFs
- VCFs must be bgzipped and tabix indexed 

## Overview
Running lumosVar involves two main steps
1. normalMetrics: analyzes a set of unmatched controls to find average read depths and position quality metrics
   - IMPORTANT - unmatched controls must be generated using the same exome capture as the tumors
2. lumosVarMain: call somatic, germline, and copy number variants

### Example dataset
http://tools.tgen.org/Files/lumosVar

### Notes on pileup engine
- lumosVar uses a custom pileup engine https://github.com/tgen/gvm
- precompiled binary is provided for the pileup engine [gvm](bin/gvm) so building is not required

## Normal Metrics 
- Inputs to normal metrics are defined in a yaml file, see [template](configTemplates/controlsConfigTemplate.yaml).  You will need to edit:
```
###input files
bamList:   BAMLIST      ###path to file contining paths to bams
regionsFile: BEDFILE    ###path to bed file defining regions targeted in exome
refGenome: REFGENOME    ###path to reference genome that was used to align bams
snpVCFpath: VCFPATH     ###path to vcfs containg population frequencies (one for each chromosome),
                        ###including part of filename before chromosome number
snpVCFname: VCFNAME     ###filename/extension of population vcf following chromosome number

###user inputs
sexList: SEXLIST        ###comma deliminated list of sex for each bam in the bamlist in the same order
gvmPath: GVMPATH        ###path to folder containing gvm executable
outfile: OUTFILE        ###path and name of output
```

- In order to correctly handle the sex chromosomes, the sex of the individuals in the bams list must be given as input.  We have provided a helper [script](scripts/guessSex.py) to determine sex from the bams if they are not known.  This script takes the same yaml file as input as the normal metrics, but only the following fields in the "input files" section are needed.
>python guessSex.py controlsConfig.yaml

The output will be written to \<BAMLIST\>.guessSex.txt.  The last line of the output contains the sexList for the yaml file

- Normal metrics is run using the [runNormalMetrics.py](scripts/runNormalMetrics.py).  It takes two input arguments, the yaml and the chromosome.  It needs to be run separately on each chromosome.
>python runNormalMetrics.py controlsConfig.yaml 21

## LumosVar Main
- Inputs to lumosVar are defined in a yaml file, see [template](configTemplates/lumosVarMainConfigTemplate.yaml).  You will need to edit:
```
#### Input Files
bamList: BAMLIST        ###path to file contining paths to bams
regionsFile: BED        ###path to bed file defining regions targeted in exome
snpVCFpath: VCFPATH     ###path to vcfs containg population frequencies (one for each chromosome),
                        ###including part of filename before chromosome number
snpVCFname: VCFNAME     ###filename/extension of population vcf following chromosome number
NormalBase: NMETRICS    ###path and filename of output form NormalMetrics step before chr number
cosmicVCF: COSMICVCF    ###path to cancer mutation count VCF
refGenome: REFGENOME    ###path to reference genome that was used to align bams

#### User Inputs
outName: OUTNAME        ###path and filename base for output files
outMat: OUTMAT          ###path and filename for "mat" file (matlab data) export
gvmPath: GVMPATH        ###path to folder containing gvm executable
workingDirectory: WORK  ###poth to folder containing files in the "work" directory in the github repo
NormalSample: NINDEX    ###position in bamList of normal sample, 0 indicates tumor only
priorF: PRIORF          ###vector of expected tumor fractions with one value per bam
                        ###for example [0.1;0.7] for a pair of bams with low and high expected tumor content
numCPU: CORES           ### number of parallel processors
```

To run lumosVar
```
./lumosVarMain lumosVarConfig.yaml
```

## LumosVar Output
- .lumosVarSNV.vcf - somatic and germline SNV/indel calls
- .lumosVarSeg.vcf - copy number calls by segment
- .exonData.tsv - copy number calls by regions in bed file
- .cloneSummary.tsv - summary of clonal variant groups
- .cloneSummary.pdf - graphical summary of clonal variant groups
- .groupLinePlots.pdf - line plots by clonal variant group
- .vafPlot.pdf - plots of variant allele fractions
- .cnaPlot.pdf - plot of copy number states
- \<SAMPLENAME\>.qualMetrics.tsv - quality metrics for candidate variant position
- .mat - matlab workspace export
