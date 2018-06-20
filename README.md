# multiSampleCaller
Calls somatic SNVs, indels, and allelic copy number jointly across multiple samples from the same patient.  These can be standard tumor/normal pair, longitudinal samples, primary/met, etc.  Can also be used for tumor only calling, ideally with a high tumor content and a low tumor content sample.

## Prerequisites
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

### SNP VCFs
-VCFs can be downloaded from 1000 genomes, Exac, or similar population genotyping projects
-VCFs must have the population allele frequency as "AF" in the INFO field
-There must be one VCF per chromsome
-VCFs must be bgzipped and tabix indexed 


## Overview
Running lumosVar involves two main steps
1. normalMetrics: analyzes a set of unmatched controls to find average read depths and position quality metrics
   - IMPORTANT - unmatched controls must be generated using the same exome capture as the tumors
2. lumosVarMain: call somatic, germline, and copy number variants

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

- Normal metrics is run using the [runNormalMetrics.py](scripts/runNormalMetrics.py)
>python runNormalMetrics.py controlsConfig.yaml
