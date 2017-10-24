# multiSampleCaller
Calls somatic SNVs, indels, and allelic copy number jointly across multiple samples from the same patient.  These can be standard tumor/normal pair, longitudinal samples, primary/met, etc.  Can also be used for tumor only calling, ideally with a high tumor content and a low tumor content sample.

PREREQUISITES
samtools and htslib (tested with 1.2, 1.3, and 1.4)
http://www.htslib.org/download/
Matlab Runtime (MCR 9.0)
http://www.mathworks.com/products/compiler/mcr/

BAM PREPERATION 
Bams should be created using bwa-mem
http://bio-bwa.sourceforge.net/
It is recommended (but not required) to run abra reassembly to improve indel calling and allele frequency estimates
https://github.com/mozack/abra

OVERVIEW
Running lumosVar involves three main steps:
1. lumosVarNormalMetrics: analyzes a set of unmatched controls to find average read depths and position quality metrics
2. gvm: reads data from bam files and creates table of read counts and metrics for each canidate variant position
3. lumosVarMain: call somatic, germline, and copy number variants

USAGE 
1.  Analyze unmatched controls:

./lumosVarNormalMetrics CONFIGNAME CHR

example:

./lumosVarNormalMetrics controlsConfig.yaml 1

This step is run seperately for each chromosome (See exampleScripts/runPrintNormalMetrics.sh for how to run all chromosomes in parallel on an HPC system).  CHR should always be an integer, and the sex chromosomes are assigned numbers starting at one plus the largest autosome.  The number of autosomes and symbols for the sex chromosomes are specified in the yaml config file.   When running, the files parsePileupData.packed.pl and indexControlMetrics.sh must be in your working directory specified in the config files. You must also have a text file that lists the paths to the control bams you wish to analyze. You specify the path to your bam list file as well as your samtools and htslib installations, your bed file and your output path in the yaml config file. See configTemplates/controlsConfig.yaml analyzeControls will generate a set of bgzip tabix indexed files (one per chromosome) with quality metrics to be used by the main caller.

2. Preprocces bam files:

Preprocessing is done by a seperate tool called gvm, which is available from:

https://github.com/tgen/gvm
 
Detailed instructions for running are provided in the gvm repository.  It takes the same config file as input as lumosVarMain.  It will produce two tab delimited text file for each chromosome, a "pos" file that contains metrics on each canidate variant position, and an "exon" file that contains metrics for each region in the BED file.

3.  Main Caller

./lumosVarMain CONFIGNAME

example

./lumosVarMain tumorConfig.yaml

This runs the main caller.  If the file specified by the outMat path in the yaml config file exists, it will read the data from the mat file.  If starting from a mat file, it will also use the sampleNames, bamList, NormalSample, and priorF from the mat file and ignore the values for those parameters specified in the yaml.  If the mat file does not exist, it will read in the data from the chromosome text files created in step two that have file name prefix matching the outName specified in the yaml file. When running, cghcbshybridnu.mat must be in your working directory.  The caller will produce an SNV and indel VCF (.tumorOnly.all.vcf), a copy number and LOH VCF (.cna.seg.vcf), a clone summary table (.cloneSummary.csv), and a summary figure (.png)


