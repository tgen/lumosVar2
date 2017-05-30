# multiSampleCaller
Calls somatic SNVs, indels, and allelic copy number jointly across multiple samples from the same patient.  These can be standard tumor/normal pair, longitudinal samples, primary/met, etc.  Can also be used for tumor only calling, ideally with a high tumor content and a low tumor content sample.

PREREQUISITES
samtools and htslib (tested with 1.2)
http://www.htslib.org/download/
Matlab Runtime (MCR 9.0)
http://www.mathworks.com/products/compiler/mcr/

BAM PREPERATION 
Bams should be created using bwa-mem
http://bio-bwa.sourceforge.net/
It is recommended (but not required) to run abra reassembly to improve indel calling and allele frequency estimates
https://github.com/mozack/abra

USAGE OVERVIEW
Running lumosVar involves three main steps:
1. printNomralMetrics: analyzes a set of unmatched controls to find average read depths and position quality metrics
2. preprocess: reads data from bam files and creates table of read counts and metrics for each canidate variant position
3. main: call somatic, germline, and copy number variants

USAGE 
First a set of unmatched controls must be analyzed to find average read depths and position quality metrics:
module load MCR/9.0
./analyzeControls controlsConfig.yaml
When running, the files parsePileupData.packed.pl and indexControlMetrics.sh must be in your working directory specified in the config files. You must also have a text file that lists the paths to the control bams you wish to analyze. You specify the path to your bam list file as well as your samtools and htslib installations, your bed file and your output path in the yaml config file. See configTemplates/controlsConfig.yaml analyzeControls will generate a set of bgzip tabix indexed files (one per chromosome) with quality metrics to be used by the main caller.

To call a tumor bam:
module load MCR/9.0
./callTumorOnly tumorConfig.yaml
When running, the files parsePileupData.packed.pl and cghcbshybridnu.mat must be in the working directory specified in the config files. You specify the path to the bam, as well as the control metrics files, your samtools and htslib installations, your bed file, and your output path in the yaml config file. See configTemplates/tumorConfig.yaml. The caller will produce an SNV and indel VCF (.tumorOnly.all.vcf), a copy number and LOH VCF (.cna.seg.vcf), a clone summary table (.cloneSummary.csv), and a summary figure (.png)
