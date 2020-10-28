function chrTable=chr2idx(inputParam)
%preprocessTumorOnly - creates data structures for tumor only calling
%calls parsePileupData.packed.pl to parse samtools output
%
% Syntax:  [T, E]=preprocessTumorOnly(inputParam,paramFile)
%
% Inputs:
%   inputParam - data structure with the following fields: regionsFile,
%       numCPU, outname, blockSize, snpVCFpath, snpVCFname,
%       workingDirectory, tabixPath, NormalBase
%   
% Outputs:
%   T - table of data by position
%   E - table of data by exon
%
% Other m-files required: none
% Other requirements: parsePileupData.packed.pl, samtools, htslib
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016

%------------- BEGIN CODE --------------

%[status,out]=system('printenv')
%profile('-memory','on');
%profile on;

sexChr=regexprep(inputParam.sexChr,'''','');
sexChr=regexp(sexChr,',','split');
if max(cellfun('length',sexChr))==0
    chrList=cellstr(num2str(inputParam.autosomes','%-d'))
else
    chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];
end
if isfield(inputParam,'contigPrefix')
    chrList=strcat(inputParam.contigPrefix,chrList)
end

[~,out]=system(['samtools view -H `head -n1 ' inputParam.bamList ' ` | grep ^@SQ']);
 
contigs=strtok(extractAfter(strsplit(out,'\n'),'SN:')); 

[lia,locb]=ismember(chrList,contigs);

chrTable=table();
chrTable.chrName=chrList(lia);
chrTable.chrIdx=locb(lia);

if sum(~lia)>0
    message=['warning chr ' strjoin(chrList(~lia),',') ' not found in bam']
end
