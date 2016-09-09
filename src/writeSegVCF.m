function writeSegVCF(segsTable,inputParam)
%writeSegVCF - writes VCF for copy number alterations
%
% Syntax: writeSegVCF(segsTable,inputParam)
%
% Inputs:
%   segsTable: matrix of segment data with columns:
%       1-'Chr',2-'StartPos',3-'EndPos',4-'segmentMean Tumor/Normal Log Ratio',
%       5-'N',6-'M',7-'F',8-'W',9-'log2FC'
%  inputParam: structure with all parameters   
%
% Outputs:
%    writes a VCF file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, callCNA, fitCNA

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016
%------------- BEGIN CODE --------------

fout=fopen([inputParam.outName '.cna.seg.vcf'],'w');

%%% print VCF header
fprintf(fout,'%s','##fileformat=VCFv4.2\n');
fprintf(fout,['##fileData=' datestr(clock) '\n']);
inputFields=fieldnames(inputParam);
for i=1:length(inputFields)
    if(isnumeric(inputParam.(inputFields{i})))
        fprintf(fout,['##INPUT=<' inputFields{i} '=' mat2str(inputParam.(inputFields{i})') '>\n']);
    else
        fprintf(fout,['##INPUT=<' inputFields{i} '=' inputParam.(inputFields{i}) '>\n']);
    end
end
fprintf(fout,['##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy Number">\n']);
fprintf(fout,['##INFO=<ID=MACN,Number=1,Type=Integer,Description="Min Allele Copy Number">\n']);
fprintf(fout,['##INFO=<ID=LOG2FC,Number=1,Type=Float,Description="log2 fold change">\n']);
fprintf(fout,['##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of Segment">\n']);
fprintf(fout,['##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural Variant">\n']);
fprintf(fout,['##INFO=<ID=END,Number=1,Type=String,Description="End of Variant">\n']);
fprintf(fout,['##INFO=<ID=CNF,Number=1,Type=Float,Description="Fraction containg Copy Number Alteration">\n']);

%%% determine alteration type
type(segsTable.N>2 & segsTable.M>0,:)={'DUP'};
type(segsTable.N>2 & segsTable.M==0,:)={'DUPLOH'};
type(segsTable.N==2 & segsTable.M==1,:)={'NONE'};
type(segsTable.N==2 & segsTable.M==0,:)={'LOH'};
type(segsTable.N<2,:)={'DEL'};

%%% construct info field
Info=cellstr(strcat('CN=',num2str(segsTable.N,'%-.0f'),';MACN=',num2str(segsTable.M,'%-.0f')));
Info=strcat(Info,';SVLEN=',num2str(segsTable.EndPos-segsTable.StartPos,'%-.0f'));
Info=strcat(Info,';SVTYPE=',type,';END=',num2str(segsTable.EndPos,'%-.0f'));

for i=1:size(segsTable.F,2)
    formatStr(segsTable.N~=2 | segsTable.M~=1,i)=cellstr(strcat(num2str(segsTable.F(segsTable.N~=2 | segsTable.M~=1,i),'%-.3f'),':',num2str(segsTable.log2FC(segsTable.N~=2 | segsTable.M~=1,i),'%-.2f')));
    formatStr(segsTable.N==2 & segsTable.M==1,i)=cellstr(strcat('NA:',num2str(segsTable.log2FC(segsTable.N==2 & segsTable.M==1,i),'%-.2f')));
end
    
formatFields=repmat({'CNF:LOG2FC'},size(segsTable,1),1);

%%% write output
outData=[num2cell(segsTable.Chr) num2cell(segsTable.StartPos) num2cell(segsTable.EndPos) cellstr(char(ones(size(segsTable,1),1)*78)) strcat('<',type,'>') num2cell(abs(mean(segsTable.log2FC,2))) cellstr(char(ones(size(segsTable,1),1)*46)) Info formatFields formatStr];
headers={'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILT', 'INFO', 'FORMAT'};
headers=[headers regexp(inputParam.sampleNames,',','split')];
for i=1:length(headers)
    fprintf(fout,'%s\t',headers{i});
end
for i=1:size(outData,1)
    fprintf(fout,strcat('\n%d\t%d\t%d\t%s\t%s\t%f\t%s\t%s\t%s',repmat('\t%s',1,size(segsTable.F,2))),outData{i,:});
end

