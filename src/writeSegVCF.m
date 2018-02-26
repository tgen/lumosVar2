function writeSegVCF(segsTable,exonRD,CNAscale,Tcell,hetPos,inputParam)
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

fout=fopen([inputParam.outName '.lumosVarSeg.vcf'],'w');

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
fprintf(fout,['##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of Segment">\n']);
fprintf(fout,['##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural Variant">\n']);
fprintf(fout,['##INFO=<ID=END,Number=1,Type=String,Description="End of Variant">\n']);
fprintf(fout,['##INFO=<ID=EXONCOUNT,Number=1,Type=Float,Description="Number of Exons in Segment">\n']);
fprintf(fout,['##INFO=<ID=HETCOUNT,Number=1,Type=Float,Description="Number of Heterozygous Variants in Segment">\n']);

fprintf(fout,['##FORMAT=<ID=CNF,Number=1,Type=Float,Description="Fraction containg Copy Number Alteration">\n']);
fprintf(fout,['##FORMAT=<ID=LOG2FC,Number=1,Type=Float,Description="log2 fold change">\n']);
fprintf(fout,['##FORMAT=<ID=MEANBAF,Number=1,Type=Float,Description="mean B-allele frequency in segment">\n']);

segsTableCond=segsTable(1,:);
r=1;
for i=2:size(segsTable,1)
    if(segsTable.N(i)==2 && segsTable.M(i)==1 && sum(segsTableCond{r,[1 4 5]}==segsTable{i,[1 4 5]})==3)
        segsTableCond.EndPos(r)=segsTable.EndPos(i);   
    elseif(sum(segsTableCond{r,[1 4 5 7]}==segsTable{i,[1 4 5 7]})==4)
        segsTableCond.EndPos(r)=segsTable.EndPos(i);
    else
        r=r+1;
        segsTableCond(r,:)=segsTable(i,:);
    end
end

for i=1:length(exonRD)
    meanTumorRDexon(:,i)=getMeanInRegions(exonRD{i}(:,1:2),exonRD{i}(:,4),segsTableCond{:,1:3});
    meanNormalRDexon(:,i)=getMeanInRegions(exonRD{i}(:,1:2),exonRD{i}(:,5),segsTableCond{:,1:3});
    meanBAF(:,i)=getMeanInRegions([Tcell{i}.Chr(hetPos) Tcell{i}.Pos(hetPos)],Tcell{i}.BcountsComb(hetPos)./Tcell{i}.ReadDepthPass(hetPos),segsTableCond{:,1:3});
end

%%% calculate log2FC
for i=1:length(exonRD)
    log2FC(:,i)=log2((CNAscale(i)./2).*meanTumorRDexon(:,i)./(meanNormalRDexon(:,i)));
end

segsTableCond.log2FC=log2FC;
segsTableCond.meanBAF=meanBAF;

idx=getPosInRegions(exonRD{1}(:,1:2),segsTableCond{:,1:3});
exonCounts=hist(idx,1:size(segsTableCond,1));
idxHet=getPosInRegions([Tcell{1}.Chr(hetPos) Tcell{1}.Pos(hetPos)],segsTableCond{:,1:3});
hetCounts=hist(idxHet,1:size(segsTableCond,1));
segsTableCond.exonCounts=exonCounts';
segsTableCond.hetCounts=hetCounts';

%%% determine alteration type
type(segsTableCond.N>2 & segsTableCond.M>0,:)={'DUP'};
type(segsTableCond.N>2 & segsTableCond.M==0,:)={'DUPLOH'};
type(segsTableCond.N==2 & segsTableCond.M==1,:)={'NONE'};
type(segsTableCond.N==2 & segsTableCond.M==0,:)={'LOH'};
type(segsTableCond.N<2,:)={'DEL'};

%%% construct info field
Info=cellstr(strcat('CN=',num2str(segsTableCond.N,'%-.0f'),';MACN=',num2str(segsTableCond.M,'%-.0f')));
Info=strcat(Info,';SVLEN=',num2str(segsTableCond.EndPos-segsTableCond.StartPos,'%-.0f'));
Info=strcat(Info,';SVTYPE=',type,';END=',num2str(segsTableCond.EndPos,'%-.0f'));
Info=strcat(Info,';EXONCOUNT=',num2str(segsTableCond.exonCounts,'%-.0f'),';HETCOUNT=',num2str(segsTableCond.hetCounts,'%-.0f'));

for i=1:size(segsTableCond.F,2)
    formatStr(segsTableCond.N~=2 | segsTableCond.M~=1,i)=cellstr(strcat(num2str(segsTableCond.F(segsTableCond.N~=2 | segsTableCond.M~=1,i),'%-.3f'),':',num2str(segsTableCond.log2FC(segsTableCond.N~=2 | segsTableCond.M~=1,i),'%-.2f'),':',num2str(segsTableCond.meanBAF(segsTableCond.N~=2 | segsTableCond.M~=1,i),'%-.2f')));
    formatStr(segsTableCond.N==2 & segsTableCond.M==1,i)=cellstr(strcat('NA:',num2str(segsTableCond.log2FC(segsTableCond.N==2 & segsTableCond.M==1,i),'%-.2f'),':',num2str(segsTableCond.meanBAF(segsTableCond.N==2 & segsTableCond.M==1,i),'%-.2f')));
end
    
formatFields=repmat({'CNF:LOG2FC:MEANBAF'},size(segsTableCond,1),1);

%%% write output
sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];

outData=[chrList(segsTableCond.Chr) num2cell(segsTableCond.StartPos) num2cell(segsTableCond.EndPos) cellstr(char(ones(size(segsTableCond,1),1)*78)) strcat('<',regexprep(type,'DUPLOH','DUP'),'>') num2cell(abs(mean(segsTableCond.log2FC,2))) cellstr(char(ones(size(segsTableCond,1),1)*46)) Info formatFields formatStr];
headers={'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILT', 'INFO', 'FORMAT'};
headers=[headers regexp(inputParam.sampleNames,',','split')];
for i=1:length(headers)
    fprintf(fout,'%s\t',headers{i});
end
for i=1:size(outData,1)
    fprintf(fout,strcat('\n%s\t%d\t%d\t%s\t%s\t%f\t%s\t%s\t%s',repmat('\t%s',1,size(segsTableCond.F,2))),outData{i,:});
end
fclose(fout);


sampleNames=char(regexp(inputParam.sampleNames,',','split')');
sampleNamesShort=cellstr(sampleNames(:,1:min(namelengthmax-10,size(sampleNames,2))));
sampleNamesShort=regexprep(cellstr(sampleNamesShort),'-','_');
exonsOut=table();
exonsOut.Chr=exonRD{1}(:,1);
exonsOut.StartPos=exonRD{1}(:,2);
exonsOut.EndPos=exonRD{1}(:,3);
for i=1:length(exonRD)
    exonsOut.(strcat('meanDP_',sampleNamesShort{i}))=exonRD{i}(:,4);
    exonsOut.(strcat('log2FC_',sampleNamesShort{i}))=log2(exonRD{i}(:,4)./exonRD{i}(:,5));
end
exonsOut.meanDP_controls=exonRD{1}(:,5);
exonsOut.meanControlQC=-10*log10(exonRD{1}(:,6));
exonsOut.N=segsTableCond.N(idx);
exonsOut.M=segsTableCond.M(idx);
exonsOut.cloneId=segsTableCond.cnaIdx(idx);

writetable(exonsOut,[inputParam.outName '.exonData.tsv'],'FileType','text');
