function [Filter,somaticDetected]=writeJointVCF(Tcell,P,fIn,cloneId,F,inputParam)
%writeVCF - writes VCF for SNV and indel calls
%
% Syntax: writeVCF(T,pSomatic,posterior,pArtifact,pGermline,pHom,cloneId,f,W,inputParam,pDataSomatic,pDataHet,pDataHom,pCNA,F)
%
% Inputs:
%   T: table of position data with the following columns: 'Chr','Pos',
%       'ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ',
%       'AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ',
%       'BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount',
%       'ControlRD','PosMapQC','perReadPass','abFrac'
%   pSomatic: posterior probability somatic variant
%   pTrust: posterior probability call should be trusted
%   pArtifact: posterior probability artificat
%   pGermline: posterior probability germline heterozygous
%   pHom: posterior probability homozygous AA
%   cloneId: most likley clone assuming somatic
%   f: sample fraction of each clone
%   W: w parameter for each clone
%   inputParam: structure with all parameters
%   pDataSomatic: likliehood somatic
%   pDataHet: likliehood germline het
%   pCNA: posterior probability of copy number alteration
%   F: sample fraction containing copy number variant
%   
% Outputs:
%    writes a VCF file
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TumorOnlyWrapper, callSNV, qualDiscrim

% Author: Rebecca F. Halperin, PhD
% Translational Genomics Research Institute
% email: rhalperin@tgen.org
% Website: https://github.com/tgen
% Last revision: 3-June-2016
%------------- BEGIN CODE --------------

f=zeros(length(Tcell),inputParam.numClones);
tIdx=setdiff(1:length(Tcell),inputParam.NormalSample);
f(tIdx,1:end)=reshape(fIn,[],inputParam.numClones);

fout=fopen([inputParam.outName '.lumosVarSNV.vcf'],'w');

%%%print VCF header
fprintf(fout,'##fileformat=VCFv4.2\n');
fprintf(fout,['##fileData=' datestr(clock) '\n']);
inputFields=fieldnames(inputParam);
for i=1:length(inputFields)
    if(isnumeric(inputParam.(inputFields{i})))
        fprintf(fout,['##INPUT=<' inputFields{i} '=' mat2str(inputParam.(inputFields{i})') '>\n']);
    else
        fprintf(fout,['##INPUT=<' inputFields{i} '=' inputParam.(inputFields{i}) '>\n']);
    end
end
for i=1:size(f,2)
    outString=['##CloneID=' num2str(i)];
    for j=1:size(f,1)
        outString=[outString ',f' num2str(j) '=' num2str(f(j,i))];
    end
    outString=[outString ',PassCount=' num2str(sum(P.Somatic(:,1)>inputParam.pSomaticThresh & max(P.trust,[],2)>inputParam.pGoodThresh & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)<inputParam.maxSomPopFreq & cloneId(:,1)==i))];
    fprintf(fout,[outString '\n']);
end

fprintf(fout,['##INFO=<ID=A,Number=1,Type=String,Description="Major Allele Base">\n']);
fprintf(fout,['##INFO=<ID=B,Number=1,Type=String,Description="Minor Allele Base">\n']);
fprintf(fout,['##INFO=<ID=JPT,Number=1,Type=Float,Description="Phred Scaled Joint Posterior Probability the Call can be Trusted">\n']);
fprintf(fout,['##INFO=<ID=JPA,Number=1,Type=Float,Description="Phred Scaled Joint Posterior Probability the Position is an Artifact">\n']);
fprintf(fout,['##INFO=<ID=JPS,Number=1,Type=Float,Description="Joint Posterior Probability of Somatic Mutation">\n']);
fprintf(fout,['##INFO=<ID=JPGAB,Number=1,Type=Float,Description="Joint Posterior Probability of No Somatic Mutation and Position is Germline AB">\n']);
fprintf(fout,['##INFO=<ID=JPGAA,Number=1,Type=Float,Description="Joint Posterior Probability of No Somatic Mutation and Position is Germline AA">\n']);
fprintf(fout,['##INFO=<ID=JPGND,Number=1,Type=Float,Description="Joint Posterior Probability of Variant Present in Germline Not Following Diploid Model">\n']);
fprintf(fout,['##INFO=<ID=ApopAF,Number=1,Type=Float,Description="population frequency of A allele">\n']);
fprintf(fout,['##INFO=<ID=BpopAF,Number=1,Type=Float,Description="population frequency of B allele">\n']);
fprintf(fout,['##INFO=<ID=CosmicCount,Number=1,Type=Float,Description="Count of mutations at position in COSMIC>\n']);
fprintf(fout,['##INFO=<ID=CloneId,Number=1,Type=Integer,Description="CloneId">\n']);
fprintf(fout,['##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy Number">\n']);
fprintf(fout,['##INFO=<ID=MACN,Number=1,Type=Integer,Description="Min Allele Copy Number">\n']);

%fprintf(fout,['##INFO=<ID=PCNA,Number=1,Type=Float,Description="Posterior Probability Copy Number is Correct">\n']);
%for i=1:size(F{1},2)
%    fprintf(fout,'%s\n',char(strcat('##INFO=<ID=',F{1}.Properties.VariableNames(i),',Number=1,Type=Float,Description="',F{1}.Properties.VariableDescriptions(i))));
%end
fprintf(fout,['##FILTER=<ID=SomaticPASS,Description="JPS>pSomaticThresh and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=SomaticLowQC,Description="JPS>0.5 and artifact filters">\n']);
fprintf(fout,['##FILTER=<ID=SomaticPairPASS,Description="PPS>pSomaticThresh and JPS<0.5 pass filters">\n']);
fprintf(fout,['##FILTER=<ID=SomaticPairLowQC,Description="PPS>0.5 and JPS<0.5 and artifact filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineHetPASS,Description="JPGAB>pGermlineThresh and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineHetLowQC,Description="JPGAB>0.5 and artifact filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineHomPASS,Description="JPGAA>pGermlineThresh and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineHomLowQC,Description="JPGAA>0.5 and artifact filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineShiftPASS,Description="JPND>pGermlineThresh and pass filters">\n']);
fprintf(fout,['##FILTER=<ID=GermlineShiftLowQC,Description="JPND>0.5 and artifact filters">\n']);
fprintf(fout,['##FILTER=<ID=NoCall,Description="PG<0.5 and PS<0.5 and pass filters>\n']);
fprintf(fout,['##FILTER=<ID=REJECT,Description="Below at least one QC filter">\n']);

%GT:DP:AD:FT:PPS:PT:PA:PLS:PL:PLND:VSF:CNF
fprintf(fout,['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n']);
fprintf(fout,['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of Reads with MQ and BQ greater than min">\n']);
fprintf(fout,['##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depths">\n']);
fprintf(fout,['##FORMAT=<ID=FT,Number=.,Type=String,Description="Filters">\n']);
fprintf(fout,['##FORMAT=<ID=PPS,Number=1,Type=Float,Description="Phred Scaled Posterior Probability Position is Somatic Based on Paired Analysis with Control">\n']);
fprintf(fout,['##FORMAT=<ID=PT,Number=1,Type=Float,Description="Phred Scaled Posterior Probability Position is Trusted">\n']);
fprintf(fout,['##FORMAT=<ID=PA,Number=1,Type=Float,Description="Phred Scaled Posterior Probability Position is Artifact">\n']);
fprintf(fout,['##FORMAT=<ID=PLS,Number=1,Type=Float,Description="Phred Scaled Likelihood Somatic">\n']);
fprintf(fout,['##FORMAT=<ID=PL,Number=G,Type=Float,Description="Phred Scaled Likelihood of Germline Genotypes">\n']);
fprintf(fout,['##FORMAT=<ID=PLND,Number=1,Type=Float,Description="Phred Scaled Likelihood of NonDiploid">\n']);
fprintf(fout,['##FORMAT=<ID=VSF,Number=1,Type=Float,Description="Fraction of Sample Containing Variant">\n']);
fprintf(fout,['##FORMAT=<ID=CNF,Number=1,Type=Float,Description="Fraction containg Copy Number Alteration">\n']);

%%% get reference and alternate bases
T=Tcell{1};
RefNT=int2ntIndels(T.RefComb);
RefNT(strcmp(RefNT,'<LongIndel>'),:)={'N'};
AltNT(T.RefComb~=T.Acomb & T.RefComb~=T.Bcomb,:)=cellstr(strcat(int2ntIndels(T.Acomb(T.RefComb~=T.Acomb & T.RefComb~=T.Bcomb)),',', int2ntIndels(T.Bcomb(T.RefComb~=T.Acomb & T.RefComb~=T.Bcomb))));
AltNT(T.RefComb==T.Acomb & P.Hom(:,1)<=0.5,:)=cellstr(int2ntIndels(T.Bcomb(T.RefComb==T.Acomb & P.Hom(:,1)<=0.5)));
AltNT(T.RefComb==T.Bcomb,:)=cellstr(int2ntIndels(T.Acomb(T.RefComb==T.Bcomb)));
AltNT(T.RefComb==T.Acomb & P.Hom(:,1)>0.5,:)={'.'};
AltNT(strcmp(RefNT,'<LongIndel>'),:)={'<DEL>'};
AltNT(strcmp(AltNT,'<LongIndel>'),:)={'<INS>'};
RefNT(~cellfun('isempty',strfind(AltNT,'N')) & T.RefComb>4,:)={'N'};
AltNT(~cellfun('isempty',strfind(AltNT,'N')) & T.RefComb>4,:)={'<DEL>'};


%%% calculate quality
Qual=-10*log10(1-mean([mean(P.trust,2) mean(1-P.artifact,2) max([P.Het P.Hom P.Somatic P.NonDip],[],2)],2));
%Qual=zeros(size(T.Pos));
%Qual(P.Somatic(:,n)>0.5,:)=min(min(-10*log10(1-P.Somatic(P.Somatic(:,n)>0.5,n)),-10*log10(1-P.trust(P.Somatic(:,n)>0.5,n))),-10*log10(P.artifact(P.Somatic(:,n)>0.5,n)));
%Qual(P.Het(:,n)>0.5,:)=min(min(-10*log10(1-P.Het(P.Het(:,n)>0.5,n)),-10*log10(1-P.trust(P.Het(:,n)>0.5,n))),-10*log10(P.artifact(P.Het(:,n)>0.5,n)));
%Qual(P.Hom(:,n)>0.5,:)=min(min(-10*log10(1-P.Hom(P.Hom(:,n)>0.5,n)),-10*log10(1-P.trust(P.Hom(:,n)>0.5,n))),-10*log10(P.artifact(P.Hom(:,n)>0.5,n)));
%Qual(P.NonDip(:,n)>0.5,:)=min(min(-10*log10(1-P.NonDip(P.NonDip(:,n)>0.5,n)),-10*log10(1-P.trust(P.NonDip(:,n)>0.5,n))),-10*log10(P.artifact(P.NonDip(:,n)>0.5,n)));

%%% assign filters
passPos=max(P.trust,[],2)>inputParam.pGoodThresh & T.RefComb>0 & T.Acomb>0 & T.Bcomb>0;
Filter(P.Somatic(:,1)>0.5,:)={'SomaticLowQC'};
Filter(P.Somatic(:,1)>inputParam.pSomaticThresh & passPos,:)={'SomaticPASS'};
Filter(max([P.Somatic(:,1) P.SomaticPair],[],2)>0.5 & min([T.ApopAF T.BpopAF],[],2)>=inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
Filter(P.Het(:,1)>0.5,:)={'GermlineHetLowQC'};
Filter(P.Het(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHetPASS'};
Filter(P.Hom(:,1)>0.5,:)={'GermlineHomLowQC'};
Filter(P.Hom(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHomPASS'};
Filter(P.NonDip(:,1)>0.5,:)={'GermlineShiftLowQC'};
Filter(P.NonDip(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineShiftPASS'};
Filter(P.Somatic(:,1)<0.5 & P.Het(:,1)<0.5 & P.Hom(:,1)<0.5 & P.NonDip(:,1)<0.5,:)={'NoCall'};
Filter(max(P.artifact,[],2)>inputParam.pGoodThresh,:)={'REJECT'};
idxSom=strncmp(Filter,'Somatic',7);
if inputParam.NormalSample>0
    for i=1:tIdx
        Filter(P.SomaticPair(:,tIdx(i))>0.5 & ~idxSom & max(P.artifact(:,[tIdx(i) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)={'SomaticPairLowQC'};
        Filter(max(P.SomaticPair(:,tIdx(i)),[],2)>inputParam.pSomaticThresh & ~idxSom & min(P.trust(:,[tIdx(i) inputParam.NormalSample]),[],2)>inputParam.pGoodThresh,:)={'SomaticPairPASS'};
    end
end

%%% construct info fields
%Info=cellstr(strcat('DP=',num2str(T.ReadDepth,'%-.0f'),';DPQC=',num2str(T.ReadDepthPass,'%-.0f')));
%Info(T.RefComb==T.Acomb,:)=strcat(Info(T.RefComb==T.Acomb,:),';AF=',num2str(T.BcountsComb(T.RefComb==T.Acomb,:)./T.ReadDepthPass(T.RefComb==T.Acomb,:),'%-.3f'));
%Info(T.RefComb==T.Bcomb,:)=strcat(Info(T.RefComb==T.Bcomb,:),';AF=',num2str(T.AcountsComb(T.RefComb==T.Bcomb,:)./T.ReadDepthPass(T.RefComb==T.Bcomb,:),'%-.3f'));
Info=cellstr(strcat('JPT=',num2str(-10*log10(1-mean(P.trust,2)),'%-.1f'),';JPA=',num2str(-10*log10(1-mean(P.artifact,2)),'%-.1f'),';JPS=',num2str(P.Somatic(:,1),'%-.5f'),';JPGAB=',num2str(P.Het(:,1),'%-.5f'),';JPGAA=',num2str(P.Hom(:,1),'%-.5f'),';JPND=',num2str(P.NonDip(:,1),'%-.5f')));
%Info=strcat(Info,';LS=',num2str(P.DataSomatic(:,n),'%-.5f'),';LGAB=',num2str(P.DataHet(:,n),'%-.5f'),';LGAA=',num2str(P.DataHom(:,n),'%-.5f'),';LND=',num2str(P.DataNonDip(:,n),'%-.5f'));
Info=strcat(Info,';A=',int2ntIndels(T.Acomb),';B=',int2ntIndels(T.Bcomb));
%Info=strcat(Info,';ACountF=',num2str(T.ACountF,'%-.0f'),';ACountR=',num2str(T.ACountR,'%-.0f'),';BCountF=',num2str(T.BCountF,'%-.0f'),';BCountR=',num2str(T.BCountR,'%-.0f'));
Info=strcat(Info,';ApopAF=',num2str(T.ApopAF,'%-.5f'),';BpopAF=',num2str(T.BpopAF,'%-.5f'),';CosmicCount=',num2str(T.CosmicCount,'%-.0f'));
Info=strcat(Info,';CloneId=',num2str(cloneId(:,1)));
Info(T.NumCopies~=2 | T.MinAlCopies~=1)=strcat(Info(T.NumCopies~=2 | T.MinAlCopies~=1),';CN=',num2str(T.NumCopies(T.NumCopies~=2 | T.MinAlCopies~=1)),';MACN=',num2str(T.MinAlCopies(T.NumCopies~=2 | T.MinAlCopies~=1)));
%Info=strcat(Info,';PCNA=',num2str(pCNA,'%-.5f'));
%for i=1:size(F,2)
%    Info=strcat(Info,';',F.Properties.VariableNames(i),'=',num2str(F{:,i},'%-.3f'));
%end
copyList=unique([T.NumCopies T.MinAlCopies],'rows');
tumorGT=cell(1,size(T,1));
for i=1:size(copyList,1)
    idx=T.NumCopies==copyList(i,1) & T.MinAlCopies==copyList(i,2);
    tumorGT(P.Hom(:,1)>0.5 & T.Acomb==T.Ref & idx)={strjoin(cellstr(num2str(zeros(copyList(i,1),1))),'/')};
    tumorGT(P.Hom(:,1)>0.5 & T.Bcomb==T.Ref & idx)={strjoin(cellstr(num2str(ones(copyList(i,1),1))),'/')};
    tumorGT(P.Het(:,1)>0.5 & T.Acomb==T.Ref & idx)={strjoin(cellstr(num2str([zeros(copyList(i,1)-copyList(i,2),1); ones(copyList(i,2),1)])),'/')};
    tumorGT(P.Het(:,1)>0.5 & T.Bcomb==T.Ref & idx)={strjoin(cellstr(num2str([zeros(copyList(i,2),1); ones(copyList(i,1)-copyList(i,2),1)])),'/')};
    tumorGT(P.Somatic(:,1)>0.5 & T.cnaF==f(1,cloneId(:,1))' & idx)={strjoin(cellstr(num2str([zeros(copyList(i,2),1); ones(copyList(i,1)-copyList(i,2),1)])),'/')};
end
tumorGT(P.Somatic(:,1)>0.5 & T.cnaF~=f(1,cloneId(:,1))')={'0/1'};
tumorGT(cellfun('isempty',tumorGT))={'.'};

germGT=cell(1,size(T,1));
germGT(P.Hom(:,1)>0.5 & T.Acomb==T.Ref)={'0/0'};
germGT(P.Hom(:,1)>0.5 & T.Acomb~=T.Ref)={'1/1'};
germGT(P.Het(:,1)>0.5)={'0/1'};
germGT(P.Somatic(:,1)>0.5)={'0/0'};
germGT(cellfun('isempty',germGT))={'.'};



for j=1:length(Tcell)
    T=Tcell{j};
    bIdx=T.ApopAFcomb>=T.BpopAFcomb;
    AF(bIdx,j)=T.BcountsComb(bIdx)./T.ReadDepthPass(bIdx);
    AF(~bIdx,j)=T.AcountsComb(~bIdx)./T.ReadDepthPass(~bIdx);
    matchIdx=f(j,cloneId(:,j))'==T.cnaF;
    sampleFrac(matchIdx,j)=(2*AF(matchIdx,j))./(T.NumCopies(matchIdx)-T.MinAlCopies(matchIdx)+2*AF(matchIdx,j)-AF(matchIdx,j).*T.NumCopies(matchIdx));
    sampleFrac(~matchIdx,j)=AF(~matchIdx,j)./(T.cnaF(~matchIdx).*T.NumCopies(~matchIdx)+2.*(1-T.cnaF(~matchIdx)));
end
    

somaticDetected=zeros(size(AF));
if inputParam.NormalSample<1
    formatStr(:,1)=strcat(germGT',':.:.:.:.:.:.:.:.:.:.:.');
    n=2;
else
    n=1;
end
for i=1:length(Tcell)
    T=Tcell{i};
    if inputParam.NormalSample==i
        formatStr(:,n)=germGT';
    else
%        n
 %       size(formatStr)
 %       size(tumorGT)
        formatStr(:,n)=tumorGT';
    end
    formatStr(:,n)=strcat(formatStr(:,n),':',num2str(T.ReadDepthPass,'%-.0f'));
    aIdx=T.RefComb==T.Acomb;
    formatStr(aIdx,n)=strcat(formatStr(aIdx,n),':',num2str(T.AcountsComb(aIdx),'%-.0f'),',',num2str(T.BcountsComb(aIdx),'%-.0f'));
    bIdx=T.RefComb==T.Bcomb;
    formatStr(bIdx,n)=strcat(formatStr(bIdx,n),':',num2str(T.BcountsComb(bIdx),'%-.0f'),',',num2str(T.AcountsComb(bIdx),'%-.0f'));
    formatStr(~aIdx & ~bIdx,n)=strcat(formatStr(~aIdx & ~bIdx,n),':NA,',num2str(T.AcountsComb(~aIdx & ~bIdx),'%-.0f'),',',num2str(T.BcountsComb(~aIdx & ~bIdx),'%-.0f'));
    filtStr(P.trust(:,i)>=inputParam.pGoodThresh)={'PASS'};
    filtStr(P.trust(:,i)<inputParam.pGoodThresh & P.artifact(:,i)<inputParam.pGoodThresh)={'LowQC'};
    filtStr(P.artifact(:,i)>=inputParam.pGoodThresh)={'REJECT'};
    if inputParam.NormalSample<1
        filtStr(P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb))=strcat(filtStr(P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb)),';SomaticDetected');
        somaticDetected(P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb),i)=1;
        filtStr(P.Somatic(:,i)>0.5 & (P.DataSomatic(:,i)<=P.DataHom(:,i) | (T.BcountsComb==0 & T.A==T.RefComb)))=strcat(filtStr(P.Somatic(:,i)>0.5 & (P.DataSomatic(:,i)<=P.DataHom(:,i) | (T.BcountsComb==0 & T.A==T.RefComb))),';SomaticNotDetected');
        formatStr(:,n)=strcat(formatStr(:,n),':',filtStr',':NA');
    else
        filtStr(P.SomaticPair(:,i)>0.5)=strcat(filtStr(P.SomaticPair(:,i)>0.5),';SomaticDetected');
        somaticDetected(P.SomaticPair(:,i)>0.5,i)=1;
        filtStr(P.Somatic(:,i)>0.5 & P.SomaticPair(:,i)<=0.5)=strcat(filtStr(P.Somatic(:,i)>0.5 & P.SomaticPair(:,i)<=0.5),';SomaticNotDetected');
        formatStr(:,n)=strcat(formatStr(:,n),':',filtStr',':',num2str(P.SomaticPair(:,i),'%-.3f'));
    end
    formatStr(:,n)=strcat(formatStr(:,n),':',num2str(-10*log10(1-P.trust(:,i)),'%-.0f'),':',num2str(-10*log10(1-P.artifact(:,i)),'%-.0f'),':',num2str(P.DataSomatic(:,i),'%-.0f'));
    formatStr(aIdx,n)=strcat(formatStr(aIdx,n),':',num2str(-10*log10(P.DataHom(aIdx,i)),'%-.0f'),',',num2str(-10*log10(P.DataHet(aIdx,i)),'%-.0f'),',NA');
    formatStr(bIdx,n)=strcat(formatStr(bIdx,n),':NA,',num2str(-10*log10(P.DataHet(bIdx,i)),'%-.0f'),',',num2str(-10*log10(P.DataHom(bIdx,i)),'%-.0f'));
    formatStr(~aIdx & ~bIdx,n)=strcat(formatStr(~aIdx & ~bIdx,n),':.');
    if inputParam.NormalSample>0
        formatStr(:,n)=strcat(formatStr(:,n),':',num2str(-10*log10(P.DataNonDip(:,i)),'%-.0f'));
    else
        formatStr(:,n)=strcat(formatStr(:,n),':NA');
    end
    formatStr(P.Somatic(:,i)>0.5,n)=strcat(formatStr(P.Somatic(:,i)>0.5,n),':',num2str(sampleFrac(P.Somatic(:,i)>0.5,i),'%-.3f'));
    formatStr(P.Somatic(:,i)<=0.5,n)=strcat(formatStr(P.Somatic(:,i)<0.5,n),':NA');
    formatStr(T.NumCopies==2 & T.MinAlCopies==1,n)=strcat(formatStr(T.NumCopies==2 & T.MinAlCopies==1,n),':NA');
    formatStr(T.NumCopies~=2 | T.MinAlCopies~=1,n)=strcat(formatStr(T.NumCopies~=2 | T.MinAlCopies~=1,n),':',num2str(T.cnaF(T.NumCopies~=2 | T.MinAlCopies~=1),'%-.3f'));
    n=n+1;
end


formatFields=repmat({'GT:DP:AD:FT:PPS:PT:PA:PLS:PL:PLND:VSF:CNF'},size(P,1),1);

%%% print output
outData=[num2cell(T.Chr) num2cell(T.Pos) cellstr(char(ones(size(T.Pos))*46)) cellstr(RefNT) squeeze(AltNT) num2cell(Qual) Filter Info formatFields formatStr];
headers={'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILT', 'INFO', 'FORMAT'};
if inputParam.NormalSample<1
    headers=[headers 'InferredGermline' regexp(inputParam.sampleNames,',','split')];
else
    headers=[headers regexp(inputParam.sampleNames,',','split')];
end
for i=1:length(headers)
    fprintf(fout,'%s\t',headers{i});
end

for i=1:size(outData,1)
    fprintf(fout,strcat('\n%d\t%d\t%s\t%s\t%s\t%f\t%s\t%s\t%s',repmat('\t%s',1,n)),outData{i,:});
end

fclose(fout);
