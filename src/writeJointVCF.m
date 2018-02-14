function sampleFrac=writeJointVCF(Tcell,P,fIn,cloneId,alleleId,Filter,somaticDetected,trustScore,artifactScore,inputParam)
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
%RefNT(strcmp(RefNT,'<LongIndel>'),:)={'N'};
delPos=cellfun('length',RefNT)>1;
longDelPos=strncmp(RefNT,'LEN',3);
svLen=zeros(height(T),1);
svLen(delPos)=1-cellfun('length',RefNT(delPos));
svLen(longDelPos)=-1*str2double(regexprep(RefNT(longDelPos),'LEN_',''));
endPos=T.Pos-svLen;
RefNT(longDelPos)={'N'};
AltNT=repmat({'.'},height(T),1);
AltNT(T.RefComb~=T.Acomb & T.RefComb~=T.Bcomb,:)=cellstr(strcat(int2ntIndels(T.Acomb(T.RefComb~=T.Acomb & T.RefComb~=T.Bcomb)),',', int2ntIndels(T.Bcomb(T.RefComb~=T.Acomb & T.RefComb~=T.Bcomb))));
AltNT(T.RefComb==T.Acomb & P.Hom(:,1)<=0.5,:)=cellstr(int2ntIndels(T.Bcomb(T.RefComb==T.Acomb & P.Hom(:,1)<=0.5)));
AltNT(T.RefComb==T.Bcomb,:)=cellstr(int2ntIndels(T.Acomb(T.RefComb==T.Bcomb)));
AltNT(T.RefComb==T.Acomb & P.Hom(:,1)>0.5,:)={'.'};
insPos=cellfun('length',AltNT)>1;
AltNT(longDelPos,:)={'<DEL>'};
svLen(insPos)=cellfun(@(x) max(cellfun('length',strsplit(x,','))),AltNT(insPos))-1;
longInsPos=strncmp(AltNT,'LEN',3);
svLen(longInsPos)=str2double(regexprep(AltNT(longInsPos),'LEN_',''));
AltNT(longInsPos)={'<INS>'};
%RefNT(~cellfun('isempty',strfind(AltNT,'N')) & T.RefComb>4,:)={'N'};
%AltNT(~cellfun('isempty',strfind(AltNT,'N')) & T.RefComb>4,:)={'<DEL>'};


%%% calculate quality
Qual=-10*log10(1-mean([trustScore artifactScore max([P.Het P.Hom P.Somatic P.NonDip],[],2)],2));
%Qual=zeros(size(T.Pos));
%Qual(P.Somatic(:,n)>0.5,:)=min(min(-10*log10(1-P.Somatic(P.Somatic(:,n)>0.5,n)),-10*log10(1-P.trust(P.Somatic(:,n)>0.5,n))),-10*log10(P.artifact(P.Somatic(:,n)>0.5,n)));
%Qual(P.Het(:,n)>0.5,:)=min(min(-10*log10(1-P.Het(P.Het(:,n)>0.5,n)),-10*log10(1-P.trust(P.Het(:,n)>0.5,n))),-10*log10(P.artifact(P.Het(:,n)>0.5,n)));
%Qual(P.Hom(:,n)>0.5,:)=min(min(-10*log10(1-P.Hom(P.Hom(:,n)>0.5,n)),-10*log10(1-P.trust(P.Hom(:,n)>0.5,n))),-10*log10(P.artifact(P.Hom(:,n)>0.5,n)));
%Qual(P.NonDip(:,n)>0.5,:)=min(min(-10*log10(1-P.NonDip(P.NonDip(:,n)>0.5,n)),-10*log10(1-P.trust(P.NonDip(:,n)>0.5,n))),-10*log10(P.artifact(P.NonDip(:,n)>0.5,n)));

%%% assign filters
% passPos=max(P.trust,[],2)>inputParam.pGoodThresh & max(P.artifact,[],2)<inputParam.pGoodThresh & Tcell{1}.RefComb>0 & Tcell{1}.Acomb>0 & Tcell{1}.Bcomb>0;
% Filter=cell(size(passPos));
% Filter(P.Somatic(:,1)>0.5,:)={'SomaticLowQC'};
% Filter(P.Somatic(:,1)>inputParam.pSomaticThresh & passPos,:)={'SomaticPASS'};
% Filter(P.Somatic(:,1)>0.5 & min([Tcell{1}.ApopAF Tcell{1}.BpopAF],[],2)>inputParam.maxSomPopFreq,:)={'SomaticDBsnp'};
% Filter(P.Het(:,1)>0.5,:)={'GermlineHetLowQC'};
% Filter(P.Het(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHetPASS'};
% Filter(P.Hom(:,1)>0.5,:)={'GermlineHomLowQC'};
% Filter(P.Hom(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineHomPASS'};
% Filter(P.NonDip(:,1)>0.5,:)={'GermlineShiftLowQC'};
% Filter(P.NonDip(:,1)>inputParam.pGermlineThresh & passPos,:)={'GermlineShiftPASS'};
% Filter(P.Somatic(:,1)<0.5 & P.Het(:,1)<0.5 & P.Hom(:,1)<0.5 & P.NonDip(:,1)<0.5,:)={'NoCall'};
% Filter(min(P.artifact,[],2)>inputParam.pGoodThresh,:)={'REJECT'};
% %tIdx=setdiff(1:length(Tcell),inputParam.NormalSample)
% if inputParam.NormalSample>0
%     idxSom=strncmp(Filter,'Somatic',7);
%     idxSomPass=strncmp(Filter,'SomaticPASS',11);
%     for i=1:length(tIdx)
%         Filter(P.SomaticPair(:,tIdx(i))>0.5 & ~idxSom & min(P.artifact(:,[tIdx(i) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)={'SomaticPairLowQC'};
%         Filter(max(P.SomaticPair(:,tIdx(i)),[],2)>inputParam.pSomaticThresh & ~idxSomPass & max(P.trust(:,[tIdx(i) inputParam.NormalSample]),[],2)>inputParam.pGoodThresh &  max(P.artifact(:,[tIdx(i) inputParam.NormalSample]),[],2)<inputParam.pGoodThresh,:)={'SomaticPairPASS'};
%     end
% end

AltNT(strncmp('SomaticPair',Filter,10) & T.RefComb==T.Acomb,:)=int2ntIndels(T.Bcomb(strncmp('SomaticPair',Filter,10) & T.RefComb==T.Acomb,:));
AltNT(strncmp('SomaticPair',Filter,10) & T.RefComb~=T.Acomb,:)=int2ntIndels(T.Acomb(strncmp('SomaticPair',Filter,10) & T.RefComb~=T.Acomb,:));

%%% construct info fields
%Info=cellstr(strcat('DP=',num2str(T.ReadDepth,'%-.0f'),';DPQC=',num2str(T.ReadDepthPass,'%-.0f')));
%Info(T.RefComb==T.Acomb,:)=strcat(Info(T.RefComb==T.Acomb,:),';AF=',num2str(T.BcountsComb(T.RefComb==T.Acomb,:)./T.ReadDepthPass(T.RefComb==T.Acomb,:),'%-.3f'));
%Info(T.RefComb==T.Bcomb,:)=strcat(Info(T.RefComb==T.Bcomb,:),';AF=',num2str(T.AcountsComb(T.RefComb==T.Bcomb,:)./T.ReadDepthPass(T.RefComb==T.Bcomb,:),'%-.3f'));
%Info=cellstr(strcat('JPT=',num2str(real(-10*log10(1-mean(P.trust,2))),'%-.1f'),';JPA=',num2str(real(-10*log10(1-mean(P.artifact,2))),'%-.1f'),';JPS=',num2str(P.Somatic(:,1),'%-.5f'),';JPGAB=',num2str(P.Het(:,1),'%-.5f'),';JPGAA=',num2str(P.Hom(:,1),'%-.5f'),';JPND=',num2str(P.NonDip(:,1),'%-.5f')));
Info=cellstr(strcat('JPT=',strsplit(sprintf('%-.1f\n',real(-10*log10(1-trustScore))))',';JPA=',strsplit(sprintf('%-.1f\n',real(-10*log10(1-artifactScore))))',';JPS=',strsplit(sprintf('%-.5f\n',P.Somatic(:,1)))',';JPGAB=',strsplit(sprintf('%-.5f\n',P.Het(:,1)))',';JPGAA=',strsplit(sprintf('%-.5f\n',P.Hom(:,1)))',';JPND=',strsplit(sprintf('%-.5f\n',P.NonDip(:,1)))'));
%Info=strcat(Info,';LS=',num2str(P.DataSomatic(:,n),'%-.5f'),';LGAB=',num2str(P.DataHet(:,n),'%-.5f'),';LGAA=',num2str(P.DataHom(:,n),'%-.5f'),';LND=',num2str(P.DataNonDip(:,n),'%-.5f'));
Info=strcat(Info,';A=',[int2ntIndels(T.Acomb);{''}],';B=',[int2ntIndels(T.Bcomb); {''}]);
%Info=strcat(Info,';ACountF=',num2str(T.ACountF,'%-.0f'),';ACountR=',num2str(T.ACountR,'%-.0f'),';BCountF=',num2str(T.BCountF,'%-.0f'),';BCountR=',num2str(T.BCountR,'%-.0f'));
Info=strcat(Info,';ApopAF=',strsplit(sprintf('%-.5f\n',T.ApopAF))',';BpopAF=',strsplit(sprintf('%-.5f\n',T.BpopAF))',';CosmicCount=',strsplit(sprintf('%-.0f\n',T.CosmicCount))');
Info=strcat(Info,';CloneId=',strsplit(sprintf('%d\n',cloneId(:,1)))');
Info([delPos; true])=strcat(Info([delPos; true]),';END=',strsplit(sprintf('%d\n',endPos(delPos)))');
Info([insPos; true])=strcat(Info([insPos; true]),';SVLEN=',strsplit(sprintf('%d\n',svLen(insPos)))');
cnPos=T.NumCopies~=2 | T.MinAlCopies~=1;
Info([cnPos; true])=strcat(Info([cnPos; true]),';CN=',strsplit(sprintf('%d\n',T.NumCopies(cnPos)))',';MACN=',strsplit(sprintf('%d\n',T.MinAlCopies(cnPos)))');
Info=Info(1:height(T));
%Info=strcat(Info,';PCNA=',num2str(pCNA,'%-.5f'));
%for i=1:size(F,2)
%    Info=strcat(Info,';',F.Properties.VariableNames(i),'=',num2str(F{:,i},'%-.3f'));
%end

gt=repmat('.',height(T),2);
if inputParam.NormalSample>0
    bIdx=Tcell{inputParam.NormalSample}.AcountsComb>=Tcell{inputParam.NormalSample}.BcountsComb;
    aIdx=Tcell{inputParam.NormalSample}.AcountsComb<Tcell{inputParam.NormalSample}.BcountsComb;    
else
    bIdx=Tcell{1}.ApopAFcomb>=Tcell{1}.BpopAFcomb;
    aIdx=Tcell{1}.ApopAFcomb<Tcell{1}.BpopAFcomb;
end
gt(bIdx & T.Acomb==T.RefComb,1)='0';
gt(bIdx & T.Acomb~=T.RefComb,1)='1';
gt(aIdx & T.Bcomb==T.RefComb,1)='0';
gt(aIdx & T.Bcomb~=T.RefComb & T.Acomb==T.RefComb,1)='1';
gt(aIdx & T.Acomb~=T.RefComb & T.Bcomb~=T.RefComb,1)='2';
gt(bIdx & T.Acomb==T.RefComb,2)='1';
gt(bIdx & T.Bcomb==T.RefComb,2)='0';
gt(bIdx & T.Acomb~=T.RefComb & T.Bcomb~=T.RefComb,2)='2';
gt(aIdx & T.Acomb==T.RefComb,2)='0';
gt(aIdx & T.Bcomb==T.RefComb,2)='1';
gt(aIdx & T.Acomb~=T.RefComb & T.Bcomb~=T.RefComb,2)='1';

copyList=unique([T.NumCopies T.MinAlCopies],'rows');
tumorGT=cell(1,size(T,1));
T=Tcell{tIdx(1)};
for i=1:size(copyList,1)
    idx=T.NumCopies==copyList(i,1) & T.MinAlCopies==copyList(i,2);
    currSomIdx=strncmp(Filter,'Somatic',7) & T.cnaF==f(tIdx(1),cloneId(:,1))' & idx;
    if copyList(i,1)>1
        tumorGT(P.Hom(:,1)>0.5 & idx)=cellstr(repmat(gt(P.Hom(:,1)>0.5 & idx,1),1,copyList(i,1)));
        tumorGT(currSomIdx & alleleId==1)=cellstr(sort([repmat(gt(currSomIdx & alleleId==1,1),1,copyList(i,1)-copyList(i,2)) repmat(gt(currSomIdx & alleleId==1,2),1,copyList(i,2))],2));
        tumorGT(currSomIdx & alleleId==2)=cellstr(sort([repmat(gt(currSomIdx & alleleId==2,1),1,copyList(i,2)) repmat(gt(currSomIdx & alleleId==2,2),1,copyList(i,1)-copyList(i,2))],2));
    elseif copyList(i,1)==1
        tumorGT(P.Hom(:,1)>0.5 & idx)=cellstr(gt(P.Hom(:,1)>0.5 & idx,1));
        tumorGT(currSomIdx)=cellstr(gt(currSomIdx,2));
    else
        tumorGT(idx)={'.'};
    end            
    %tumorGT(P.Hom(:,1)>0.5 & T.Acomb==T.Ref & idx)={strjoin(cellstr(num2str(zeros(copyList(i,1),1))),'/')};
    %tumorGT(P.Hom(:,1)>0.5 & T.Bcomb==T.Ref & idx)={strjoin(cellstr(num2str(ones(copyList(i,1),1))),'/')};
    tumorGT(P.Het(:,1)>0.5 & T.Acomb==T.Ref & idx)={strjoin(cellstr(num2str([zeros(copyList(i,1)-copyList(i,2),1); ones(copyList(i,2),1)])),'')};
    tumorGT(P.Het(:,1)>0.5 & T.Bcomb==T.Ref & idx)={strjoin(cellstr(num2str([zeros(copyList(i,2),1); ones(copyList(i,1)-copyList(i,2),1)])),'')};
    tumorGT(P.Het(:,1)>0.5 & T.Acomb~=T.Ref & T.Bcomb~=T.Ref & idx)={strjoin(cellstr(num2str([ones(copyList(i,1)-copyList(i,2),1); 2*ones(copyList(i,2),1)])),'')};
%     tumorGT(strncmp(Filter,'Somatic',7) & T.cnaF==f(1,cloneId(:,1))' & idx)={strjoin(cellstr(num2str([zeros(copyList(i,2),1); ones(copyList(i,1)-copyList(i,2),1)])),'/')};
end
currSomIdx=strncmp(Filter,'Somatic',7) & T.cnaF~=f(tIdx(1),cloneId(:,1))';
tumorGT(currSomIdx)=cellstr(sort(gt(currSomIdx,:),2));
tumorGT(cellfun('isempty',tumorGT))={'.'};
tumorGT=regexprep(tumorGT,'([0-9])','$1\');
tumorGT=regexprep(tumorGT,'\\$','');


germGT=cell(1,size(T,1));
germGT(P.Hom(:,1)>0.5 | P.Somatic(:,1)>0.5)=cellstr(repmat(gt(P.Hom(:,1)>0.5 | P.Somatic(:,1)>0.5,1),1,2));
germGT(P.Het(:,1)>0.5)=cellstr(sort(gt(P.Het(:,1)>0.5,:),2));
% germGT(P.Hom(:,1)>0.5 & T.Acomb==T.Ref)={'0/0'};
% germGT(P.Hom(:,1)>0.5 & T.Acomb~=T.Ref)={'1/1'};
% germGT(P.Het(:,1)>0.5)={'0/1'};
% germGT(P.Somatic(:,1)>0.5)={'0/0'};
germGT(cellfun('isempty',germGT))={'.'};
germGT=regexprep(germGT,'([0-9])','$1\');
germGT=regexprep(germGT,'\\$','');





if inputParam.NormalSample>0
    bIdx=Tcell{inputParam.NormalSample}.AcountsComb>=Tcell{inputParam.NormalSample}.BcountsComb;
else
    bIdx=T.ApopAFcomb>=T.BpopAFcomb;
end

AF=NaN(height(Tcell{1}),length(Tcell));
sampleFrac=NaN(height(Tcell{1}),length(Tcell));
for j=1:length(Tcell)
    T=Tcell{j};
    AF(bIdx,j)=T.BcountsComb(bIdx)./T.ReadDepthPass(bIdx);
    AF(~bIdx,j)=T.AcountsComb(~bIdx)./T.ReadDepthPass(~bIdx);
    matchIdx=f(j,cloneId(:,j))'==T.cnaF;
    sfMajor=AF(:,j).*(T.cnaF.*T.NumCopies+2.*(1-T.cnaF))./(T.NumCopies-T.MinAlCopies);
    sfMinor=AF(:,j).*(T.cnaF.*T.NumCopies+2.*(1-T.cnaF));
    %pairIdx=strncmp(Filter,'SomaticPair',10);
    sampleFrac(matchIdx & alleleId==2,j)=sfMajor(matchIdx & alleleId==2);
    sampleFrac(~matchIdx | alleleId==1,j)=sfMinor(~matchIdx | alleleId==1);
    %sampleFrac(pairIdx,j)=min(sfMajor(pairIdx),sfMinor(pairIdx));
end


%somaticDetected=zeros(size(AF));
if inputParam.NormalSample<1
    formatStr(:,1)=[strcat(germGT',':.:.:.:.:.:.:.:.:.:.:.'); {''}];
    n=2;
else
    n=1;
end
for i=1:length(Tcell)
    T=Tcell{i};
    if inputParam.NormalSample==i
        formatStr(:,n)=[germGT'; {''}];
    else
%        n
 %       size(formatStr)
 %       size(tumorGT)
        formatStr(:,n)=[tumorGT'; {''}];
    end
    formatStr(:,n)=strcat(formatStr(:,n),':',strsplit(sprintf('%-.0f\n',T.ReadDepthPass))');
    aIdx=T.RefComb==T.Acomb;
    formatStr([aIdx; true],n)=strcat(formatStr([aIdx; true],n),':',strsplit(sprintf('%-.0f\n',T.AcountsComb(aIdx)))',',',strsplit(sprintf('%-.0f\n',T.BcountsComb(aIdx)))');
    bIdx=T.RefComb==T.Bcomb;
    formatStr([bIdx; true],n)=strcat(formatStr([bIdx; true],n),':',strsplit(sprintf('%-.0f\n',T.BcountsComb(bIdx)))',',',strsplit(sprintf('%-.0f\n',T.AcountsComb(bIdx)))');
    formatStr([~aIdx & ~bIdx; true],n)=strcat(formatStr([~aIdx & ~bIdx; true],n),':NA,',strsplit(sprintf('%-.0f\n',T.AcountsComb(~aIdx & ~bIdx)))',',',strsplit(sprintf('%-.0f\n',T.BcountsComb(~aIdx & ~bIdx)))');
    filtStr=repmat({'REJECT'},height(T),1);
    filtStr(P.trust(:,i)>=inputParam.pGoodThresh)={'PASS'};
    filtStr(P.trust(:,i)<inputParam.pGoodThresh & P.artifact(:,i)<inputParam.pGoodThresh)={'LowQC'};
    if inputParam.NormalSample<1
        filtStr(somaticDetected(:,i)==1)=strcat(filtStr(somaticDetected(:,i)==1),';SomaticDetected');
        %somaticDetected(P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb),i)=1;
        filtStr(P.Somatic(:,i)>0.5 & ~somaticDetected(:,i))=strcat(filtStr(P.Somatic(:,i)>0.5 & ~somaticDetected(:,i)),';SomaticNotDetected');
        formatStr(:,n)=strcat(formatStr(:,n),':',[filtStr; {''}],':NA');
    else
        %somaticDetected(P.SomaticPair(:,i)>0.5 | (P.Somatic(:,i)>0.5 & P.DataSomatic(:,i)>P.DataHom(:,i) & (T.BcountsComb>0 | T.A~=T.RefComb)),i)=1;
        filtStr(somaticDetected(:,i)==1)=strcat(filtStr(somaticDetected(:,i)==1),';SomaticDetected');
        filtStr(P.Somatic(:,i)>0.5 & ~somaticDetected(:,i))=strcat(filtStr(P.Somatic(:,i)>0.5 & ~somaticDetected(:,i)),';SomaticNotDetected');
        formatStr(:,n)=strcat(formatStr(:,n),':',[filtStr; {''}],':',strsplit(sprintf('%-.3f\n',P.SomaticPair(:,i)))');
    end
    formatStr(:,n)=strcat(formatStr(:,n),':',strsplit(sprintf('%-.0f\n',-10*log10(1-P.trust(:,i))))',':',strsplit(sprintf('%-.0f\n',-10*log10(1-P.artifact(:,i))))',':',strsplit(sprintf('%-.0f\n',P.DataSomatic(:,i)))');
    formatStr([aIdx; true],n)=strcat(formatStr([aIdx; true],n),':',strsplit(sprintf('%-.0f\n',-10*log10(P.DataHom(aIdx,i))))',',',strsplit(sprintf('%-.0f\n',-10*log10(P.DataHet(aIdx,i))))',',NA');
    formatStr([bIdx; true],n)=strcat(formatStr([bIdx; true],n),':NA,',strsplit(sprintf('%-.0f\n',-10*log10(P.DataHet(bIdx,i))))',',',strsplit(sprintf('%-.0f\n',-10*log10(P.DataHom(bIdx,i))))');
    formatStr(~aIdx & ~bIdx,n)=strcat(formatStr(~aIdx & ~bIdx,n),':.');
    if inputParam.NormalSample>0
        formatStr(:,n)=strcat(formatStr(:,n),':',strsplit(sprintf('%-.0f\n',-10*log10(P.DataNonDip(:,i))))');
    else
        formatStr(:,n)=strcat(formatStr(:,n),':NA');
    end
    formatStr([strncmp(Filter,'Somatic',7); true],n)=strcat(formatStr([strncmp(Filter,'Somatic',7); true],n),':',strsplit(sprintf('%-.3f\n',sampleFrac(strncmp(Filter,'Somatic',7),i)))');
    formatStr(~strncmp(Filter,'Somatic',7),n)=strcat(formatStr(~strncmp(Filter,'Somatic',7),n),':NA');
    formatStr(T.NumCopies==2 & T.MinAlCopies==1,n)=strcat(formatStr(T.NumCopies==2 & T.MinAlCopies==1,n),':NA');
    formatStr([T.NumCopies~=2 | T.MinAlCopies~=1; true],n)=strcat(formatStr([T.NumCopies~=2 | T.MinAlCopies~=1; true],n),':',strsplit(sprintf('%-.3f\n',T.cnaF(T.NumCopies~=2 | T.MinAlCopies~=1)))');
    n=n+1;
end
formatStr=formatStr(1:height(T),:);


formatFields=repmat({'GT:DP:AD:FT:PPS:PT:PA:PLS:PL:PLND:VSF:CNF'},size(P,1),1);
sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];

%%% print output
outData=[chrList(T.Chr) num2cell(T.Pos) cellstr(char(ones(size(T.Pos))*46)) cellstr(RefNT) squeeze(AltNT) num2cell(Qual) Filter Info formatFields formatStr];
headers={'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'};
if inputParam.NormalSample<1
    headers=[headers 'InferredGermline' regexp(inputParam.sampleNames,',','split')];
else
    headers=[headers regexp(inputParam.sampleNames,',','split')];
end
for i=1:length(headers)
    fprintf(fout,'%s\t',headers{i});
end

for i=1:size(outData,1)
    fprintf(fout,strcat('\n%s\t%d\t%s\t%s\t%s\t%f\t%s\t%s\t%s',repmat('\t%s',1,n)),outData{i,:});
end

fclose(fout);
