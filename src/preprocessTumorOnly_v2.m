function [T, E]=preprocessTumorOnly_v2(inputParam,paramFile)
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

sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr']

regTable=readtable(inputParam.regionsFile,'FileType','text','Delimiter','\t','ReadVariableNames',false);
[lia,locb]=ismember(regTable{:,1},chrList);
Regions=[locb(lia) regTable{lia,2:3}];

%%%import bed fileregTable=readtable(regionsFile,'FileType','text','Delimiter','\t','ReadVariableNames',false);
[lia,locb]=ismember(regTable{:,1},chrList);
Regions=[locb(lia) regTable{lia,2:3}];

fid=fopen(inputParam.bamList);
bamList=textscan(fid,'%s');
sampleCount=length(bamList{1});
fclose(fid);

%%% process by chromosome

parfor chrIdx=1:length(chrList)
    outPosFile=[inputParam.outName '_' chrList{chrIdx} '_pos.txt'];
    outExonFile=[inputParam.outName '_' chrList{chrIdx} '_exon.txt'];
    if exist(outPosFile,'file') && exist(outExonFile,'file')
        fpos=fopen(outPosFile,'a+');
        fexon=fopen(outExonFile,'a+');
        [q,w] = system(['tail -n 1 ' outExonFile]);
        data=str2double(regexp(w,'\t','split'));
        currRegion=double(Regions(Regions(:,1)==chrIdx,:));
        if(size(data,2)==8*sampleCount)
            idx=find(data(2)==currRegion(:,2));
            if(idx==size(currRegion,1))
                continue;
            else
                currRegion=currRegion(idx+1:end,:);
            end
        end
    else
        fpos=fopen(outPosFile,'w');
        fexon=fopen(outExonFile,'w');
        currRegion=double(Regions(Regions(:,1)==chrIdx,:));
    end
    outname=[inputParam.outName '_' chrList{chrIdx} '_log.txt'];
    fout=fopen(outname,'w');
    largeIdx=find(currRegion(:,3)-currRegion(:,2)>inputParam.blockSize);
    %%%%breakup regions larger than blockSize
    if ~isempty(largeIdx)
        subCount=ceil((currRegion(largeIdx,3)-currRegion(largeIdx,2))./inputParam.blockSize);
        newSize=(currRegion(largeIdx,3)-currRegion(largeIdx,2))./subCount;
        newRegions=nan(sum(subCount),3);
        newRegions(1:subCount(1),:)=[double(chrIdx)*ones(subCount(1),1) round([currRegion(largeIdx(1),2):newSize(1):currRegion(largeIdx(1),3)-newSize(1)]') round([currRegion(largeIdx(1),2)+newSize(1)-1:newSize(1):currRegion(largeIdx(1),3)]')];
        for i=2:length(largeIdx)
            newRegions(sum(subCount(1:i-1))+1:sum(subCount(1:i)),:)=[double(chrIdx)*ones(subCount(i),1) round([currRegion(largeIdx(i),2):newSize(i):currRegion(largeIdx(i),3)-newSize(i)]') round([currRegion(largeIdx(i),2)+newSize(i)-1:newSize(i):currRegion(largeIdx(i),3)]')];
        end
        currRegion=[currRegion(currRegion(:,3)-currRegion(:,2)<=inputParam.blockSize,:); newRegions];
        currRegion=sortrows(currRegion,2);
    end
    %%%% make sure snpVCF is valid
    snpVCF=[inputParam.snpVCFpath chrList{chrIdx} inputParam.snpVCFname];
    if(~exist(snpVCF,'file'))
        fprintf(fout,'%s\n',['snpVCF file not found in: ' snpVCF]);
        continue;
    else
        fprintf(fout,'%s\n',['snpVCF file found: ' snpVCF]);
    end
    %%% make sure normal data is valid
    NormalPath=[inputParam.NormalBase chrList{chrIdx} '.txt.gz'];
    if(~exist(NormalPath,'file'))
        fprintf(fout,'%s\n',['NormalData file not found in: ' NormalPath]);
        continue;
    else
        fprintf(fout,'%s\n',['NormalData file found: ' NormalPath]);
    end
    %fid=fopen(inputParam.bamList);
    %bamList=textscan(fid,'%s');
    %sampleCount=length(bamList{1});
    %fclose(fid);
    %%%get data by region
    startIdx=1;
    %endIdx=1;
    dataIdx=1;
    %tempData=zeros(round(sum(currRegion(:,3)-currRegion(:,2))./100),27);
    %tempExonRD=zeros(size(currRegion,1),8,sampleCount);
    while(startIdx<=size(currRegion,1))
        %%% define block
        endIdx=find(cumsum(currRegion(startIdx:end,3)-currRegion(startIdx:end,2))>inputParam.blockSize,1)+startIdx-1;
        if isempty(endIdx)
            endIdx=size(currRegion,1);
        end
        block=[chrList{chrIdx} ':' num2str(currRegion(startIdx,2)) '-' num2str(currRegion(endIdx,3))]
        fprintf(fout,'\n%s',['Analyzing ' block]);
        %%% get tumor data
        cd(inputParam.workingDirectory);
        output=strsplit(perl('parsePileupData.pl',paramFile,block,num2str(sampleCount)),'\n');
        idx1=~cellfun('isempty',regexp(output,'^\d'));
        TumorData=str2num(char(output(idx1)'));
        fprintf(fout,'\n%s',['TumorData has length:' num2str(size(TumorData,1))]);
        idx2=~cellfun('isempty',regexp(output,'^\@'));
        temp=char(output(idx2)');
        fprintf(fout,'\n%s',output{~idx1 & ~idx2});
        readDepth=str2num(temp(:,2:end));
        fprintf(fout,'\n%s',['readDepth has length:' num2str(size(readDepth,1))]);
        %%% get normal data
        cd(inputParam.tabixPath);
        [status,output]=system(['./tabix ' NormalPath ' ' block]);
        NormalData=str2num(char(output));
        fprintf(fout,'\n%s',['NormalData has length:' num2str(size(NormalData,1))]);
        cd(inputParam.workingDirectory);
        currRegion(startIdx:endIdx,:);
        if isempty(readDepth) || isempty(NormalData)
            startIdx=endIdx+1;
            continue
        end
        %%% find exon means
        ids=unique(readDepth(:,1))
        tempExonRD=zeros(length(startIdx:endIdx),8,sampleCount);
        for i=1:length(ids)
            currReadDepth=readDepth(readDepth(:,1)==ids(i),2:end);
            tempExonRD(:,:,i)=[currRegion(startIdx:endIdx,1:3) getMeanInRegions(currReadDepth(:,1:2),currReadDepth(:,3),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,3),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,4),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,5),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,6),currRegion(startIdx:endIdx,:))];
        end
        fprintf(fexon,[strjoin(cellstr(repmat('%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f',sampleCount,1)),'\\t') '\n'],reshape(tempExonRD,[],8*sampleCount)');
        startIdx=endIdx+1;
        %%%% combine tumor and normal data
        if isempty(TumorData)
            continue
        end
        [Lia,Locb]=ismember(TumorData(:,3),NormalData(:,2));
        tempData=[TumorData(Lia,:) NormalData(Locb(Lia),3:end)];
        size(tempData);
        fprintf(fpos,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n',tempData');
        dataIdx=dataIdx+sum(Lia);
        fprintf(fout,'\n%s',[' found ' num2str(sum(Lia)) ' canidate positions']);
        waitbar(endIdx./length(currRegion));       
    end   
    %dlmwrite([inputParam.outName '_' chrList{chrIdx} '_pos.txt'],tempData(tempData(:,1)>0,:),'precision',9);
    %dlmwrite([inputParam.outName '_' chrList{chrIdx} '_exon.txt'],reshape(tempExonRD(tempExonRD(:,1)>0,:,:),[],8*sampleCount),'precision',9); 
    %AllData{chrIdx}=tempData(tempData(:,1)>0,:);
    %AllExonData{chrIdx}=tempExonRD(tempExonRD(:,1)>0,:,:); 
end
message='finished processing chromosomes'

for i=1:length(chrList)
    outPosFile=[inputParam.outName '_' chrList{i} '_pos.txt'];
    if(system(['[ -s ' outPosFile ' ]'])==0)    
        AllData{i}=dlmread(outPosFile);
    end
    outExonFile=[inputParam.outName '_' chrList{i} '_exon.txt'];
    AllExonData{i}=reshape(dlmread(outExonFile),[],8,sampleCount);
end

%%%% create position data table
ColHeaders={'Chr','Pos','ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ','AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ','BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount','ControlRD','PosMapQC','perReadPass','abFrac'};
matLen=cellfun(@numel,AllData)./27;
dataMat=zeros(sum(matLen),length(ColHeaders)+1);
currIdx=1;
for i=1:length(AllData)
    dataMat(currIdx:currIdx+matLen(i)-1,:)=AllData{i};
    currIdx=currIdx+matLen(i);
end
if isempty(dataMat)
    T={};
else
    ids=unique(dataMat(:,1));
    for i=1:length(ids)
        T{i}=array2table(dataMat(dataMat(:,1)==ids(i),2:end),'VariableNames',ColHeaders);
    end
end
clear dataMat AllData;
message='finished combining tumor data'

%%% create exon data table
exonColHeaders={'Chr','StartPos','EndPos','TumorRD','NormalRD', 'MapQC', 'perReadPass', 'abFrac'};
exonMatLen=cellfun(@numel,AllExonData)./(8*sampleCount);
exonRD=zeros(sum(exonMatLen),sum(length(exonColHeaders)),sampleCount);
currIdx=1;
for i=1:length(AllExonData)
    exonRD(currIdx:currIdx+exonMatLen(i)-1,:,:)=AllExonData{i};
    currIdx=currIdx+exonMatLen(i);
end
for i=1:size(exonRD,3)
    E{i}=array2table(exonRD(:,:,i),'VariableNames',exonColHeaders);
end

%message='finished combining exon data'
%return
