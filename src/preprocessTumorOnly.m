function [T, E]=preprocessTumorOnly(inputParam,paramFile)
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

%%%import bed file
fid=fopen(inputParam.regionsFile);
Regions=cell2mat(textscan(fid,'%d%d%d %*[^\n]'));
fclose(fid);

fid=fopen(inputParam.bamList);
bamList=textscan(fid,'%s');
sampleCount=length(bamList{1});
fclose(fid);

%%% process by chromosome
chrList=[1:22];
parfor chrIdx=1:length(chrList)
    if exist([inputParam.outName '_' num2str(chrList(chrIdx)) '_pos.txt'],'file') & exist([inputParam.outName '_' num2str(chrList(chrIdx)) '_exon.txt'],'file')
        %fid=fopen([inputParam.outName '_' num2str(chrList(chrIdx)) '_pos.txt']);
        %tempData=textscan(fid,'%u64');
        %fclose(fid);
        tempData=dlmread([inputParam.outName '_' num2str(chrList(chrIdx)) '_pos.txt']);
        tempExonRD=dlmread([inputParam.outName '_' num2str(chrList(chrIdx)) '_exon.txt']);
        AllData{chrIdx}=tempData;
        AllExonData{chrIdx}=reshape(tempExonRD,[],8,sampleCount); 
        continue;
    end
    outname=[inputParam.outName '_' num2str(chrList(chrIdx)) '_log.txt'];
    fout=fopen(outname,'w');
    currRegion=double(Regions(Regions(:,1)==chrList(chrIdx),:));
    largeIdx=find(currRegion(:,3)-currRegion(:,2)>inputParam.blockSize);
    %%%%breakup regions larger than blockSize
    if ~isempty(largeIdx)
        subCount=ceil((currRegion(largeIdx,3)-currRegion(largeIdx,2))./inputParam.blockSize);
        newSize=(currRegion(largeIdx,3)-currRegion(largeIdx,2))./subCount;
        newRegions=nan(sum(subCount),3);
        newRegions(1:subCount(1),:)=[double(chrList(chrIdx))*ones(subCount(1),1) round([currRegion(largeIdx(1),2):newSize(1):currRegion(largeIdx(1),3)-newSize(1)]') round([currRegion(largeIdx(1),2)+newSize(1)-1:newSize(1):currRegion(largeIdx(1),3)]')];
        for i=2:length(largeIdx)
            newRegions(sum(subCount(1:i-1))+1:sum(subCount(1:i)),:)=[double(chrList(chrIdx))*ones(subCount(i),1) round([currRegion(largeIdx(i),2):newSize(i):currRegion(largeIdx(i),3)-newSize(i)]') round([currRegion(largeIdx(i),2)+newSize(i)-1:newSize(i):currRegion(largeIdx(i),3)]')];
        end
        currRegion=[currRegion(currRegion(:,3)-currRegion(:,2)<=inputParam.blockSize,:); newRegions];
        currRegion=sortrows(currRegion,2);
    end
    %%%% make sure snpVCF is valid
    snpVCF=[inputParam.snpVCFpath num2str(chrList(chrIdx)) inputParam.snpVCFname];
    if(~exist(snpVCF,'file'))
        fprintf(fout,'%s\n',['snpVCF file not found in: ' snpVCF]);
        continue;
    else
        fprintf(fout,'%s\n',['snpVCF file found: ' snpVCF]);
    end
    %%% make sure normal data is valid
    NormalPath=[inputParam.NormalBase num2str(chrList(chrIdx)) '.txt.gz'];
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
    endIdx=1;
    dataIdx=1;
    tempData=zeros(round(sum(currRegion(:,3)-currRegion(:,2))./100),27);
    tempExonRD=zeros(size(currRegion,1),8,sampleCount);
    while(startIdx<=size(currRegion,1))
        %%% define block
        endIdx=find(cumsum(currRegion(startIdx:end,3)-currRegion(startIdx:end,2))>inputParam.blockSize,1)+startIdx-1;
        if isempty(endIdx)
            endIdx=size(currRegion,1);
        end
        block=[num2str(chrList(chrIdx)) ':' num2str(currRegion(startIdx,2)) '-' num2str(currRegion(endIdx,3))]
        fprintf(fout,'\n%s',['Analyzing ' block]);
        %%% get tumor data
        cd(inputParam.workingDirectory);
        output=strsplit(perl('parsePileupData.pl',paramFile,block,num2str(sampleCount)),'\n');
        %fprintf(fout,'\n%s',output{:});
        idx=~cellfun('isempty',regexp(output,'^\d'));
        TumorData=str2num(char(output(idx)'));
        fprintf(fout,'\n%s',['TumorData has length:' num2str(size(TumorData,1))]);
        idx=~cellfun('isempty',regexp(output,'^\@'));
        temp=char(output(idx)');
        readDepth=str2num(temp(:,2:end));
        fprintf(fout,'\n%s',['readDepth has length:' num2str(size(readDepth,1))]);
        %%% get normal data
        cd(inputParam.tabixPath);
        [status,output]=system(['./tabix ' NormalPath ' ' block]);
        NormalData=str2num(char(output));
        fprintf(fout,'\n%s',['NormalData has length:' num2str(size(NormalData,1))]);
        cd(inputParam.workingDirectory);
        currRegion(startIdx:endIdx,:);
        if isempty(readDepth)
            startIdx=endIdx+1;
            continue
        end
        %%% find exon means
        ids=unique(readDepth(:,1))
        for i=1:length(ids)
            currReadDepth=readDepth(readDepth(:,1)==ids(i),2:end);
            tempExonRD(startIdx:endIdx,:,i)=[currRegion(startIdx:endIdx,1:3) getMeanInRegions(currReadDepth(:,1:2),currReadDepth(:,3),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,3),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,4),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,5),currRegion(startIdx:endIdx,:)) getMeanInRegionsExcludeNaN(NormalData(:,1:2),NormalData(:,6),currRegion(startIdx:endIdx,:))];
        end
        startIdx=endIdx+1;
        %%%% combine tumor and normal data
        if isempty(TumorData)
            continue
        end
        [Lia,Locb]=ismember(TumorData(:,3),NormalData(:,2));
        tempData(dataIdx:dataIdx+sum(Lia)-1,:)=[TumorData(Lia,:) NormalData(Locb(Lia),3:end)];
        dataIdx=dataIdx+sum(Lia);
        fprintf(fout,'\n%s',[' found ' num2str(sum(Lia)) ' canidate positions']);
        waitbar(endIdx./length(currRegion));       
    end   
    dlmwrite([inputParam.outName '_' num2str(chrList(chrIdx)) '_pos.txt'],tempData(tempData(:,1)>0,:),'precision',9);
    dlmwrite([inputParam.outName '_' num2str(chrList(chrIdx)) '_exon.txt'],reshape(tempExonRD(tempExonRD(:,1)>0,:,:),[],8*sampleCount),'precision',9); 
    AllData{chrIdx}=tempData(tempData(:,1)>0,:);
    AllExonData{chrIdx}=tempExonRD(tempExonRD(:,1)>0,:,:); 
end

message='finished processing chromosomes'

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
