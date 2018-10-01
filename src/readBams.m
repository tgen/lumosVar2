function [T, E]=readBams(inputParam,paramFile)
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

sexChr=regexp(inputParam.sexChr,',','split');
if max(cellfun('length',(regexp(inputParam.sexChr,',','split'))))==0
    chrList=cellstr(num2str(inputParam.autosomes','%-d'));
else
    chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr'];
end

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
failList=[];
parfor chrIdx=1:length(chrList)
    outPosFile=[inputParam.outName '_' chrList{chrIdx} '_pos.txt'];
    outExonFile=[inputParam.outName '_' chrList{chrIdx} '_exon.txt'];
    if exist(outPosFile,'file') && exist(outExonFile,'file')
        ferror=fopen([inputParam.outName '_' chrList{chrIdx} '.log.txt'],'a+');
        [q,w] = system(['tail -n 1 ' outExonFile]);
        data=str2double(regexp(w,'\t','split'));
        currRegion=double(Regions(Regions(:,1)==chrIdx,:));
        if(size(data,2)==8*sampleCount)
           idx=find(data(2)==currRegion(:,2));
           if(idx==size(currRegion,1))
                fprintf(ferror,'%s\n',['bam data already read for ' chrList{chrIdx}]);
                continue;
            end
        end
        fprintf(ferror,'%s\n',['restarting gvm on ' chrList{chrIdx}]);
    else
        ferror=fopen([inputParam.outName '_' chrList{chrIdx} '.log.txt'],'w');
    end
    [status,out]=system([inputParam.gvmPath ' --conf ' paramFile ' --chr ' chrList{chrIdx}]);
    fprintf(ferror,'%s\n',out);
    if(status==0)
        fprintf(ferror,'%s\n',['gvm completed on ' chrList{chrIdx}]);
    else
        fprintf(ferror,'%s\n',['gvm failed with status ' num2str(status) ' on chr ' chrList{chrIdx}]);
        failList=[failList chrIdx];
    end
end
if ~isempty(failList)
    message=['gvm failed on ' num2str(failList)]
    return;
else
    message='finished processing chromosomes'
end

for i=1:length(chrList)
    outPosFile=[inputParam.outName '_' chrList{i} '_pos.txt'];
    if(system(['[ -s ' outPosFile ' ]'])==0)    
        AllData{i}=dlmread(outPosFile);
    end
    outExonFile=[inputParam.outName '_' chrList{i} '_exon.txt'];
    AllExonData{i}=reshape(dlmread(outExonFile),[],8,sampleCount);
end
%profile off;
%profsave;
%profile resume;

%%%% create position data table
ColHeaders={'Chr','Pos','ReadDepth','ReadDepthPass','Ref','A','ACountF','ACountR','AmeanBQ','AmeanMQ','AmeanPMM','AmeanReadPos','B','BCountF','BCountR','BmeanBQ','BmeanMQ','BmeanPMM','BmeanReadPos','ApopAF','BpopAF','CosmicCount','ControlRD','PosMapQC','perReadPass','abFrac'};
matLen=cellfun(@numel,AllData)./27;
dataMat=zeros(sum(matLen),length(ColHeaders)+1);
currIdx=1;
for i=1:length(AllData)
    dataMat(currIdx:currIdx+matLen(i)-1,:)=AllData{i};
    currIdx=currIdx+matLen(i);
end
clear AllData;
%profile off;
%profsave;
%profile resume;
pack;
if isempty(dataMat)
    T={};                           
else
    ids=unique(dataMat(:,1));
    for i=1:length(ids)
        T{i}=array2table(dataMat(dataMat(:,1)==ids(i),2:end),'VariableNames',ColHeaders);
    end
end
clear dataMat;
message='finished combining tumor data'
%profile off;
%profsave;
%profile resume;

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

%profile off;
%profsave;
%profile resume;
%message='finished combining exon data'
%return
