
function exitCode=printNormalMetrics(configFile)
%printNormalMetrics - gets mean read depths and 
%position quality scores for a set of normal bams
%calls parsePileupData.packed.pl to parse samtools output
%writes a bgziped tabix index table for each chromosome
%prerequisite for running TumorOnlyWrapper.m
%
% Syntax:  printNormalMetrics(configFile)
%
% Inputs:
%   configFile - contains input parameters in yaml format
%   
% Outputs:
%   writes one table per chromosome.  Each table is bgziped and tabix indexed. 
%       tables have columns: 1-'Chr',2-'Pos',3-'ControlRD', 4-'PosMapQC', 
%       5-'perReadPass',6-'abFrac'
%
% Other m-files required: calculateNormalMetrics.m
% Other requirements: parsePileupData.packed.pl, indexControlMetrics.sh,
%   samtools, htslib
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

%%% read in parameters
inputParam=readInputs(configFile);
regionsFile=inputParam.regionsFile;
outfile=inputParam.outfile;
blockSize=inputParam.blockSize;
priorMapError=inputParam.priorMapError;

%%% read in bedfile
regTable=readtable(regionsFile,'FileType','text','Delimiter','\t');
size(regTable)
if ~isnumeric(regTable{:,1})
    chr=cellfun(@str2num,regTable{:,1},'UniformOutput',0);
    size(chr)
    pos=~cellfun(@isempty,chr);
    sum(pos)
    Regions=[cell2mat(chr(pos)) regTable{pos,2:3}];
else
    Regions=regTable{:,1:3};
end


% fid=fopen(regionsFile);
% Regions=cell2mat(textscan(fid,'%d%d%d %*[^\n]'));
% fclose(fid);

%%% start parrellel pool
delete(gcp);
parpool;

%%% analyze by chromoxome
chrList=[1:22];
parfor chrIdx=1:length(chrList)
    %%% open output files
    if(exist([outfile '_' num2str(chrList(chrIdx)) '.txt.gz'],'file'))
        continue;
    end
    if(exist([outfile '_' num2str(chrList(chrIdx)) '.txt'],'file'))
        fout=fopen([outfile '_' num2str(chrList(chrIdx)) '.txt'],'a+');
        ferror=fopen([outfile '_' num2str(chrList(chrIdx)) '.errorLog.txt'],'a+');
        [q,w] = system(['tail -n 1 ' outfile '_' num2str(chrList(chrIdx)) '.txt']);
        data=str2double(regexp(w,'\t','split'));
        currRegion=double(Regions(Regions(:,1)==chrList(chrIdx),:));
        if(size(data,2)==6)
            idx=getPosInRegions([chrList(chrIdx) data(2)],currRegion);
            currRegion=currRegion(idx:end,:);
            currRegion(1,2)=data(2);
        end
    else
        fout=fopen([outfile '_' num2str(chrList(chrIdx)) '.txt'],'w');
        ferror=fopen([outfile '_' num2str(chrList(chrIdx)) '.errorLog.txt'],'w');
        currRegion=double(Regions(Regions(:,1)==chrList(chrIdx),:));
    end
    largeIdx=find(currRegion(:,3)-currRegion(:,2)>blockSize);
    subCount=ceil((currRegion(largeIdx,3)-currRegion(largeIdx,2))./blockSize);
    newSize=(currRegion(largeIdx,3)-currRegion(largeIdx,2))./subCount;
    newRegions=nan(sum(subCount),3);
    newRegions(1:subCount(1),:)=[double(chrList(chrIdx))*ones(subCount(1),1) round([currRegion(largeIdx(1),2):newSize(1):currRegion(largeIdx(1),3)-newSize(1)]') round([currRegion(largeIdx(1),2)+newSize(1)-1:newSize(1):currRegion(largeIdx(1),3)]')];
    for i=2:length(largeIdx)
        newRegions(sum(subCount(1:i-1))+1:sum(subCount(1:i)),:)=[double(chrList(chrIdx))*ones(subCount(i),1) round([currRegion(largeIdx(i),2):newSize(i):currRegion(largeIdx(i),3)-newSize(i)]') round([currRegion(largeIdx(i),2)+newSize(i)-1:newSize(i):currRegion(largeIdx(i),3)]')];
    end
    currRegion=[currRegion(currRegion(:,3)-currRegion(:,2)<=blockSize,:); newRegions];
    currRegion=sortrows(currRegion,2);
    
    startIdx=1;
    endIdx=1;
    %%% split up regions greater than block size
    
    %%% analyze by region
    while(endIdx<length(currRegion))
        endIdx=find(cumsum(currRegion(startIdx:end,3)-currRegion(startIdx:end,2))>blockSize,1)+startIdx-1;
        if isempty(endIdx)
            endIdx=length(currRegion);
        end
        %%% bet pileup data
        block=[num2str(chrList(chrIdx)) ':' num2str(currRegion(startIdx,2)) '-' num2str(currRegion(endIdx,3))]
        output=strsplit(perl('parsePileupData.packed.pl',configFile,block,'0'),'\n');
        idx=~cellfun('isempty',regexp(output,'^\d'));
        NormalData=str2num(char(output(idx)'));
        if(isempty(NormalData))
            fprintf(ferror,'%s\n',['Normal Data not found in ' block]);
            for i=1:length(output)
                fprintf(ferror,'%s\n',output{i});
            end
            startIdx=endIdx+1;
            continue;
        end
        %%% calculate normal metrics
        output={};
        sampleIds=unique(NormalData(:,1));
        allPos=[currRegion(startIdx,2):currRegion(endIdx,3)]';
        idx=getPosInRegions([double(chrList(chrIdx)).*ones(currRegion(endIdx,3)-currRegion(startIdx,2)+1,1) allPos],currRegion(startIdx:endIdx,:));
        outPos=allPos(~isnan(idx));
        normalMetrics=nan(length(outPos),6,length(sampleIds));
        for i=1:length(sampleIds)
            sampleData=NormalData(NormalData(:,1)==sampleIds(i),:);
            [Lia,Locb]=ismember(outPos,sampleData(:,3));
            currData=[double(chrList(chrIdx)).*ones(size(outPos)) outPos nan(length(outPos),size(NormalData,2)-3)];
            currData(Lia,:)=sampleData(Locb(Lia),2:end);
            normalMetrics(:,:,sampleIds(i))=calculateNormalMetrics(currData,priorMapError);
        end
        %%% find mean across samples and print to file
        meanNormalMetrics=nanmean(normalMetrics,3);
        fprintf(fout,'%d\t%d\t%f\t%f\t%f\t%f\n',meanNormalMetrics');
        startIdx=endIdx+1;
    end
    fclose(fout);
    %%% bgzip and tabix index output file
    system(['sh indexControlMetrics.sh ' outfile '_' num2str(chrList(chrIdx)) '.txt ' inputParam.tabixPath]);
end

exitCode=0;
