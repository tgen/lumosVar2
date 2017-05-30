
function printNormalMetrics_byChr_v2(configFile,step)
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

chrIdx=str2num(step);
step 
chrIdx
sexChr=regexp(inputParam.sexChr,',','split');
chrList=[cellstr(num2str(inputParam.autosomes','%-d')); sexChr']
[lia,locb]=ismember(chrList(chrIdx),sexChr);
sexList=regexp(inputParam.sexList,',','split');
ploidyList=nan(size(sexList));
if lia
    mIdx=strcmp('M',sexList);
    ploidyList(mIdx)=inputParam.M(locb);
    fIdx=strcmp('F',sexList);
    ploidyList(fIdx)=inputParam.F(locb);
else
    ploidyList=2*ones(size(sexList))
end


%%% read in bedfile
regTable=readtable(regionsFile,'FileType','text','Delimiter','\t','ReadVariableNames',false);
[lia,locb]=ismember(regTable{:,1},chrList);
Regions=[locb(lia) regTable{lia,2:3}];
% size(regTable)
% if ~isnumeric(regTable{:,1})
%     chr=cellfun(@str2num,regTable{:,1},'UniformOutput',0);
%     size(chr)
%     pos=~cellfun(@isempty,chr);
%     sum(pos)
%     Regions=[cell2mat(chr(pos)) regTable{pos,2:3}];
% else
%     Regions=regTable{:,1:3};
% end



%for chrIdx=1:length(chrList)
%%% open output files
if(~exist([outfile '_' chrList{chrIdx} '.txt.gz'],'file'))
    if(exist([outfile '_' chrList{chrIdx} '.txt'],'file'))
        fout=fopen([outfile '_' chrList{chrIdx} '.txt'],'a+');
        ferror=fopen([outfile '_' chrList{chrIdx} '.errorLog.txt'],'a+');
        [q,w] = system(['tail -n 1 ' outfile '_' chrList{chrIdx} '.txt']);
        data=str2double(regexp(w,'\t','split'));
        currRegion=double(Regions(Regions(:,1)==chrIdx,:));
        size(currRegion)
        if(size(data,2)==6)
            idx=getPosInRegions([chrIdx data(2)],currRegion);
            currRegion(idx,:);
            currRegion=currRegion(idx:end,:);
            currRegion(1,2)=data(2);
        end
    else
        fout=fopen([outfile '_' chrList{chrIdx} '.txt'],'w');
        ferror=fopen([outfile '_' chrList{chrIdx} '.errorLog.txt'],'w');
        currRegion=double(Regions(Regions(:,1)==chrIdx,:));
    end
    largeIdx=find(currRegion(:,3)-currRegion(:,2)>blockSize);
    subCount=ceil((currRegion(largeIdx,3)-currRegion(largeIdx,2))./blockSize);
    newSize=(currRegion(largeIdx,3)-currRegion(largeIdx,2))./subCount;
    newRegions=nan(sum(subCount),3);
    if ~isempty(subCount)
        newRegions(1:subCount(1),:)=[chrIdx*ones(subCount(1),1) round([currRegion(largeIdx(1),2):newSize(1):currRegion(largeIdx(1),3)-newSize(1)]') round([currRegion(largeIdx(1),2)+newSize(1)-1:newSize(1):currRegion(largeIdx(1),3)]')];
        for i=2:length(largeIdx)
            newRegions(sum(subCount(1:i-1))+1:sum(subCount(1:i)),:)=[chrIdx*ones(subCount(i),1) round([currRegion(largeIdx(i),2):newSize(i):currRegion(largeIdx(i),3)-newSize(i)]') round([currRegion(largeIdx(i),2)+newSize(i)-1:newSize(i):currRegion(largeIdx(i),3)]')];
        end
        currRegion=[currRegion(currRegion(:,3)-currRegion(:,2)<=blockSize,:); newRegions];
        currRegion=sortrows(currRegion,2);
    end
    
    startIdx=1;
    endIdx=1;
    %%% split up regions greater than block size
    
    %%% analyze by region
    while(endIdx<size(currRegion,1))
        %currRegion
        [startIdx endIdx]
        endIdx=find(cumsum(currRegion(startIdx:end,3)-currRegion(startIdx:end,2))>blockSize,1)+startIdx-1;
        if isempty(endIdx)
            endIdx=size(currRegion,1);
        end
        %%% bet pileup data
        block=[chrList{chrIdx} ':' num2str(currRegion(startIdx,2)) '-' num2str(currRegion(endIdx,3))]
        output=strsplit(perl('parsePileupData.packed.pl',configFile,block,'0'),'\n');
        idx=~cellfun('isempty',regexp(output,'^\d'));
        NormalData=str2num(char(regexprep(output(idx),chrList{chrIdx},num2str(chrIdx))'));
        fprintf(ferror,'%s\n',output{~idx});
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
        idx=getPosInRegions([chrIdx.*ones(currRegion(endIdx,3)-currRegion(startIdx,2)+1,1) allPos],currRegion(startIdx:endIdx,:));
        outPos=allPos(~isnan(idx));
        normalMetrics=nan(length(outPos),6,length(sampleIds));
        for i=1:length(sampleIds)
            sampleData=NormalData(NormalData(:,1)==sampleIds(i),:);
            [Lia,Locb]=ismember(outPos,sampleData(:,3));
            currData=[chrIdx.*ones(size(outPos)) outPos nan(length(outPos),size(NormalData,2)-3)];
            currData(Lia,:)=sampleData(Locb(Lia),2:end);
            normalMetrics(:,:,sampleIds(i))=calculateNormalMetrics_v2(currData,priorMapError,ploidyList(i));
        end
        %%% find mean across samples and print to file
        meanNormalMetrics=nanmean(normalMetrics,3);
        fprintf(fout,'%d\t%d\t%f\t%f\t%f\t%f\n',meanNormalMetrics');
        startIdx=endIdx+1;
    end
    fclose(fout);
    %%% bgzip and tabix index output file
    system(['sh indexControlMetrics.sh ' outfile '_' chrList{chrIdx} '.txt ' inputParam.tabixPath]);
end

exit(0);
