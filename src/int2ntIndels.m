function nt=int2ntIndels(indelInt)
% converts integer representation of nucleotides to strings, assuming a
% base 5 representation of the nucleotide and numbers greater than 5^13
% represent the length of the indel rather than the sequence
%
% Syntax:  nt=int2ntIndels(indelInt)
% 
% Inputs:
%   indelInt - integer representing nucleotide
%   
% Outputs:
%   nt - string representation of nucleotide
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: LumosVarMain, writeJointVCF

%%%get length of long indels
indelStr(indelInt>5^13)=cellstr(strcat('L',num2str(indelInt(indelInt>5^13)-5^13)));

%%%convert indels to base 5
indelStr(indelInt<=5^13 & indelInt>0)=regexprep(cellstr(dec2base(indelInt(indelInt<=5^13 & indelInt>0),5)),'^0+','');

%%% place holder for problematic positions
indelStr(indelInt<=0)={'X'};

%%%find positions that should have nt sequences
pos=cellfun('isempty',strfind(indelStr,'e')) & cellfun('isempty',strfind(indelStr,'L')) & cellfun('isempty',strfind(indelStr,'X'));

%%%make character array of sequences
indelChar=char(indelStr');
ntChar='';
for i=1:size(indelChar,2)
    if(sum(~isspace(indelChar(:,i)) & pos'>0))
        ntChar(~isspace(indelChar(:,i)) & pos',i)=int2nt(str2double(cellstr((indelChar(~isspace(indelChar(:,i)) & pos',i)))));
    end
end

%%%get nt values
nt=regexprep(cellstr(ntChar),'[\x0]','');
nt(~pos)=indelStr(~pos);
nt(indelInt<=0)={'N'};

return;