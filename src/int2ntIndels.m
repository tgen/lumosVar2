function nt=int2ntIndels(indelInt)

indelStr(indelInt>2^52)={'INF'};
%indelStr(indelInt<=2^52)=num2str(dec2base(indelInt(indelInt<=2^52),5),'%64d');
indelStr(indelInt<=2^52 & indelInt>0)=regexprep(cellstr(dec2base(indelInt(indelInt<=2^52 & indelInt>0),5)),'^0+','');
indelStr(indelInt<=0)={'X'};

pos=cellfun('isempty',strfind(indelStr,'e')) & cellfun('isempty',strfind(indelStr,'I')) & cellfun('isempty',strfind(indelStr,'X'));
indelChar=char(indelStr');

ntChar='';
for i=1:size(indelChar,2)
    ntChar(~isspace(indelChar(:,i)) & pos',i)=int2nt(str2num(indelChar(~isspace(indelChar(:,i)) & pos',i)));
end

nt=regexprep(cellstr(ntChar),'[\x0]','');
nt(~pos)={'*LongIndel*'};
nt(indelInt<=0)={'*'};

return;