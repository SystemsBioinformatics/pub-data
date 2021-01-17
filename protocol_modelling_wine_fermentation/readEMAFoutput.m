function [required, alternatives] = readEMAFoutput(fileName)

[~,s] = xlsread(fileName);

required = split(regexprep(s{4,1},'REQUIRED,',''),',');
for i =  6:size(s,1)
    alternatives{i-5} = regexprep(s{i,1},['ALT' num2str(i-5) ','],'');
end
alternatives = alternatives';
end