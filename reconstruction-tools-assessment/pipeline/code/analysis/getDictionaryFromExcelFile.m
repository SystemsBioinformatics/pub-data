function [k,v,kt,vt] = getDictionaryFromExcelFile(fileName)

[~,s] = xlsread(fileName);
kt = [{'language_key'};s(:,1)];
k = [{'keys'};s(:,2)];
vt = [{'language_values'};s(:,3)];
v = [{'values'};s(:,4)];

end