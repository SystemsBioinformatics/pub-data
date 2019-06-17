function getTotalNumberOfNewRelationships

global rootFolder

[n1,s1] = xlsread(fullfile(rootFolder, 'results', 'lpl', 'newRelationships.xls'),'new');

[n2,s2] = xlsread(fullfile(rootFolder, 'results', 'bpe', 'newRelationships.xls'),'new');

sT = [s1;s2];
for i = length(sT):-1:2
    pos1 = find(strcmp(sT(1:i-1,1),sT{i,1}));
    pos2 = find(strcmp(sT(1:i-1,2),sT{i,2}));
    
    if ~isempty(pos1) && ~isempty(pos2)
        if ~isempty(intersect(pos1,pos2))
            sT(i,:) = [];
        end
    end
end

labels = {'ID in manual model', 'ID in automatic model', 'Equation in manual model', 'Equation in automatic model', 'Translated equation in automatic model'};
sT = [labels; sT];
xlswrite('summaryNewRelations',sT)

end