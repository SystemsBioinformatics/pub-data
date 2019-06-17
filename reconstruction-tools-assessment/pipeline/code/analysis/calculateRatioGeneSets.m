function [ratio, tableGeneComparison]= calculateRatioGeneSets(models, names, geneMatrix, geneIDs, geneExcelFileName, species)

coverage = zeros(length(models)-1,1);
additional = zeros(length(models)-1,1);
nOriginalModel = length(models{1}.genes);
total = zeros(length(models)-1,1);
for i = 1:length(models)-1
    total(i) = length(find(geneMatrix(1:end,i+1)));
    additional(i) = length(find(geneMatrix(nOriginalModel+1:end,i+1)));
    coverage(i) = length(find(geneMatrix(1:nOriginalModel,i+1)));
end
coverage_p = 100*coverage/nOriginalModel;
additional_p = 100*additional/nOriginalModel;
ratio = coverage_p./additional_p;
% ratio = coverage_p./(additional_p.*(100-coverage_p));
names = names(2:end);
if size(names,2)>size(names,1)
    names = names';
end
tableGeneComparison = [[{'model'} {'gene coverage'} {'additional genes'} {'genes coverage(%)'} {'additional genes (%)'} {'ratio coverage/additional'}]; ...
    [names, num2cell(coverage), num2cell(additional), num2cell(coverage_p), num2cell(additional_p), num2cell(ratio)]];

if exist([geneExcelFileName '_' species '.xls'], 'file') ==2
    delete([geneExcelFileName '_' species '.xls'])
end
xlswrite([geneExcelFileName '_' species], geneIDs, 'genes');
xlswrite([geneExcelFileName '_' species], tableGeneComparison, 'summary');

end