function [ratio, correctedRatio, tableMetComparison] = calculateRatioMetSets(models, names, metMatrix, metIDs, posExcludedMets, metExcelFileName, species)

%calculate coverage
coverage = zeros(length(models)-1,1);
additional = zeros(length(models)-1,1);
coverage_p = zeros(length(models)-1,1);
coverage2_p = zeros(length(models)-1,1);
additional_p = zeros(length(models)-1,1);
additional2_p = zeros(length(models)-1,1);

nOriginalModel = length(models{1}.mets);
n_specific = length(posExcludedMets);
n_corrected = nOriginalModel - n_specific;

for i = 1:length(models)
    additional(i) = length(find(metMatrix(nOriginalModel+1:end,i)));
    coverage(i) = length(find(metMatrix(1:nOriginalModel,i)));
    additional_p(i) = 100*length(find(metMatrix(nOriginalModel+1:end,i)))/nOriginalModel;
    coverage_p(i) = 100*length(find(metMatrix(1:nOriginalModel,i)))/nOriginalModel;
    
    additional2_p(i) = 100*length(find(metMatrix(nOriginalModel+1:end,i)))/n_corrected;
    coverage2_p(i) = 100*length(find(metMatrix(1:nOriginalModel,i)))/(n_corrected);
end
ratio = coverage_p./additional_p;
correctedRatio = coverage2_p./additional2_p;

if size(names,2)>size(names,1)
    names = names';
end

tableLabels = names;

tableMetComparison = [[{'model'} {'met coverage'} {'additional mets'} {'met coverage(%)'} {'met coverage(%)'} {'additional mets(%)'} {'ratio coverage/additional'} {'ratio coverage/additional'}]; ...
    [tableLabels, num2cell(coverage), num2cell(additional), num2cell(coverage_p), num2cell(coverage2_p), num2cell(additional_p), num2cell(ratio), num2cell(correctedRatio)]];

if exist([metExcelFileName '_' species '.xlsx'], 'file') ==2
    delete([metExcelFileName '_' species '.xlsx'])
end
xlswrite([metExcelFileName '_' species], metIDs,'metabolites');
xlswrite([metExcelFileName '_' species], tableMetComparison, 'summary');

end
