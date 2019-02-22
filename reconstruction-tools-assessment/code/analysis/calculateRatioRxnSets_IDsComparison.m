function [ratio, correctedRatio, tableRxnComparison] = ...
    calculateRatioRxnSets_IDsComparison(models, names, rxnMatrix, rxnIDs, posExcludedRxns, rxnExcelFileName, species)

if size(names,2)>size(names,1)
    names = names';
end

%calculate coverage
coverage = zeros(length(models)-1,1);
additional = zeros(length(models)-1,1);
sumOfRxns = zeros(length(models)-1,1);
total = zeros(length(models)-1,1);
coverage_p = zeros(length(models)-1,1);
coverage2_p = zeros(length(models)-1,1);
additional_p = zeros(length(models)-1,1);
additional2_p = zeros(length(models)-1,1);

nOriginalModel = length(models{1}.rxns);
n_specific = length(posExcludedRxns);
n_corrected = nOriginalModel - n_specific;

for i = 1:length(models)-1
    additional(i) = length(find(rxnMatrix(nOriginalModel+1:end,i+1)));
    coverage(i) = length(find(rxnMatrix(1:nOriginalModel,i+1)));
    sumOfRxns(i) = coverage(i)+additional(i);
    total(i) = length(models{i+1}.rxns);
    additional_p(i) = 100*length(find(rxnMatrix(nOriginalModel+1:end,i+1)))/nOriginalModel;
    additional2_p(i) = 100*length(find(rxnMatrix(nOriginalModel+1:end,i+1)))/n_corrected;
    
    coverage_p(i) = 100*length(find(rxnMatrix(1:nOriginalModel,i+1)))/nOriginalModel;
    coverage2_p(i) = 100*length(find(rxnMatrix(1:nOriginalModel,i+1)))/(n_corrected);
    
end
ratio = coverage_p./additional_p;
correctedRatio = coverage2_p./additional2_p;

tableLabels = [names(2:end)];
tableRxnComparison = [[{'model'} {'rxn coverage'} {'additional rxns'} {'sum rxns'} {'total rxns model'} {'rxn coverage (%)'} {'rxn coverage, corrected(%)'} {'additional rxns (%)'} {'additional rxns, corrected (%)'} {'ratio coverage/additional'} {'ratio coverage/additional, corrected'}]; ...
    [tableLabels, num2cell(coverage), num2cell(additional), num2cell(sumOfRxns), num2cell(total), num2cell(coverage_p), num2cell(coverage2_p), num2cell(additional_p), num2cell(additional2_p), num2cell(ratio) num2cell(correctedRatio)]];

if exist([rxnExcelFileName '_' species '.xlsx'], 'file') ==2
    delete([rxnExcelFileName '_' species '.xlsx'])
end
xlswrite([rxnExcelFileName '_' species], rxnIDs,'reactions');
xlswrite([rxnExcelFileName '_' species], tableRxnComparison, 'summary');
end