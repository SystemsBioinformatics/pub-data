function [ratio, correctedRatio,correctedRatio2,correctedRatio3, tableRxnComparison] = ...
    calculateRatioRxnSets(models, names, rxnMatrix, rxnIDs, rxnMatrix2, rxnIDs2, Eqs2, posExcludedRxns, rxnExcelFileName, species)

if size(names,2)>size(names,1)
    names = names';
end

%calculate coverage
coverage = zeros(length(models)-1,1);
partialCoverage = zeros(length(models)-1,1);
totalCoverage = zeros(length(models)-1,1);
sumOfRxns = zeros(length(models)-1,1);
total = zeros(length(models)-1,1);

additional = zeros(length(models)-1,1);
partialCoverage_p = zeros(length(models)-1,1);
partialCoverage_p2 = zeros(length(models)-1,1);
coverage_p = zeros(length(models)-1,1);
coverage2_p = zeros(length(models)-1,1);
coverage3_p = zeros(length(models)-1,1);
coverage4_p = zeros(length(models)-1,1);
additional_p = zeros(length(models)-1,1);
additional2_p = zeros(length(models)-1,1);


nOriginalModel = length(models{1}.rxns);
n_specific = length(posExcludedRxns);
n_corrected = nOriginalModel - n_specific;

for i = 1:length(models)-1
    total(i) = length(models{i+1}.rxns);
    additional(i) = length(find(rxnMatrix(nOriginalModel+1:end,i+1)));
    coverage(i) = length(find(rxnMatrix(1:nOriginalModel,i+1)));
    partialCoverage(i) = length(find(rxnMatrix2(1:nOriginalModel,i+1)));
    totalCoverage(i) = coverage(i) + partialCoverage(i);
    sumOfRxns(i) = totalCoverage(i)+additional(i);
    
    additional_p(i) = 100*additional(i)/nOriginalModel;
    additional2_p(i) = 100*additional(i)/n_corrected;   
    partialCoverage_p(i) = 100*partialCoverage(i)/nOriginalModel;
    partialCoverage_p2(i) = 100*partialCoverage(i)/(n_corrected);
    coverage_p(i) = 100*coverage(i)/nOriginalModel;
    coverage2_p(i) = 100*coverage(i)/(n_corrected);
    coverage3_p(i) = 100*totalCoverage(i)/(nOriginalModel);
    coverage4_p(i) = 100*totalCoverage(i)/(n_corrected);
end
ratio = coverage_p./additional_p;
correctedRatio = coverage2_p./additional2_p;
correctedRatio2 = coverage3_p./additional_p;
correctedRatio3 = coverage4_p./additional2_p;

tableLabels = names(2:end);

tableRxnComparison = [[{'model'} {'full rxn coverage'} {'partial rxn coverage'} {'total rxn coverage'}, {'additional rxns'} {'sum rxns'} {'total rxns model'} {'full rxns coverage(%)'} {'full rxns coverage, corrected(%)'}  {'partial rxns coverage(%)'} {'partial rxns coverage, corrected(%)'} {'total rxn coverage(%)'} {'total rxn coverage, corrected(%)'} {'additional rxns (%)'} {'additional rxns, corrected (%)'} {'ratio full coverage/additional'} {'ratio full coverage/additional, corrected'} {'ratio total coverage/additional'} {'ratio total coverage/additional, corrected'}]; ...
    [tableLabels, num2cell(coverage), num2cell(partialCoverage), num2cell(totalCoverage) , num2cell(additional), num2cell(sumOfRxns), num2cell(total), num2cell(coverage_p), num2cell(coverage2_p),  num2cell(partialCoverage_p), num2cell(partialCoverage_p2), num2cell(coverage3_p), num2cell(coverage4_p), num2cell(additional_p), num2cell(additional2_p), num2cell(ratio), num2cell(correctedRatio) num2cell(correctedRatio2) num2cell(correctedRatio3)]];

if exist([rxnExcelFileName '_' species '.xlsx'], 'file') ==2
    delete([rxnExcelFileName '_' species '.xlsx'])
end
xlswrite([rxnExcelFileName '_' species], rxnIDs,'reactions');
xlswrite([rxnExcelFileName '_' species], rxnIDs2,'partial reactions');
xlswrite([rxnExcelFileName '_' species], Eqs2,'partial equations');
xlswrite([rxnExcelFileName '_' species], tableRxnComparison, 'summary');
end