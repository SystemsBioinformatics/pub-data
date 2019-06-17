function JC_sim_r_all = getRxnSimilarityForModels(models, ...
    species, abbr, excelFileName, option)

if nargin < 5
    option = 1;
end

modelRef = models{1};
model1 = modelRef;
JC_sim_r_all = zeros(length(models),1);
for i = 1:length(models)
    fprintf('calculating similary model %2.0f...' ,i);
    model2 = models{i};
    
    JC_sim_r = getJaccardSimilarityRxns(model1, model2, option);
    JC_sim_r_all(i) = JC_sim_r;
    fprintf(' done: progress %2.0f %%\n', 100*i/length(models));
end
info = [abbr', num2cell(JC_sim_r_all)];
xlswrite([excelFileName '_' species],info,'rxns')

end