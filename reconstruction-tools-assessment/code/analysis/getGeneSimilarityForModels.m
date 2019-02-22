function JC_sim_g_all = getGeneSimilarityForModels(models, ...
    species, abbr, excelFileName)

modelRef = models{1};

JC_sim_g_all = zeros(length(models),1);
for i = 1:length(models)
    fprintf('calculating similary model %2.0f...' ,i);
    model1 = modelRef;
    model2 = models{i};
    
    JC_sim_g = getJaccardSimilarityGenes(model1, model2);
    JC_sim_g_all(i) = JC_sim_g;
    fprintf(' done: progress %2.0f %%\n', 100*i/length(models));
end
info = [abbr', num2cell(JC_sim_g_all)];
xlswrite([excelFileName '_' species],info,'genes')

end