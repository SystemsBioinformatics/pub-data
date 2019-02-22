function [allMetrics, table] = getTopologicalAnalysis(models, labels, species, excelFileName)

allMetrics = cell(size(models));

for i = 1:length(models)
    fprintf('calculating topological metrics: progress %2.0f %%\n', 100*i/length(models));
    model_i = models{i};
    n_mets = length(model_i.mets);
    
    if ~isfield(model_i, 'rules') && isfield(model_i, 'grRules')
        model_i = createRulesFromgrRules(model_i);
    end
    
    [inDegree, outDegree] = calculateInAndOutDegree(model_i);
    orphansRxns_wrtGenes = findOrphanRxns(model_i);
    deadEnds = detectDeadEnds(model_i);
    [metConn, rxnParticipation] = networkTopology(model_i);
   
    orphansRxns_wrtMets = find(arrayfun(@(x) all(ismember(find(model_i.S(:, x)), find(metConn==1))), 1:length(model_i.rxns)));  
    
    metrics.n_deadEnds = length(deadEnds);
    metrics.n_orphansRxns_wrtGenes = length(orphansRxns_wrtGenes);
    metrics.n_orphansRxns_wrtMets = length(orphansRxns_wrtMets);
    
    metrics.ratio_deadEnds = length(deadEnds)/n_mets;
    metrics.ratio_orphansRxns_wrtGenes = length(orphansRxns_wrtGenes)/n_mets;
    metrics.ratio_orphansRxns_wrtMets = length(orphansRxns_wrtMets)/n_mets;
    
    metrics.n10_metConn = prctile(sort(metConn), 10);
    metrics.n50__metConn = prctile(sort(metConn), 50);
    metrics.n90__metConn = prctile(sort(metConn), 90);
    metrics.n95__metConn = prctile(sort(metConn), 95);
    metrics.n99__metConn = prctile(sort(metConn), 99);
    
    metrics.n10_inDegree = prctile(sort(inDegree), 10);
    metrics.n50__inDegree = prctile(sort(inDegree), 50);
    metrics.n90__inDegree = prctile(sort(inDegree), 90);
    metrics.n95__inDegree = prctile(sort(inDegree), 95);
    metrics.n99__inDegree = prctile(sort(inDegree), 99);
    
    metrics.n10_outDegree = prctile(sort(outDegree), 10);
    metrics.n50__outDegree = prctile(sort(outDegree), 50);
    metrics.n90__outDegree = prctile(sort(outDegree), 90);
    metrics.n95__outDegree = prctile(sort(outDegree), 95);
    metrics.n99__outDegree = prctile(sort(outDegree), 99);
    
    metrics.n10_rxnParticipation = prctile(sort(rxnParticipation), 10);
    metrics.n50__rxnParticipation = prctile(sort(rxnParticipation), 50);
    metrics.n90__rxnParticipation = prctile(sort(rxnParticipation), 90);
    metrics.n95__rxnParticipation = prctile(sort(rxnParticipation), 95);
    metrics.n99__rxnParticipation = prctile(sort(rxnParticipation), 99);
    
    allMetrics{i} = metrics;
    
end

names = fieldnames(allMetrics{1});
table = cell(length(names)+1, length(models)+1);
for i = 1:length(names)
    table{i+1,1} = names{i};
end

for i = 1:length(labels)
    table{1,i+1} = labels{i};
end

for i = 1:length(names)
    for j = 1:length(models)
        table{i+1,j+1} = getfield(allMetrics{j}, names{i});
    end
end

if exist([excelFileName '_' species '.xlsx'], 'file') ==2
    delete([excelFileName '_' species '.xlsx'])
end
xlswrite([excelFileName '_' species], table);


end