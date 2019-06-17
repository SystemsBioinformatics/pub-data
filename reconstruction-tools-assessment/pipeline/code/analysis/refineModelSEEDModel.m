function model = refineModelSEEDModel(model, inputFastaFile)
[compSymbols, compNames] = getCompSymbols;
[genes, ~, ~, locus_tag] = readFastaProteinsGeneral(inputFastaFile);
genes = regexprep(genes,'lcl\|','');
[~, ind1, ind2] = intersect(genes, model.genes);
model.genes(ind2) = locus_tag(ind1);
% model.genes = setdiff(model.genes,{'Unknown','lcl'}); %ojo, las rules est?n malas
model = creategrRulesField(model);
model = removeGenes(model, {'Unknown'});
model = createRulesFromgrRules(model);

model.rxns = regexprep(model.rxns, {'R_','_e0','_c0'}, {'','',''});
model.mets = regexprep(model.mets,strcat('\[(', compSymbols, ')0\]$'),'_$1');

[~, s] = xlsread(which('compounds.xlsx'));
ids = s(2:end,1);
formulas = s(2:end,4);
formulasMS = cell(length(model.mets),1);
for i = 1:length(model.mets)
    met_wo_comps = regexprep(model.mets{i}, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
    pos = find(strcmp(ids, met_wo_comps));
    if ~isempty(pos) && length(pos)==1
        formulasMS{i} = formulas{pos};
    end
end
model.metFormulas = formulasMS;
model.comps{1} = 'c';
model.comps{2} = 'e';
model.compNames{1} = 'cytosol';
model.compNames{2} = 'extracelular space';
for j = 1:length(model.metFormulas)
    if isempty(model.metFormulas{j}) || ~isempty(strfind(model.metFormulas{j}, 'R')) || ~isempty(strfind(model.metFormulas{j}, 'X'))
        model.metCharges(j) = nan;
    end
end

end