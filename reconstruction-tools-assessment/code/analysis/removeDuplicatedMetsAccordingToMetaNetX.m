function [modelWithOutDuplicatedMets, removedMets, modelWithOutDuplicatedMets_translated] = removeDuplicatedMetsAccordingToMetaNetX(model, idsType, metOtherIDs, metMNXIDs)
compSymbols = getCompSymbols;
mets = model.mets;
mets_wo_comps = regexprep(mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
metComp = cellfun(@(x) regexp(x,'.*_(.*)$','tokens'), model.mets, 'UniformOutput', 0);
for i = 1:length(metComp); metComp{i} = metComp{i}{1}{1}; end;

removedMets = {};
modelWithOutDuplicatedMets = model;

translations = cell(size(mets));
for i = 1:length(mets)
    pos = find(strcmp(metOtherIDs, [idsType ':' mets_wo_comps{i}]));
    if ~isempty(pos)
        translations{i} = [metMNXIDs{pos} '_' metComp{i}];
    end
end

[~, indices, kept] = getDuplicated(translations);
if ~isempty(indices)
    modelWithOutDuplicatedMets = model;
    removedMets = mets(indices);
    keptMets = mets(kept);
    for i = 1:length(removedMets)
        modelWithOutDuplicatedMets = changeMetIdentifier(modelWithOutDuplicatedMets, removedMets{i}, keptMets{i});
    end
    notEmptyNotDuplicated = setdiff(1:length(translations), union(find(cellfun(@isempty, translations)), indices));
else
    notEmptyNotDuplicated = setdiff(1:length(translations), find(cellfun(@isempty, translations)));
    
end

keys = mets(notEmptyNotDuplicated);
values = translations(notEmptyNotDuplicated);
pos = cell2mat(arrayfun(@(x)find(strcmp(x,modelWithOutDuplicatedMets.mets)),keys,'UniformOutput',false))';
modelWithOutDuplicatedMets_translated = modelWithOutDuplicatedMets;
modelWithOutDuplicatedMets_translated.mets(pos) = values;

end