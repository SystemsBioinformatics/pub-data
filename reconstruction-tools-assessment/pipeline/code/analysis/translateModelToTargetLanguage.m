function [model, removedEmptyRxnsAfterTranslation, keptRxns] = translateModelToTargetLanguage(model, inputLanguage, outputLanguage, modelRef, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, keepOriginalRxnNames)
[compSymbols, ~] = getCompSymbols;
removedEmptyRxnsAfterTranslation = {};
keptRxns = {};

if nargin<9
    keepOriginalRxnNames = 1;
end

load('D:\Dropbox\Databases\MNX\rxnLanguages');
load('D:\Dropbox\Databases\MNX\metLanguages');
possibleLanguages = intersect(metLanguages, rxnLanguages);
if ~ismember(inputLanguage, possibleLanguages) || ~ismember(outputLanguage, possibleLanguages)
    return;
end

if nargin<5
    load('D:\Dropbox\Databases\MNX\rxnOtherIDs');
    load('D:\Dropbox\Databases\MNX\rxnMNXIDs');
    load('D:\Dropbox\Databases\MNX\metOtherIDs');
    load('D:\Dropbox\Databases\MNX\metMNXIDs');
end

mets_wo_comps = regexprep(model.mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
translations = getMetIDsInTargetLanguageFromInputLanguage(mets_wo_comps,inputLanguage,outputLanguage, metOtherIDs, metMNXIDs);
comps = getMetCompartments(model);
mets = model.mets;
rmMets = {};
for i = 1:length(mets)
    if ~isempty(translations{i})
        inModel = ismember(strcat(translations{i}, '_', comps{i}),modelRef.mets);
        if any(inModel)
            specific_translations = translations{i}(inModel);
            if length(specific_translations)==1
                newName = strcat(specific_translations{1}, '_', comps{i});
                model = changeMetIdentifier(model, mets{i}, newName);
            else
                rmMets = [rmMets; mets(i)];
            end
        end
    end
end

i = 0;
while ~isempty(rmMets)
    i = i + 1;
    if i ==1; dioVueltaSinCambios = 1; end

    pos = find(strcmp(mets, rmMets(i)));
    inModel = ismember(strcat(translations{pos}, '_', comps{pos}),modelRef.mets);
    specific_translations = translations{pos}(inModel);
    rxns1 = getRxnEquationFromMets(model, rmMets(i),0);
    hits =  cellfun(@(y) length(find(cellfun(@(x) reactionFormulaInModel(modelRef, x, 1), regexprep(rxns1{1},[' ' rmMets{i} ' '], [' ' y '_c '])))), specific_translations);
    [maximum, maxPos] = max(hits);
    if length(find(hits==maximum)) == 1
        newName = strcat(specific_translations{maxPos}, '_', comps{pos});
        model = changeMetIdentifier(model, rmMets{i}, newName);
        if i ==1
            rmMets = rmMets(2:end);
        elseif i == length(rmMets)
            rmMets = rmMets(1:end-1);
        else
            rmMets = [rmMets(1:i-1) rmMets(i+1:end)];
        end
        i = 0;
        dioVueltaSinCambios = 0;
    end
    %     rxns2 =
    if (i == length(rmMets) && dioVueltaSinCambios) || isempty(rmMets)
        break;
    end
end

[model, removedRxns, keptRxns] = removeDuplicatedRxnsFromModel(model,0);

if ~keepOriginalRxnNames
    if ~isempty(removedRxns)
        removedEmptyRxnsAfterTranslation = removedRxns;
        disp('')
    end
    eqs = getRxn_cobraFormat(model, 1:length(model.rxns));
    for i = 1:length(model.rxns)
        [reactionFormulaWasFound, posRxn] = reactionFormulaInModel(modelRef, eqs{i} , 1, 0);
        if reactionFormulaWasFound
            model.rxns{i} = modelRef.rxns{posRxn};
        end
    end
end

end