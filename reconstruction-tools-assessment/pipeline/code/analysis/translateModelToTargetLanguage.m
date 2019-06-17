function [modelOut, removedEmptyRxnsAfterTranslation, keptRxns, n_translated,...
    p_translated, dictionary, metsTranslatedIntoTheSameID, requireManualCheck] = ...
    translateModelToTargetLanguage(model, inputLanguage, outputLanguage, modelRef, ...
    metOtherIDs, metMNXIDs, keepOriginalRxnNames, keepDuplicatedRxns, relaxCompartments, useInternalIDs)

global rootFolder

[compSymbols, ~] = getCompSymbols;
removedEmptyRxnsAfterTranslation = {};
keptRxns = {};
metsTranslatedIntoTheSameID = {};
requireManualCheck = cell(1000,4); for i = 1:length(requireManualCheck); requireManualCheck{i} = ''; end;
n_requireManualCheck = 0;
n_translated = 0;
p_translated = 0;
dictionary = {};

if nargin<7; keepOriginalRxnNames = 1; end
if nargin<8; keepDuplicatedRxns = 0; end;
if nargin<9; relaxCompartments = 1; end;
if nargin<10; useInternalIDs = 0; end;

load(fullfile(rootFolder, 'MNX', 'rxnLanguages'));
load(fullfile(rootFolder, 'MNX', 'metLanguages'));
possibleLanguages = union(intersect(metLanguages, rxnLanguages),'MNX');
if ~ismember(inputLanguage, possibleLanguages) || ~ismember(outputLanguage, possibleLanguages)
    return;
end
if strcmp(inputLanguage, outputLanguage); 
    modelOut = model;
    
    return; 
end

if nargin<5
    load(fullfile(rootFolder, 'MNX', 'rxnOtherIDs'));
    load(fullfile(rootFolder, 'MNX', 'rxnMNXIDs'));
    load(fullfile(rootFolder, 'MNX', 'metOtherIDs'));
    load(fullfile(rootFolder, 'MNX', 'metMNXIDs'));
end

dictionary = cell(5000,4);

modelOut = model;

modelUsesSB = 0;
if ~isempty(strfind(modelOut.mets{1},'['))
    modelUsesSB = 1;
end
mets_wo_comps = regexprep(modelOut.mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
mets_wo_comps_Ref = regexprep(modelRef.mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
translations = getMetIDsInTargetLanguageFromInputLanguage(mets_wo_comps,inputLanguage,outputLanguage, metOtherIDs, metMNXIDs);
if useInternalIDs 
    posEmpty = find(cellfun(@isempty, translations));
    if isfield(model, 'metMetaNetXID')
        posMNXNotEmpy = find(cellfun(@isempty, model.metMetaNetXID)==0);
    else
        posMNXNotEmpy = [];
    end
    
    if (strcmp(outputLanguage,'bigg') && isfield(model,'metBiGGID')) ||...
            (strcmp(outputLanguage,'kegg') && isfield(model,'metKEGGID')) ||...
            (strcmp(outputLanguage,'seed') && isfield(model,'metSEEDID')) ||...
            (strcmp(outputLanguage,'metacyc') && isfield(model,'metMetaCycID'))
        switch outputLanguage
            case 'bigg'
                posNotEmpySpecificOutputLanguage = find(cellfun(@isempty, model.metBiGGID)==0);
            case 'kegg'
                posNotEmpySpecificOutputLanguage = find(cellfun(@isempty, model.metKEGGID)==0);
            case 'seed'
                posNotEmpySpecificOutputLanguage = find(cellfun(@isempty, model.metSEEDID)==0);
            case 'metacyc'
                posNotEmpySpecificOutputLanguage = find(cellfun(@isempty, model.metMetaCycID)==0);
        end
    else
        posNotEmpySpecificOutputLanguage = [];
    end
    posInt = intersect(posEmpty,posMNXNotEmpy);
    if ~isempty(posInt)
        MNXids = model.metMetaNetXID(posInt);
        MNXtranslations = getMetIDsInTargetLanguageFromInputLanguage(MNXids,'MNX',outputLanguage, metOtherIDs, metMNXIDs);
        for i = 1:length(MNXtranslations)
            if ~isempty(MNXtranslations{i})
                MNXtranslations{i} = strsplit(MNXtranslations{i}{1},';');
            end
        end
        translations(posInt) = MNXtranslations;
    end
    if ~isempty(setdiff(posNotEmpySpecificOutputLanguage, posMNXNotEmpy))
        posInt = intersect(posEmpty,posNotEmpySpecificOutputLanguage);
        if ~isempty(posInt)
            switch outputLanguage
                case 'bigg'
                    specificLanguageIDs = model.metBiGGID(posInt);
                case 'kegg'
                    specificLanguageIDs = model.metKEGGID(posInt);
                case 'seed'
                    specificLanguageIDs = model.metSEEDID(posInt);
                case 'metacyc'
                    specificLanguageIDs = model.metMetaCycID(posInt);
            end
            MNXtranslations = getMetIDsInTargetLanguageFromInputLanguage(specificLanguageIDs,inputLanguage,outputLanguage, metOtherIDs, metMNXIDs);
            for i = 1:length(MNXtranslations)
                if ~isempty(MNXtranslations{i})
                    MNXtranslations{i} = strsplit(MNXtranslations{i}{1},';');
                end
            end
            translations(posInt) = MNXtranslations;
        end
    end
    
end

comps = getMetCompartmentsFromModel(modelOut);
mets = modelOut.mets;
metNames = modelOut.metNames;
if isfield(modelOut, 'metSIDs')
    metSIDs = modelOut.metSIDs;
end
n_translated = 0;
rmMets = {};
for i = 1:length(mets)
    if ~isempty(translations{i})
        if relaxCompartments
            inModel = ismember(translations{i},mets_wo_comps_Ref);
        else
            inModel = ismember(strcat(translations{i}, '_', comps{i}),modelRef.mets);
        end
        if any(inModel)
            specific_translations = translations{i}(inModel);
            if length(specific_translations)==1
                n_translated = n_translated +1;
                dictionary{n_translated,1} = mets{i};
                if  modelUsesSB
                    newName = strcat(specific_translations{1}, '[', comps{i},']');
                else
                    newName = strcat(specific_translations{1}, '_', comps{i});
                end
                dictionary{n_translated,2} = newName;
                dictionary{n_translated,3} = metNames{i};
                if isfield(modelOut, 'metSIDs')
                    dictionary{n_translated,4} = metSIDs{i};
                end
                modelOut = changeMetIdentifier(modelOut, mets{i}, newName);
                if length(modelOut.mets)<length(model.mets)
                    disp('')
                end
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
    if relaxCompartments
        inModel = ismember(translations{pos},mets_wo_comps_Ref);
    else
        inModel = ismember(strcat(translations{pos}, '_', comps{pos}),modelRef.mets);
    end
    specific_translations = translations{pos}(inModel);
    rxns1 = getRxnEquationFromMets(modelOut, rmMets(i),0);
    compartment = getCompartmentsFromMetList(rmMets(i));
    hits =  cellfun(@(y) length(find(cellfun(@(x) reactionFormulaInModel(modelRef, x, 1), regexprep(rxns1{1},{[' ' rmMets{i} '$'],[' ' rmMets{i} ' ']}, {[' ' y '_' compartment{1}],[' ' y '_' compartment{1} ' ']})))), specific_translations);
    [maximum, maxPos] = max(hits);
    if length(find(hits==maximum)) == 1
        dictionary{n_translated,1} = rmMets{i};
        n_translated = n_translated +1;
        if  modelUsesSB
            newName = strcat(specific_translations{maxPos}, '[', comps{pos},']');
        else
            newName = strcat(specific_translations{maxPos}, '_', comps{pos});
        end
        dictionary{n_translated,2} = newName;
        dictionary{n_translated,3} = metNames{i};
        if isfield(modelOut, 'metSIDs')
            dictionary{n_translated,4} = metSIDs{i};
        end
        modelOut = changeMetIdentifier(modelOut, rmMets{i}, newName);
        if i ==1
            rmMets = rmMets(2:end);
        elseif i == length(rmMets)
            rmMets = rmMets(1:end-1);
        else
            rmMets = [rmMets(1:i-1); rmMets(i+1:end)];
        end
        i = 0;
        dioVueltaSinCambios = 0;
    else
        if ~ismember(rmMets{i},requireManualCheck(:,1))
            n_requireManualCheck = n_requireManualCheck+1;
            a = rmMets(i);
            b = {strjoin(specific_translations,',')};
            c = model.metNames(pos);
            if isfield(modelOut, 'metSIDs')
                d = model.metSIDs(pos);
            else
                d = {''};
            end
            info = [a, b, c, d];
            requireManualCheck(n_requireManualCheck,:) = info;
        end
    end
    %     rxns2 =
    if (i == length(rmMets) && dioVueltaSinCambios) || isempty(rmMets)
        break;
    end
end
requireManualCheck = requireManualCheck(1:n_requireManualCheck,:);

posTranslated = find(ismember(requireManualCheck(:,1),modelOut.mets)==0);
requireManualCheck(posTranslated,:) = [];

if ~keepDuplicatedRxns
    [modelOut, removedRxns, keptRxns] = removeDuplicatedRxnsFromModel(modelOut,0);
end

if ~keepOriginalRxnNames
    if ~keepDuplicatedRxns && ~isempty(removedRxns)
        removedEmptyRxnsAfterTranslation = removedRxns;
        disp('')
    end
    eqs = getRxn_cobraFormat(modelOut, 1:length(modelOut.rxns));
    for i = 1:length(modelOut.rxns)
        [reactionFormulaWasFound, posRxn] = reactionFormulaInModel(modelRef, eqs{i} , 1, 0);
        if reactionFormulaWasFound
            modelOut.rxns{i} = modelRef.rxns{posRxn};
        end
    end
end
p_translated = n_translated/length(mets);
dictionary = dictionary(1:n_translated,:);

%find merged models
[duplicated, indices, kept] = getDuplicated(dictionary(:,2));
if ~isempty(duplicated)
    metsTranslatedIntoTheSameID = cell(length(duplicated),2);
    for i = 1:length(duplicated)
        metsTranslatedIntoTheSameID{i,1} = dictionary(union(indices(i),kept(i)),1);
        metsTranslatedIntoTheSameID{i,2} = duplicated{i};
        pos=cell2mat(arrayfun(@(x)find(strcmp(x,model.mets)),metsTranslatedIntoTheSameID{i,1},'UniformOutput',false))'; 
        metsTranslatedIntoTheSameID{i,3} = model.metNames(pos);
        if isfield(model, 'metSIDs')
            metsTranslatedIntoTheSameID{i,4} = metSIDs(pos);
        end
    end
end

end