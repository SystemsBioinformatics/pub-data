function [dictionaryIsBeneficial, beneficial_pairs, n_rxns_improved, dictionaryIsBeneficial_perPair] = canBenefitFromDictionary(model, modelRef, dictionary)

dictionary = dictionary(2:end,:);
dictionaryIsBeneficial_perPair = zeros(size(dictionary,1),1);
dictionaryIsBeneficial = 0;
beneficial_pairs = {};
n_rxns_improved = 0;

for i = 1:size(dictionary,1)
    pos1 = find(strcmp(model.mets, dictionary{i,1}));
    if ~isempty(pos1)
        posRxns = find(model.S(pos1,:));
        if ~isempty(posRxns)
            eqs1 = getRxn_cobraFormat(model, posRxns);
            rxnFormulaWasFound = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqs1,'UniformOutput',0);
            n_rxns_recovered_wo_dict = length(find(cell2mat(rxnFormulaWasFound)));
            
            %
            model2 = model;
            model2.mets{pos1} = dictionary{i,2};
            eqs2 = getRxn_cobraFormat(model2, posRxns);
            rxnFormulaWasFound2 = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqs2,'UniformOutput',0);
            n_rxns_recovered_w_dict = length(find(cell2mat(rxnFormulaWasFound2)));
            if n_rxns_recovered_w_dict>n_rxns_recovered_wo_dict
                dictionaryIsBeneficial_perPair(i) = 1;
            elseif n_rxns_recovered_w_dict<n_rxns_recovered_wo_dict
                dictionaryIsBeneficial_perPair(i) = -1;
            end
        end
    end
end

posAllRxnsInvolved = [];
metKeys= dictionary(:,1);
metValues= dictionary(:,2);
metKeysInModel = metKeys(ismember(metKeys, model.mets));
metValuesInModel = metValues(ismember(metKeys, model.mets));

for i = 1:length(metKeysInModel)
    pos1 = find(strcmp(model.mets, metKeysInModel{i}));
    if ~isempty(pos1)
        posRxns = find(model.S(pos1,:));
        if size(posRxns,2)>size(posRxns,1)
           posRxns = posRxns'; 
        end
        if ~isempty(posRxns)
            posAllRxnsInvolved = [posAllRxnsInvolved; posRxns];
        end
    end
end

if isempty(posAllRxnsInvolved)
   return; 
end

model2 = model;
eqsInit = getRxn_cobraFormat(model2, posAllRxnsInvolved);
rxnFormulaWasFound = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqsInit,'UniformOutput',0);
n_rxns_recoved_inital = length(find(cell2mat(rxnFormulaWasFound))); 
keepMetKeys = zeros(size(metKeysInModel));
for i = 1:length(metKeysInModel)
    modelAux = model2;
    
    pos1 = find(strcmp(model.mets,  metKeysInModel{i}));
    if ~isempty(pos1)
        eqs1 = getRxn_cobraFormat(model2, posAllRxnsInvolved);
        rxnFormulaWasFound = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqs1,'UniformOutput',0);
        n_rxns_recovered_wo_dict = length(find(cell2mat(rxnFormulaWasFound)));
        model2.mets{pos1} = metValuesInModel{i};
        eqs2 = getRxn_cobraFormat(model2, posAllRxnsInvolved);
        rxnFormulaWasFound2 = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqs2,'UniformOutput',0);
        n_rxns_recovered_w_dict = length(find(cell2mat(rxnFormulaWasFound2)));
        if n_rxns_recovered_w_dict>=n_rxns_recovered_wo_dict
            keepMetKeys(i) = 1;
        else
            model2 = modelAux;
        end  
    end
end

eqsFinal = getRxn_cobraFormat(model2, posAllRxnsInvolved);
rxnFormulaWasFound = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqsFinal,'UniformOutput',0);
n_rxns_recoved_final = length(find(cell2mat(rxnFormulaWasFound))); 
if n_rxns_recoved_final-n_rxns_recoved_inital>0
    dictionaryIsBeneficial = 1;
    n_rxns_improved = n_rxns_recoved_final-n_rxns_recoved_inital;
end

posConditional = intersect(find(dictionaryIsBeneficial_perPair ==0), find(ismember(metKeys, model.mets)));
conditional = dictionary(posConditional,:);
% keepConditional = zeros(size(conditional,1),1);
for i = 1:size(conditional,1)
    %apply all mets except 
    posOthers = setdiff(find(ismember(metKeys, model.mets)), posConditional(i));
    model2 = model;
    for j = 1:length(posOthers)
        pos1 = find(strcmp(model.mets,  dictionary{posOthers(j),1}));
        model2.mets{pos1} = dictionary{posOthers(j),2};
    end
    eqsBefore = getRxn_cobraFormat(model2, posAllRxnsInvolved);
    rxnFormulaWasFound = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqsBefore,'UniformOutput',0);
    n_rxns_recoved_before = length(find(cell2mat(rxnFormulaWasFound)));
    
    %apply last change
    pos1 = find(strcmp(model.mets,  dictionary{posConditional(i),1}));
    model2.mets{pos1} = dictionary{posConditional(i),2};
    
    eqsAfter = getRxn_cobraFormat(model2, posAllRxnsInvolved);
    rxnFormulaWasFound = cellfun(@(x) reactionFormulaInModel(modelRef, x, 1, 1), eqsAfter,'UniformOutput',0);
    n_rxns_recoved_after = length(find(cell2mat(rxnFormulaWasFound)));
    
    if n_rxns_recoved_after<=n_rxns_recoved_before
        pos = find(strcmp(metKeysInModel, conditional{i}));
        keepMetKeys(pos) = 0;
    end
end

beneficial_pairs = [metKeysInModel(keepMetKeys==1), metValuesInModel(keepMetKeys==1)];

end