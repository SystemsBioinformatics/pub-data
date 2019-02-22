function [model, removedRxns, keptRxns] = removeDuplicatedRxnsFromModel(model, makeReport, fileName)
removedRxns = {};
keptRxns = {};
if makeReport
    Info = [{'removed reactions'},{'rxn equation'}, {'kept reactions'},{'rxn equation'}];
    Info2 = [{'action'}, {'removed rxn ID'}, {'removed rxn rule'}, {'move rule to rxn ID'}, {'kept rxn rule'}];
end

posEliminadas = zeros(length(model.rxns),1);
cont = 0;
%elminate empty reactions
for i=1:length(model.rxns)
    if ~any(full(model.S(:,i)))
        cont = cont + 1;
        posEliminadas(cont) = i;
    end
end
if cont>0
    posEliminadas = posEliminadas(1:cont);
    rxnsEliminar = model.rxns(unique(posEliminadas));
    model = removeRxns(model, rxnsEliminar);
    Info = [Info; [rxnsEliminar repmat({''},length(rxnsEliminar),1)]];
end

redundantSets = cell(length(model.rxns),1);
cont = 0;
duplicated = [];
for i = 1:length(model.rxns)
    if ismember(i, duplicated); continue; end

    [redundant, pos] = isRxnDuplicated(model, i, 1);
    if redundant
        cont = cont + 1;
        redundantSets{cont} = [i;pos];
        duplicated = union(duplicated, [i;pos]);
    end
end
redundantSets = redundantSets(1:cont);

posEliminadas = [];
posMantenidas = [];
noEliminadas = [];
for i = 1:length(redundantSets)
    set_i = redundantSets{i};
    rev = arrayfun(@(x) model.lb(x) < 0 && model.ub(x) > 0, set_i);
    if any(rev)
        posRev = set_i(find(rev));
        posRev = posRev(1);
        posEliminadas_i = setdiff(set_i, posRev);
        posEliminadas = [posEliminadas; posEliminadas_i];
        posMantenidas = [posMantenidas; repmat(posRev,length(posEliminadas_i),1)];
        noEliminadas = [noEliminadas; posRev];
        if makeReport
            for j = 1:length(posEliminadas_i)
                if ~isempty(model.grRules{posEliminadas_i(j)})
                    Info2 = [Info2; [{'rule transfer'}, model.rxns(posEliminadas_i(j)), model.grRules(posEliminadas_i(j)),  model.rxns(posRev), model.grRules(posRev)]];
                end
            end
        end
        
    else
        posEliminadas_i =setdiff(set_i, set_i(1));
        posEliminadas = [posEliminadas; posEliminadas_i];
        posMantenidas = [posMantenidas; repmat(set_i(1),length(posEliminadas_i),1)];
        noEliminadas = [noEliminadas; set_i(1)];
        if makeReport
            for j = 1:length(posEliminadas_i)
                if ~isempty(model.grRules{posEliminadas_i(j)})
                    Info2 = [Info2; [{'rule transfer'}, model.rxns(posEliminadas_i(j)), model.grRules(posEliminadas_i(j)),  model.rxns(set_i(1)), model.grRules(set_i(1))]];
                end
            end
        end
    end
end
if ~isempty(posEliminadas)
    eq1 = getRxn_cobraFormat(model, posEliminadas);
    eq2 = getRxn_cobraFormat(model, posMantenidas);
    if makeReport
        Info = [Info; [model.rxns(posEliminadas), eq1, model.rxns(posMantenidas), eq2]];
    end
    removedRxns = model.rxns(posEliminadas);
    keptRxns = model.rxns(posMantenidas);
    model = removeReactions(model, removedRxns);
end

if makeReport 
    if exist([fileName '.xls'], 'file')==2
        delete([fileName '.xls'])
    end
    xlswrite(fileName, Info, 'removed Rxns')
    xlswrite(fileName, Info2,'rules')
end