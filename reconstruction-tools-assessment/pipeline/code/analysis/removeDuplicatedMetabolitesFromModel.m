function [model, duplicatedMets, removedRxns] = removeDuplicatedMetabolitesFromModel(model)

duplicatedMets = {};
removedRxns = {};

[mets, pos, posAll] = getDuplicatedMetabolites(model);
if ~isempty(mets)
    posAllRemovedMets = [];
    for i = 1:length(mets)
        posRemovedMets = setdiff(posAll{i}, pos(i));
        posAllRemovedMets = union(posAllRemovedMets, posRemovedMets);
        for j = 1:length(posRemovedMets)
            rxnsToBeModified = getRxnsFromMets(model, posRemovedMets(j));
            posKeptMet = pos(i);
            model = changeDuplicatedMetsInRxns(model, rxnsToBeModified, posKeptMet, posRemovedMets(j));
        end
    end
    duplicatedMets = model.mets(posAllRemovedMets);
    model.mets(posAllRemovedMets) = strcat(model.mets(posAllRemovedMets),'_deleted');
    rxns_before = model.rxns;
    model = removeMetabolites(model, model.mets(posAllRemovedMets), 1);
    rxns_after = model.rxns;
    removedRxns = setdiff(rxns_before, rxns_after);
end

end