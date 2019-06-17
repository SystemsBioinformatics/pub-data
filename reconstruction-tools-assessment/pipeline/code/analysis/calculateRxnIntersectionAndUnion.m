function [n_intersection, n_union] = calculateRxnIntersectionAndUnion(model1, model2, option)

if nargin < 3
    option = 1;
end

rxns1 = model1.rxns;
rxns2 = model2.rxns;
n_rxns1 = length(rxns1);
n_rxns2 = length(rxns2);

switch option
    case 1
        
        n_intersection = length(intersect(rxns1,rxns2));
        n_union = length(union(rxns1,rxns2));
    case 2
        if n_rxns1<n_rxns2
            eqs = getRxn_cobraFormat(model1);
            reactionFormulaWasFound = ...
            cellfun(@(x) reactionFormulaInModel(model2, x, 1, 1), eqs,'UniformOutput',0);
        else
            eqs = getRxn_cobraFormat(model2);
            reactionFormulaWasFound = ...
            cellfun(@(x) reactionFormulaInModel(model1, x, 1, 1), eqs,'UniformOutput',0);
        end
        n_intersection = length(find(cell2mat(reactionFormulaWasFound)));
        n_union = n_intersection + (n_rxns1-n_intersection) + (n_rxns2-n_intersection);
end        
        

end