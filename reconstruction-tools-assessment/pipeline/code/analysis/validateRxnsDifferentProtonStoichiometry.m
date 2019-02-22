function valid = validateRxnsDifferentProtonStoichiometry(rxnEq1, rxnEq2)
valid = 1;
% disp(rxnEq1)
[metaboliteList1, stoichCoeffList1, ~] = parseRxnFormula(rxnEq1);
posHs1 = find(~cellfun(@isempty, regexp(metaboliteList1, '^h_')));
comps1 = getCompartmentsFromMetList(metaboliteList1(posHs1));
unComps1 = unique(comps1);

[metaboliteList2, stoichCoeffList2, ~] = parseRxnFormula(rxnEq2);
posHs2 = find(~cellfun(@isempty, regexp(metaboliteList2, '^h_')));
comps2 = getCompartmentsFromMetList(metaboliteList2(posHs2));
unComps2 = unique(comps2);

if length(unComps1)~=length(unComps2) && (length(unComps2)>=2 || length(unComps1)>=2)
    valid = 0;
end

end