function [model,mediaExchangeRxns,nutrients] = setMediaFromExcelFileWithRxnIDs(model, fileName,sheet)


[~,s] = xlsread(fileName,sheet);
mediaExchangeRxns = {};
for i =1:size(s,1)
    mediaExchangeRxns = union(mediaExchangeRxns,split(s(i,2),';'));
end
mediaExchangeRxns = setdiff(mediaExchangeRxns,'');

valueMediaExchangeRxns = 10*ones(size(mediaExchangeRxns));

if isempty(find(cellfun(@isempty, regexp(model.rxns,'^EX_.*_e$'))==0))
    pos_exchange = find(findExcRxns(model));
    exchange = model.rxns(pos_exchange);
    isNutrient = arrayfun(@(x) full(model.S(find(model.S(:,pos_exchange(x))), pos_exchange(x)))>0, 1:length(pos_exchange));
    nutrients = exchange(isNutrient);
    nutrients = nutrients(endsWith(nutrients,'_ex_'));
elseif ~isempty(find(cellfun(@isempty, regexp(model.rxns,'^EX_.*_e$'))==0))
    pos_exchange = findExcRxnsWithIDs(model);
    nutrients = model.rxns(pos_exchange);
end

model = setMediaFromRxns(model,mediaExchangeRxns,valueMediaExchangeRxns,nutrients);

end