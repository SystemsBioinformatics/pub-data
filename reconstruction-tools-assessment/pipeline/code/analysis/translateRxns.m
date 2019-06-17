function [modelOut, dict, n, p, requiereManualCheck] = translateRxns(model, relaxProtons, referenceDatabase)

if allMetsInCOBRAFormat(model.mets)
    model = transformModelToCBMPYFormat(model);
end

requiereManualCheck = {};
modelOut = model; 
n = 0;
dict = cell(length(model.rxns),1);

eqs = getRxn_cobraFormat(modelOut, 1:length(modelOut.rxns));
for i = 1:length(modelOut.rxns)

    [reactionFormulaWasFound, posRxn] = reactionFormulaInModel(referenceDatabase, eqs{i} , 1, relaxProtons);
    if reactionFormulaWasFound
        if length(posRxn) ==1
            if ~strcmp(modelOut.rxns{i},referenceDatabase.rxns{posRxn})
                if ~ismember(referenceDatabase.rxns{posRxn},modelOut.rxns)
                    n = n+1;
                    dict{n,1} = modelOut.rxns{i};
                    dict{n,2} = referenceDatabase.rxns{posRxn};
                    modelOut.rxns{i} = referenceDatabase.rxns{posRxn};
                    if length(unique(modelOut.rxns))<length(unique(model.rxns))
                        disp('')
                    end
                else
                    a = modelOut.rxns(i);
                    b = referenceDatabase.rxns{posRxn};
                    info = [a, b];
                    requiereManualCheck = [requiereManualCheck; info];
                end
            end
        else
            if ~ismember(modelOut.rxns{i},referenceDatabase.rxns(posRxn))
                a = modelOut.rxns(i);
                b = {strjoin(referenceDatabase.rxns(posRxn),',')};
                info = [a, b];
                requiereManualCheck = [requiereManualCheck; info];
            end
        end
    end
end
p = n/length(model.rxns);
dict = dict(1:n,:);
end