function createInputsForEMAF(model, biomassRxnID, baseDir, modelName, constraints_ids, constraints_lb, constraints_ub,posEX)

cd(baseDir);
fidg = fopen(['EMAF_' modelName '.bat'],'w+');

model = changeObjective(model, biomassRxnID);
model = transformModelToCOBRAFormat(model);
if nargin < 8 || isempty(posEX)
    posEX = findExcRxnsWithIDs(model);
end
EXRxns = model.rxns(posEX);
% model = changeRxnBounds(model, EXRxns, -1000, 'l');
% model = changeRxnBounds(model, EXRxns, 1000, 'u');
fba = optimizeCbModel(model);
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model,'sRxns',EXRxns);


if all(arrayfun(@(x) full(model.S(find(model.S(:,x)), x))<0, posEX))
    %it is assumed that reactions in the backward sense will be in the form '
    %-> 1 fru[e]'. For that all the exchange reactions in the original model
    %should be in the form "1 fru[e] -> "
    EXRxns_b = strcat(EXRxns,'_b');
    is = ismember(EXRxns_b, modelIrrev.rxns);
    posExchangeBackwards = getPosOfElementsInArray(EXRxns_b(is), modelIrrev.rxns);
    if ~isempty(find(is==0))
        posExchangeBackwards = [posExchangeBackwards, getPosOfElementsInArray(strcat(EXRxns(~is),'_r'), modelIrrev.rxns)];
    end
    posEXi = intersect(findExcRxnsWithIDs(modelIrrev), posExchangeBackwards);
    EXiRxns = modelIrrev.rxns(posEXi);
elseif all(arrayfun(@(x) full(model.S(find(model.S(:,x)), x))<0, posEX)==0)
    
    is = ismember(EXRxns, modelIrrev.rxns);
    if ~isempty(find(is==0))
        pos = find(is==0);
        for i =1:length(pos)
            EXRxns{pos(i)} = [EXRxns{pos(i)},'_r'];
        end
    end
    EXiRxns = EXRxns;

end
if exist([baseDir filesep 'emaf'], 'dir')~=7; mkdir([baseDir filesep 'emaf']); end
if exist([baseDir filesep 'emaf' filesep 'lpx'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'lpx']); end
if exist([baseDir filesep 'emaf' filesep 'models'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'models']); end
if exist([baseDir filesep 'emaf' filesep 'minimizationsets'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'minimizationsets']); end
if exist([baseDir filesep 'emaf' filesep 'constraints'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'constraints']); end

cd([baseDir filesep 'emaf' filesep 'models'])
newName = [modelName '_irrev.xml'];
if exist(newName,'file')~=2
    writeCbModel(modelIrrev, 'fileName', 'model', 'format', 'sbml');
    removeSquareBracketsFromBIGGModel('model.xml', newName);
    delete('model.xml');
end

cd([baseDir filesep 'emaf' filesep 'minimizationsets'])
if isempty(find(~cellfun(@isempty, regexp(EXiRxns, '^R_'))))
    EXiRxns = strcat('R_', EXiRxns);
end
data = ['Reaction ID'; EXiRxns];
writeCellArrayInCSVFormat(data, 'all.csv');

cd([baseDir filesep 'emaf' filesep 'constraints'])
if isempty(find(~cellfun(@isempty, regexp(constraints_ids, '^R_'))))
    constraints_ids = strcat('R_', constraints_ids);
end
info = [constraints_ids, num2cell(constraints_lb), num2cell(constraints_ub)];
for i = 1:size(info,1);  for j = 2:size(info,2); info{i,j} = num2str(info{i,j}); end; end;
data = [[{'Reaction ID'},{'LB'},{'UB'}]; info];
writeCellArrayInCSVFormat(data, 'basic.csv');

setModelNameForEMAF(newName, [baseDir filesep 'emaf']);

filePath = regexprep([baseDir filesep 'emaf'],'\\','\\\\');
fprintf(fidg, ['cd /d "' filePath '"\n']);
filePath = regexprep([baseDir filesep 'emaf' filesep 'runMedia3.py'],'\\','\\\\');
fprintf(fidg, ['python "' filePath '"\n']);
filePath = regexprep([baseDir filesep 'emaf' filesep 'pushRunMedia3.py'],'\\','\\\\');
fprintf(fidg, ['python "' filePath '"\n']);

fclose(fidg);

cd(baseDir)
end