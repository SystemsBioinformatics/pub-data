function identifyNewRelationshipForModel(model, database, modelName, language, newFileName)

global rootFolder
rootFolder = fileparts(which('initSystemsBioinformaticsToolbox'));
% load(fullfile(rootFolder, 'pipeline', 'reconstructions', 'manually_curated_models', 'BPE','iBP1870'))
% load(fullfile(rootFolder, 'BIGG', 'bigg_85.mat'))
% modelName = 'iBP1870';
% newFileName = [modelName '_additionalMapping'];
% language = 'sysbio';
% database = bigg;

posProtonsInModel = find(~cellfun(@isempty, regexp(database.mets, '^cpd00067_.*$|h_.*$|^PROTON_.*$|^C00080_.*$')));
databaseWOP = removeMetabolites(database, database.mets(posProtonsInModel),0);


if exist(['ids_curation_model_' modelName, '.xls'], 'file')==2
    iteration = xlsread(['ids_curation_model_' modelName], 'lastIteration');
    
    [~,s] = xlsread(['ids_curation_model_' modelName], 'accepted');
    if size(s,1)>1
        keys = s(2:end,1);
        values = s(2:end,2);
        
        for j = 1:length(keys)
            posMet1 = find(strcmp(model.mets,keys(j)));
            if ~isempty(posMet1)
                model = changeMetIdentifier(model, model.mets{posMet1}, values{j});
            end
            %                     model.mets{posMet1} = values{j};
            
            posMet1 = find(strcmp(model.mets,keys(j)));
            if ~isempty(posMet1)
                model = changeMetIdentifier(model, model.mets{posMet1}, values{j});
            end
        end
        %         [SysBio_keys_general, SysBio_values_general] = addMetabolitesToLocalDictionary(keys, values);
    end
else
    iteration = 0;
end

while 1
    iteration = iteration+1;
    if iteration ==1 && exist(['ids_curation_model_' modelName, '.xls'], 'file')~=2
        xlswrite(['ids_curation_model_' modelName], {'keys', 'values'} ,'rejected');
        xlswrite(['ids_curation_model_' modelName], {'keys', 'values'} ,'accepted');
        xlswrite(['ids_curation_model_' modelName],0,'lastIteration')
    end
    
    
    %filter out all the reactions with full coverage of metabolites
    allMetsInDatabase = zeros(size(model.rxns));
    for j = 1:length(model.rxns)
        allMetsInDatabase(j) = length(intersect(model.mets(model.S(:,j)~=0), database.mets)) == length(find(model.S(:,j)));
    end
    threeOrMoreMetabolites = zeros(size(model.rxns));
    for j = 1:length(model.rxns)
        threeOrMoreMetabolites(j) = length(model.mets(model.S(:,j)~=0))>=3;
    end
    allMetsInDatabaseButOne = zeros(size(model.rxns));
    for j = 1:length(model.rxns)
        allMetsInDatabaseButOne(j) = length(intersect(model.mets(model.S(:,j)~=0), database.mets)) == length(find(model.S(:,j)))-1;
    end
    
    posCandidateRxns = intersect(find(allMetsInDatabaseButOne), find(threeOrMoreMetabolites));
    if isempty(posCandidateRxns)
        break;
    end
    
    posI = arrayfun(@(x) findEquivalentRxnsWithMissingMetRelationship(model, x, database, databaseWOP), posCandidateRxns, 'UniformOutput', 0);

%     eqs = getRxn_cobraFormat(model, posCandidateRxns);     
%     posI = cellfun(@(x) findEquivalentRxnsWithMissingMetRelationship(model, x, database, databaseWOP), eqs, 'UniformOutput', 0);
    notEmpty = find(~cellfun(@isempty, posI));
    if isempty(notEmpty)
        break;
    end
    uniqueRxn = notEmpty;
    
    info = cell(length(uniqueRxn)+1,9);
    info(1,:) = {'met1','metName1','rxnID1','eq1','met2','metName2','rxnID2','eq2','acceptedByUser'};
    cont = 1;
    for j = 1:length(uniqueRxn)
        posRxn1 = posCandidateRxns(uniqueRxn(j));
        posRxn2 = posI{uniqueRxn(j)};
        eq1 = getRxn_cobraFormat(model, posRxn1);
        eq2 = getRxn_cobraFormat(database, posRxn2);
        rxnID1 = model.rxns(posRxn1);
        rxnID2 = database.rxns(posRxn2);
        for k = 1:length(posRxn2)
            mets1 = model.mets(model.S(:,posRxn1)~=0);
            mets2 = databaseWOP.mets(databaseWOP.S(:, posRxn2(k))~=0);
            common = intersect(mets1, mets2);
            met1 = setdiff(mets1, common);
            pos_met1 = find(strcmp(model.mets,met1));
            met1Name = model.metNames(pos_met1);
            met2 = setdiff(mets2, common);
            pos_met2 = find(strcmp(databaseWOP.mets,met2));
            met2Name = databaseWOP.metNames(pos_met2);
            info_j = [met1 met1Name rxnID1 eq1 met2 met2Name, rxnID2(k) eq2(k)];
            info(cont+1,1:8) = info_j;
            cont = cont + 1;
        end
    end
    
    %reduce list to unique mets
    info2 = info;
    for j = size(info2,1):-1:3
        id1 = info2(j,1);
        id2 = info2(j,5);
        
        is1 = ismember(id1, info2(2:j-1,1));
        is2 = ismember(id2, info2(2:j-1,5));
        if is1 && is2
            if ~isempty(intersect(find(strcmp(info2(2:j-1,1), id1)), find(strcmp(info2(2:j-1,5), id2))))
                info2(j,:) = [];
            end
        end
    end
    [~, ind] = sort(info2(2:end,1));
    info2 = [info2(1,:); info2(ind+1,:)];
    
    %sort the inf
    pairs = zeros(length(info(2:end,1)),1);
    for j = 2:length(info2(1:end,1))
        pos1 = find(strcmp(info(2:end,1), info2(j,1)));
        pos2 = find(strcmp(info(2:end,5), info2(j,5)));
        pairs(intersect(pos1,pos2)) = j-1;
    end
    [~, ind ] = sort(pairs);
    info_sorted = [info(1,:); info(ind+1,:)];
    
    %eliminate those already rejected
    [~,s]=xlsread(['ids_curation_model_' modelName], 'rejected');
    if size(s,1)>1
        keys_to_eliminate = s(2:end,1);
        values_to_eliminate = s(2:end,2);
        for j = 1:length(keys_to_eliminate)
            pos1 = find(strcmp(info2(2:end,1),keys_to_eliminate{j}));
            if ~isempty(pos1)
                [is, pos2] = ismember(values_to_eliminate{j},info2(pos1+1,5));
                if is
                    info2(pos1(pos2)+1,:) = [];
                end
            end
            
            pos1 = find(strcmp(info_sorted(2:end,1),keys_to_eliminate{j}));
            if ~isempty(pos1)
                [is, pos2] = ismember(values_to_eliminate{j},info_sorted(pos1+1,5));
                if is
                    info_sorted(pos1(pos2)+1,:) = [];
                end
            end
            
        end
    end
    
    if size(info2,1)==1
        break;
    end
    
    xlswrite(['ids_curation_model_' modelName], info_sorted, ['a_' num2str(iteration)])
    xlswrite(['ids_curation_model_' modelName], info2, num2str(iteration))
    
    pause()
    %we read the curation from the user
    [n,s] = xlsread(['ids_curation_model_' modelName], num2str(iteration));
    %accepted
    
    if isempty(n)
        xlswrite(excelFilePath,num2str(iteration),'lastIteration')
        break;
    end
    
    
    same = find(n(:,1));
    if ~isempty(same)
        mets1 = s(same+1,1);
        isInBigg = ismember(mets1,database.mets);
        mets1_dict = mets1;
        mets1_dict(isInBigg==1) = strcat('bigg:',  mets1_dict(isInBigg==1));
        mets1_dict(isInBigg==0) = strcat([language ':'],  mets1_dict(isInBigg==0));
        mets1_dict = removeCompartmentFromMets(mets1_dict);
        
        mets2 = s(same+1,5);
        mets2_dict = mets2;
        mets2_dict = strcat('bigg:',mets2_dict);
        mets2_dict = removeCompartmentFromMets(mets2_dict);
        for j = 1:length(mets1)
            posMet1 = find(strcmp(model.mets,mets1(j)));
            if ~isempty(posMet1)
                model = changeMetIdentifier(model, model.mets{posMet1}, mets2{j});
            end
            
            posMet1 = find(strcmp(model.mets,mets1(j)));
            if ~isempty(posMet1)
                model = changeMetIdentifier(model, model.mets{posMet1}, mets2{j});
            end
        end
        
%         [SysBio_keys_general, SysBio_values_general] = addMetabolitesToLocalDictionary(mets1, mets2);
        excelFilePath = ['ids_curation_model_' modelName];
        excelFileSheet = 'accepted';
        addMetabolitesToDictionaryInExcel(excelFilePath, excelFileSheet, mets1, mets2)
    end
    
    %rejected
    diff = find(n(:,1)==0);
    if ~isempty(diff)
        mets1d = s(diff+1,1);
        mets2d = s(diff+1,5);
        %$PENDING leer antes de grabar
        excelFilePath = ['ids_curation_model_' modelName];
        excelFileSheet = 'rejected';
        addMetabolitesToDictionaryInExcel(excelFilePath, excelFileSheet, mets1d, mets2d)
%         [SysBioRej_keys, SysBioRej_values] = addMetabolitesToLocalDictionaryRejections(mets1d, mets2d);
    end
    
    xlswrite(excelFilePath,num2str(iteration),'lastIteration')
    
end
save(newFileName,'model')
end