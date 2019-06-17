function identifyNewRelationships(models, idsType, metOtherIDs, metMNXIDs, bigg, species, modelsFileName)

global rootFolder

if exist(fullfile(rootFolder, 'MNX', 'SysBio_keys.mat'),'file')==2
    load(fullfile(rootFolder, 'MNX', 'SysBio_keys.mat'))
    load(fullfile(rootFolder, 'MNX', 'SysBio_values.mat'))
else
    SysBio_keys = {'keys'};
    SysBio_values = {'values'};
    save(fullfile(rootFolder, 'MNX', 'SysBio_keys.mat'),'SysBio_keys');
    save(fullfile(rootFolder, 'MNX', 'SysBio_values.mat'),'SysBio_values');
end

if exist(fullfile(rootFolder, 'MNX', 'SysBioRej_keys.mat'),'file')==2
    load(fullfile(rootFolder, 'MNX', 'SysBioRej_keys.mat'))
    load(fullfile(rootFolder, 'MNX', 'SysBioRej_values.mat'))
else
    SysBioRej_keys = {'keys'};
    SysBioRej_values = {'values'};
    save(fullfile(rootFolder, 'MNX', 'SysBioRej_keys.mat'),'SysBioRej_keys');
    save(fullfile(rootFolder, 'MNX', 'SysBioRej_values.mat'),'SysBioRej_values');
end

for i = 2:length(models)
    fprintf('analyzing models: progress %2.0f %%\n', 100*(i-1)/(length(models)-1));

    if ~strcmp(idsType(i), idsType{1})
        manualModel = models{1};
        model_i = models{i};
        posProtonsInModel = find(~cellfun(@isempty, regexp(model_i.mets, '^h_.*$')));
        model_i = removeMetabolites(model_i, model_i.mets(posProtonsInModel),0);
        
        posProtonsInModel = find(~cellfun(@isempty, regexp(manualModel.mets, '^h_.*$')));
        manualModelWOP = removeMetabolites(manualModel, manualModel.mets(posProtonsInModel),0);
        
        %test if the model can befenit from the dictionary
        dictionary = [SysBio_keys, SysBio_values];
        [dictionaryIsBeneficial, beneficial_pairs, n_rxns_improved, dictionaryIsBeneficial_perPair] = canBenefitFromDictionary(models{i}, models{1}, dictionary);
        
        if dictionaryIsBeneficial
            %apply changes
            keys = beneficial_pairs(:,1);
            values = beneficial_pairs(:,2);
            for j = 1:size(beneficial_pairs,1)
                
                posMet1 = find(strcmp(model_i.mets,keys(j)));
                if ~isempty(posMet1)
                    model_i = changeMetIdentifier(model_i, model_i.mets{posMet1}, values{j});
                end
                
                posMet1 = find(strcmp(models{i}.mets,keys(j)));
                if ~isempty(posMet1)
                    models{i} = changeMetIdentifier(models{i}, models{i}.mets{posMet1}, values{j});
                end
            end

            excelFilePath = ['ids_curation_model_' num2str(i)];
            excelFileSheet = 'accepted';
            if exist(['ids_curation_model_' num2str(i), '.xls'], 'file')~=2
                xlswrite(['ids_curation_model_' num2str(i)], {'keys', 'values'} ,'rejected');
                xlswrite(['ids_curation_model_' num2str(i)], {'keys', 'values'} ,'accepted');
                xlswrite(['ids_curation_model_' num2str(i)],0,'lastIteration')
            end
            addMetabolitesToDictionaryInExcel(excelFilePath, excelFileSheet, keys, values)
        end
        
        if exist(['ids_curation_model_' num2str(i), '.xls'], 'file')==2
            iteration = xlsread(['ids_curation_model_' num2str(i)], 'lastIteration');
            
            [~,s] = xlsread(['ids_curation_model_' num2str(i)], 'accepted');
            if size(s,1)>1
                keys = s(2:end,1);
                values = s(2:end,2);
                
                for j = 1:length(keys)
                    posMet1 = find(strcmp(model_i.mets,keys(j)));
                    if ~isempty(posMet1)
                        model_i = changeMetIdentifier(model_i, model_i.mets{posMet1}, values{j});
                    end
                    %                     model_i.mets{posMet1} = values{j};
                    
                    posMet1 = find(strcmp(models{i}.mets,keys(j)));
                    if ~isempty(posMet1)
                        models{i} = changeMetIdentifier(models{i}, models{i}.mets{posMet1}, values{j});
                    end
                end
                [SysBio_keys, SysBio_values] = addMetabolitesToLocalDictionary(keys, values);
            end
        else
            iteration = 0;
        end
        
        while 1
            iteration = iteration+1;
            if iteration ==1 && exist(['ids_curation_model_' num2str(i), '.xls'], 'file')~=2
                xlswrite(['ids_curation_model_' num2str(i)], {'keys', 'values'} ,'rejected');
                xlswrite(['ids_curation_model_' num2str(i)], {'keys', 'values'} ,'accepted');
                xlswrite(['ids_curation_model_' num2str(i)],0,'lastIteration')
            end
            
            %filter out all the reactions with full coverage of metabolites
            allMetsInDatabase = zeros(size(model_i.rxns));
            for j = 1:length(model_i.rxns)
                allMetsInDatabase(j) = length(intersect(model_i.mets(model_i.S(:,j)~=0), manualModel.mets)) == length(find(model_i.S(:,j)));
            end
            threeOrMoreMetabolites = zeros(size(model_i.rxns));
            for j = 1:length(model_i.rxns)
                threeOrMoreMetabolites(j) = length(model_i.mets(model_i.S(:,j)~=0))>=3;
            end
            allMetsInDatabaseButOne = zeros(size(model_i.rxns));
            for j = 1:length(model_i.rxns)
                allMetsInDatabaseButOne(j) = length(intersect(model_i.mets(model_i.S(:,j)~=0), manualModel.mets)) == length(find(model_i.S(:,j)))-1;
            end
            
            posCandidateRxns = intersect(find(allMetsInDatabaseButOne), find(threeOrMoreMetabolites));
            
            if isempty(posCandidateRxns)
                break;
            end
            
            posI = arrayfun(@(x) findEquivalentRxnsWithMissingMetRelationship(models{i}, x, manualModel, manualModelWOP), posCandidateRxns, 'UniformOutput', 0);
            
%             eqs = getRxn_cobraFormat(models{i}, posCandidateRxns);
%             posI = cellfun(@(x) findEquivalentRxnsWithMissingMetRelationship(models{i}, x, manualModel, manualModelWOP), eqs, 'UniformOutput', 0);
            
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
                eq1 = getRxn_cobraFormat(models{i}, posRxn1);
                eq2 = getRxn_cobraFormat(manualModel, posRxn2);
                rxnID1 = models{i}.rxns(posRxn1);
                rxnID2 = manualModel.rxns(posRxn2);
                for k = 1:length(posRxn2)
                    mets1 = model_i.mets(model_i.S(:,posRxn1)~=0);
                    mets2 = manualModelWOP.mets(manualModelWOP.S(:, posRxn2(k))~=0);
                    common = intersect(mets1, mets2);
                    met1 = setdiff(mets1, common);
                    pos_met1 = find(strcmp(model_i.mets,met1));
                    met1Name = model_i.metNames(pos_met1);
                    met2 = setdiff(mets2, common);
                    pos_met2 = find(strcmp(manualModelWOP.mets,met2));
                    met2Name = manualModelWOP.metNames(pos_met2);
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
            [~,s]=xlsread(['ids_curation_model_' num2str(i)], 'rejected');
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
            
            if length(SysBioRej_keys)>1
                keys_to_eliminate = SysBioRej_keys(2:end);
                values_to_eliminate = SysBioRej_values(2:end);
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
            
            xlswrite(['ids_curation_model_' num2str(i)], info_sorted, ['a_' num2str(iteration)])
            xlswrite(['ids_curation_model_' num2str(i)], info2, num2str(iteration))

            pause()
            %we read the curation from the user
            [n,s] = xlsread(['ids_curation_model_' num2str(i)], num2str(iteration));
            %accepted
            
            if isempty(n)
                xlswrite(excelFilePath,num2str(iteration),'lastIteration')
                break;
            end
            
            same = find(n(:,1));
            if ~isempty(same)
                mets1 = s(same+1,1);
                isInBigg = ismember(mets1,bigg.mets);
                mets1_dict = mets1;
                mets1_dict(isInBigg==1) = strcat('bigg:',  mets1_dict(isInBigg==1));
                mets1_dict(isInBigg==0) = strcat([idsType{i} ':'],  mets1_dict(isInBigg==0));
                mets1_dict = removeCompartmentFromMets(mets1_dict);
                
                mets2 = s(same+1,5);
                mets2_dict = mets2;
                mets2_dict = strcat('bigg:',mets2_dict);
                mets2_dict = removeCompartmentFromMets(mets2_dict);
                for j = 1:length(mets1)
                    posMet1 = find(strcmp(model_i.mets,mets1(j)));
                    if ~isempty(posMet1)
                        model_i = changeMetIdentifier(model_i, model_i.mets{posMet1}, mets2{j});
                    end
                    
                    posMet1 = find(strcmp(models{i}.mets,mets1(j)));
                    if ~isempty(posMet1)
                        models{i} = changeMetIdentifier(models{i}, models{i}.mets{posMet1}, mets2{j});
                    end
                end
                
                [SysBio_keys, SysBio_values] = addMetabolitesToLocalDictionary(mets1, mets2);
                excelFilePath = ['ids_curation_model_' num2str(i)];
                excelFileSheet = 'accepted';
                addMetabolitesToDictionaryInExcel(excelFilePath, excelFileSheet, mets1, mets2)
            end
            
            %rejected
            diff = find(n(:,1)==0);
            if ~isempty(diff)
                mets1d = s(diff+1,1);
                mets2d = s(diff+1,5);
                %$PENDING leer antes de grabar
                excelFilePath = ['ids_curation_model_' num2str(i)];
                excelFileSheet = 'rejected';
                addMetabolitesToDictionaryInExcel(excelFilePath, excelFileSheet, mets1d, mets2d)
                [SysBioRej_keys, SysBioRej_values] = addMetabolitesToLocalDictionaryRejections(mets1d, mets2d);
            end

            xlswrite(excelFilePath,num2str(iteration),'lastIteration')
            
        end
    end
    
end
models_additional_ids = models;
save([modelsFileName '_' species], 'models_additional_ids')

end