function [metMatrix, metIDs, rxnMatrix, rxnIDs, rxnMatrix2, rxnIDs2, Eqs2, rxnMatrix3, rxnIDs3] = ...
    compareNetworks3(models, idsType, names, fromBeginning, species, metOtherIDs, metMNXIDs, ...
    rxnOtherIDs, rxnMNXIDs, posExcludedMets, posExcludedRxns, metMatrixFileName, metIDsFileName, ...
    rxnMatrixFileName, rxnIDsFileName, rxnMatrixFileName2, rxnIDsFileName2, ...
    eqsFileName, metExcelFileName, rxnExcelFileName)

global rootFolder

if size(names,2)>size(names,1)
    names = names';
end

% initialize rxn table
rxnIDs = cell(20000,length(models)+1);
rxnIDs2 = cell(20000,length(models)+1);
Eqs2 = cell(20000,length(models)+1);
metIDs = cell(20000,length(models));

for i = 1:length(models)+1
    emptys =  find(cellfun(@isempty,rxnIDs(:,i))==1);
    for j = 1:length(emptys); rxnIDs{emptys(j),i} = ''; end;
end

for i = 1:length(models)+1
    emptys =  find(cellfun(@isempty,rxnIDs2(:,i))==1);
    for j = 1:length(emptys); rxnIDs2{emptys(j),i} = ''; end;
end

for i = 1:length(models)+1
    emptys =  find(cellfun(@isempty,Eqs2(:,i))==1);
    for j = 1:length(emptys); Eqs2{emptys(j),i} = ''; end;
end

for i = 1:length(models)
    emptys =  find(cellfun(@isempty,metIDs(:,i))==1);
    for j = 1:length(emptys); metIDs{emptys(j),i} = ''; end;
end

%% metabolite comparison

if fromBeginning
    if exist(fullfile(rootFolder, 'results', species, [metIDsFileName '_' species '.mat']), 'file')==2;
        delete(fullfile(rootFolder, 'results', species, [metIDsFileName '_' species '.mat']))
    end
    if exist(fullfile(rootFolder, 'results', species, [metMatrixFileName '_' species '.mat']), 'file')==2;
        delete(fullfile(rootFolder, 'results', species, [metMatrixFileName '_' species '.mat']))
    end
    
    % fill first column of rxn comparison table
    metIDs(1:length(models{1}.mets), 1) = models{1}.mets;
    
    fprintf('comparing metabolites in models: progress %2.0f %%\n', 0);
    for i = 2:length(models)
        met_i = models{i}.mets;
        diff = met_i;
        for j =1:i-1
            mets_j = metIDs(:,j);
            [int, ~, ind2] = intersect(diff, mets_j);
            if ~isempty(int)
                diff = setdiff(diff, int);
                metIDs(ind2, i) = metIDs(ind2, j);
            end
            if j ==i-1
                lastPosition = 0;
                for k = 1:i-1
                    lastPosition_k = find(cellfun(@isempty,metIDs(:,k))==0);
                    lastPosition_k = lastPosition_k(end);
                    if lastPosition_k>lastPosition
                        lastPosition = lastPosition_k;
                    end
                end
                if ~isempty(diff)
                    metIDs(lastPosition+1:lastPosition+length(diff),i) = diff;
                end
            end
        end
        fprintf('comparing metabolites in models: progress %2.0f %%\n', 100*(i-1)/length(models));
    end
    
    
    maximo = 0;
    for i = 1:length(models)
        indGene = find(cellfun(@isempty,metIDs(:,i))==0);
        finalGene = indGene(end);
        
        if finalGene>maximo
            maximo = finalGene;
        end
    end
    
    metIDs = metIDs(1:maximo,:);
    metMatrix = metIDs;
    
    for i = 1:length(models)
        pos_empty = find(cellfun(@isempty, metIDs(:,i)));
        pos_full =  find(cellfun(@isempty, metIDs(:,i)) ==0);
        for j = 1:length(pos_empty)
            metMatrix{pos_empty(j),i} = 0;
        end
        for j = 1:length(pos_full)
            metMatrix{pos_full(j),i} = 1;
        end
    end
    
    metMatrix = cell2mat(metMatrix);
    save([metMatrixFileName '_' species], 'metMatrix')
    save([metIDsFileName '_' species], 'metIDs')
    
else
    load([metMatrixFileName '_' species])
    load([metIDsFileName '_' species])
end

mapa = [0.95 1 1
    0.95 1 1
    51/255 153/255 255/255];

ColumnLabelsValue = names;
h = HeatMap(metMatrix,'ColumnLabels', ColumnLabelsValue','Colormap',mapa);
addYLabel(h, 'Metabolites', 'FontSize', 12)
addXLabel(h, 'Tools', 'FontSize', 12)
addTitle(h, 'Comparison of metabolites', 'FontSize', 14)
f = plot(h);
set(gcf,'renderer','opengl','renderermode','manual')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [metExcelFileName '_' species, '.pdf']);
close(gcf);
close(gcf);
close all hidden

%% reaction comparison
tic
if exist(fullfile(rootFolder, 'results', species, [rxnMatrixFileName '_' species '.mat']), 'file')==2 && ...
        exist(fullfile(rootFolder, 'results', species, [rxnIDsFileName '_' species '.mat']), 'file')==2 &&...
        exist(fullfile(rootFolder, 'results', species, [rxnMatrixFileName2 '_' species '.mat']), 'file')==2 &&...
        exist(fullfile(rootFolder, 'results', species, [rxnIDsFileName2 '_' species '.mat']), 'file')==2 &&...
        exist(fullfile(rootFolder, 'results', species, [eqsFileName '_' species '.mat']), 'file')==2 &&...
        ~fromBeginning
    load(fullfile(rootFolder, 'results', species, [rxnMatrixFileName '_' species '.mat']));
    load(fullfile(rootFolder, 'results', species, [rxnIDsFileName '_' species '.mat']))
    load(fullfile(rootFolder, 'results', species, [rxnMatrixFileName2 '_' species '.mat']));
    load(fullfile(rootFolder, 'results', species, [rxnIDsFileName2 '_' species '.mat']))
    load(fullfile(rootFolder, 'results', species, [eqsFileName '_' species '.mat']))
    
elseif fromBeginning
    if exist(fullfile(rootFolder, 'results', species, [rxnMatrixFileName '_' species '.mat']), 'file')==2
        delete(fullfile(rootFolder, 'results', species, [rxnMatrixFileName '_' species '.mat']));
    end
    if exist(fullfile(rootFolder, 'results', species, [rxnIDsFileName '_' species '.mat']), 'file')==2
        delete(fullfile(rootFolder, 'results', species, [rxnIDsFileName '_' species '.mat']));
    end
    if exist(fullfile(rootFolder, 'results', species, [rxnMatrixFileName2 '_' species '.mat']), 'file')==2
        delete(fullfile(rootFolder, 'results', species, [rxnMatrixFileName2 '_' species '.mat']));
    end
    if exist(fullfile(rootFolder, 'results', species, [rxnIDsFileName2 '_' species '.mat']), 'file')==2
        delete(fullfile(rootFolder, 'results', species, [rxnIDsFileName2 '_' species '.mat']));
    end
    
    % fill first column of rxn comparison table
    [inter, ~,  posInMNX]= intersect(strcat([idsType{1} ':'],models{1}.rxns), rxnOtherIDs);
    MNXIDs_1 = rxnMNXIDs(posInMNX);
    pos = cell2mat(arrayfun(@(x)find(strcmp(x,strcat([idsType{1} ':'],models{1}.rxns))),inter,'UniformOutput',false))';
    rxnIDs(1:length(models{1}.rxns), 2) = models{1}.rxns;
    rxnIDs(pos, 1) = MNXIDs_1;
    
    rxnIDs2(1:length(models{1}.rxns), 2) = models{1}.rxns;
    rxnIDs2(pos, 1) = MNXIDs_1;
    
    for i = 1:length(models{1}.rxns)
        Eqs2(i,2) = getRxn_cobraFormat(models{1}, i);
    end
    
    % construct the rest of rxn comparison table
    fprintf('comparing rxns in models: progress %2.0f %%\n', 0);
    for i = 2:length(models)
        model_i = models{i};
        rxns_i = model_i.rxns;
        diff = rxns_i;
        rxns_original = rxnIDs(:,2);
        [rxnFormulaWasFound, posRxn, rxnOppositeDirection, matchWithoutRemovingProtons, ...
            additionalMatchsIfRemovingProtons, posRxnsFullMatch, posRxnsPartialMatch, validPartialMatchs] = ...
            cellfun(@(x) reactionFormulaInModel(models{1}, x, 1, 1), getRxn_cobraFormat(model_i, rxns_i),'UniformOutput',0);
        posExactMatch = intersect(find(cell2mat(rxnFormulaWasFound)), find(cell2mat(matchWithoutRemovingProtons)));
        int = rxns_i(posExactMatch);
        ind2 = cell2mat(posRxnsFullMatch(posExactMatch));
        if ~isempty(int)
            diff = setdiff(diff, int);
            rxnIDs(ind2, i+1) = int;
        end
        posMatchWithProtonFlexibility = find(~cellfun(@isempty, posRxnsPartialMatch));
        int2 = rxns_i(posMatchWithProtonFlexibility);
        
        
        if ~isempty(int2)
            valids = cell(1000,1);
            contValids = 0;
            for j = 1:length(int2)
                ind2_j = posRxnsPartialMatch{posMatchWithProtonFlexibility(j)};
                valid = validPartialMatchs{posMatchWithProtonFlexibility(j)};
                valid = valid(arrayfun(@(x)find(posRxn{posMatchWithProtonFlexibility(j)} ==x),ind2_j));
                for k = 1:length(ind2_j)
                    if valid(k)
                        Eqs2(ind2_j(k), i+1) = getRxn_cobraFormat(model_i, int2(j));
                        rxnIDs2(ind2_j(k), i+1) = int2(j);
                        contValids = contValids + 1;
                        valids{contValids} = int2{j};
                    end
                end
            end
            valids = valids(1:contValids);
            diff = setdiff(diff, valids);
        end
        
        lastPosition = find(cellfun(@isempty,rxnIDs(:,2))==0);
        lastPosition = lastPosition(end);
        if ~isempty(diff)
            rxnIDs(lastPosition+1:lastPosition+length(diff),i+1) = diff;
        end
        fprintf('comparing rxns in models: progress %2.0f %%\n', 100*(i-1)/(length(models)-1));
        
        %         setdiff(models{i}.rxns, rxnIDs(:,i+1))
    end
    
    final = 0;
    for i = 1:length(models)+1
        ind = find(cellfun(@isempty,rxnIDs(:,i))==0);
        if ~isempty(ind)
            if ind(end)>final
                final = ind(end);
            end
        end
    end
    rxnIDs = rxnIDs(1:final,:);
    
    final = 0;
    for i = 1:length(models)+1
        ind = find(cellfun(@isempty,rxnIDs2(:,i))==0);
        if ~isempty(ind)
            if ind(end)>final
                final = ind(end);
            end
        end
    end
    rxnIDs2 = rxnIDs2(1:final,:);
    
    rxnMatrix = rxnIDs;
    for i = 1:length(models)+1
        pos_empty = find(cellfun(@isempty, rxnIDs(:,i)));
        pos_full =  find(cellfun(@isempty, rxnIDs(:,i)) ==0);
        for j = 1:length(pos_empty)
            rxnMatrix{pos_empty(j),i} = 0;
        end
        for j = 1:length(pos_full)
            rxnMatrix{pos_full(j),i} = 1;
        end
    end
    rxnMatrix = cell2mat(rxnMatrix);
    aux = rxnMatrix(:,1);
    rxnMatrix = [rxnMatrix(:,2:end), aux];
    
    %verify that matrices are consistent
    
    for i = 1:size(rxnMatrix,1)
        for j = 1:size(rxnMatrix,2)-1
            if rxnMatrix(i,j)==1 && isempty(rxnIDs{i,j+1})
                fprintf('error tipo 1: i: %1.0f j: %1.0f\n', i, j)
            elseif rxnMatrix(i,j)==0 && ~isempty(rxnIDs{i,j+1})
                fprintf('error tipo 2: i: %1.0f j: %1.0f\n', i, j)
            end
        end
    end
    
    rxnMatrix2 = rxnIDs2;
    for i = 1:length(models)+1
        pos_empty = find(cellfun(@isempty, rxnIDs2(:,i)));
        pos_full =  find(cellfun(@isempty, rxnIDs2(:,i)) ==0);
        for j = 1:length(pos_empty)
            rxnMatrix2{pos_empty(j),i} = 0;
        end
        for j = 1:length(pos_full)
            rxnMatrix2{pos_full(j),i} = 1;
        end
    end
    rxnMatrix2 = cell2mat(rxnMatrix2);
    aux = rxnMatrix2(:,1);
    rxnMatrix2 = [rxnMatrix2(:,2:end), aux];
    
    %verify that matrices are consistent
    
    for i = 1:size(rxnMatrix2,1)
        for j = 1:size(rxnMatrix2,2)-1
            if rxnMatrix2(i,j)==1 && isempty(rxnIDs2{i,j+1})
                fprintf('error tipo 1: i: %1.0f j: %1.0f\n', i, j)
            elseif rxnMatrix2(i,j)==0 && ~isempty(rxnIDs2{i,j+1})
                fprintf('error tipo 2: i: %1.0f j: %1.0f\n', i, j)
            end
        end
    end
    
    %     n_rxns = length(models{1}.rxns);
    %     rxnMatrix3 = rxnMatrix;
    %     rxnIDs3 = rxnIDs;
    %     for i = 1:n_rxns
    %         for j = 3:length(models)+1
    %             if rxnMatrix2(i,j-1)==1
    %                 rxnMatrix3(i,j-1)=1;
    %                 rxnIDs3(i,j)=rxnIDs2(i,j);
    %             end
    %         end
    %     end
    %
    %     n_rxns = length(models{1}.rxns);
    %     rxnMatrix3 = rxnMatrix;
    %     rxnIDs3 = rxnIDs;
    %     for j = 3:length(models)+1
    %         last = find(~cellfun(@isempty, rxnIDs3(:,j)));
    %         last = last(end);
    %         auxPile = rxnIDs3(n_rxns+1:last,j);
    %         newLast = last;
    %         for i = 1:length(models{1}.rxns)
    %             if rxnMatrix2(i,j-1)==1 && rxnMatrix(i,j-1)==0
    %                 disp(j)
    %                 disp(i)
    %                 disp('here')
    %                 rxnIDs3(i,j)=rxnIDs2(i,j);
    %                 pos = find(strcmp(auxPile, rxnIDs2{i,j}));
    %                 disp(pos)
    %                 auxPile(pos) = [];
    %             end
    %         end
    %     end
    
    save([rxnMatrixFileName '_' species], 'rxnMatrix')
    save([rxnIDsFileName '_' species], 'rxnIDs')
    save([rxnMatrixFileName2 '_' species], 'rxnMatrix2')
    save([rxnIDsFileName2 '_' species], 'rxnIDs2')
    save([eqsFileName '_' species], 'Eqs2')
    %     save([rxnMatrixFileName3 '_' species], 'rxnMatrix3')
    %     save([rxnIDsFileName3 '_' species], 'rxnIDs3')
end
toc

ColumnLabelsValue = [ names;'MetaNetX'];

mapa = [0.95 1 1
    0.95 1 1
    51/255 153/255 255/255];

h = HeatMap(rxnMatrix,'ColumnLabels', ColumnLabelsValue','Colormap',mapa);
addYLabel(h, 'Reactions', 'FontSize', 12)
addXLabel(h, 'Tools', 'FontSize', 12)
addTitle(h, 'Comparison of reactions', 'FontSize', 14)
f = plot(h);
set(gcf,'renderer','opengl','renderermode','manual')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [rxnExcelFileName '_' species, '.pdf']);
close(gcf);
close(gcf);
close all hidden

end