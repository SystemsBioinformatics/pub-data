function [rxnMatrix, rxnIDs1, rxnEqs1, rxnMatrix2, rxnIDs2, rxnEqs2, rxnMatrix3, rxnIDs3, rxnEqs3] = ...
    compareRxnNetworksWithEquations(baseDirectory, models, idsType, names, fromBeginning, species, ...
    rxnOtherIDs, rxnMNXIDs, rxnMatrixFileName1, rxnIDsFileName1, eqsFileName1, rxnMatrixFileName2, rxnIDsFileName2, ...
    eqsFileName2, rxnMatrixFileName3, rxnIDsFileName3,  eqsFileName3, generatePlot, plotFileName, plotFileName2)

% This function compares a set of networks, in terms of reactions.
% Reactions are compared though equations. 
%
% USAGE:
%
%     [rxnMatrix, rxnIDs1, rxnMatrix2, rxnIDs2, rxnEqs2, rxnMatrix3, rxnIDs3] = ...
%     compareRxnNetworksWithEquations(baseDirectory, models, idsType, names, fromBeginning, species, ...
%     rxnOtherIDs, rxnMNXIDs, rxnMatrixFileName1, rxnIDsFileName1, rxnMatrixFileName2, rxnIDsFileName2, ...
%     rxnMatrixFileName3, rxnIDsFileName3, eqsFileName2, plotFileName, plotFileName2)
%
% INPUT:
%    models:         cell array of models. 
%
% OUTPUTS:
%    rxnMatrix:      matrix with 1
%
% EXAMPLE:
%
% .. Authors:
%       - Sebastián Mendoza 17/01/2019


if exist([baseDirectory filesep species filesep rxnMatrixFileName1 '_' species '.mat'], 'file')==2 && ...
        exist([baseDirectory filesep species filesep rxnIDsFileName1 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep rxnMatrixFileName2 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep rxnIDsFileName2 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep rxnMatrixFileName3 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep rxnIDsFileName3 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep eqsFileName1 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep eqsFileName2 '_' species '.mat'], 'file')==2 &&...
        exist([baseDirectory filesep species filesep eqsFileName3 '_' species '.mat'], 'file')==2 &&...
        ~fromBeginning
    
    load([baseDirectory filesep species filesep rxnMatrixFileName1 '_' species '.mat']);
    load([baseDirectory filesep species filesep rxnIDsFileName1 '_' species '.mat'])
    load([baseDirectory filesep species filesep rxnMatrixFileName2 '_' species '.mat']);
    load([baseDirectory filesep species filesep rxnIDsFileName2 '_' species '.mat'])
    load([baseDirectory filesep species filesep rxnMatrixFileName3 '_' species '.mat']);
    load([baseDirectory filesep species filesep rxnIDsFileName3 '_' species '.mat'])
    load([baseDirectory filesep species filesep eqsFileName1 '_' species '.mat'])
    load([baseDirectory filesep species filesep eqsFileName2 '_' species '.mat'])
    load([baseDirectory filesep species filesep eqsFileName3 '_' species '.mat'])
    
elseif fromBeginning
    if exist([baseDirectory filesep species filesep rxnMatrixFileName1 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep rxnMatrixFileName1 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep rxnIDsFileName1 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep rxnIDsFileName1 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep rxnMatrixFileName2 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep rxnMatrixFileName2 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep rxnIDsFileName2 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep rxnIDsFileName2 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep rxnMatrixFileName3 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep rxnMatrixFileName3 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep rxnIDsFileName3 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep rxnIDsFileName3 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep eqsFileName1 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep eqsFileName1 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep eqsFileName2 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep eqsFileName2 '_' species '.mat']);
    end
    if exist([baseDirectory filesep species filesep eqsFileName3 '_' species '.mat'], 'file')==2
        delete([baseDirectory filesep species filesep eqsFileName3 '_' species '.mat']);
    end
    
    
    % initialize rxn table
    rxnIDs1 = cell(5000,length(models)+1);
    rxnIDs2 = cell(5000,length(models)+1);
    rxnEqs1 = cell(5000,length(models)+1);
    rxnEqs2 = cell(5000,length(models)+1);

    for i = 1:length(models)+1
        emptys =  find(cellfun(@isempty,rxnIDs1(:,i))==1);
        for j = 1:length(emptys); rxnIDs1{emptys(j),i} = ''; end;
    end
    
    for i = 1:length(models)+1
        emptys =  find(cellfun(@isempty,rxnIDs2(:,i))==1);
        for j = 1:length(emptys); rxnIDs2{emptys(j),i} = ''; end;
    end
    
    for i = 1:length(models)+1
        emptys =  find(cellfun(@isempty,rxnEqs2(:,i))==1);
        for j = 1:length(emptys); rxnEqs2{emptys(j),i} = ''; end;
    end
    
    % fill first column of rxn comparison table
    [inter, ~,  posInMNX]= intersect(strcat([idsType{1} ':'],models{1}.rxns), rxnOtherIDs);
    MNXIDs_1 = rxnMNXIDs(posInMNX);
    pos = cell2mat(arrayfun(@(x)find(strcmp(x,strcat([idsType{1} ':'],models{1}.rxns))),inter,'UniformOutput',false))';
    rxnIDs1(1:length(models{1}.rxns), 2) = models{1}.rxns;
    rxnIDs1(pos, 1) = MNXIDs_1;
    
    rxnIDs2(1:length(models{1}.rxns), 2) = models{1}.rxns;
    rxnIDs2(pos, 1) = MNXIDs_1;
    
    for i = 1:length(models{1}.rxns)
        rxnEqs1(i,2) = getRxn_cobraFormat(models{1}, i);
        rxnEqs2(i,2) = rxnEqs1(i,2);
    end
    
    % construct the rest of rxn comparison table
    fprintf('comparing rxns in models: progress %2.0f %%\n', 0);
    for i = 2:length(models)
        model_i = models{i};
        rxns_i = model_i.rxns;
        diff = rxns_i;
        [rxnFormulaWasFound, posRxn, ~, matchWithoutRemovingProtons, ...
            ~, posRxnsFullMatch, posRxnsPartialMatch, validPartialMatchs] = ...
            cellfun(@(x) reactionFormulaInModel(models{1}, x, 1, 1), getRxn_cobraFormat(model_i, rxns_i),'UniformOutput',0);
        posExactMatch = intersect(find(cell2mat(rxnFormulaWasFound)), find(cell2mat(matchWithoutRemovingProtons)));
        int = rxns_i(posExactMatch);
        ind2 = cell2mat(posRxnsFullMatch(posExactMatch));
        if ~isempty(int)
            diff = setdiff(diff, int);
            rxnIDs1(ind2, i+1) = int;
            for j = 1:length(int)
                rxnEqs1(ind2, i+1) = getRxn_cobraFormat(model_i, int(j));
            end
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
                        rxnEqs2(ind2_j(k), i+1) = getRxn_cobraFormat(model_i, int2(j));
                        rxnIDs2(ind2_j(k), i+1) = int2(j);
                        contValids = contValids + 1;
                        valids{contValids} = int2{j};
                    end
                end
            end
            valids = valids(1:contValids);
            diff = setdiff(diff, valids);
        end
        
        lastPosition = find(cellfun(@isempty,rxnIDs1(:,2))==0);
        lastPosition = lastPosition(end);
        if ~isempty(diff)
            rxnIDs1(lastPosition+1:lastPosition+length(diff),i+1) = diff;
        end
        fprintf('comparing rxns in models: progress %2.0f %%\n', 100*(i-1)/(length(models)-1));
    end
    
    final = 0;
    for i = 1:length(models)+1
        ind = find(cellfun(@isempty,rxnIDs1(:,i))==0);
        if ~isempty(ind)
            if ind(end)>final
                final = ind(end);
            end
        end
    end
    rxnIDs1 = rxnIDs1(1:final,:);
    
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
    
    rxnMatrix = rxnIDs1;
    for i = 1:length(models)+1
        pos_empty = find(cellfun(@isempty, rxnIDs1(:,i)));
        pos_full =  find(cellfun(@isempty, rxnIDs1(:,i)) ==0);
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
            if rxnMatrix(i,j)==1 && isempty(rxnIDs1{i,j+1})
                fprintf('error tipo 1: i: %1.0f j: %1.0f\n', i, j)
            elseif rxnMatrix(i,j)==0 && ~isempty(rxnIDs1{i,j+1})
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
                fprintf('error 1: i: %1.0f j: %1.0f\n', i, j)
            elseif rxnMatrix2(i,j)==0 && ~isempty(rxnIDs2{i,j+1})
                fprintf('error 2: i: %1.0f j: %1.0f\n', i, j)
            end
        end
    end
    
    n_rxns = length(models{1}.rxns);
    rxnMatrix3 = rxnMatrix;
    rxnIDs3 = rxnIDs1;
    rxnEqs3 = rxnEqs1;
    for i = 1:n_rxns
        for j = 3:length(models)+1
            if rxnMatrix2(i,j-1)==1
                rxnMatrix3(i,j-1)=1;
                rxnIDs3(i,j)=rxnIDs2(i,j);
                rxnEqs3(i,j)=rxnEqs2(i,j);
            end
        end
    end
    
    save([rxnMatrixFileName1 '_' species], 'rxnMatrix')
    save([rxnIDsFileName1 '_' species], 'rxnIDs1')
    save([eqsFileName1 '_' species], 'rxnEqs1')
    save([rxnMatrixFileName2 '_' species], 'rxnMatrix2')
    save([rxnIDsFileName2 '_' species], 'rxnIDs2')
    save([eqsFileName2 '_' species], 'rxnEqs2')
    save([rxnMatrixFileName3 '_' species], 'rxnMatrix3')
    save([rxnIDsFileName3 '_' species], 'rxnIDs3')
    save([eqsFileName3 '_' species], 'rxnEqs3')
end

if generatePlot
    if size(names,2)>size(names,1); names = names'; end
    ColumnLabelsValue = [ names;'MetaNetX'];
    
    map = [0.95 1 1
        0.95 1 1
        51/255 153/255 255/255];
    
    h = HeatMap(rxnMatrix,'ColumnLabels', ColumnLabelsValue','Colormap',map);
    addYLabel(h, 'Reactions', 'FontSize', 12)
    addXLabel(h, 'Tools', 'FontSize', 12)
    addTitle(h, 'Comparison of reactions (only full match)', 'FontSize', 14)
    f = plot(h);
    set(gcf,'renderer','opengl','renderermode','manual')
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 33 15]);
    set(gcf,'PaperOrientation','landscape')
    set(gcf,'PaperPosition', [-1 1 30 19])
    saveas(gcf, [plotFileName '_' species, '.pdf']);
    close(gcf);
    close(gcf);
    close all hidden
    
    h = HeatMap(rxnMatrix3,'ColumnLabels', ColumnLabelsValue','Colormap',map);
    addYLabel(h, 'Reactions', 'FontSize', 12)
    addXLabel(h, 'Tools', 'FontSize', 12)
    addTitle(h, 'Comparison of reactions (full and partial match)', 'FontSize', 14)
    f = plot(h);
    set(gcf,'renderer','opengl','renderermode','manual')
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 33 15]);
    set(gcf,'PaperOrientation','landscape')
    set(gcf,'PaperPosition', [-1 1 30 19])
    saveas(gcf, [plotFileName2 '_' species, '.pdf']);
    close(gcf);
    close(gcf);
    close all hidden
end

end