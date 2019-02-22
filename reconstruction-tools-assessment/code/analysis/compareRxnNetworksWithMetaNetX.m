function [rxnMatrix, rxnIDs] = compareRxnNetworksWithMetaNetX(models, idsType, names,...
    fromBeginning, species, rxnOtherIDs, rxnMNXIDs, rxnMatrixFileName, rxnIDsFileName, pdfFileName)
% This function compares a set of networks, in terms of metabolites and
% reactions. The comparison is done by string comparison. It is assummed
% that the reference network is written in BIGG language. In the case of
% the networks which are not written in BIGG, MetaNetX is used to map
% metabolites and reactiones
%
% USAGE:
%
%    [metMatrix, rxnMatrix, metIDs, rxnIDs, tableMetComparison, tableRxnComparison] = ...
%     compareNetworks2(models, idsType, names, fromBeginning, species, metOtherIDs, metMNXIDs...
%     , rxnOtherIDs, rxnMNXIDs, posExcludedRxns, posExcludedMets, ...
%     metMatrixFileName, metIDsFileName, rxnMatrixFileName, rxnIDsFileName
%
% INPUT:
%    models:         cell array of models. 
%
% OUTPUTS:
%    metMatrix:      matrix with 1
%
% EXAMPLE:
%
%    formula = '0.01 cdpdag-SC[m] + 0.01 pg-SC[m]  -> 0.01 clpn-SC[m] + cmp[m] + h[m]'
%
%    [metaboliteList, stoichCoeffList, revFlag] = parseRxnFormula(formula)
%
%    %metaboliteList = 'cdpdag-SC[m]'    'pg-SC[m]'    'clpn-SC[m]'    'cmp[m]'    'h[m]'
%    %stoichCoeffList = -0.01 -0.01 0.01 1 1
%    %revFlag = false
%
% .. Authors:
%       - Sebastián Mendoza 17/01/2019

% initialize rxn table
rxnIDs = cell(10000,length(models)+1);
for i = 1:length(models)+1
    emptys =  find(cellfun(@isempty,rxnIDs(:,i))==1);
    for j = 1:length(emptys); rxnIDs{emptys(j),i} = ''; end;
end

if exist(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnMatrixFileName '_' species '.mat'], 'file')==2 && ...
        exist(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnIDsFileName '_' species '.mat'], 'file')==2 ...
        && ~fromBeginning
    load(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnMatrixFileName '_' species '.mat']);
    load(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnIDsFileName '_' species '.mat'])
elseif fromBeginning
    if exist(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnMatrixFileName '_' species '.mat'], 'file')==2 
        delete(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnMatrixFileName '_' species '.mat']);        
    end
    if exist(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnIDsFileName '_' species '.mat'], 'file')==2 
        delete(['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species filesep rxnIDsFileName '_' species '.mat']);
    end
    
    % fill first column of rxn comparison table
    [inter, ~,  posInMNX]= intersect(strcat([idsType{1} ':'],models{1}.rxns), rxnOtherIDs);
    MNXIDs_1 = rxnMNXIDs(posInMNX);
    pos = cell2mat(arrayfun(@(x)find(strcmp(x,strcat([idsType{1} ':'],models{1}.rxns))),inter,'UniformOutput',false))';
    rxnIDs(1:length(models{1}.rxns), 2) = models{1}.rxns;
    rxnIDs(pos, 1) = MNXIDs_1;
    
    % construct the rest of rxn comparison table
    fprintf('comparing rxns in models: progress %2.0f %%\n', 0);
    for i = 2:length(models)
        if strcmp('bigg', idsType{i})
            rxns_i = models{i}.rxns;
            diff = rxns_i;
            for j =2:i
                rxns_j = rxnIDs(:,j);
                [int, ~, ind2] = intersect(diff, rxns_j);
                if ~isempty(int)
                    diff = setdiff(diff, int);
                    rxnIDs(ind2, i+1) = rxnIDs(ind2, j);
                end
                if j ==i
                    lastPosition = 0;
                    for k = 2:i
                        lastPosition_k = find(cellfun(@isempty,rxnIDs(:,k))==0);
                        lastPosition_k = lastPosition_k(end);
                        if lastPosition_k>lastPosition
                            lastPosition = lastPosition_k;
                        end
                    end
                    if ~isempty(diff)
                        rxnIDs(lastPosition+1:lastPosition+length(diff),i+1) = diff;
                        rxnDiff = strcat([idsType{i} ':'],diff);
                        translation = cell(size(diff));
                        for k = 1:length(diff)                 
                            pos = find(strcmp(rxnOtherIDs, rxnDiff{k}));
                            if ~isempty(pos)
                                %has translation
                                translation_rxn_k = rxnMNXIDs(pos);
                                translation(k) = translation_rxn_k;
                            else
                                translation{k} = '';
                            end
                        end
                        rxnIDs(lastPosition+1:lastPosition+length(diff),1) = translation;

                    end
                end
            end
            
        else
            
            pile_withoutTranslation = cell(2000,1);
            pile_withTranslation = cell(2000,1);
            translation = cell(2000,1);
            cont_1 = 0; 
            cont_2 = 0;
            
            for j = 1:length(models{i}.rxns)
                
                if strcmp(models{i}.rxns{j},'CARBODEHYDRAT-RXN')
                    disp('')
                end
                
                pos = find(strcmp(rxnOtherIDs, [idsType{i} ':' models{i}.rxns{j}]));
                
                if ~isempty(pos)
                    %has translation
                    if length(pos) ==1
                        translation_j = rxnMNXIDs(pos);
                        pos2 = find(strcmp(rxnIDs(:,1),translation_j ));
                        if ~isempty(pos2)
                            %is already in table
                            if length(pos2) ==1
                                if isempty(rxnIDs{pos2,i+1})
                                    rxnIDs{pos2,i+1} = models{i}.rxns{j};
                                else
                                    %if there is already one reaction in
                                    %that space, then j is a duplicated
                                    %reaction in the model i. Then, it is
                                    %added to the same spot
                                    rxnIDs{pos2,i+1} = strcat(rxnIDs{pos2,i+1}, ';', models{i}.rxns{j});
                                end
                            else
                                cont_1 = cont_1 + 1;
                                pile_withoutTranslation{cont_1} = models{i}.rxns{j};
                            end
                        else
                            %is not in table
                            cont_2 = cont_2 + 1;
                            pile_withTranslation{cont_2} = models{i}.rxns{j};
                            translation(cont_2) = rxnMNXIDs(pos);
                            
                        end                       
                    else
                        cont_1 = cont_1 + 1;
                        pile_withoutTranslation{cont_1} = models{i}.rxns{j};
                    end
                    
                else
                    %hasn't translation
                    cont_1 = cont_1 + 1;
                    pile_withoutTranslation{cont_1} = models{i}.rxns{j}; 
                end              
            end
            
            ind = find(cellfun(@isempty,rxnIDs(:,i))==0);
            final = ind(end);
            
            translation = translation(1:cont_2);
            pile_withTranslation = pile_withTranslation(1:cont_2);
            pile_withoutTranslation = pile_withoutTranslation(1:cont_1);
            
            rxnIDs(final+1:final+length(translation),1) = translation;
            rxnIDs(final+1:final+length(pile_withTranslation),i+1) = pile_withTranslation;
            
            ind2 = find(cellfun(@isempty,rxnIDs(:,1))==0);
            final2 = ind2(end);
            rxnIDs(final2+1:final2+length(pile_withoutTranslation),i+1) = pile_withoutTranslation;        
        end  
        fprintf('comparing rxns in models: progress %2.0f %%\n', 100*(i-1)/(length(models)-1));
    end
    
    final = 0;
    for i = 1:length(models)+1
        ind = find(cellfun(@isempty,rxnIDs(:,i))==0);
        if ind(end)>final
           final = ind(end); 
        end
    end
    rxnIDs = rxnIDs(1:final,:);
    
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
    
    save([rxnMatrixFileName '_' species], 'rxnMatrix')
    save([rxnIDsFileName '_' species], 'rxnIDs')
end


mapa = [0.95 1 1
    0.95 1 1
    51/255 153/255 255/255];

if size(names,2)>size(names,1)
    names = names';
end

ColumnLabelsValue = [ names;'MetaNetX'];
h = HeatMap(rxnMatrix,'ColumnLabels', ColumnLabelsValue','Colormap',mapa);
addYLabel(h, 'Reactions', 'FontSize', 12)
addXLabel(h, 'Tools', 'FontSize', 12)
addTitle(h, 'Comparison of reactions', 'FontSize', 14)
plot(h);
set(gcf,'renderer','opengl','renderermode','manual')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [pdfFileName '_' species, '.pdf']);
close(gcf);
close(gcf);
close all hidden

end