function [metMatrix, metIDs] = compareMetNetworksWithStrings(baseDirectory, models, names,...
    fromBeginning, species, metMatrixFileName, metIDsFileName, generatePlot, plotFileName)
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
% .. Authors:
%       - Sebastián Mendoza 17/01/2019

if fromBeginning
    if exist([baseDirectory filesep species filesep metIDsFileName '_' species], 'file')==2;
        delete([baseDirectory filesep species filesep metIDsFileName '_' species])
    end
    if exist([baseDirectory filesep species filesep metMatrixFileName '_' species], 'file')==2;
        delete([baseDirectory filesep species filesep metMatrixFileName '_' species])
    end
    
    metIDs = cell(10000,length(models));
    
    for i = 1:length(models)
        emptys =  find(cellfun(@isempty,metIDs(:,i))==1);
        for j = 1:length(emptys); metIDs{emptys(j),i} = ''; end;
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

if generatePlot
    mapa = [0.95 1 1
        0.95 1 1
        51/255 153/255 255/255];
    
    if size(names,2)>size(names,1)
        names = names';
    end
    
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
    saveas(gcf, [plotFileName '_' species, '.pdf']);
    close(gcf);
    close(gcf);
    close all hidden
end

end