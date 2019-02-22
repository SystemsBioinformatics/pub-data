function [geneMatrix, geneIDs] = geneComparison(baseDirectory, models, names, ...
    fromBeginning, species, generatePlot, plotFileName)

% initialize gene table

if fromBeginning
    if exist([baseDirectory filesep species filesep 'geneIDs_' species], 'file')==2;
        delete([baseDirectory filesep species filesep 'geneIDs_' species])
    end
    if exist([baseDirectory filesep species filesep 'geneMatrix_' species], 'file')==2;
        delete([baseDirectory filesep species filesep 'geneMatrix_' species])
    end
    geneIDs = cell(5000,length(models));
    for i = 1:length(models)
        emptys =  find(cellfun(@isempty,geneIDs(:,i))==1);
        for j = 1:length(emptys); geneIDs{emptys(j),i} = ''; end;
    end
    
    % fill first column of rxn comparison table
    geneIDs(1:length(models{1}.genes), 1) = models{1}.genes;
    
    fprintf('comparing genes in models: progress %2.0f %%\n', 0);
    for i = 2:length(models)
        genes_i = models{i}.genes;
        diff = genes_i;
        for j =1:i-1
            genes_j = geneIDs(:,j);
            [int, ~, ind2] = intersect(diff, genes_j);
            if ~isempty(int)
                diff = setdiff(diff, int);
                geneIDs(ind2, i) = geneIDs(ind2, j);
            end
            if j ==i-1
                lastPosition = find(cellfun(@isempty,geneIDs(:,j))==0);
                lastPosition = lastPosition(end);
                if ~isempty(diff)
                    geneIDs(lastPosition+1:lastPosition+length(diff),i) = diff;
                end
            end
        end
        fprintf('comparing genes in models: progress %2.0f %%\n', 100*(i-1)/length(models));
    end
    
    maximum = 0;
    for i = 1:length(models)
        indGene = find(cellfun(@isempty,geneIDs(:,i))==0);
        finalGene = indGene(end);
        
        if finalGene>maximum
            maximum = finalGene;
        end
    end
    
    geneIDs = geneIDs(1:maximum,:);
    geneMatrix = geneIDs;
    
    for i = 1:length(models)
        pos_empty = find(cellfun(@isempty, geneIDs(:,i)));
        pos_full =  find(cellfun(@isempty, geneIDs(:,i)) ==0);
        for j = 1:length(pos_empty)
            geneMatrix{pos_empty(j),i} = 0;
        end
        for j = 1:length(pos_full)
            geneMatrix{pos_full(j),i} = 1;
        end
    end
    
    geneMatrix = cell2mat(geneMatrix);
    save([baseDirectory filesep species filesep 'geneMatrix_' species], 'geneMatrix')
    save([baseDirectory filesep species filesep 'geneIDs_' species], 'geneIDs')
    
else
    load([baseDirectory filesep species filesep 'geneMatrix_' species] )
    load([baseDirectory filesep species filesep 'geneIDs_' species])
end

if generatePlot
    ColumnLabelsValue = names;
    mapa = [0.95 1 1
        0.95 1 1
        51/255 153/255 255/255];
    
    h = HeatMap(geneMatrix,'ColumnLabels', ColumnLabelsValue,'Colormap',mapa);
    addYLabel(h, 'Genes', 'FontSize', 12)
    addXLabel(h, 'Tools', 'FontSize', 12)
    addTitle(h, 'Comparison of genes', 'FontSize', 14)
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