function classifyNotCovered(models, fastRxnMatrix5,shortAbbr)

n_reconstructions = length(models)-1;
ylabels = 'percentage of reactions';
titleName = 'distribution of not recovered reactions';

p = [];
p1 = [];
p2 = [];
h = [];
h3 = [];
h4 = [];
h5 = [];
h6 = [];
posNotCytosol = getPosOfElementsInArray(findTransRxns(transformModelToCOBRAFormat(models{1})), models{1}.rxns);
posNotCytosol = union(posNotCytosol, find(findExcRxns(models{1})));
for i = 2:length(models); 
    %not found in draft
    pos = find(fastRxnMatrix5(1:length(models{1}.rxns),i)==0);  
    %associated genes
    iGenes = getInvolvedGenes(models{1},pos); 
    pEmpty = 100*length(find(cellfun(@isempty, iGenes)))/length(iGenes);
    %how many, of the reactions that were not recovered in the manual model, are not gene-associated
    p = [p; pEmpty];
    p1 = [p1; 100*length(intersect(posNotCytosol,pos(find(cellfun(@isempty, iGenes)))))/length(iGenes)];
    p2 =  [p2; p(i-1)-p1(i-1)];
    %of the ones that are gene-associated
    notEmpty = find(cellfun(@isempty, iGenes)==0);
    unGenes ={};
    for j = 1:length(notEmpty)
        notEmpty_j = iGenes{notEmpty(j)};
        unGenes = union(unGenes, notEmpty_j);
    end
    howMuchInDraft = 100*length(intersect(unGenes, models{i}.genes))/length(unGenes);
    %how much of those genes are in the draft model
    h = [h; howMuchInDraft];
    
    %how much percentage of genes of each reaction (with GPRs) missing are present in the draft
    h2 = arrayfun(@(y) 100*length(find(cellfun(@(x) ismember(x,models{i}.genes), iGenes{notEmpty(y)})))/length(iGenes{notEmpty(y)}), 1:length(notEmpty));
    %how many have at least one gene missing
    h3 = [h3; 100*length(intersect(find(h2<100), find(h2>=50)))/length(iGenes)];
    h4 = [h4; 100*length(find(h2<50))/length(iGenes)];
    
    h5 = [h5; 100*length(intersect(posNotCytosol,pos(notEmpty(find(h2==100)))))/length(iGenes)];
    h6 = [h6; 100-p(i-1)-h3(i-1)-h4(i-1)-h5(i-1)];
end;

numbers = [p1';p2';h4';h3';h5';h6'];

c = shortAbbr(2:end);
f = figure;
h = bar(numbers','stacked');
limits = [0 n_reconstructions+1];
ylim([0 100]);
xlim(limits);
% set(h,{'FaceColor'},{[0 0.3 0.6]; [0 0.5 0.6];[0 0.7 0.6 ];[0 0.9 0.6 ];[0.9 0.9 0 ]});
ylabel(ylabels, 'FontSize', 10);
ax=gca;
ax.XTick = 1:1:length(c);
set(gca, 'FontSize', 10)
% set(gca, 'XTickLabel', c, 'FontSize', 8)
set(gca,'Xtick',[]);
legend({'[p] or [e], no GPR', '[c], no GPR', '>= 50% of genes missing', '< 50% of genes missing', '[p] or [e], gene associated', 'others'},'Location','southoutside','Orientation', 'horizontal', 'FontSize', 10)

ylimit = 100;
altura = ylimit*80/3500;
% htxt = text(xb(1,:)-0.3,repmat(-80,1,length(xb)), c, 'FontSize', FontSize+1); %for 3500 reactions
% htxt = text(xb(1,:)-0.3,repmat(-36,1,length(xb)), c, 'FontSize',FontSize+1); %for 1600 reactions
htxt = text([1:n_reconstructions]-0.3,repmat(-altura,1,length(1:[n_reconstructions])), c, 'FontSize',8);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [titleName, '.pdf']);
close(gcf);


p = [];
p1 = [];
p2 = [];
h = [];
h3 = [];
h4 = [];
h5 = [];
h6 = [];
posNotCytosol = getPosOfElementsInArray(findTransRxns(transformModelToCOBRAFormat(models{1})), models{1}.rxns);
posNotCytosol = union(posNotCytosol, find(findExcRxns(models{1})));
for i = 2:length(models); 
    %not found in draft
    pos = find(fastRxnMatrix5(1:length(models{1}.rxns),i)==0);  
    %associated genes
    iGenes = getInvolvedGenes(models{1},pos); 
    pEmpty = length(find(cellfun(@isempty, iGenes)));
    %how many, of the reactions that were not recovered in the manual model, are not gene-associated
    p = [p; pEmpty];
    p1 = [p1; length(intersect(posNotCytosol,pos(find(cellfun(@isempty, iGenes)))))];
    p2 =  [p2; p(i-1)-p1(i-1)];
    %of the ones that are gene-associated
    notEmpty = find(cellfun(@isempty, iGenes)==0);
    unGenes ={};
    for j = 1:length(notEmpty)
        notEmpty_j = iGenes{notEmpty(j)};
        unGenes = union(unGenes, notEmpty_j);
    end
    howMuchInDraft = 100*length(intersect(unGenes, models{i}.genes))/length(unGenes);
    %how much of those genes are in the draft model
    h = [h; howMuchInDraft];
    
    %how much percentage of genes of each reaction (with GPRs) missing are present in the draft
    h2 = arrayfun(@(y) 100*length(find(cellfun(@(x) ismember(x,models{i}.genes), iGenes{notEmpty(y)})))/length(iGenes{notEmpty(y)}), 1:length(notEmpty));
    %how many have at least one gene missing
    h3 = [h3; length(intersect(find(h2<100), find(h2>=50)))];
    h4 = [h4; length(find(h2<50))];
    
    h5 = [h5; length(intersect(posNotCytosol,pos(notEmpty(find(h2==100)))))];
    h6 = [h6; length(iGenes)-p(i-1)-h3(i-1)-h4(i-1)-h5(i-1)];
end;

numbers = [p1';p2';h4';h3';h5';h6'];

n1 = max(sum(numbers,1)) - rem(max(sum(numbers,1)),200) + 200;
n2 = max(sum(numbers,1)) - rem(max(sum(numbers,1)),500) + 500;
if n1 < n2 
    ylimit = n1;
else
    ylimit = n2;
end

c = shortAbbr(2:end);
f = figure;
h = bar(numbers','stacked');
limits = [0 n_reconstructions+1];
ylim([0 ylimit]);
xlim(limits);
% set(h,{'FaceColor'},{[0 0.3 0.6]; [0 0.5 0.6];[0 0.7 0.6 ];[0 0.9 0.6 ];[0.9 0.9 0 ]});
ylabel(ylabels, 'FontSize', 10);
ax=gca;
ax.XTick = 1:1:length(c);
set(gca, 'FontSize', 10)
% set(gca, 'XTickLabel', c, 'FontSize', 8)
set(gca,'Xtick',[]);
legend({'[p] or [e], no GPR', '[c], no GPR', '>= 50% of genes missing', '< 50% of genes missing', '[p] or [e], gene associated', 'others'},'Location','southoutside','Orientation', 'horizontal', 'FontSize', 10)

altura = ylimit*80/3500;
% htxt = text(xb(1,:)-0.3,repmat(-80,1,length(xb)), c, 'FontSize', FontSize+1); %for 3500 reactions
% htxt = text(xb(1,:)-0.3,repmat(-36,1,length(xb)), c, 'FontSize',FontSize+1); %for 1600 reactions
htxt = text([1:n_reconstructions]-0.3,repmat(-altura,1,length(1:[n_reconstructions])), c, 'FontSize',8);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [titleName, '_b.pdf']);
close(gcf);

end