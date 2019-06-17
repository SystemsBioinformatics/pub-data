function plotSimilarity(names,ratios,JDs,titleName)

unique_names = unique(names);
colors = {'b','r','g','c','m','y','k'};
f = figure;
for i = 1:length(unique_names)
    pos = getPosOfElementsInArray(unique_names(i),names);
    xs = JDs(pos);
    ys = ratios(pos);
    
    hold on
    plot(xs, ys,[colors{i} '.'],'MarkerSize',30)
end
ylimit = max(ratios)*1.2;
ylim([0 ylimit]);
xlim([0 1]);
legend(unique_names,'Location', 'SouthOutside','Orientation','horizontal')
ylabel('Ratio')
xlabel('Jaccard Distance')

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [titleName, '.pdf']);

close(gcf);

end