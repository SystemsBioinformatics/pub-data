function getReconstructionComparisonPlotFromMatrix(matrix, matrix2, models, xlabels, names, idsType, ylabels,...
    titleName, vector, excludedElements, correctByCoverage, species)

if nargin < 7
    vector = [11, 120, 90, 10, 5, 5, 10];
end
if nargin < 9
    correctByCoverage = 0;
end

percentageCutoff = vector(1);
addToStringOver10 = vector(2);
addToStringBelow10 = vector(3);
padval1= vector(4);
padval2= vector(5);
padval3= vector(6);
addLine= vector(7);

n_reconstructions = length(models)-1;
n_elements = length(find(matrix(:,1)));
n_excluded = length(excludedElements);
if n_excluded==0 
    n_grupos = 3;
else
    n_grupos = 4;
end
if ~isempty(matrix2)
    n_grupos = n_grupos + 1;
end
y = zeros(size(matrix,2)-1, n_grupos);
percentages = zeros(n_grupos,size(matrix,2)-1);
numbers = zeros(n_grupos,size(matrix,2)-1);

for i = 2:size(matrix,2)
    n_coverage = length(find(matrix(1:n_elements,i)));
    if ~isempty(matrix2)
        n_partialCoverage = length(find(matrix2(1:n_elements,i)));
    else 
        n_partialCoverage = 0;
    end
    n_notCovered = (n_elements-n_excluded) - n_coverage - n_partialCoverage;
    n_additional = length(find(matrix(n_elements+1:end,i)));
    if n_excluded ==0
        if ~isempty(matrix2)
            y_i = [n_coverage, n_partialCoverage, n_notCovered, n_additional];
        else
            y_i = [n_coverage, n_notCovered, n_additional];
        end
    else
        if ~isempty(matrix2)
            y_i = [n_coverage, n_partialCoverage, n_notCovered ,n_excluded, n_additional];
        else
            y_i = [n_coverage, n_notCovered ,n_excluded, n_additional];
        end
    end
    y(i-1,:) = y_i;
    if n_excluded ==0
        percentages(:,i-1) = [100*n_coverage/(n_elements-n_excluded); 100*n_notCovered/(n_elements-n_excluded); 100*n_additional/(n_elements-n_excluded)];
        numbers(:,i-1) =  [n_coverage; n_notCovered; n_additional];
        if ~isempty(matrix2)
            percentages(:,i-1) = [100*n_coverage/(n_elements-n_excluded); 100*n_partialCoverage/(n_elements-n_excluded); 100*n_notCovered/(n_elements-n_excluded); 100*n_additional/(n_elements-n_excluded)];
            numbers(:,i-1) =  [n_coverage; n_partialCoverage; n_notCovered; n_additional];
        else
            percentages(:,i-1) = [100*n_coverage/(n_elements-n_excluded); 100*n_notCovered/(n_elements-n_excluded); 100*n_additional/(n_elements-n_excluded)];
            numbers(:,i-1) =  [n_coverage; n_notCovered; n_additional];
        end
    else 
        if ~isempty(matrix2)
            percentages(:,i-1) = [100*n_coverage/(n_elements-n_excluded); 100*n_partialCoverage/(n_elements-n_excluded); 100*n_notCovered/(n_elements-n_excluded); 100*n_excluded/n_elements; 100*n_additional/(n_elements-n_excluded)];
            numbers(:,i-1) = [n_coverage; n_partialCoverage; n_notCovered; n_excluded; n_additional];
        else
            percentages(:,i-1) = [100*n_coverage/(n_elements-n_excluded); 100*n_notCovered/(n_elements-n_excluded); 100*n_excluded/n_elements; 100*n_additional/(n_elements-n_excluded)];
            numbers(:,i-1) = [n_coverage; n_notCovered; n_excluded; n_additional];
        end
    end
end
n1 = max(sum(numbers,1)) - rem(max(sum(numbers,1)),200) + 200;
n2 = max(sum(numbers,1)) - rem(max(sum(numbers,1)),500) + 500;
if n1 < n2 
    ylimit = n1;
else
    ylimit = n2;
end
if correctByCoverage
    if ~isempty(strfind(ylabels,'metabolites'))
        [~,s] = xlsread(['MNX_coverage_' species '.xls'],'Mets');
        [~,s2] = xlsread(['metabolites_comparison_' species '.xls'],'metabolites');
        for i = 2:length(models)
            in_model_i = setdiff(unique(s2(1:n_elements,i+1)),'');
            not_in_model_i = setdiff(unique(s2(n_elements+1:end,i+1)),'');
            pos_end = find(~cellfun(@isempty, s(1:end,2*i)));
            pos_end = pos_end(end);
            if pos_end == 1
                not_in_MNX = {};
            else
                not_in_MNX = s(2:pos_end,2*i);
            end
            pos_end = find(~cellfun(@isempty, s(1:end,2*i-1)));
            pos_end = pos_end(end);
            if pos_end == 1
                in_MNX = {};
            else
                in_MNX = s(2:pos_end,2*i-1);
            end
            not_in_model_not_in_MNX = intersect(not_in_MNX, not_in_model_i);
            not_in_model_in_MNX = intersect(in_MNX, not_in_model_i);
            in_model_in_MNX = intersect(in_MNX, in_model_i);
            percentage_in_model_of_elements_in_MNX = length(in_model_in_MNX)/length(in_MNX);
            to_add = ceil(length(not_in_MNX)*percentage_in_model_of_elements_in_MNX);
            numbers(1,i-1) = numbers(1,i-1) + to_add;
            numbers(2,i-1) = (n_elements-n_excluded) - numbers(1,i-1); 
            numbers(end,i-1) = numbers(end,i-1) - to_add;
            if n_excluded ==0
                percentages(:,i-1) = (100*numbers(:,i-1))./[(n_elements-n_excluded), (n_elements-n_excluded), (n_elements-n_excluded)]';
            else
                percentages(:,i-1) = (100*numbers(:,i-1))./[(n_elements-n_excluded), (n_elements-n_excluded), (n_elements), (n_elements-n_excluded)]';
            end
        end
    else
    end
end

c = xlabels(2:end);
f = figure;
h = bar(numbers','stacked');
limits = [0 n_reconstructions+1];
ylim([0 ylimit]);
xlim(limits);
if n_grupos ==5
    set(h,{'FaceColor'},{[0 0.3 0.6]; [0 0.5 0.6];[0 0.7 0.6 ];[0 0.9 0.6 ];[0.9 0.9 0 ]});
elseif n_grupos ==4
    set(h,{'FaceColor'},{[0 0.3 0.6];[0 0.45 0.6 ];[0 0.6 0.6 ];[0.9 0.9 0 ]});
elseif n_grupos ==3
    set(h,{'FaceColor'},{[0 0.3 0.6];[0 0.45 0.6 ];[0.9 0.9 0 ]});
end
ylabel(ylabels, 'FontSize', 10);
ax=gca;
ax.XTick = 1:1:length(c);
set(gca, 'FontSize', 10)
% set(gca, 'XTickLabel', c, 'FontSize', 8)
set(gca,'Xtick',[]);
% set(gca, 'YTickLabel', c, 'FontSize', 14)
% title(titleName, 'FontSize', 14);
yb = cat(1, h.YData);
for j = n_grupos:-1:2
    for i = 1:size(yb,2)
        yb(j,i) = sum(yb(1:j,i));
    end
end

xb = bsxfun(@plus, h(1).XData, [h.XOffset]');
if n_grupos==4
    yb = [yb(1:2,:);yb(4,:)];
    xb = [xb(1:2,:);xb(4,:)];
    percentages = [percentages(1:2,:);percentages(4,:)];
elseif n_grupos==5
    yb = [yb(1:3,:);yb(5,:)];
    xb = [xb(1:3,:);xb(5,:)];
    percentages = [percentages(1:3,:);percentages(5,:)];
end

percentages_filtered = zeros(size(percentages,1),size(percentages,2));
xb_filtered = zeros(size(percentages,1),size(percentages,2));
yb_filtered = zeros(size(percentages,1),size(percentages,2));
whiteBars = zeros(size(percentages,1),size(percentages,2));
% greens = zeros(size(percentages,1),size(percentages,2));
cont = 0;

xlines = [];
ylines = [];
whiteLines = [];
for j = 1:size(percentages,2)
    for i = 1:size(percentages,1)
        cont = cont+1;
        if percentages(i,j) > percentageCutoff
            yb_filtered(cont) = yb(i,j);
        elseif percentages(i,j) < 9.5
            xlines = [xlines; j j];
            if i==1
                ylines = [ylines; yb(i,j)/2, yb(i,j)+addLine];
                whiteLines = [whiteLines 1];
            elseif i==2 && n_grupos ==5
                ylines = [ylines; (yb(i,j)+yb(i-1,j))/2, yb(i,j)+addLine];
                whiteLines = [whiteLines 1];
            else
                ylines = [ylines; (yb(i,j)+yb(i-1,j))/2, yb(i,j)+addLine];
                whiteLines = [whiteLines 0];
            end
            yb_filtered(cont) = yb(i,j)+addToStringBelow10;
        else
            xlines = [xlines; j j];
            if i==1
                ylines = [ylines; yb(i,j)/2, yb(i,j)+addLine];
                whiteLines = [whiteLines 1];
            elseif i==2 && n_grupos ==5
                ylines = [ylines; (yb(i,j)+yb(i-1,j))/2, yb(i,j)+addLine];
                whiteLines = [whiteLines 1];
            else
                ylines = [ylines; (yb(i,j)+yb(i-1,j))/2, yb(i,j)+addLine];
                whiteLines = [whiteLines 0];
            end
            yb_filtered(cont) = yb(i,j)+addToStringOver10;
        end
        if i ==1
            yb_filtered(cont) = yb_filtered(cont) - padval1;
        elseif i ==2
            yb_filtered(cont) = yb_filtered(cont) - padval2;
        else
            yb_filtered(cont) = yb_filtered(cont) - padval3;
        end
        percentages_filtered(cont) = percentages(i,j);
        xb_filtered(cont) = xb(i,j);
        if i==1
            whiteBars(cont) = 1;
        elseif i ==2 && percentages(i,j) > percentageCutoff
            whiteBars(cont) = 1;
        elseif n_grupos ==5 && i ==2
            whiteBars(cont) = 1;
        elseif n_grupos ==5 && i ==3
            whiteBars(cont) = 1;
        end
        
    end
end
percentages_filtered = percentages_filtered(1:cont);
xb_filtered = xb_filtered(1:cont);
yb_filtered = yb_filtered(1:cont);
whiteBars = whiteBars(1:cont);

hold on;
for i = 1:size(xlines,1)
    if ~whiteLines(i)
        line(xlines(i,:), ylines(i,:), 'Color','black','LineWidth',0.6);
    else
        line(xlines(i,:), ylines(i,:), 'Color','white','LineWidth',0.6);
    end
end

FontSize =9; %before 9
if n_grupos == 4
    FontSize = 8;  
end    
% htxt = text(xb_filtered(:),yb_filtered(:), strcat(cellstr(num2str(percentages_filtered(:),'%3.1f\n')), ' %'), 'rotation', 90, 'horiz', 'right','FontSize', FontSize);
htxt = text(xb_filtered(:),yb_filtered(:), cellstr(num2str(percentages_filtered(:),'%3.1f\n')), 'rotation', 90, 'horiz', 'right','FontSize', FontSize);
set(htxt(find(whiteBars)), 'color', 'w'); % for legibility
% set(htxt(find(greens)), 'color', 'w'); % for legibility

if n_grupos ==3
    legend({'In model', 'Not in model', 'Additional'},'Location','southoutside', 'Orientation', 'horizontal', 'FontSize', 10)
elseif n_grupos ==4
    legend({'In model', 'Not in model', 'Organism-specific', 'Additional'},'Location','southoutside','Orientation', 'horizontal', 'FontSize', 10)
elseif  n_grupos ==5
    legend({'In model, full match', 'In model, partial match', 'Not in model', 'Organism-specific', 'Additional'},'Location','southoutside','Orientation', 'horizontal', 'FontSize', 10)
end

altura = ylimit*80/3500;
% htxt = text(xb(1,:)-0.3,repmat(-80,1,length(xb)), c, 'FontSize', FontSize+1); %for 3500 reactions
% htxt = text(xb(1,:)-0.3,repmat(-36,1,length(xb)), c, 'FontSize',FontSize+1); %for 1600 reactions
htxt = text(xb(1,:)-0.3,repmat(-altura,1,length(xb)), c, 'FontSize',8);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, [titleName, '.pdf']);
% saveas(gcf, [titleName, '.tif']);
% 
% 
% set(gcf,'color','w');
% imagewd = getframe(gcf);
% imwrite(imagewd.cdata, [titleName, '2.tif'], 'Resolution', 600)
close(gcf);

%determine groups
groups = {};
groups_tools = {};
curr_group = 1;
groups{1} = 2;
groups_tools{1} = names{groups{curr_group}(1)};

for i = 3:length(names)
    if ~isempty(find(strcmp(groups_tools, names{i}))) 
        pos = find(strcmp(groups_tools, names{i}));
        groups{pos} = union(groups{pos}, i);
    else
        curr_group = curr_group + 1;
        groups{curr_group} = i;
        groups_tools{curr_group} = names{i};
    end
end
%calculate variation in coverage, additional reactions, per tool 
coverage_av_per_tool = cell(size(groups));
additional_av_per_tool = cell(size(groups));
coverage_std_per_tool= cell(size(groups));
additional_std_per_tool = cell(size(groups));
for i = 1:length(groups)
    coverage_av_per_tool{i}= mean(percentages(1,groups{i}-1));
    additional_av_per_tool{i} = mean(percentages(3,groups{i}-1));
    coverage_std_per_tool{i}= std(percentages(1,groups{i}-1));
    additional_std_per_tool{i} = std(percentages(3,groups{i}-1));
end


info1 = [{'tools', 'coverage_av_per_tool','additional_av_per_tool','coverage_std_per_tool','additional_std_per_tool',};
    groups_tools' coverage_av_per_tool' additional_av_per_tool' coverage_std_per_tool' additional_std_per_tool'];

%determine groups
groups = {};
groups_tools = {};
curr_group = 1;
groups{1} = 2;
groups_tools{1} = names{groups{curr_group}(1)};
groups_languages{1} = idsType{groups{curr_group}(1)};

for i = 3:length(names)
    if ~isempty(find(strcmp(groups_tools, names{i}))) && ~isempty(find(strcmp(groups_languages, idsType{i}))) && ~isempty(intersect(find(strcmp(groups_tools, names{i})),find(strcmp(groups_languages, idsType{i})))) 
        pos = intersect(find(strcmp(groups_tools, names{i})),find(strcmp(groups_languages, idsType{i})));
        groups{pos} = union(groups{pos}, i);
    else
        curr_group = curr_group + 1;
        groups{curr_group} = i;
        groups_tools{curr_group} = names{i};
        groups_languages{curr_group} = idsType{i};
    end
end
%calculate variation in coverage, additional reactions, per tool 
coverage_av_per_tool_per_language = cell(size(groups));
additional_av_per_tool_per_language = cell(size(groups));
coverage_std_per_tool_per_language= cell(size(groups));
additional_std_per_tool_per_language = cell(size(groups));
for i = 1:length(groups)
    coverage_av_per_tool_per_language{i}= mean(percentages(1,groups{i}-1));
    additional_av_per_tool_per_language{i} = mean(percentages(3,groups{i}-1));
    coverage_std_per_tool_per_language{i}= std(percentages(1,groups{i}-1));
    additional_std_per_tool_per_language{i} = std(percentages(3,groups{i}-1));
end

info2 = [{'tools', 'language', 'coverage_av_per_tool_per_language','additional_av_per_tool_per_language','coverage_std_per_tool_per_language','additional_std_per_tool_per_language',};
    groups_tools', groups_languages', coverage_av_per_tool_per_language', additional_av_per_tool_per_language' coverage_std_per_tool_per_language' additional_std_per_tool_per_language'];
xlswrite([titleName '_statistics_per_tool_' species], info1, 'per tool');
xlswrite([titleName '_statistics_per_tool_' species], info2, 'per language');

end