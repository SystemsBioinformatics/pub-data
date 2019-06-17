function [pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency)

global rootFolder

[~,s] = xlsread(fullfile(rootFolder, 'pipeline', 'code', 'analysis', 'KEGG_pathways.xlsx'));
general = cell(12,1);
specific = cell(12,1);
codes_specific = cell(12,1);

cont = 0;
list = {};
codes = {};
for i = 1:length(s)
    if isempty(strfind(s{i},'0'))
        if ~isempty(list)
            specific{cont} = list;
            codes_specific{cont} = codes;
        end
        cont = cont+1;
        general{cont} = s{i};
        list = {};
        codes = {};
    else
        pos = strfind(s{i}, '  ');
        code = s{i}(1:pos-1);
        pathway = s{i}(pos+2:end);
        list = [list; pathway];
        codes = [codes; code];
        if i ==length(s)
            specific{cont} = list;
            codes_specific{cont} = codes;
        end
    end
end

pathways2 = pathways;
frequency2 = frequency;

for i = 1:length(general)
    general_frequence = 0;
    specific_pathways_i = specific{i};
    for j = 1:length(specific_pathways_i)
        [is, pos] = ismember(specific_pathways_i{j},pathways2);
        if is
            general_frequence = general_frequence + frequency2(pos);
            pathways2(pos) = [];
            frequency2(pos) = [];
        end
    end
    
    if general_frequence>0
        pos = 0;
        allPos = find(frequency2<=general_frequence);
        if ~isempty(allPos)
            pos = allPos(end);
        end
        pos = pos+1;
        if pos==1
            pathways2 = [general{i};pathways2];
            frequency2 = [general_frequence;frequency2];
        elseif pos>length(pathways2)
            pathways2 = [pathways2;general{i}];
            frequency2 = [frequency2;general_frequence];
        else
            pathways2 = [pathways2(1:pos-1);general{i};pathways2(pos:end)];
            frequency2 = [frequency2(1:pos-1);general_frequence;frequency2(pos:end)];
        end
    end
end

end