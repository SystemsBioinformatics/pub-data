function [metMXNlikeIDs,metOtherIDs] = getMNX_like_dictionaryForLocalDictionary(keys,values,language_keys,language_values)

metOtherIDs = cell(1000,1); for i = 1:length(metOtherIDs); metOtherIDs{i} = ''; end;
metMXNlikeIDs = cell(1000,1); for i = 1:length(metMXNlikeIDs); metMXNlikeIDs{i} = ''; end;

keys = removeCompartmentFromMets(keys);
values = removeCompartmentFromMets(values);

base = 'MNXlike';

n_elements = 0;
n_MNXlikeIDs = 0;
for i = 1:length(keys)
    if ismember([language_keys{i} ':' keys{i}],metOtherIDs) && ~ismember([language_values{i} ':' values{i}],metOtherIDs)
        pos = find(strcmp(metOtherIDs,[language_keys{i} ':' keys{i}]));
        n_elements = n_elements+1;
        metOtherIDs{n_elements} = [language_values{i} ':' values{i}];
        metMXNlikeIDs{n_elements} = metMXNlikeIDs{pos};
    elseif ismember([language_values{i} ':' values{i}],metOtherIDs) && ~ismember([language_keys{i} ':' keys{i}],metOtherIDs)
        pos = find(strcmp(metOtherIDs,[language_values{i} ':' values{i}]));
        n_elements = n_elements+1;
        metOtherIDs{n_elements} = [language_keys{i} ':' keys{i}];
        metMXNlikeIDs{n_elements} = metMXNlikeIDs{pos};
    elseif ismember([language_values{i} ':' values{i}],metOtherIDs) && ismember([language_keys{i} ':' keys{i}],metOtherIDs)
        disp(['pair ' [language_keys{i} ':' keys{i}] ,' ' [language_values{i} ':' values{i}] ' already exist in the dictionary'])
        
    else
        n_MNXlikeIDs = n_MNXlikeIDs+1;
        n_elements = n_elements+1;
        metOtherIDs{n_elements} = [language_keys{i} ':' keys{i}];
        metMXNlikeIDs{n_elements} = createNewID(base, 5, n_MNXlikeIDs);
        n_elements = n_elements+1;
        metOtherIDs{n_elements} = [language_values{i} ':' values{i}];
        metMXNlikeIDs{n_elements} = createNewID(base, 5, n_MNXlikeIDs);
    end
end

metOtherIDs = metOtherIDs(1:n_elements);
metMXNlikeIDs = metMXNlikeIDs(1:n_elements);

[metMXNlikeIDs, posSorted] = sort(metMXNlikeIDs);
metOtherIDs = metOtherIDs(posSorted);
end