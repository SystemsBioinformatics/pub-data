function [finalScores, fieldCompleteness, cellCompleteness, tableCompleteness] = assessModelCompleteness(models, names, species)


fields = {'rxns';
    'rxnNames';
    'rules';
    'subSystems';
    'mets'
    'metNames'
    'metFormulas'
    'metCharges'
    'genes'};
fieldCompleteness = zeros(length(fields),length(models));
cellCompleteness = zeros(length(fields),length(models));

for i = 1:length(models)
    if isfield(models{i}, 'subSystems')
        if iscell(models{i}.subSystems)
            subs = models{i}.subSystems;
            if iscell(subs{1})
                for j = 1:length(subs)
                    subs{j} = subs{j}{1};
                end
                models{i}.subSystems = subs;
            end
        end
    end
    if isfield(models{i}, 'metCharges') && isfield(models{i}, 'metFormulas')
        for j = 1:length(models{i}.metFormulas)
            if isempty(models{i}.metFormulas{j}) || ~isempty(strfind(models{i}.metFormulas{j}, 'R')) || ~isempty(strfind(models{i}.metFormulas{j}, 'X'))
                models{i}.metCharges(j) = nan;
            end
        end
    elseif isfield(models{i}, 'metFormulas') && ~isfield(models{i}, 'metCharges')
        
    end
    if isfield(models{i}, 'grRules') && ~isfield(models{i}, 'rules')
        models{i} = createRulesFromgrRules(models{i});
    end
end

fprintf('assessing cell completeness: progress %2.0f %%\n', 0); 
for i = 1:length(models)
    model_i = models{i};
    for j = 1:length(fields)
        if isfield(models{i},fields{j})
            field = eval(['model_i.' fields{j}]);
            fieldCompleteness(j,i) = 1;
            if iscell(field)
                cellCompleteness(j,i) = length(find(cellfun(@isempty, field)==0))/length(field);
            elseif isnumeric( eval(['model_i.' fields{j}]))
                cellCompleteness(j,i) = length(find(arrayfun(@isnan, field)==0))/length(field);
            end
            
        else
            fieldCompleteness(j,i) = 0;
            cellCompleteness(j,i) = 0;
        end
    end
    fprintf('assessing cell completeness: progress %2.0f %%\n', 100*i/length(models)); 
end

finalScores = mean(cellCompleteness,1);

if size(finalScores,2)>size(finalScores,1)
    finalScores = finalScores';
end
fieldCompleteness = fieldCompleteness';
cellCompleteness = cellCompleteness';

if size(names,2)>size(names,1)
    names = names';
end
tableCompleteness = [[{'model'}, {'final_score'} strcat('is_present: ', fields') strcat('completeness: ', fields')]; [names, num2cell(finalScores), num2cell(fieldCompleteness), num2cell(cellCompleteness)]];


if exist(['completeness_assessment_' species '.xlsx'], 'file') ==2
    delete(['completeness_assessment_' species '.xlsx'])
end
xlswrite(['completeness_assessment_' species], tableCompleteness);

end