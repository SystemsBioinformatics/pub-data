function exportMultipleSolutionsToJson(model,solutions, fileName)

fid = fopen(fileName,'w+');

fstring = '[';
for i = 1:size(solutions,2)
    lstring = regexprep(strjoin(strcat('"', model.rxns, '":', arrayfun(@(x) num2str(x), solutions(:,i), 'UniformOutput', false)),', '), ':',': ');
    string = ['{' lstring '}'];
    if i==1
        fstring = [fstring string];
    else
        fstring = [fstring ', ' string];
    end
end
fstring = [fstring, ']'];
fprintf(fid, fstring);

fclose(fid);

end