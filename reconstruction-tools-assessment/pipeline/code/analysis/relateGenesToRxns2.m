function relateGenesToRxns2(models, abbr, geneMatrix, rxnMatrix, rxnIDs, rxnMatrix2, rxnIDs2, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, species, excelFileName)

global rootFolder

cd([rootFolder filesep 'results' filesep species])
abbr = abbr(2:end);
manualModel = models{1};
n_genes = length(manualModel.genes); 
n_rxns = length(manualModel.rxns);

rxnMatrix3 = rxnMatrix;
rxnIDs3 = rxnIDs;
for i = 1:n_rxns
    for j = 3:length(models)+1
        if rxnMatrix2(i,j-1)==1
            rxnMatrix3(i,j-1)=1;
            rxnIDs3(i,j)=rxnIDs2(i,j);
        end
    end
end
  
gene_rxns_Matrix = cell(size(geneMatrix));
gene_rxns_Matrix_specific = cell(size(geneMatrix));

n_NotRecovered = cell(size(geneMatrix));
n_Additional = cell(size(geneMatrix));
n_AdditionalInModel = cell(size(geneMatrix));
n_AdditionalNotInModel = cell(size(geneMatrix));

percentagesNotRecovered = cell(size(geneMatrix));
percentagesAdditional = cell(size(geneMatrix));
percentagesAdditionalInModel = cell(size(geneMatrix));
percentagesAdditionalNotInModel = cell(size(geneMatrix));

p_recoveredRxns_perGene_perModel = cell(length(models)-1,1);
p_additionalRxns_perGene_perModel = cell(length(models)-1,1);
p_additionalRxnsInModel_perGene_perModel = cell(length(models)-1,1);
p_additionalRxnsNotInModel_perGene_perModel = cell(length(models)-1,1);

n_recoveredRxns_perGene_perModel = cell(length(models)-1,1);
n_additionalRxns_perGene_perModel = cell(length(models)-1,1);
n_additionalRxnsInModel_perGene_perModel = cell(length(models)-1,1);
n_additionalRxnsNotInModel_perGene_perModel = cell(length(models)-1,1);
cases_perModel = cell(length(models)-1,1);
tabulatedCases_perModel = cell(length(models)-1,1);
percentagesTabulationCases_perModel = zeros(length(models)-1,8);

base = pwd;
if ~exist([base filesep 'analysis_gene_rxns'],'dir')
    mkdir([base filesep 'analysis_gene_rxns'])
end
cd([base filesep 'analysis_gene_rxns']);
base2 = pwd;

for i = 2:length(models) 
    cd(base2)
    if ~exist([base2 filesep 'model' num2str(i)],'dir')
        mkdir([base2 filesep 'model' num2str(i)])
    end
    cd([base2 filesep 'model' num2str(i)]);
    if ~exist([base2 filesep 'model' num2str(i) filesep 'toCheck'],'dir')
        mkdir([base2 filesep 'model' num2str(i) filesep 'toCheck'])
    end
    toCheck = [base2 filesep 'model' num2str(i) filesep 'toCheck'];
    
    pos_genes_i = find(geneMatrix(1:n_genes,i));
    genes_i = manualModel.genes(find(geneMatrix(1:n_genes,i)));
    
    n_recoveredRxns_perGene = zeros(size(genes_i));
    n_additionalRxns_perGene = zeros(size(genes_i));
    n_additionalRxnsInModel_perGene = zeros(size(genes_i));
    n_additionalRxnsNotInModel_perGene = zeros(size(genes_i));
    
    p_recoveredRxns_perGene = zeros(size(genes_i));
    p_additionalRxns_perGene = zeros(size(genes_i));
    p_additionalRxnsInModel_perGene = zeros(size(genes_i));
    p_additionalRxnsNotInModel_perGene = zeros(size(genes_i));
    

    model_i = models{i};
    if isfield(model_i, 'grRules') && ~isfield(model_i, 'rules')
        model_i = createRulesFromgrRules(model_i);
    end
    model_i = creategrRulesField(model_i);

    if i ==1
       disp('') 
    end
    
    cases = zeros(n_genes,1);
    for j = 1:n_genes
        
        % 8 cases
%         recovered  additional_in_model   additional_not_in_model
%       all recovered, at least one add, at least one add. 
% 1          0 0 0 
% 2          0 0 1
% 3          0 1 0 
% 4          1 0 0 
% 5          0 1 1
% 6          1 0 1
% 7          1 1 0
% 8         1 1 1
        fprintf(['i:' num2str(i) ' j:' num2str(j) ' \n'])
        relatedRxns_manual = manualModel.rxns(find(~cellfun(@isempty, strfind(manualModel.grRules, manualModel.genes{j}))));
        relatedRxns_i = model_i.rxns(find(~cellfun(@isempty, strfind(model_i.grRules, manualModel.genes{j}))));

        if isempty(relatedRxns_manual)
           disp('') 
        end
        
        logical_rxnWasRecovered = zeros(size(relatedRxns_manual));
        logical_rxnIsAdditional = ones(size(relatedRxns_i));
        pos_match = zeros(size(relatedRxns_manual));
        n_recovered = 0;
        for k = 1:length(relatedRxns_manual)
            pos = find(strcmp(rxnIDs3(1:n_rxns,2) ,relatedRxns_manual{k}));
            rxn = rxnIDs3(pos, i+1);
            if ~isempty(rxn) 
                [is, posMatch] = ismember(rxn,relatedRxns_i);
                if is
                    n_recovered = n_recovered+1;
                    logical_rxnWasRecovered(k) = 1;
                    pos_match(k) = posMatch;
                    logical_rxnIsAdditional(posMatch) = 0;
                end
            end
        end
        n_NotRecovered{j,i} = length(find(logical_rxnWasRecovered==0));
        
        n_additional = length(find(logical_rxnIsAdditional));
        st1 = '';
        st2 = '';
        st3 = '';
        if n_additional > 0 || ~isempty(find(logical_rxnWasRecovered==0, 1))
            fileID = fopen(['gene' num2str(j) '.txt'],'w+');
        end
        
        if ~isempty(find(logical_rxnWasRecovered==0, 1))
            notRecovered = relatedRxns_manual(logical_rxnWasRecovered==0);
            eqs_notRecovered = getRxn_cobraFormat(manualModel, notRecovered);
            eqs_notRecovered2 = getRxn_cobraFormat(manualModel, notRecovered, 1);

            st_1 = strcat(notRecovered, ':' ,eqs_notRecovered);
            st1 = strjoin (st_1,', ');
            st1 = strcat('notRecovered: ', st1);
            percentagesNotRecovered{j,i} = 100*n_recovered/length(relatedRxns_manual);
            if ismember(j, pos_genes_i)
                p_recoveredRxns_perGene(pos_genes_i==j) = percentagesNotRecovered{j,i};
                n_recoveredRxns_perGene(pos_genes_i==j) = n_recovered;
            end
            
            fprintf(fileID,'notRecovered:\n');
            for k = 1:length(notRecovered)
                fprintf(fileID, [strcat(notRecovered{k}, ':' ,eqs_notRecovered{k}, ' (', eqs_notRecovered2{k} ,')') '\n']);
            end
            fprintf(fileID,'\n\n');
        end
        if n_additional > 0
            n_Additional{j,i} = n_additional;
            additionalRxns = relatedRxns_i(logical_rxnIsAdditional==1);
            eqs_additional = getRxn_cobraFormat(model_i, additionalRxns);
            eqs_additional2 = getRxn_cobraFormat(model_i, additionalRxns,1);
%             inManualModel = ismember(additionalRxns, rxnIDs3(1:n_rxns,2));
            inManualModel = cellfun(@(x) reactionFormulaInModel(models{1},x, 1,1), eqs_additional);
            percentagesAdditional{j,i} = 100*n_additional/length(relatedRxns_manual);
            if ismember(j, pos_genes_i)
                p_additionalRxns_perGene(pos_genes_i==j) = percentagesAdditional{j,i};
                n_additionalRxns_perGene(pos_genes_i==j) = n_additional;
            end
            
            rxnsInManualModel = additionalRxns(inManualModel==1);
            rxnsNotInManualModel = additionalRxns(inManualModel==0);
            
            rxnEqsInManualModel = eqs_additional(inManualModel==1);
            rxnEqsInManualModel2 = eqs_additional2(inManualModel==1);
            
            rxnEqsNotInManualModel = eqs_additional(inManualModel==0);
            rxnEqsNotInManualModel2 = eqs_additional2(inManualModel==0);
            
            if ~isempty(rxnsNotInManualModel)
                st_1 = strcat(rxnsNotInManualModel, ':', rxnEqsNotInManualModel);
                st2 = strjoin (st_1,', ');
                st2 = strcat('. AdditionalNotInModel: ', st2);
                n_AdditionalNotInModel{j,i} = length(rxnsNotInManualModel);
                percentagesAdditionalNotInModel{j,i} = 100*length(n_AdditionalNotInModel{j,i})/length(relatedRxns_manual);
                                
                if ismember(j, pos_genes_i)
                    p_additionalRxnsNotInModel_perGene(pos_genes_i==j) = percentagesAdditionalNotInModel{j,i};
                    n_additionalRxnsNotInModel_perGene(pos_genes_i==j) = n_AdditionalNotInModel{j,i};
                end
                
                fprintf(fileID,'AdditionalNotInModel:\n');
                for k = 1:length(rxnsNotInManualModel)
                    fprintf(fileID, [strcat(rxnsNotInManualModel{k}, ':', rxnEqsNotInManualModel{k}, ' (', rxnEqsNotInManualModel2{k} ,')') '\n']);
                end
                fprintf(fileID,'\n\n');
                
                [metIDs, allMetIDs, allMetNames, allMetPos] = getMetaboliteIDsFromRxns(model_i, rxnsNotInManualModel);
                %AQUI QUEDE
                [diff, pos] = setdiff(allMetIDs, models{1}.mets);
                if ~isempty(pos)
                    st2 = strcat(st2, strcat('. AdditionalMetabolites: ', strjoin(strcat(allMetIDs(pos), ' (', allMetNames(pos), ')'), ', ')));
                    fprintf(fileID,'AdditionalMetabolites:\n');
                    for k = 1:length(pos)
                        fprintf(fileID, [strcat(allMetIDs{pos(k)}, ' (', allMetNames{pos(k)}, ')') '\n']);
                    end
                    fprintf(fileID,'\n\n');
                    
                end
                
                
                if ~isempty(find(logical_rxnWasRecovered==0, 1))
                    copyfile(['gene' num2str(j) '.txt'],[toCheck filesep 'gene' num2str(j) '.txt']);
                end
                
            end
            if ~isempty(rxnsInManualModel)
                st_1 = strcat(rxnsInManualModel, ':', rxnEqsInManualModel);
                st3 = strjoin (st_1,', ');
                st3 = strcat('. AdditionalInModel: ', st3);
                n_AdditionalInModel{j,i} = length(rxnsInManualModel);
                percentagesAdditionalInModel{j,i} = 100*length(n_AdditionalInModel{j,i})/length(relatedRxns_manual);
                
                if ismember(j, pos_genes_i)
                    p_additionalRxnsInModel_perGene(pos_genes_i==j) = percentagesAdditionalInModel{j,i};
                    n_additionalRxnsInModel_perGene(pos_genes_i==j) = n_AdditionalInModel{j,i};
                end
                
                fprintf(fileID,'AdditionalInModel:\n');
                for k = 1:length(rxnsInManualModel)
                    fprintf(fileID, [strcat(rxnsInManualModel{k}, ':', rxnEqsInManualModel{k}, ' (', rxnEqsInManualModel2{k} ,')') '\n']);
                end
            end

        end
        
        if n_NotRecovered{j,i}>0 
            if n_additional==0
                cases(j) = 1;
            elseif isempty(n_AdditionalInModel{j,i}) && n_AdditionalNotInModel{j,i}>0
                cases(j) = 2;
            elseif n_AdditionalInModel{j,i}>0 && isempty(n_AdditionalNotInModel{j,i})
                cases(j) = 3;
            else
                cases(j) = 5;
            end
        else
            if n_additional==0
                cases(j) = 4;
            elseif isempty(n_AdditionalInModel{j,i}) && n_AdditionalNotInModel{j,i}>0
                cases(j) = 6;
            elseif n_AdditionalInModel{j,i}>0 && isempty(n_AdditionalNotInModel{j,i})
                cases(j) = 7;
            else
                cases(j) = 8;
            end 
        end
        
        if n_additional > 0 || ~isempty(find(logical_rxnWasRecovered==0, 1))
            fclose(fileID);
        end
        
        gene_rxns_Matrix{j,i} = strcat(st1, st2, st3);
        if ismember(j, pos_genes_i)
            gene_rxns_Matrix_specific{j,i} = strcat(st1, st2, st3);
        end
    end
    cases_perModel{i-1} = cases;
    tabulatedCases_perModel{i-1} = tabulateWithFixedLabels(cases, [1:8]');
    percentagesTabulationCases_perModel(i-1,:) = tabulatedCases_perModel{i-1}(:,3)';
    
    p_recoveredRxns_perGene_perModel{i-1} = p_recoveredRxns_perGene;
    p_additionalRxns_perGene_perModel{i-1} = p_additionalRxns_perGene;
    p_additionalRxnsInModel_perGene_perModel{i-1} = p_additionalRxnsInModel_perGene;
    p_additionalRxnsNotInModel_perGene_perModel{i-1} = p_additionalRxnsNotInModel_perGene;
    
    n_recoveredRxns_perGene_perModel{i-1} = n_recoveredRxns_perGene;
    n_additionalRxns_perGene_perModel{i-1} = n_additionalRxns_perGene;
    n_additionalRxnsInModel_perGene_perModel{i-1} = n_additionalRxnsInModel_perGene;
    n_additionalRxnsNotInModel_perGene_perModel{i-1} = n_additionalRxnsNotInModel_perGene;
    
end

p_averageRecoveredRxns_perGene_perModel = cellfun(@(x) mean(x),p_recoveredRxns_perGene_perModel);
p_averageAdditionalRxns_perGene_perModel = cellfun(@(x) mean(x),p_additionalRxns_perGene_perModel);
p_averageAdditionalRxnsInModel_perGene_perModel = cellfun(@(x) mean(x),p_additionalRxnsInModel_perGene_perModel);
p_averageAdditionalRxnsNotInModel_perGene_perModel = cellfun(@(x) mean(x),p_additionalRxnsNotInModel_perGene_perModel);

n_averageRecoveredRxns_perGene_perModel = cellfun(@(x) mean(x),n_recoveredRxns_perGene_perModel);
n_averageAdditionalRxns_perGene_perModel = cellfun(@(x) mean(x),n_additionalRxns_perGene_perModel);
n_averageAdditionalRxnsInModel_perGene_perModel = cellfun(@(x) mean(x),n_additionalRxnsInModel_perGene_perModel);
n_averageAdditionalRxnsNotInModel_perGene_perModel = cellfun(@(x) mean(x),n_additionalRxnsNotInModel_perGene_perModel);

gene_rxns_Matrix(1:n_genes,1) = models{1}.genes;
gene_rxns_Matrix = [[{' '} abbr];gene_rxns_Matrix];

gene_rxns_Matrix_specific(1:n_genes,1) = models{1}.genes;
gene_rxns_Matrix_specific = [[{' '} abbr];gene_rxns_Matrix_specific];

cd(base);
save(['gene_rxns_Matrix_' species], 'gene_rxns_Matrix');
% save(['percentagesNotRecovered_' species], 'percentagesNotRecovered');
% save(['percentagesAdditional_' species], 'percentagesAdditional');

if exist([excelFileName '_' species '.xlsx'],'file')==2
    delete([excelFileName '_' species '.xlsx']);
end

labels = {'n_recovered_av', 'n_additional_av', 'n_additionalInModel_av', 'n_additionalNotInModel_av', 'p_recovered_av', 'p_additional_av', 'p_additionalInModel_av', 'p_additionalNotInModel_av'};
info = [n_averageRecoveredRxns_perGene_perModel, n_averageAdditionalRxns_perGene_perModel, n_averageAdditionalRxnsInModel_perGene_perModel, n_averageAdditionalRxnsNotInModel_perGene_perModel,...
    p_averageRecoveredRxns_perGene_perModel, p_averageAdditionalRxns_perGene_perModel, p_averageAdditionalRxnsInModel_perGene_perModel, p_averageAdditionalRxnsNotInModel_perGene_perModel];
info = [labels; num2cell(info)];
info = [[{'models'};abbr'], info];

xlswrite([excelFileName '_' species '.xlsx'],info,'summary');
% xlswrite([excelFileName '_' species '.xlsx'],averageRecoveredRxns_perGene_perModel,'avRec');
% xlswrite([excelFileName '_' species '.xlsx'],averageAdditionalRxns_perGene_perModel,'avAdd');
xlswrite([excelFileName '_' species '.xlsx'],gene_rxns_Matrix,'gene_rxns_Matrix');
xlswrite([excelFileName '_' species '.xlsx'],gene_rxns_Matrix_specific,'gene_rxns_Matrix_sp');
xlswrite([excelFileName '_' species '.xlsx'],percentagesNotRecovered,'percentagesNotRecovered');
xlswrite([excelFileName '_' species '.xlsx'],percentagesAdditional,'percentagesAdditional');
xlswrite([excelFileName '_' species '.xlsx'],[[{''};abbr'],num2cell([1:8;percentagesTabulationCases_perModel])],'percentagesCasesPerModel');

end