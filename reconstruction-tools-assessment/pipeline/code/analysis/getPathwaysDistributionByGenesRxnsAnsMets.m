function getPathwaysDistributionByGenesRxnsAnsMets(models, geneIDs, geneMatrix, oldRxnIDs, oldRxnMatrix, rxnIDs, rxnMatrix,...
    rxnIDs2, rxnMatrix2, metIDs, metMatrix,  names, species, bigg)
global rootFolder
if exist(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],'file')==2
    delete(['pathwayDistributionGeneRxnsMets_' species '.xlsx']);
end

posFirstRaven = find(~cellfun(@isempty, strfind(names,'RAVEN')));
posFirstRaven = posFirstRaven(1);

fileName = fullfile(rootFolder, 'pipeline', 'reconstructions', 'Inputs', upper(species), 'annotation.gb');
[Genes, Secuencias, oldGenes, Enzimas, EC]=...
    TransformarFASTA_ParaPantographGeneral(fileName, '', 1, 0);

n_genes = length(models{1}.genes);
%how many genes were recovered by all the networks
pos1 = find(sum(geneMatrix(1:n_genes, 2:length(models)),2)==length(models)-1);
n_g1 = length(pos1);

%export the list of genes
fileid = fopen('geneList.txt','w+');
for i = 1:n_g1;  fprintf(fileid, '%s\n',geneIDs{pos1(i),1}); end
fclose(fileid);

gpfg_filepath = which('getPathwaysFromGenes.py'); 
geneList_path = fullfile(rootFolder, 'results', species, 'geneList.txt');
workingFolder = fullfile(rootFolder, 'results', species);
system(['python "' gpfg_filepath '" "' species '" "' geneList_path '" "' workingFolder '"'])

filename = 'pathwaysDistributionFromGenes.csv';
delimiterIn = '\t';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
pathways = A.textdata;
frequency = A.data;
[pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency);
[f2, pos] = sort(frequency2, 'descend');
p2 = pathways2(pos);
allData_g1 = [p2 num2cell(f2)];
[is, pos] = ismember('Global and overview maps',allData_g1(:,1)); 
if is; allData_g1(pos,:) = []; f2(pos) = []; end;
allData_g1 = [allData_g1, num2cell(100*f2(1:end)/sum(f2(1:end)))];
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],allData_g1,'allData_g1');

%how many genes were not recovered by any the networks
pos3 = find(sum(geneMatrix(1:n_genes, 2:length(models)),2)==0);
n_g3 = length(pos3);
%which are the functions associated
g3 = models{1}.genes(pos3);
g3_r = cellfun(@(x) strjoin(getRxn_cobraFormat(models{1}, find(~cellfun(@isempty, strfind(models{1}.grRules, x)))), ';'), g3, 'UniformOutput', 0);
g3_rn = cellfun(@(x) strjoin(getRxn_cobraFormat(models{1}, find(~cellfun(@isempty, strfind(models{1}.grRules, x))),1), ';'), g3, 'UniformOutput', 0);

posAux = cell2mat(arrayfun(@(x)find(strcmp(x,Genes)),g3,'UniformOutput',false))';
functions3 = Enzimas(posAux);
f3 = cell(size(g3)); for i = 1:length(g3); pos = find(strcmp(Genes, g3(i))); if ~isempty(pos); f3{i} = Enzimas{pos}; end; end;
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],functions3,'functions3');
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],[g3 f3 g3_r g3_rn],'functions3_r');

fileid = fopen('geneList.txt','w+');
for i = 1:n_g3;  fprintf(fileid, '%s\n',geneIDs{pos3(i),1}); end
fclose(fileid);

system(['python "' gpfg_filepath '" "' species '" "' geneList_path '" "' workingFolder '"'])

filename = 'pathwaysDistributionFromGenes.csv';
delimiterIn = '\t';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
pathways = A.textdata;
frequency = A.data;
[pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency);
[f2, pos] = sort(frequency2, 'descend');
p2 = pathways2(pos);
allData_g3 = [p2 num2cell(f2)];
[is, pos] = ismember('Global and overview maps',allData_g3(:,1)); 
if is; allData_g3(pos,:) = []; f2(pos) = []; end;
allData_g3 = [allData_g3, num2cell(100*f2(1:end)/sum(f2(1:end)))];
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],allData_g3,'allData_g3');


%how many additional genes were predicted by all the networks
pos2 = find(sum(geneMatrix(n_genes+1:end, 2:length(models)),2)==length(models)-1);
n_g2 = length(pos2);
%which are the functions associated
g2 = geneIDs(n_genes+pos2,2);
models{4} = creategrRulesField(models{4});
g2_r = cellfun(@(x) strjoin(getRxn_cobraFormat(models{4}, find(~cellfun(@isempty, strfind(models{4}.grRules, x)))), ';'), g2, 'UniformOutput', 0);
g2_rn = cellfun(@(x) strjoin(getRxn_cobraFormat(models{4}, find(~cellfun(@isempty, strfind(models{4}.grRules, x))),1), ';'), g2, 'UniformOutput', 0);

posAux = cell2mat(arrayfun(@(x)find(strcmp(x,Genes)),g2,'UniformOutput',false))';
functions2 = Enzimas(posAux);
f2 = cell(size(g2)); for i = 1:length(g2); pos = find(strcmp(Genes, g2(i))); if ~isempty(pos); f3{i} = Enzimas{pos}; end; end;

xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],functions2,'functions2');
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],[g2 f2 g2_r g2_rn],'functions2_r');
info = [[{'genes in manual model and in all drafts'};{'genes not in manual model but in all drafts'};{'genes in manual model and in any draft'}], num2cell([n_g1; n_g2; n_g3])];
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],info,'summary_genes');


%% reactions
n_rxns = length(models{1}.rxns);

rxnMatrix3 = rxnMatrix;
rxnIDs3 = rxnIDs;
for i = 1:n_rxns
    for j = 2:length(models)
        if rxnMatrix2(i,j)==1
            rxnMatrix3(i,j)=1;
            rxnIDs3(i,j+1)=rxnIDs2(i,j+1);
        end
    end
end

%which reactions were found by all
pos1 = find(sum(rxnMatrix3(1:n_rxns, 2:length(models)),2)==length(models)-1);
n_r1 = length(pos1);
r1 = rxnIDs3(pos1,posFirstRaven+1);

%export the list of genes
fileid = fopen('rxnList.txt','w+');
for i = 1:n_r1; fprintf(fileid, '%s\n',r1{i}); end
fclose(fileid);

gpfr_filepath = which('getPathwaysFromRxns.py'); 
rxnList_path = fullfile(rootFolder, 'results', species, 'rxnList.txt');
workingFolder = fullfile(rootFolder, 'results', species);
system(['python "' gpfr_filepath '" "' species '" "' rxnList_path '" "' workingFolder '"'])

filename = 'pathwaysDistributionFromRxns.csv';
delimiterIn = '\t';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
pathways = A.textdata;
frequency = A.data;
[pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency);
[f2, pos] = sort(frequency2, 'descend');
p2 = pathways2(pos);
allData_r1 = [p2 num2cell(f2)];
[is, pos] = ismember('Global and overview maps',allData_r1(:,1)); 
if is; allData_r1(pos,:) = []; f2(pos) = []; end;
allData_r1 = [allData_r1, num2cell(100*f2(1:end)/sum(f2(1:end)))];
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],allData_r1,'allData_r1');

%which reactions were not found by any
pos2 = find(sum(rxnMatrix3(1:n_rxns, 2:length(models)),2)==0);
n_r2 = length(pos2);
r2 = models{1}.rxns(pos2);
noGPR_pos = find(cellfun(@isempty, models{1}.grRules(pos2)));
noGPR_rxns = r2(noGPR_pos);
GPR_rxns = setdiff(r2,noGPR_rxns);

r1_MO_specific = intersect(findOrganismSpecificRxnsAndMets(models{1}, 1, ['_' upper(species(1:2))], ''), r2);
r1_not_MO_specific = setdiff(r2, r1_MO_specific);

r1_not_MO_specific_exchange = r1_not_MO_specific(find(~cellfun(@isempty, strfind(r1_not_MO_specific, 'EX'))));
r1_not_MO_specific_exchange_in_bigg = intersect(r1_not_MO_specific_exchange, bigg.rxns);
r1_not_MO_specific_exchange_not_in_bigg = setdiff(r1_not_MO_specific_exchange, bigg.rxns);
r1_not_MO_specific_exchange_w_gpr = intersect(r1_not_MO_specific_exchange, GPR_rxns);
r1_not_MO_specific_exchange_wo_gpr = setdiff(r1_not_MO_specific_exchange, GPR_rxns);

r1_not_MO_specific_exchange_in_bigg_w_gpr = intersect(r1_not_MO_specific_exchange_in_bigg, r1_not_MO_specific_exchange_w_gpr);
r1_not_MO_specific_exchange_in_bigg_wo_gpr = intersect(r1_not_MO_specific_exchange_in_bigg, r1_not_MO_specific_exchange_wo_gpr);
r1_not_MO_specific_exchange_not_in_bigg_w_gpr = intersect(r1_not_MO_specific_exchange_not_in_bigg, r1_not_MO_specific_exchange_w_gpr);
r1_not_MO_specific_exchange_not_in_bigg_wo_gpr = intersect(r1_not_MO_specific_exchange_not_in_bigg, r1_not_MO_specific_exchange_wo_gpr);

r1_not_MO_specific_transport = intersect(models{1}.rxns(find(getTransportRxns(models{1}))), r1_not_MO_specific);
r1_not_MO_specific_transport_in_bigg = intersect(r1_not_MO_specific_transport, bigg.rxns);
r1_not_MO_specific_transport_not_in_bigg = setdiff(r1_not_MO_specific_transport, bigg.rxns);
r1_not_MO_specific_transport_w_gpr = intersect(r1_not_MO_specific_transport, GPR_rxns);
r1_not_MO_specific_transport_wo_gpr = setdiff(r1_not_MO_specific_transport, GPR_rxns);

r1_not_MO_specific_transport_in_bigg_w_gpr = intersect(r1_not_MO_specific_transport_in_bigg, r1_not_MO_specific_transport_w_gpr);
r1_not_MO_specific_transport_in_bigg_wo_gpr = intersect(r1_not_MO_specific_transport_in_bigg, r1_not_MO_specific_transport_wo_gpr);
r1_not_MO_specific_transport_not_in_bigg_w_gpr = intersect(r1_not_MO_specific_transport_not_in_bigg, r1_not_MO_specific_transport_w_gpr);
r1_not_MO_specific_transport_not_in_bigg_wo_gpr = intersect(r1_not_MO_specific_transport_not_in_bigg, r1_not_MO_specific_transport_wo_gpr);

r1_not_MO_specific_remaining = setdiff(r1_not_MO_specific, union(r1_not_MO_specific_exchange, r1_not_MO_specific_transport));
r1_not_MO_specific_remaining_in_bigg = intersect(r1_not_MO_specific_remaining, bigg.rxns);
r1_not_MO_specific_remaining_not_in_bigg = setdiff(r1_not_MO_specific_remaining, bigg.rxns);
r1_not_MO_specific_remaining_w_gpr = intersect(r1_not_MO_specific_remaining, GPR_rxns);
r1_not_MO_specific_remaining_wo_gpr = setdiff(r1_not_MO_specific_remaining, GPR_rxns);

r1_not_MO_specific_remaining_in_bigg_w_gpr = intersect(r1_not_MO_specific_remaining_in_bigg, r1_not_MO_specific_remaining_w_gpr);
r1_not_MO_specific_remaining_in_bigg_wo_gpr = intersect(r1_not_MO_specific_remaining_in_bigg, r1_not_MO_specific_remaining_wo_gpr);
r1_not_MO_specific_remaining_not_in_bigg_w_gpr = intersect(r1_not_MO_specific_remaining_not_in_bigg, r1_not_MO_specific_remaining_w_gpr);
r1_not_MO_specific_remaining_not_in_bigg_wo_gpr = intersect(r1_not_MO_specific_remaining_not_in_bigg, r1_not_MO_specific_remaining_wo_gpr);


numbers = [n_r2; length(r1_MO_specific); length(r1_not_MO_specific);...
    length(r1_not_MO_specific_exchange); length(r1_not_MO_specific_exchange_in_bigg); length(r1_not_MO_specific_exchange_in_bigg_w_gpr); ...
    length(r1_not_MO_specific_exchange_in_bigg_wo_gpr); length(r1_not_MO_specific_exchange_not_in_bigg);length(r1_not_MO_specific_exchange_not_in_bigg_w_gpr);...
    length(r1_not_MO_specific_exchange_not_in_bigg_wo_gpr);...
    length(r1_not_MO_specific_transport); length(r1_not_MO_specific_transport_in_bigg); length(r1_not_MO_specific_transport_in_bigg_w_gpr);...
    length(r1_not_MO_specific_transport_in_bigg_wo_gpr); length(r1_not_MO_specific_transport_not_in_bigg);...
    length(r1_not_MO_specific_transport_not_in_bigg_w_gpr); length(r1_not_MO_specific_transport_not_in_bigg_wo_gpr);...
    length(r1_not_MO_specific_remaining); length(r1_not_MO_specific_remaining_in_bigg); length(r1_not_MO_specific_remaining_in_bigg_w_gpr);...
    length(r1_not_MO_specific_remaining_in_bigg_wo_gpr); length(r1_not_MO_specific_remaining_not_in_bigg);...
    length(r1_not_MO_specific_remaining_not_in_bigg_w_gpr); length(r1_not_MO_specific_remaining_not_in_bigg_wo_gpr)];

r1_in_bigg = intersect(r2,bigg.rxns);
r1_exchange = r1_in_bigg(find(~cellfun(@isempty, strfind(r1_in_bigg, 'EX'))));
r1_transport = intersect(models{1}.rxns(find(getTransportRxns(models{1}))), r1_in_bigg);
[r1_remaining, posrem] = setdiff(r1_in_bigg,union(r1_exchange, r1_transport));
pathwaysRemaining = models{1}.subSystems(posrem);

r1_not_in_bigg = setdiff(r2,bigg.rxns);
r1_not_MO_specific_not_in_BIGG = setdiff(r1_not_MO_specific,bigg.rxns);

%which reactions were found by all but they are not in the manually curated
%model

% load(['rxnMatrix_' species '.mat']);
% load(['rxnIDs_' species '.mat']);

rxnMatrix3 = oldRxnMatrix;
rxnIDs3 = oldRxnIDs;
for i = 1:n_rxns
    for j = 2:length(models)
        if rxnMatrix2(i,j)==1
            rxnMatrix3(i,j)=1;
            rxnIDs3(i,j+1)=rxnIDs2(i,j+1);
        end
    end
end

pos3 = find(sum(rxnMatrix3(n_rxns+1:end, 2:length(models)),2)==length(models)-1);
if ~isempty(pos3)
    n_r3 = length(pos3);
    r3 = rxnIDs3(n_rxns+pos3,posFirstRaven+1);
    
    fileid = fopen('rxnList.txt','w+');
    for i = 1:n_r3; fprintf(fileid, '%s\n',r3{i}); end
    fclose(fileid);
    
    system(['python "' gpfr_filepath '" "' species '" "' rxnList_path '" "' workingFolder '"'])
    
    filename = 'pathwaysDistributionFromRxns.csv';
    delimiterIn = '\t';
    headerlinesIn = 0;
    A = importdata(filename,delimiterIn,headerlinesIn);
    pathways = A.textdata;
    frequency = A.data;
    [pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency);
    [f2, pos] = sort(frequency2, 'descend');
    p2 = pathways2(pos);
    allData_r3 = [p2 num2cell(f2)];
    [is, pos] = ismember('Global and overview maps',allData_r3(:,1));
    if is; allData_r3(pos,:) = []; f2(pos) = []; end;
    allData_r3 = [allData_r3, num2cell(100*f2(1:end)/sum(f2(1:end)))];
    xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],allData_r3,'allData_r3');
    
    info = [[{'reactions in manual model and in all drafts'};{'reactions not in manual model but in all drafts'};{'reactions in manual model and in any draft'}], num2cell([n_r1; n_r3; n_r2])];
    xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],info,'summary_reactions');
end
%% metabolites
n_mets = length(models{1}.mets);
% load(['metMatrix_' species '.mat']);
% load(['metIDs_' species '.mat']);

%which metabolites were found by all
pos1 = find(sum(metMatrix(1:n_mets, 2:length(models)),2)==length(models)-1);
n_m1 = length(pos1);
m1 = metIDs(pos1,posFirstRaven+1);
% m1 = models{posFirstRaven}.mets(pos1);

fileid = fopen('metList.txt','w+');
for i = 1:n_m1; fprintf(fileid, '%s\n',regexprep(m1{i},'_c','')); end
fclose(fileid);

gpfm_filepath = which('getPathwaysFromMets.py'); 
metList_path = fullfile(rootFolder, 'results', species, 'metList.txt');
workingFolder = fullfile(rootFolder, 'results', species);
system(['python "' gpfm_filepath '" "' species '" "' metList_path '" "' workingFolder '"'])

filename = 'pathwaysDistributionFromMets.csv';
delimiterIn = '\t';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
pathways = A.textdata;
frequency = A.data;
[pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency);
[f2, pos] = sort(frequency2, 'descend');
p2 = pathways2(pos);
allData_m1 = [p2 num2cell(f2)];
[is, pos] = ismember('Global and overview maps',allData_m1(:,1)); 
if is; allData_m1(pos,:) = []; f2(pos) = []; end;
allData_m1 = [allData_m1, num2cell(100*f2(1:end)/sum(f2(1:end)))];
xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],allData_m1,'allData_m1');

%which metabolites were found by none
pos2 = find(sum(metMatrix(1:n_mets, 2:length(models)),2)==0);
n_m2 = length(pos2);
m2 = models{1}.mets(pos2);

m2_MO_specific = m2(find(~cellfun(@isempty, strfind(m2, upper(species(1:2))))));
m2_in_bigg = intersect(m2,bigg.mets);
m2_not_in_bigg = setdiff(m2,bigg.mets);
m2_not_MO_specific_not_in_BIGG = setdiff(setdiff(m2,m2(find(~cellfun(@isempty, strfind(m2, upper(species)))))),bigg.mets);

%which metabolites were found by all but they are not in the original model
pos3 = find(sum(metMatrix(n_mets+1:end, 2:length(models)),2)==length(models)-1);
if ~isempty(pos3)
    n_m3 = length(pos3);
    m3 = metIDs(n_mets+pos3, posFirstRaven+1);
    
    fileid = fopen('metList.txt','w+');
    for i = 1:n_m3; fprintf(fileid, '%s\n',regexprep(m3{i},'_c','')); end
    fclose(fileid);
    
    system(['python "' gpfm_filepath '" "' species '" "' metList_path '" "' workingFolder '"'])
    
    filename = 'pathwaysDistributionFromMets.csv';
    delimiterIn = '\t';
    headerlinesIn = 0;
    A = importdata(filename,delimiterIn,headerlinesIn);
    pathways = A.textdata;
    frequency = A.data;
    [pathways2, frequency2] = reducePathwaysDistributionFromSpecificToGeneral(pathways, frequency);
    [f2, pos] = sort(frequency2, 'descend');
    p2 = pathways2(pos);
    allData_m3 = [p2 num2cell(f2)];
    [is, pos] = ismember('Global and overview maps',allData_m3(:,1));
    if is; allData_m3(pos,:) = []; f2(pos) = []; end;
    allData_m3 = [allData_m3, num2cell(100*f2(1:end)/sum(f2(1:end)))];
    xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],allData_m3,'allData_m3');
    
    info = [[{'metabolites in manual model and in all drafts'};{'metabolites not in manual model but in all drafts'};{'metabolites in manual model and in any draft'}], num2cell([n_m1; n_m3; n_m2])];
    xlswrite(['pathwayDistributionGeneRxnsMets_' species '.xlsx'],info,'summary_metabolites');
end
end
