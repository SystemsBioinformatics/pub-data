function quantifyCoverageInDatabases(model_original, model_v1, model_v2, modelLanguage, species)


userModel_original = model_original;
userModel_v1 = model_v1;
userModel_v2 = model_v2;

global rootFolder
rootFolder = fileparts(which('initSystemsBioinformaticsToolbox'));
cd(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models', upper(species)))

load(fullfile(rootFolder, 'BIGG', 'bigg_85_more_fixediNF517.mat'))
bigg = bigg_more_fixediNF517;
load(fullfile(rootFolder, 'KEGG', 'keggRxns_v_87.1.mat'))
kegg = model;
% kegg.mets = strcat(kegg.mets,'_c');
% load(fullfile(rootFolder, 'KEGG', 'keggMets_v_87.1.mat'))
% posMetsInKEGGModel = find(ismember(removeCompartmentFromMets(kegg.mets), model.mets)==1);
% kegg.metNames = kegg.mets;
% metswc = removeCompartmentFromMets(kegg.mets(posMetsInKEGGModel));
% for i = 1:length(posMetsInKEGGModel); 
%     kegg.metNames{posMetsInKEGGModel(i)} = model.metNames{find(strcmp(model.mets,metswc{i}))};
% end;
% model = kegg;
% save(fullfile(rootFolder, 'KEGG', 'keggRxns_v_87.1.mat'), 'model');
load(fullfile(rootFolder, 'MetaCyc', 'metacyc_22.0.mat'))
metacyc = transformModelToCBMPYFormat(metacyc);
load(fullfile(rootFolder, 'ModelSEED', 'ModelSEED.mat'))
seed = model;
seed = transformModelToCBMPYFormat(seed);
load(fullfile(rootFolder, 'MNX', 'rxnOtherIDs.mat'))
load(fullfile(rootFolder, 'MNX', 'rxnMNXIDs.mat'))
load(fullfile(rootFolder, 'MNX', 'metOtherIDs.mat'))
load(fullfile(rootFolder, 'MNX', 'metMNXIDs.mat'))
load(fullfile(rootFolder, 'MNX', 'MNX.mat'))
MNX = model;

load(fullfile(rootFolder, 'MNX', 'SysBio_keys.mat'))
load(fullfile(rootFolder, 'MNX', 'SysBio_values.mat'))
load(fullfile(rootFolder, 'MNX', 'SysBio_language_keys.mat'))
load(fullfile(rootFolder, 'MNX', 'SysBio_language_values.mat'))

% load(fullfile(rootFolder, 'MNX', 'MNX.mat'))
mets_bigg_wo_comps = removeCompartmentFromMets(bigg.mets);
mets_kegg_wo_comps = regexprep(kegg.mets,'_[cpe]$','');
mets_metacyc_wo_comps = removeCompartmentFromMets(metacyc.mets);
mets_seed_wo_comps = removeCompartmentFromMets(seed.mets);

isMetBIGG =~cellfun(@isempty, strfind(metOtherIDs, 'bigg:'));
isMetKEGG = ~cellfun(@isempty, strfind(metOtherIDs, 'kegg:'));
isMetMetaCyc = ~cellfun(@isempty, strfind(metOtherIDs, 'metacyc:'));
isMetModelSEED = ~cellfun(@isempty, strfind(metOtherIDs, 'seed:'));

met_biggIDs = metOtherIDs(isMetBIGG);
met_keggIDs = metOtherIDs(isMetKEGG);
met_metacycIDs = metOtherIDs(isMetMetaCyc);
met_seedIDs = metOtherIDs(isMetModelSEED);

metMNX_biggIDs = metMNXIDs(isMetBIGG);
metMNX_keggIDs = metMNXIDs(isMetKEGG);
metMNX_metacycIDs = metMNXIDs(isMetMetaCyc);
metMNX_seedIDs = metMNXIDs(isMetModelSEED);

[compSymbols, ~] = getCompSymbols;

%% model original

mets_wo_comps_original = regexprep(userModel_original.mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
mets_original = strcat([modelLanguage ':'],mets_wo_comps_original);

[~, pos, ~] = intersect(metOtherIDs,mets_original);
metMNXIDs_model_original = metMNXIDs(pos);
mets_original_isInMNX = mets_original(cellfun(@(x) ismember(x,metOtherIDs), mets_original));
MNXIDs_mets_original_isInMNX = getMetIDsInTargetLanguageFromInputLanguage(regexprep(mets_original_isInMNX,'^bigg:',''),modelLanguage,'MNX',metOtherIDs, metMNXIDs);

n_mets_original = length(userModel_original.mets);
n_met_mnx_bigg_original = length(find(cellfun(@(x) ismember(x,mets_bigg_wo_comps), removeCompartmentFromMets(userModel_original.mets))));
n_met_mnx_kegg_original = length(find(cellfun(@(x) ismember(x,metMNX_keggIDs), MNXIDs_mets_original_isInMNX)));
n_met_mnx_metacyc_original = length(find(cellfun(@(x) ismember(x,metMNX_metacycIDs), MNXIDs_mets_original_isInMNX)));
n_met_mnx_seed_original = length(find(cellfun(@(x) ismember(x,metMNX_seedIDs), MNXIDs_mets_original_isInMNX)));

model_bigg = translateModelToTargetLanguage(userModel_original, modelLanguage, 'bigg', bigg, metOtherIDs, metMNXIDs,1,0,1,1);
model_kegg = translateModelToTargetLanguage(userModel_original, modelLanguage, 'kegg', kegg, metOtherIDs, metMNXIDs,1,0,1,1);
model_metacyc = translateModelToTargetLanguage(userModel_original, modelLanguage, 'metacyc', metacyc, metOtherIDs, metMNXIDs,1,0,1,1);
model_seed = translateModelToTargetLanguage(userModel_original, modelLanguage, 'seed', seed, metOtherIDs, metMNXIDs,1,0,1,1);

rxnFormulaWasFound_bigg = cellfun(@(x) reactionFormulaInModel(bigg, x, 1, 1), getRxn_cobraFormat(model_bigg, model_bigg.rxns),'UniformOutput',0);
rxnFormulaWasFound_kegg = cellfun(@(x) reactionFormulaInModel(kegg, x, 1, 1), getRxn_cobraFormat(model_kegg, model_kegg.rxns),'UniformOutput',0);
rxnFormulaWasFound_metacyc = cellfun(@(x) reactionFormulaInModel(metacyc, x, 1, 1), getRxn_cobraFormat(model_metacyc, model_metacyc.rxns),'UniformOutput',0);
rxnFormulaWasFound_seed = cellfun(@(x) reactionFormulaInModel(seed, x, 1, 1), getRxn_cobraFormat(model_seed, model_seed.rxns),'UniformOutput',0);

n_rxn_mnx_bigg_original = length(find(cell2mat(rxnFormulaWasFound_bigg)));
n_rxn_mnx_kegg_original = length(find(cell2mat(rxnFormulaWasFound_kegg)));
n_rxn_mnx_metacyc_original = length(find(cell2mat(rxnFormulaWasFound_metacyc)));
n_rxn_mnx_seed_original = length(find(cell2mat(rxnFormulaWasFound_seed))); 

%% model version 1

mets_wo_comps_v1 = regexprep(userModel_v1.mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
mets_v1 = strcat([modelLanguage ':'],mets_wo_comps_v1);

[~, pos, ~] = intersect(metOtherIDs,mets_v1);
metMNXIDs_model_v1 = metMNXIDs(pos);

mets_v1_isInMNX = mets_v1(cellfun(@(x) ismember(x,metOtherIDs), mets_v1));
MNXIDs_mets_v1_isInMNX = getMetIDsInTargetLanguageFromInputLanguage(regexprep(mets_v1_isInMNX,'^bigg:',''),modelLanguage,'MNX',metOtherIDs, metMNXIDs);

n_mets_v1 = length(userModel_v1.mets);
n_met_mnx_bigg_v1 = length(find(cellfun(@(x) ismember(x,mets_bigg_wo_comps), removeCompartmentFromMets(userModel_v1.mets))));
n_met_mnx_kegg_v1 = length(find(cellfun(@(x) ismember(x,metMNX_keggIDs), MNXIDs_mets_v1_isInMNX)));
n_met_mnx_metacyc_v1 = length(find(cellfun(@(x) ismember(x,metMNX_metacycIDs), MNXIDs_mets_v1_isInMNX)));
n_met_mnx_seed_v1 = length(find(cellfun(@(x) ismember(x,metMNX_seedIDs), MNXIDs_mets_v1_isInMNX)));

model_bigg = userModel_v1;
model_kegg = translateModelToTargetLanguage(userModel_v1, modelLanguage, 'kegg', kegg, metOtherIDs, metMNXIDs,1,0,1,1);
model_metacyc = translateModelToTargetLanguage(userModel_v1, modelLanguage, 'metacyc', metacyc, metOtherIDs, metMNXIDs,1,0,1,1);
model_seed = translateModelToTargetLanguage(userModel_v1, modelLanguage, 'seed', seed, metOtherIDs, metMNXIDs,1,0,1,1);

rxnFormulaWasFound_bigg = cellfun(@(x) reactionFormulaInModel(bigg, x, 1, 1), getRxn_cobraFormat(model_bigg, model_bigg.rxns),'UniformOutput',0);
rxnFormulaWasFound_kegg = cellfun(@(x) reactionFormulaInModel(kegg, x, 1, 1), getRxn_cobraFormat(model_kegg, model_kegg.rxns),'UniformOutput',0);
rxnFormulaWasFound_metacyc = cellfun(@(x) reactionFormulaInModel(metacyc, x, 1, 1), getRxn_cobraFormat(model_metacyc, model_metacyc.rxns),'UniformOutput',0);
rxnFormulaWasFound_seed = cellfun(@(x) reactionFormulaInModel(seed, x, 1, 1), getRxn_cobraFormat(model_seed, model_seed.rxns),'UniformOutput',0);

n_rxn_mnx_bigg_v1 = length(find(cell2mat(rxnFormulaWasFound_bigg)));
n_rxn_mnx_kegg_v1 = length(find(cell2mat(rxnFormulaWasFound_kegg)));
n_rxn_mnx_metacyc_v1 = length(find(cell2mat(rxnFormulaWasFound_metacyc)));
n_rxn_mnx_seed_v1 = length(find(cell2mat(rxnFormulaWasFound_seed))); 

%% model version 2
mets_wo_comps_v2 = regexprep(userModel_v2.mets, strcat('_', compSymbols, '$|\[', compSymbols ,'\]$'), '');
mets_v2 = strcat([modelLanguage ':'],mets_wo_comps_v2);

[~, pos, pos2] = intersect(metOtherIDs,mets_v2);
metMNXIDs_model_v2 = metMNXIDs(pos);

mets_v2_isInMNX = mets_v2(cellfun(@(x) ismember(x,metOtherIDs), mets_v2));
MNXIDs_mets_v2_isInMNX = getMetIDsInTargetLanguageFromInputLanguage(regexprep(mets_v2_isInMNX,'^bigg:',''),modelLanguage,'MNX',metOtherIDs, metMNXIDs);


n_mets_v2 = length(userModel_v2.mets);
n_met_mnx_bigg_v2 = length(find(cellfun(@(x) ismember(x,mets_bigg_wo_comps), removeCompartmentFromMets(userModel_v2.mets))));
n_met_mnx_kegg_v2 = length(find(cellfun(@(x) ismember(x,metMNX_keggIDs), MNXIDs_mets_v2_isInMNX)));
n_met_mnx_metacyc_v2 = length(find(cellfun(@(x) ismember(x,metMNX_metacycIDs), MNXIDs_mets_v2_isInMNX)));
n_met_mnx_seed_v2 = length(find(cellfun(@(x) ismember(x,metMNX_seedIDs), MNXIDs_mets_v2_isInMNX)));

% sample_data = [n_met_mnx_bigg (n_mets-n_met_mnx_bigg);...
%     n_met_mnx_kegg (n_mets-n_met_mnx_kegg);...
%     n_met_mnx_metacyc (n_mets-n_met_mnx_metacyc);...
%     n_met_mnx_seed (n_mets-n_met_mnx_seed)];
% bar(sample_data,'stacked')

% A loop that does num2str conversion only if value is >0
% for i=1:size(sample_data,1);
%     for j=1:size(sample_data,2);
%         if j ==1; c = 'w'; else c = 'k'; end;
%         if sample_data(i,j)>0;
%         labels_stacked=num2str(sample_data(i,j),'%.0f');
%         hText = text(i, sum(sample_data(i,1:j),2), labels_stacked);
%         set(hText, 'VerticalAlignment','top', 'HorizontalAlignment', 'center','FontSize',10, 'Color',c);
%         end
%     end
% end

% a = mets;
% a(pos2) = metMNXIDs_model;
% b = intersect(metMNXIDs_model, metMNX_biggIDs);
% c = intersect(metMNXIDs_model, metMNX_keggIDs);
% d = intersect(metMNXIDs_model, metMNX_metacycIDs);
% e = intersect(metMNXIDs_model, metMNX_seedIDs);
% 
% exportListToTXT({'iLP728';'bigg';'kegg';'metacyc';'seed'},'groupNames.txt')
% exportListToTXT(a,'list1.txt')
% exportListToTXT(b,'list2.txt')
% exportListToTXT(c,'list3.txt')
% exportListToTXT(d,'list4.txt')
% exportListToTXT(e,'list5.txt')

% system(['python "' fullfile(rootFolder,'pipeline','code','analysis','plotVenn.py') '" "' [rootFolder filesep 'pipeline' filesep 'reconstructions' filesep 'manually_curated_models' filesep 'LPL'] '"'])

% model_bigg = translateModelToTargetLanguage(userModel_v2, modelLanguage, 'bigg', bigg, metOtherIDs, metMNXIDs,1,0,1);
model_bigg = userModel_v2;
[model_kegg, removedEmptyRxnsAfterTranslation_kegg, keptRxns_kegg, n_translated_kegg,...
    p_translated_kegg, dictionary_kegg, metsTranslatedIntoTheSameID_kegg, requireManualCheck_kegg] = ...
    translateModelToTargetLanguage(userModel_v2, modelLanguage, 'kegg', kegg, metOtherIDs, metMNXIDs,1,0,1,1);
[model_metacyc, removedEmptyRxnsAfterTranslation_metacyc, keptRxns_metacyc, n_translated_metacyc,...
    p_translated_metacyc, dictionary_metacyc, metsTranslatedIntoTheSameID_metacyc, requireManualCheck_metacyc] = ...
    translateModelToTargetLanguage(userModel_v2, modelLanguage, 'metacyc', metacyc, metOtherIDs, metMNXIDs,1,0,1,1);
[model_seed, removedEmptyRxnsAfterTranslation_seed, keptRxns_seed, n_translated_seed,...
    p_translated_seed, dictionary_seed, metsTranslatedIntoTheSameID_seed, requireManualCheck_seed] = ...
    translateModelToTargetLanguage(userModel_v2, modelLanguage, 'seed', seed, metOtherIDs, metMNXIDs,1,0,1,1);

rxnFormulaWasFound_bigg = cellfun(@(x) reactionFormulaInModel(bigg, x, 1, 1), getRxn_cobraFormat(model_bigg, model_bigg.rxns),'UniformOutput',0);
rxnFormulaWasFound_kegg = cellfun(@(x) reactionFormulaInModel(kegg, x, 1, 1), getRxn_cobraFormat(model_kegg, model_kegg.rxns),'UniformOutput',0);
rxnFormulaWasFound_metacyc = cellfun(@(x) reactionFormulaInModel(metacyc, x, 1, 1), getRxn_cobraFormat(model_metacyc, model_metacyc.rxns),'UniformOutput',0);
rxnFormulaWasFound_seed = cellfun(@(x) reactionFormulaInModel(seed, x, 1, 1), getRxn_cobraFormat(model_seed, model_seed.rxns),'UniformOutput',0);

n_rxn_mnx_bigg_v2 = length(find(cell2mat(rxnFormulaWasFound_bigg)));
n_rxn_mnx_kegg_v2 = length(find(cell2mat(rxnFormulaWasFound_kegg)));
n_rxn_mnx_metacyc_v2 = length(find(cell2mat(rxnFormulaWasFound_metacyc)));
n_rxn_mnx_seed_v2 = length(find(cell2mat(rxnFormulaWasFound_seed))); 

%% plot general progress mets, just manually curated modeels

databases = {'bigg', 'kegg', 'metacyc', 'seed'};

matrix = zeros(3,4,2);
n_mets = length(userModel_v2.mets);
sample_data1 = [n_met_mnx_bigg_original (n_mets-n_met_mnx_bigg_original);...
    n_met_mnx_kegg_original (n_mets-n_met_mnx_kegg_original);...
    n_met_mnx_metacyc_original (n_mets-n_met_mnx_metacyc_original);...
    n_met_mnx_seed_original (n_mets-n_met_mnx_seed_original)];
sample_data2 = [n_met_mnx_bigg_v1 (n_mets-n_met_mnx_bigg_v1);...
    n_met_mnx_kegg_v1 (n_mets-n_met_mnx_kegg_v1);...
    n_met_mnx_metacyc_v1 (n_mets-n_met_mnx_metacyc_v1);...
    n_met_mnx_seed_v1 (n_mets-n_met_mnx_seed_v1)];
sample_data3 = [n_met_mnx_bigg_v2 (n_mets-n_met_mnx_bigg_v2);...
    n_met_mnx_kegg_v2 (n_mets-n_met_mnx_kegg_v2);...
    n_met_mnx_metacyc_v2 (n_mets-n_met_mnx_metacyc_v2);...
    n_met_mnx_seed_v2 (n_mets-n_met_mnx_seed_v2)];

matrix(1,:,:) = sample_data1;
matrix(2,:,:) = sample_data2;
matrix(3,:,:) = sample_data3;

plotBarStackGroups(matrix,{'original','partially mapped','fully mapped'})
% A loop that does num2str conversion only if value is >0
for k = 1:size(matrix,1)
    for i=1:size(matrix,2);
        for j=1:size(matrix,3);
            if j ==1; c = 'w'; else c = 'k'; end;
            labels_stacked=num2str(matrix(k,i,j),'%.0f');
            hText = text((k-1)+0.75+(0.165*(i-1)), sum(matrix(k,i,1:j),3), labels_stacked); % for 4 samples
            set(hText, 'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c);
            hText2 = text((k-1)+0.75+(0.165*(i-1))+(length(databases{i})-4)*0.01, n_mets*1.04+(length(databases{i})-4)*n_mets*0.004, databases{i}); % for 4 samples
            set(hText2,'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c,'Rotation',45);
        end
    end
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, ['general_iterations_mets_manualModel' species '.pdf']);


return;
%% export inconsistencies

mnx_ids = getRxnIDsInTargetLanguageFromInputLanguage(model_bigg.rxns, 'bigg','MNX',rxnOtherIDs, rxnMNXIDs);
kegg_ids_cells = getRxnIDsInTargetLanguageFromInputLanguage(model_bigg.rxns, 'bigg','kegg',rxnOtherIDs, rxnMNXIDs);
kegg_ids = kegg_ids_cells;
for i = 1:length(kegg_ids); if ~isempty(kegg_ids{i}); kegg_ids{i} = strjoin(kegg_ids{i},','); end; end;
metacyc_ids_cells = getRxnIDsInTargetLanguageFromInputLanguage(model_bigg.rxns, 'bigg','metacyc',rxnOtherIDs, rxnMNXIDs);
metacyc_ids = metacyc_ids_cells;
for i = 1:length(metacyc_ids); if ~isempty(metacyc_ids{i}); metacyc_ids{i} = strjoin(metacyc_ids{i},','); end; end;
seed_ids_cells = getRxnIDsInTargetLanguageFromInputLanguage(model_bigg.rxns, 'bigg','seed',rxnOtherIDs, rxnMNXIDs);
seed_ids = seed_ids_cells;
for i = 1:length(seed_ids); if ~isempty(seed_ids{i}); seed_ids{i} = strjoin(seed_ids{i},','); end; end;

n_rxn_mnx_mnx2 = length(find(~cellfun(@isempty,mnx_ids)));
n_rxn_mnx_kegg2 = length(find(~cellfun(@isempty,kegg_ids)));
n_rxn_mnx_metacyc2 = length(find(~cellfun(@isempty,metacyc_ids)));
n_rxn_mnx_seed2 = length(find(~cellfun(@isempty,seed_ids)));


posNotFoundKEGG = intersect(find(cell2mat(rxnFormulaWasFound_kegg)==0), find(~cellfun(@isempty,kegg_ids(getPosOfElementsInArray(model_kegg.rxns, model_bigg.rxns)))));
posNotFoundMetaCyc = intersect(find(cell2mat(rxnFormulaWasFound_metacyc)==0), find(~cellfun(@isempty,metacyc_ids(getPosOfElementsInArray(model_metacyc.rxns, model_bigg.rxns)))));
posNotFoundSEED = intersect(find(cell2mat(rxnFormulaWasFound_seed)==0), find(~cellfun(@isempty,seed_ids(getPosOfElementsInArray(model_seed.rxns, model_bigg.rxns)))));

notFoundKEGG = model_kegg.rxns(posNotFoundKEGG);
notFoundMetaCyc = model_metacyc.rxns(posNotFoundMetaCyc);
notFoundSEED = model_seed.rxns(posNotFoundSEED);

%which of the reactions which were not found in model_kegg have a MNX ID?
info_kegg = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundKEGG)
    ids_i = kegg_ids_cells{getPosOfElementsInArray(notFoundKEGG(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_kegg(n_start,1) = notFoundKEGG(i);
    info_kegg(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundKEGG(i));
    info_kegg(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundKEGG(i),1);
    info_kegg(n_start,4) = getRxn_cobraFormat(model_kegg,notFoundKEGG(i));
    info_kegg(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundKEGG(i), model_bigg.rxns));
    info_kegg(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundKEGG(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_kegg(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},kegg.rxns)
            info_kegg(n_start:n_ends,8) = getRxn_cobraFormat(kegg, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_kegg = info_kegg(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_kegg.xlsx'),info_kegg)

info_metacyc = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundMetaCyc)
    ids_i = metacyc_ids_cells{getPosOfElementsInArray(notFoundMetaCyc(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_metacyc(n_start,1) = notFoundMetaCyc(i);
    info_metacyc(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc(i));
    info_metacyc(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc(i),1);
    info_metacyc(n_start,4) = getRxn_cobraFormat(model_metacyc,notFoundMetaCyc(i));
    info_metacyc(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundMetaCyc(i), model_bigg.rxns));
    info_metacyc(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundMetaCyc(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_metacyc(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},metacyc.rxns)
            info_metacyc(n_start:n_ends,8) = getRxn_cobraFormat(metacyc, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_metacyc = info_metacyc(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_metacyc.xlsx'),info_metacyc)

info_seed = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundSEED)
    ids_i = seed_ids_cells{getPosOfElementsInArray(notFoundSEED(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_seed(n_start,1) = notFoundSEED(i);
    info_seed(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundSEED(i));
    info_seed(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundSEED(i),1);
    info_seed(n_start,4) = getRxn_cobraFormat(model_seed,notFoundSEED(i));
    info_seed(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundSEED(i), model_bigg.rxns));
    info_seed(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundSEED(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_seed(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},seed.rxns)
            info_seed(n_start:n_ends,8) = getRxn_cobraFormat(seed, ids_i);
        end
    end
    n_start = n_ends+1;
end
info_seed = info_seed(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_seed.xlsx'),info_seed)


%% with additional identifiers
dictionary = [SysBio_keys, SysBio_values];
[model_kegg_improved, removedEmptyRxnsAfterTranslation_kegg_improved, keptRxns_kegg_improved, n_translated_kegg_improved,...
    p_translated_kegg_improved, dictionary_kegg_improved, metsTranslatedIntoTheSameID_kegg_improved, requireManualCheck_kegg_improved] = ... 
     transformModelBasedOnSynonymsFromLocalDictionary(model_kegg, kegg, dictionary, SysBio_language_keys, SysBio_language_values, 'bigg','kegg');
[model_metacyc_improved, removedEmptyRxnsAfterTranslation_metacyc_improved, keptRxns_metacyc_improved, n_translated_metacyc_improved,...
    p_translated_metacyc_improved, dictionary_metacyc_improved, metsTranslatedIntoTheSameID_metacyc_improved, requireManualCheck_metacyc_improved] = ...
 transformModelBasedOnSynonymsFromLocalDictionary(model_metacyc, metacyc, dictionary, SysBio_language_keys, SysBio_language_values, 'bigg','metacyc');
[model_seed_improved, removedEmptyRxnsAfterTranslation_seed_improved, keptRxns_seed_improved, n_translated_seed_improved,...
    p_translated_seed_improved, dictionary_seed_improved, metsTranslatedIntoTheSameID_seed_improved, requireManualCheck_seed_improved] = ...
 transformModelBasedOnSynonymsFromLocalDictionary(model_seed, seed, dictionary,SysBio_language_keys, SysBio_language_values, 'bigg','seed');

rxnFormulaWasFound_kegg_improved = cellfun(@(x) reactionFormulaInModel(kegg, x, 1, 1), getRxn_cobraFormat(model_kegg_improved, model_kegg_improved.rxns),'UniformOutput',0);
rxnFormulaWasFound_metacyc_improved = cellfun(@(x) reactionFormulaInModel(metacyc, x, 1, 1), getRxn_cobraFormat(model_metacyc_improved, model_metacyc_improved.rxns),'UniformOutput',0);
rxnFormulaWasFound_seed_improved = cellfun(@(x) reactionFormulaInModel(seed, x, 1, 1), getRxn_cobraFormat(model_seed_improved, model_seed_improved.rxns),'UniformOutput',0);

n_rxn_mnx_kegg_improved = length(find(cell2mat(rxnFormulaWasFound_kegg_improved)));
n_rxn_mnx_metacyc_improved = length(find(cell2mat(rxnFormulaWasFound_metacyc_improved)));
n_rxn_mnx_seed_improved = length(find(cell2mat(rxnFormulaWasFound_seed_improved))); 

posNotFoundKEGG_improved = intersect(find(cell2mat(rxnFormulaWasFound_kegg_improved)==0), find(~cellfun(@isempty,kegg_ids(getPosOfElementsInArray(model_kegg_improved.rxns, model_bigg.rxns)))));
posNotFoundMetaCyc_improved = intersect(find(cell2mat(rxnFormulaWasFound_metacyc_improved)==0), find(~cellfun(@isempty,metacyc_ids(getPosOfElementsInArray(model_metacyc_improved.rxns, model_bigg.rxns)))));
posNotFoundSEED_improved = intersect(find(cell2mat(rxnFormulaWasFound_seed_improved)==0), find(~cellfun(@isempty,seed_ids(getPosOfElementsInArray(model_seed_improved.rxns, model_bigg.rxns)))));

notFoundKEGG_improved = model_kegg_improved.rxns(posNotFoundKEGG_improved);
notFoundMetaCyc_improved = model_metacyc_improved.rxns(posNotFoundMetaCyc_improved);
notFoundSEED_improved = model_seed_improved.rxns(posNotFoundSEED_improved);

info_kegg_improved = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundKEGG_improved)
    disp(i)
    ids_i = kegg_ids_cells{getPosOfElementsInArray(notFoundKEGG_improved(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_kegg_improved(n_start,1) = notFoundKEGG_improved(i);
    info_kegg_improved(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundKEGG_improved(i));
    info_kegg_improved(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundKEGG_improved(i),1);
    info_kegg_improved(n_start,4) = getRxn_cobraFormat(model_kegg_improved,notFoundKEGG_improved(i));
    info_kegg_improved(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundKEGG_improved(i), model_bigg.rxns));
    info_kegg_improved(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundKEGG_improved(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_kegg_improved(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},kegg.rxns)
            info_kegg_improved(n_start:n_ends,8) = getRxn_cobraFormat(kegg, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_kegg_improved = info_kegg_improved(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_kegg_improved.xlsx'),info_kegg_improved)

info_metacyc_improved = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundMetaCyc_improved)
    ids_i = metacyc_ids_cells{getPosOfElementsInArray(notFoundMetaCyc_improved(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_metacyc_improved(n_start,1) = notFoundMetaCyc_improved(i);
    info_metacyc_improved(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc_improved(i));
    info_metacyc_improved(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc_improved(i),1);
    info_metacyc_improved(n_start,4) = getRxn_cobraFormat(model_metacyc_improved,notFoundMetaCyc_improved(i));
    info_metacyc_improved(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundMetaCyc_improved(i), model_bigg.rxns));
    info_metacyc_improved(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundMetaCyc_improved(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_metacyc_improved(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},metacyc.rxns)
            info_metacyc_improved(n_start:n_ends,8) = getRxn_cobraFormat(metacyc, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_metacyc_improved = info_metacyc_improved(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_metacyc_improved.xlsx'),info_metacyc_improved)

info_seed_improved = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundSEED_improved)
    ids_i = seed_ids_cells{getPosOfElementsInArray(notFoundSEED_improved(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_seed_improved(n_start,1) = notFoundSEED_improved(i);
    info_seed_improved(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundSEED_improved(i));
    info_seed_improved(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundSEED_improved(i),1);
    info_seed_improved(n_start,4) = getRxn_cobraFormat(model_seed_improved,notFoundSEED_improved(i));
    info_seed_improved(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundSEED_improved(i), model_bigg.rxns));
    info_seed_improved(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundSEED_improved(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_seed_improved(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},seed.rxns)
            info_seed_improved(n_start:n_ends,8) = getRxn_cobraFormat(seed, ids_i);
        end
    end
    n_start = n_ends+1;
end
info_seed_improved = info_seed_improved(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_seed_improved.xlsx'),info_seed_improved)


%% identify new relationships

% identifyNewRelationshipForModel(model_bigg,bigg,[species '_bigg'],'bigg',[species '_bigg'])
load([species '_bigg'])
model_bigg_improved2 = model; 

identifyNewRelationshipForModel(model_kegg_improved,kegg,[species '_kegg'],'kegg',[species '_kegg'])
load([species '_kegg'])
model_kegg_improved2 = model; 

identifyNewRelationshipForModel(model_seed_improved,seed,[species '_seed'],'seed',[species '_seed'])
load([species '_seed'])
model_seed_improved2 = model; 

identifyNewRelationshipForModel(model_metacyc_improved,metacyc,[species '_metacyc'],'metacyc',[species '_metacyc'])
load([species '_metacyc'])
model_metacyc_improved2 = model; 

n_mets_v2_improved2 = length(model_bigg_improved2.mets);
n_met_bigg_v2_improved2 = length(find(cellfun(@(x) ismember(x,mets_bigg_wo_comps), removeCompartmentFromMets(model_bigg_improved2.mets))));
n_met_kegg_v2_improved2 = length(find(cellfun(@(x) ismember(x,mets_kegg_wo_comps), removeCompartmentFromMets(model_kegg_improved2.mets))));
n_met_metacyc_v2_improved2 = length(find(cellfun(@(x) ismember(x,mets_metacyc_wo_comps), removeCompartmentFromMets(model_metacyc_improved2.mets))));
n_met_seed_v2_improved2 = length(find(cellfun(@(x) ismember(x,mets_seed_wo_comps), removeCompartmentFromMets(model_seed_improved2.mets))));

rxnFormulaWasFound_kegg_improved2  = cellfun(@(x) reactionFormulaInModel(kegg, x, 1, 1), getRxn_cobraFormat(model_kegg_improved2 , model_kegg_improved2 .rxns),'UniformOutput',0);
rxnFormulaWasFound_metacyc_improved2  = cellfun(@(x) reactionFormulaInModel(metacyc, x, 1, 1), getRxn_cobraFormat(model_metacyc_improved2 , model_metacyc_improved2 .rxns),'UniformOutput',0);
rxnFormulaWasFound_seed_improved2  = cellfun(@(x) reactionFormulaInModel(seed, x, 1, 1), getRxn_cobraFormat(model_seed_improved2 , model_seed_improved2 .rxns),'UniformOutput',0);

n_rxn_mnx_kegg_improved2  = length(find(cell2mat(rxnFormulaWasFound_kegg_improved2 )));
n_rxn_mnx_metacyc_improved2  = length(find(cell2mat(rxnFormulaWasFound_metacyc_improved2 )));
n_rxn_mnx_seed_improved2  = length(find(cell2mat(rxnFormulaWasFound_seed_improved2 ))); 

posNotFoundKEGG_improved2  = intersect(find(cell2mat(rxnFormulaWasFound_kegg_improved2 )==0), find(~cellfun(@isempty,kegg_ids(getPosOfElementsInArray(model_kegg_improved2.rxns, model_bigg.rxns)))));
posNotFoundMetaCyc_improved2  = intersect(find(cell2mat(rxnFormulaWasFound_metacyc_improved2 )==0), find(~cellfun(@isempty,metacyc_ids(getPosOfElementsInArray(model_metacyc_improved2.rxns, model_bigg.rxns)))));
posNotFoundSEED_improved2  = intersect(find(cell2mat(rxnFormulaWasFound_seed_improved2 )==0), find(~cellfun(@isempty,seed_ids(getPosOfElementsInArray(model_seed_improved2.rxns, model_bigg.rxns)))));

notFoundKEGG_improved2  = model_kegg_improved2.rxns(posNotFoundKEGG_improved2 );
notFoundMetaCyc_improved2  = model_metacyc_improved2.rxns(posNotFoundMetaCyc_improved2 );
notFoundSEED_improved2  = model_seed_improved2.rxns(posNotFoundSEED_improved2 );

info_kegg_improved2  = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundKEGG_improved2)
    disp(i)
    if i ==1
       disp('') 
    end
    ids_i = kegg_ids_cells{getPosOfElementsInArray(notFoundKEGG_improved2(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_kegg_improved2 (n_start,1) = notFoundKEGG_improved2 (i);
    info_kegg_improved2 (n_start,2) = getRxn_cobraFormat(model_bigg,notFoundKEGG_improved2 (i));
    info_kegg_improved2 (n_start,3) = getRxn_cobraFormat(model_bigg,notFoundKEGG_improved2 (i),1);
    info_kegg_improved2 (n_start,4) = getRxn_cobraFormat(model_kegg_improved2 ,notFoundKEGG_improved2 (i));
    info_kegg_improved2 (n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundKEGG_improved2(i), model_bigg.rxns));
    info_kegg_improved2 (n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundKEGG_improved2(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_kegg_improved2 (n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},kegg.rxns)
            info_kegg_improved2 (n_start:n_ends,8) = getRxn_cobraFormat(kegg, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_kegg_improved2  = info_kegg_improved2 (1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_kegg_improved2 .xlsx'),info_kegg_improved2 )


info_metacyc_improved2 = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundMetaCyc_improved2)
    ids_i = metacyc_ids_cells{getPosOfElementsInArray(notFoundMetaCyc_improved2(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_metacyc_improved2(n_start,1) = notFoundMetaCyc_improved2(i);
    info_metacyc_improved2(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc_improved2(i));
    info_metacyc_improved2(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc_improved2(i),1);
    info_metacyc_improved2(n_start,4) = getRxn_cobraFormat(model_metacyc_improved2,notFoundMetaCyc_improved2(i));
    info_metacyc_improved2(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundMetaCyc_improved2(i), model_bigg.rxns));
    info_metacyc_improved2(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundMetaCyc_improved2(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_metacyc_improved2(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},metacyc.rxns)
            info_metacyc_improved2(n_start:n_ends,8) = getRxn_cobraFormat(metacyc, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_metacyc_improved2 = info_metacyc_improved2(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_metacyc_improved2.xlsx'),info_metacyc_improved2)

info_seed_improved2 = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundSEED_improved2)
    ids_i = seed_ids_cells{getPosOfElementsInArray(notFoundSEED_improved2(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_seed_improved2(n_start,1) = notFoundSEED_improved2(i);
    info_seed_improved2(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundSEED_improved2(i));
    info_seed_improved2(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundSEED_improved2(i),1);
    info_seed_improved2(n_start,4) = getRxn_cobraFormat(model_seed_improved2,notFoundSEED_improved2(i));
    info_seed_improved2(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundSEED_improved2(i), model_bigg.rxns));
    info_seed_improved2(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundSEED_improved2(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_seed_improved2(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},seed.rxns)
            info_seed_improved2(n_start:n_ends,8) = getRxn_cobraFormat(seed, ids_i);
        end
    end
    n_start = n_ends+1;
end
info_seed_improved2 = info_seed_improved2(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_seed_improved2.xlsx'),info_seed_improved2)



%% manual curation 

%add to the dictionary
[k,v,kt,vt] = getDictionaryFromExcelFile(fullfile(rootFolder, 'MNX','manualDictionary.xlsx'));

[model_kegg_improved3, removedEmptyRxnsAfterTranslation_kegg_improved3, keptRxns_kegg_improved3, n_translated_kegg_improved3,...
    p_translated_kegg_improved3, dictionary_kegg_improved3, metsTranslatedIntoTheSameID_kegg_improved3, requireManualCheck_kegg_improved3] = ... 
     transformModelBasedOnSynonymsFromLocalDictionary(model_kegg_improved2, kegg, [k v], kt, vt, 'bigg','kegg');
[model_metacyc_improved3, removedEmptyRxnsAfterTranslation_metacyc_improved3, keptRxns_metacyc_improved3, n_translated_metacyc_improved3,...
    p_translated_metacyc_improved3, dictionary_metacyc_improved3, metsTranslatedIntoTheSameID_metacyc_improved3, requireManualCheck_metacyc_improved3] = ...
 transformModelBasedOnSynonymsFromLocalDictionary(model_metacyc_improved2, metacyc, [k v], kt, vt, 'bigg','metacyc');
[model_seed_improved3, removedEmptyRxnsAfterTranslation_seed_improved3, keptRxns_seed_improved3, n_translated_seed_improved3,...
    p_translated_seed_improved3, dictionary_seed_improved3, metsTranslatedIntoTheSameID_seed_improved3, requireManualCheck_seed_improved3] = ...
 transformModelBasedOnSynonymsFromLocalDictionary(model_seed_improved2, seed, [k v], kt, vt, 'bigg','seed');


n_mets_v2_improved3 = length(model_bigg_improved2.mets);
n_met_kegg_v2_improved3 = length(intersect(bigg.mets,model_kegg_improved3.mets));
n_met_metacyc_v2_improved3 = length(intersect(bigg.mets,model_metacyc_improved3.mets));
n_met_seed_v2_improved3 = length(intersect(bigg.mets,model_seed_improved3.mets));

rxnFormulaWasFound_kegg_improved3 = cellfun(@(x) reactionFormulaInModel(kegg, x, 1, 1), getRxn_cobraFormat(model_kegg_improved3, model_kegg_improved3.rxns),'UniformOutput',0);
rxnFormulaWasFound_metacyc_improved3 = cellfun(@(x) reactionFormulaInModel(metacyc, x, 1, 1), getRxn_cobraFormat(model_metacyc_improved3, model_metacyc_improved3.rxns),'UniformOutput',0);
rxnFormulaWasFound_seed_improved3 = cellfun(@(x) reactionFormulaInModel(seed, x, 1, 1), getRxn_cobraFormat(model_seed_improved3, model_seed_improved3.rxns),'UniformOutput',0);

n_rxn_mnx_kegg_improved3 = length(find(cell2mat(rxnFormulaWasFound_kegg_improved3)));
n_rxn_mnx_metacyc_improved3 = length(find(cell2mat(rxnFormulaWasFound_metacyc_improved3)));
n_rxn_mnx_seed_improved3 = length(find(cell2mat(rxnFormulaWasFound_seed_improved3))); 

posNotFoundKEGG_improved3 = intersect(find(cell2mat(rxnFormulaWasFound_kegg_improved3)==0), find(~cellfun(@isempty,kegg_ids(getPosOfElementsInArray(model_kegg_improved3.rxns, model_bigg.rxns)))));
posNotFoundMetaCyc_improved3 = intersect(find(cell2mat(rxnFormulaWasFound_metacyc_improved3)==0), find(~cellfun(@isempty,metacyc_ids(getPosOfElementsInArray(model_metacyc_improved3.rxns, model_bigg.rxns)))));
posNotFoundSEED_improved3 = intersect(find(cell2mat(rxnFormulaWasFound_seed_improved3)==0), find(~cellfun(@isempty,seed_ids(getPosOfElementsInArray(model_seed_improved3.rxns, model_bigg.rxns)))));

notFoundKEGG_improved3 = model_kegg_improved3.rxns(posNotFoundKEGG_improved3);
notFoundMetaCyc_improved3 = model_metacyc_improved3.rxns(posNotFoundMetaCyc_improved3);
notFoundSEED_improved3 = model_seed_improved3.rxns(posNotFoundSEED_improved3);

info_kegg_improved3 = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundKEGG_improved3)
    ids_i = kegg_ids_cells{getPosOfElementsInArray(notFoundKEGG_improved3(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_kegg_improved3(n_start,1) = notFoundKEGG_improved3(i);
    info_kegg_improved3(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundKEGG_improved3(i));
    info_kegg_improved3(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundKEGG_improved3(i),1);
    info_kegg_improved3(n_start,4) = getRxn_cobraFormat(model_kegg_improved3,notFoundKEGG_improved3(i));
    info_kegg_improved3(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundKEGG_improved3(i), model_bigg.rxns));
    info_kegg_improved3(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundKEGG_improved3(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_kegg_improved3(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},kegg.rxns)
            info_kegg_improved3(n_start:n_ends,8) = getRxn_cobraFormat(kegg, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_kegg_improved3 = info_kegg_improved3(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_kegg_improved3.xlsx'),info_kegg_improved3)

info_metacyc_improved3 = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundMetaCyc_improved3)
    ids_i = metacyc_ids_cells{getPosOfElementsInArray(notFoundMetaCyc_improved3(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_metacyc_improved3(n_start,1) = notFoundMetaCyc_improved3(i);
    info_metacyc_improved3(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc_improved3(i));
    info_metacyc_improved3(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundMetaCyc_improved3(i),1);
    info_metacyc_improved3(n_start,4) = getRxn_cobraFormat(model_metacyc_improved3,notFoundMetaCyc_improved3(i));
    info_metacyc_improved3(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundMetaCyc_improved3(i), model_bigg.rxns));
    info_metacyc_improved3(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundMetaCyc_improved3(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_metacyc_improved3(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},metacyc.rxns)
            info_metacyc_improved3(n_start:n_ends,8) = getRxn_cobraFormat(metacyc, ids_i{j});
        end
    end
    n_start = n_ends+1;
end
info_metacyc_improved3 = info_metacyc_improved3(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_metacyc_improved3.xlsx'),info_metacyc_improved3)

info_seed_improved3 = cell(1000,9);
n_start = 1;
for i = 1:length(notFoundSEED_improved3)
    disp(i)
    ids_i = seed_ids_cells{getPosOfElementsInArray(notFoundSEED_improved3(i), model_bigg.rxns)};
    n_ends = n_start+length(ids_i)-1;
    info_seed_improved3(n_start,1) = notFoundSEED_improved3(i);
    info_seed_improved3(n_start,2) = getRxn_cobraFormat(model_bigg,notFoundSEED_improved3(i));
    info_seed_improved3(n_start,3) = getRxn_cobraFormat(model_bigg,notFoundSEED_improved3(i),1);
    info_seed_improved3(n_start,4) = getRxn_cobraFormat(model_seed_improved3,notFoundSEED_improved3(i));
    info_seed_improved3(n_start,5) = mnx_ids(getPosOfElementsInArray(notFoundSEED_improved3(i), model_bigg.rxns));
    info_seed_improved3(n_start,6) = getRxn_cobraFormat(MNX,mnx_ids(getPosOfElementsInArray(notFoundSEED_improved3(i), model_bigg.rxns)));
    for j = 1:length(ids_i)
        info_seed_improved3(n_start:n_ends,7) = ids_i;
    end
    for j = 1:length(ids_i)
        if ismember(ids_i{j},seed.rxns)
            info_seed_improved3(n_start:n_ends,8) = getRxn_cobraFormat(seed, ids_i);
        end
    end
    n_start = n_ends+1;
end
info_seed_improved3 = info_seed_improved3(1:n_ends,:);
xlswrite(fullfile(rootFolder,'pipeline','reconstructions','manually_curated_models',upper(species),'inconsistency_seed_improved3.xlsx'),info_seed_improved3)


%% plot general progress mets

databases = {'bigg', 'kegg', 'metacyc', 'seed'};

matrix = zeros(4,4,2);
n_mets = length(userModel_v2.mets);
sample_data1 = [n_met_mnx_bigg_original (n_mets-n_met_mnx_bigg_original);...
    n_met_mnx_kegg_original (n_mets-n_met_mnx_kegg_original);...
    n_met_mnx_metacyc_original (n_mets-n_met_mnx_metacyc_original);...
    n_met_mnx_seed_original (n_mets-n_met_mnx_seed_original)];
sample_data2 = [n_met_mnx_bigg_v1 (n_mets-n_met_mnx_bigg_v1);...
    n_met_mnx_kegg_v1 (n_mets-n_met_mnx_kegg_v1);...
    n_met_mnx_metacyc_v1 (n_mets-n_met_mnx_metacyc_v1);...
    n_met_mnx_seed_v1 (n_mets-n_met_mnx_seed_v1)];
sample_data3 = [n_met_mnx_bigg_v2 (n_mets-n_met_mnx_bigg_v2);...
    n_met_mnx_kegg_v2 (n_mets-n_met_mnx_kegg_v2);...
    n_met_mnx_metacyc_v2 (n_mets-n_met_mnx_metacyc_v2);...
    n_met_mnx_seed_v2 (n_mets-n_met_mnx_seed_v2)];
sample_data4 = [n_met_bigg_v2_improved2 (n_mets-n_met_bigg_v2_improved2);...
    n_met_kegg_v2_improved2 (n_mets-n_met_kegg_v2_improved2);...
    n_met_metacyc_v2_improved2 (n_mets-n_met_metacyc_v2_improved2);...
    n_met_seed_v2_improved2 (n_mets-n_met_seed_v2_improved2)];
matrix(1,:,:) = sample_data1;
matrix(2,:,:) = sample_data2;
matrix(3,:,:) = sample_data3;
matrix(4,:,:) = sample_data4;

plotBarStackGroups(matrix,{'iter1','iter2','iter3','iter4'})
% A loop that does num2str conversion only if value is >0
for k = 1:size(matrix,1)
    for i=1:size(matrix,2);
        for j=1:size(matrix,3);
            if j ==1; c = 'w'; else c = 'k'; end;
            labels_stacked=num2str(matrix(k,i,j),'%.0f');
            hText = text((k-1)+0.75+(0.165*(i-1)), sum(matrix(k,i,1:j),3), labels_stacked); % for 4 samples
            set(hText, 'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c);
            hText2 = text((k-1)+0.75+(0.165*(i-1))+(length(databases{i})-4)*0.01, n_mets*1.04+(length(databases{i})-4)*n_mets*0.004, databases{i}); % for 4 samples
            set(hText2,'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c,'Rotation',45);
        end
    end
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, ['general_iterations_mets_' species '.pdf']);


%% plot local progress
matrix = zeros(4,4,2);
n_rxns = length(userModel_v2.rxns);
sample_data1 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_v2 (n_rxns-n_rxn_mnx_kegg_v2);...
    n_rxn_mnx_metacyc_v2 (n_rxns-n_rxn_mnx_metacyc_v2);...
    n_rxn_mnx_seed_v2 (n_rxns-n_rxn_mnx_seed_v2)];
sample_data2 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_improved (n_rxns-n_rxn_mnx_kegg_improved);...
    n_rxn_mnx_metacyc_improved (n_rxns-n_rxn_mnx_metacyc_improved);...
    n_rxn_mnx_seed_improved (n_rxns-n_rxn_mnx_seed_improved)];
sample_data3 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_improved2 (n_rxns-n_rxn_mnx_kegg_improved2);...
    n_rxn_mnx_metacyc_improved2 (n_rxns-n_rxn_mnx_metacyc_improved2);...
    n_rxn_mnx_seed_improved2 (n_rxns-n_rxn_mnx_seed_improved2)];
sample_data4 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_improved3 (n_rxns-n_rxn_mnx_kegg_improved3);...
    n_rxn_mnx_metacyc_improved3 (n_rxns-n_rxn_mnx_metacyc_improved3);...
    n_rxn_mnx_seed_improved3 (n_rxns-n_rxn_mnx_seed_improved3)];
matrix(1,:,:) = sample_data1;
matrix(2,:,:) = sample_data2;
matrix(3,:,:) = sample_data3;
matrix(4,:,:) = sample_data4;

plotBarStackGroups(matrix,{'iter1','iter2','iter3','iter4'})
% A loop that does num2str conversion only if value is >0
for k = 1:size(matrix,1)
    for i=1:size(matrix,2);
        for j=1:size(matrix,3);
            if j ==1; c = 'w'; else c = 'k'; end;
            labels_stacked=num2str(matrix(k,i,j),'%.0f');
            hText = text((k-1)+0.75+(0.165*(i-1)), sum(matrix(k,i,1:j),3), labels_stacked); % for 4 samples
            set(hText, 'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c);
            hText2 = text((k-1)+0.75+(0.165*(i-1))+(length(databases{i})-4)*0.01, n_rxns*1.04+(length(databases{i})-4)*n_rxns*0.004, databases{i});
            set(hText2,'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c,'Rotation',45);
        end
    end
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, ['local_iterations_rxns_' species '.pdf']);


%% plot general progress rxns
matrix = zeros(4,4,2);
n_rxns = length(userModel_v2.rxns);
sample_data1 = [n_rxn_mnx_bigg_original (n_rxns-n_rxn_mnx_bigg_original);...
    n_rxn_mnx_kegg_original (n_rxns-n_rxn_mnx_kegg_original);...
    n_rxn_mnx_metacyc_original (n_rxns-n_rxn_mnx_metacyc_original);...
    n_rxn_mnx_seed_original (n_rxns-n_rxn_mnx_seed_original)];
sample_data2 = [n_rxn_mnx_bigg_v1 (n_rxns-n_rxn_mnx_bigg_v1);...
    n_rxn_mnx_kegg_v1 (n_rxns-n_rxn_mnx_kegg_v1);...
    n_rxn_mnx_metacyc_v1 (n_rxns-n_rxn_mnx_metacyc_v1);...
    n_rxn_mnx_seed_v1 (n_rxns-n_rxn_mnx_seed_v1)];
sample_data3 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_v2 (n_rxns-n_rxn_mnx_kegg_v2);...
    n_rxn_mnx_metacyc_v2 (n_rxns-n_rxn_mnx_metacyc_v2);...
    n_rxn_mnx_seed_v2 (n_rxns-n_rxn_mnx_seed_v2)];
sample_data4 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_improved2 (n_rxns-n_rxn_mnx_kegg_improved2);...
    n_rxn_mnx_metacyc_improved2 (n_rxns-n_rxn_mnx_metacyc_improved2);...
    n_rxn_mnx_seed_improved2 (n_rxns-n_rxn_mnx_seed_improved2)];
matrix(1,:,:) = sample_data1;
matrix(2,:,:) = sample_data2;
matrix(3,:,:) = sample_data3;
matrix(4,:,:) = sample_data4;

plotBarStackGroups(matrix,{'iter1','iter2','iter3','iter4'})
% A loop that does num2str conversion only if value is >0
for k = 1:size(matrix,1)
    for i=1:size(matrix,2);
        for j=1:size(matrix,3);
            if j ==1; c = 'w'; else c = 'k'; end;
            labels_stacked=num2str(matrix(k,i,j),'%.0f');
            hText = text((k-1)+0.75+(0.165*(i-1)), sum(matrix(k,i,1:j),3), labels_stacked); % for 4 samples
            set(hText, 'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c);
            hText2 = text((k-1)+0.75+(0.165*(i-1))+(length(databases{i})-4)*0.01, n_rxns*1.04+(length(databases{i})-4)*n_rxns*0.004, databases{i});
            set(hText2,'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c,'Rotation',45);
        end
    end
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, ['general_iterations_rxns_' species '.pdf']);

%% plot general progress rxns, just manual model

matrix = zeros(3,4,2);
n_rxns = length(userModel_v2.rxns);
sample_data1 = [n_rxn_mnx_bigg_original (n_rxns-n_rxn_mnx_bigg_original);...
    n_rxn_mnx_kegg_original (n_rxns-n_rxn_mnx_kegg_original);...
    n_rxn_mnx_metacyc_original (n_rxns-n_rxn_mnx_metacyc_original);...
    n_rxn_mnx_seed_original (n_rxns-n_rxn_mnx_seed_original)];
sample_data2 = [n_rxn_mnx_bigg_v1 (n_rxns-n_rxn_mnx_bigg_v1);...
    n_rxn_mnx_kegg_v1 (n_rxns-n_rxn_mnx_kegg_v1);...
    n_rxn_mnx_metacyc_v1 (n_rxns-n_rxn_mnx_metacyc_v1);...
    n_rxn_mnx_seed_v1 (n_rxns-n_rxn_mnx_seed_v1)];
sample_data3 = [n_rxn_mnx_bigg_v2 (n_rxns-n_rxn_mnx_bigg_v2);...
    n_rxn_mnx_kegg_v2 (n_rxns-n_rxn_mnx_kegg_v2);...
    n_rxn_mnx_metacyc_v2 (n_rxns-n_rxn_mnx_metacyc_v2);...
    n_rxn_mnx_seed_v2 (n_rxns-n_rxn_mnx_seed_v2)];
matrix(1,:,:) = sample_data1;
matrix(2,:,:) = sample_data2;
matrix(3,:,:) = sample_data3;

plotBarStackGroups(matrix,{'original','partially mapped','fully mapped'})
% A loop that does num2str conversion only if value is >0
for k = 1:size(matrix,1)
    for i=1:size(matrix,2);
        for j=1:size(matrix,3);
            if j ==1; c = 'w'; else c = 'k'; end;
            labels_stacked=num2str(matrix(k,i,j),'%.0f');
            hText = text((k-1)+0.75+(0.165*(i-1)), sum(matrix(k,i,1:j),3), labels_stacked); % for 4 samples
            set(hText, 'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c);
            hText2 = text((k-1)+0.75+(0.165*(i-1))+(length(databases{i})-4)*0.01, n_rxns*1.04+(length(databases{i})-4)*n_rxns*0.004, databases{i});
            set(hText2,'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',10, 'Color',c,'Rotation',45);
        end
    end
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [-1 1 30 19])
saveas(gcf, ['general_iterations_manualModel_rxns_' species '.pdf']);

% plotBarStackGroups(matrix,{'iter1','iter2','iter3'})
% % A loop that does num2str conversion only if value is >0
% for k = 1:size(matrix,1)
%     for i=1:size(matrix,2);
%         for j=1:size(matrix,3);
%             if j ==1; c = 'w'; else c = 'k'; end;
%             labels_stacked=num2str(matrix(k,i,j),'%.0f');
%             hText = text((k-1)+0.75+(0.165*(i-1)), sum(matrix(k,i,1:j),3), labels_stacked);
%             set(hText, 'VerticalAlignment','top','HorizontalAlignment', 'center','FontSize',8, 'Color',c);
%         end
%     end
% end
end