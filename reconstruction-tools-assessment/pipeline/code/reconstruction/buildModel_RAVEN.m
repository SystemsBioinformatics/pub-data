function buildModel_RAVEN(species)

baseFolder = 'D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\RA\RAVEN_runs';
inputsBaseFolder = 'D:\Dropbox\Research_Projects\Review_reconstruction\reconstructions\Inputs';
inputFile = [inputsBaseFolder filesep upper(species) filesep 'protein_fasta.faa'];
n_run = 0;


%% runs with default parameters

% organismID = species
% fastaFile = inputFile
% dataDir = prok100_kegg82
% outDir = 'D:\Dropbox\Review_reconstruction\reconstructions\RA\output1
% keepUndefinedStoich not supplied (default = true)
% keepIncomplete not supplied (default = true)
% keepGeneral not supplied (default)
% cutOff not supplied (default 10^-50)
% minScoreRatioG not supplied (default)
% minScoreRatioKO not supplied (default)
% maxPhylDist not supplied (default)
% nSequences not supplied (default = inf)
% seqIdentity not supplied (default = -1)

n_run = n_run+1; %1
model=getKEGGModelForOrganism(species,inputFile,'prok100_kegg82'...
    ,[species '_run' numstr(n_run)]);
cd([baseFolder filesep species '_run' numstr(n_run)])
save([species '_' numstr(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' numstr(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %2
model=getKEGGModelForOrganism(species,inputFile,'prok90_kegg82'...
    ,[species '_run' numstr(n_run)]);
cd([baseFolder filesep species '_run' numstr(n_run)])
save([species '_' numstr(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' numstr(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %3
model=getKEGGModelForOrganism(species,inputFile,'prok50_kegg82'...
    ,[species '_run' numstr(n_run)]);
cd([baseFolder filesep species '_run' numstr(n_run)])
save([species '_' numstr(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' numstr(n_run) '_ra.xml'])
cd(baseFolder)

%% runs with default tutorial parameters

% organismID = species
% fastaFile = inputFile
% dataDir = prok100_kegg82
% outDir = 'D:\Dropbox\Review_reconstruction\reconstructions\RA\output1
% keepUndefinedStoich = false (default = true)
% keepIncomplete = false (default = true)
% keepGeneral = false (default)
% cutOff = 10^-30 (default 10^-50)
% minScoreRatioG = 0.8 (default)
% minScoreRatioKO = 0.3 (default)
% maxPhylDist = -1 (default)
% nSequences not supplied (default = inf)
% seqIdentity not supplied (default = -1)

n_run = n_run+1; %4
model=getKEGGModelForOrganism(species,inputFile,'prok100_kegg82'...
    ,[species '_run' numstr(n_run)],false,false,false,10^-30,0.8,0.3,-1);
cd([baseFolder filesep species '_run' numstr(n_run)])
save([species '_' numstr(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' numstr(n_run) '_ra.xml'])
cd(baseFolder)


n_run = n_run+1; %5
model=getKEGGModelForOrganism(species,inputFile,'prok90_kegg82'...
    ,[species '_run' numstr(n_run)],false,false,false,10^-30,0.8,0.3,-1);
cd([baseFolder filesep species '_run' numstr(n_run)])
save([species '_' numstr(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' numstr(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %6
model=getKEGGModelForOrganism(species,inputFile,'prok50_kegg82'...
    ,[species '_run' numstr(n_run)],false,false,false,10^-30,0.8,0.3,-1);
cd([baseFolder filesep species '_run' numstr(n_run)])
save([species '_' numstr(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' numstr(n_run) '_ra.xml'])
cd(baseFolder)

%% Metacyc
n_run = n_run+1; %7
if ~isdir([baseFolder filesep species '_run' num2str(n_run)])
    mkdir([baseFolder filesep species '_run' num2str(n_run)]);
end
cd([baseFolder filesep species '_run' num2str(n_run)]);
model = getMetaCycModelForOrganism([species num2str(n_run) 'MetaCyc'],inputFile,1);
save([species '_' num2str(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' num2str(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %8
if ~isdir([baseFolder filesep species '_run' num2str(n_run)])
    mkdir([baseFolder filesep species '_run' num2str(n_run)]);
end
cd([baseFolder filesep species '_run' num2str(n_run)]);
model = getMetaCycModelForOrganism([species num2str(n_run) 'MetaCyc'],inputFile,1,0,0,90,45);
save([species '_' num2str(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' num2str(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %9
if ~isdir([baseFolder filesep species '_run' num2str(n_run)])
    mkdir([baseFolder filesep species '_run' num2str(n_run)]);
end
cd([baseFolder filesep species '_run' num2str(n_run)]);
model = getMetaCycModelForOrganism([species num2str(n_run) 'MetaCyc'],inputFile,1,0,0,110,45);
save([species '_' num2str(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' num2str(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %10
if ~isdir([baseFolder filesep species '_run' num2str(n_run)])
    mkdir([baseFolder filesep species '_run' num2str(n_run)]);
end
cd([baseFolder filesep species '_run' num2str(n_run)]);
model = getMetaCycModelForOrganism([species num2str(n_run) 'MetaCyc'],inputFile,1,0,0,100,55);
save([species '_' num2str(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' num2str(n_run) '_ra.xml'])
cd(baseFolder)

n_run = n_run+1; %11
if ~isdir([baseFolder filesep species '_run' num2str(n_run)])
    mkdir([baseFolder filesep species '_run' num2str(n_run)]);
end
cd([baseFolder filesep species '_run' num2str(n_run)]);
model = getMetaCycModelForOrganism([species num2str(n_run) 'MetaCyc'],inputFile,1,0,0,100,35);
save([species '_' num2str(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' num2str(n_run) '_ra.xml'])
cd(baseFolder)

% Annotation KEGG
n_run = n_run+1; %12
if ~isdir([baseFolder filesep species '_run' num2str(n_run)])
    mkdir([baseFolder filesep species '_run' num2str(n_run)]);
end
cd([baseFolder filesep species '_run' num2str(n_run)]);
model = getKEGGModelForOrganism(species,'','',[baseFolder filesep species '_run' num2str(n_run)],0,0);
save([species '_' num2str(n_run) '_ra.mat'],'model')
writeCbModel(model, 'format','sbml', 'filename', [species '_' num2str(n_run) '_ra.xml'])

% n_run = n_run+1; %13
% modelIDs = {'iNF517'};
% load('D:\Dropbox\Databases\BIGG\iNF517.mat');
% models{1} = eval(modelIDs{1});
% models{1}.id = modelIDs{1};
% refFastaFiles = {'D:\Dropbox\Databases\AuremeInputs\iNF517\FAA_Model.faa'};
% blastStructure = getBlast(species,inputFile,modelIDs,refFastaFiles);
% draftModel = getModelFromHomology(models, blastStructure, species, 1);

cd(baseFolder)

end