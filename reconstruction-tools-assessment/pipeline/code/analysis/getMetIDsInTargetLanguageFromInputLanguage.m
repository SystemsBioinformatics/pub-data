function metIDs_outputLanguage = getMetIDsInTargetLanguageFromInputLanguage(metIDs_inputLanguage,inputLanguage,outputLanguage, metOtherIDs, metMNXIDs)

global rootFolder

% load('D:\Dropbox\Databases\MNX\metLanguages');
if nargin < 4
    load(fullfile(rootFolder, 'MNX', 'metOtherIDs'));
    load(fullfile(rootFolder, 'MNX', 'metMNXIDs'));
end

metIDs_outputLanguage = cell(size(metIDs_inputLanguage));

if strcmp(inputLanguage,'MNX')
    for i = 1:length(metIDs_inputLanguage)
        MNXID = metIDs_inputLanguage{i};
        pos2 = find(strcmp(metMNXIDs,MNXID));
        otherIDs = metOtherIDs(pos2);
        
        posOutputID = find(~cellfun(@isempty, strfind(otherIDs,outputLanguage)));
        if ~isempty(posOutputID)
            metIDs_outputLanguage(i) = {regexprep(otherIDs(posOutputID),[outputLanguage ':'],'')};
        end
    end
    
elseif strcmp(outputLanguage,'MNX')
    for i = 1:length(metIDs_inputLanguage)
        pos = find(strcmp(metOtherIDs,[inputLanguage ':' metIDs_inputLanguage{i}]));
        if ~isempty(pos)
            metIDs_outputLanguage(i) = metMNXIDs(pos);
        end
    end
else
    for i = 1:length(metIDs_inputLanguage)
        pos = find(strcmp(metOtherIDs,[inputLanguage ':' metIDs_inputLanguage{i}]));
        if ~isempty(pos)
            MNXID = metMNXIDs(pos);
            pos2 = find(strcmp(metMNXIDs,MNXID));
            
            otherIDs = metOtherIDs(pos2);
            
            posOutputID = find(~cellfun(@isempty, strfind(otherIDs,outputLanguage)));
            if ~isempty(posOutputID)
                metIDs_outputLanguage(i) = {regexprep(otherIDs(posOutputID),[outputLanguage ':'],'')};
            end
            
        end
    end
end

for i = 1:length(metIDs_outputLanguage); if isempty(metIDs_outputLanguage{i}); metIDs_outputLanguage{i} = ''; end; end;