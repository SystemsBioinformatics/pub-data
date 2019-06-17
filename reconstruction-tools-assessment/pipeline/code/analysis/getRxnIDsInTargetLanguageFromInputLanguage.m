function [rxnIDs_outputLanguage, MNXIDs] = getRxnIDsInTargetLanguageFromInputLanguage(rxnIDs_inputLanguage,inputLanguage,outputLanguage, rxnOtherIDs, rxnMNXIDs)

% load('D:\Dropbox\Databases\MNX\rxnLanguages');
% load('D:\Dropbox\Databases\MNX\rxnOtherIDs');
% load('D:\Dropbox\Databases\MNX\rxnMNXIDs');

rxnIDs_outputLanguage = cell(size(rxnIDs_inputLanguage));
MNXIDs = cell(size(rxnIDs_inputLanguage));

if strcmp(inputLanguage,'MNX')
    for i = 1:length(rxnIDs_inputLanguage)
        MNXID = rxnIDs_inputLanguage{i};
        pos2 = find(strcmp(rxnMNXIDs,MNXID));
        otherIDs = rxnOtherIDs(pos2);
        
        posOutputID = find(~cellfun(@isempty, strfind(otherIDs,outputLanguage)));
        if ~isempty(posOutputID)
            rxnIDs_outputLanguage(i) = {regexprep(otherIDs(posOutputID),[outputLanguage ':'],'')};
        end
    end
elseif strcmp(outputLanguage,'MNX')
    for i = 1:length(rxnIDs_inputLanguage)
        pos = find(strcmp(rxnOtherIDs,[inputLanguage ':' rxnIDs_inputLanguage{i}]));
        if ~isempty(pos)
            rxnIDs_outputLanguage(i) = rxnMNXIDs(pos);
        end
    end
else
    
    for i = 1:length(rxnIDs_inputLanguage)
        pos = find(strcmp(rxnOtherIDs,[inputLanguage ':' rxnIDs_inputLanguage{i}]));
        if ~isempty(pos)
            MNXID = rxnMNXIDs(pos);
            pos2 = find(strcmp(rxnMNXIDs,MNXID));
            MNXIDs(i) = MNXID;
            
            otherIDs = rxnOtherIDs(pos2);
            
            posOutputID = find(~cellfun(@isempty, strfind(otherIDs,outputLanguage)));
            if ~isempty(posOutputID)
                rxnIDs_outputLanguage(i) = {regexprep(otherIDs(posOutputID),[outputLanguage ':'],'')};
            end            
        end
    end    
end

end