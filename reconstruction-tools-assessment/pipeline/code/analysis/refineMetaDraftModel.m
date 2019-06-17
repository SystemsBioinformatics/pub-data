function model = refineMetaDraftModel(model, bigg)

if allMetsInCOBRAFormat(model.mets)
   model = transformModelToCBMPYFormat(model); 
end
% model.mets = regexprep(model.mets,{'(.*)\[(.*)\]$','Cytosol','Extra_organism','C_c','C_e'},{'$1_$2','c','e','c','e'}); % metadraft v.0.8.1
% the following commented code is for MetaDraft v.0.8.1
% for i = 1:length(model.genes) 
%     pos = strfind(model.genes{i},'_');
%     if ~isempty(pos)
%         pos = pos(pos>6);
%         if ~isempty(pos)
%             pos = pos(1);
%             model.genes{i} =  model.genes{i}(1:pos-1);
%         end
%     end
% end

model = addChargesFromBigg(model, bigg);

end