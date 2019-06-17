function model = refineAuReMeModel(model, bigg)
if allMetsInCOBRAFormat(model.mets)
    model = transformModelToCBMPYFormat(model);
end
model = addFormulasAndChargesFromBigg(model, bigg);
end