function JC_sim_m = getJaccardSimilarityMets(model1, model2)

mets1 = model1.mets;
mets2 = model2.mets;
n_intersection = length(intersect(mets1,mets2));
n_union = length(union(mets1,mets2));
JC_sim_m = 1 - n_intersection/n_union;

end