function JC_sim_r = getJaccardSimilarityRxns(model1, model2)

rxns1 = model1.rxns;
rxns2 = model2.rxns;
n_intersection = length(intersect(rxns1,rxns2));
n_union = length(union(rxns1,rxns2));
JC_sim_r = 1 - n_intersection/n_union;

end