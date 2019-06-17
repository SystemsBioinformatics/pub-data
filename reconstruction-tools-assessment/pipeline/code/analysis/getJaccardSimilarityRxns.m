function JC_sim_r = getJaccardSimilarityRxns(model1, model2, option)

[n_intersection, n_union] = calculateRxnIntersectionAndUnion(model1, model2, option);
JC_sim_r = 1 - n_intersection/n_union;

end