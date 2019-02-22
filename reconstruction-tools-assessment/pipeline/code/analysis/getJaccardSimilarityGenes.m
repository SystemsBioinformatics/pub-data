function JC_sim_g = getJaccardSimilarityGenes(model1, model2)

genes1 = model1.genes;
genes2 = model2.genes;
n_intersection = length(intersect(genes1,genes2));
n_union = length(union(genes1,genes2));
JC_sim_g = 1 - n_intersection/n_union;

end