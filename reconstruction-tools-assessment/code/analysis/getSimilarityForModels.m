function [JC_sim_r_all, JC_sim_m_all, JC_sim_g_all] = getSimilarityForModels(models, ...
    species, abbr, excelFileName)

JC_sim_r_all = getRxnSimilarityForModels(models, species, abbr, excelFileName);
JC_sim_m_all = getMetSimilarityForModels(models, species, abbr, excelFileName);
JC_sim_g_all = getGeneSimilarityForModels(models, species, abbr, excelFileName);


end