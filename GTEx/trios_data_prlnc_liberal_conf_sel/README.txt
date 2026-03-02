=======================================================
MRGN Results for Liberal Confounder Selection (NO FDR)
=======================================================

This folder contains the results from applying MRGN to trios with liberally selected confounders for top5 tissues by sample size (GMAC tissues)
The confounder selection procedure is the GMAC regression based procedure but no FDR correction is applied to the selected confounders

significant confounders are selected using the alpha < 0.01 and alpha < 0.05 cutoff

folders:

/list_output/ - the entire output as returned by get.conf.trios.R saved as .RData file
/infer_trio_results/ - The entire output from analysis of each trio with its selected confs: format is 14 (vector returned by infer.trio) x Number of Trios "list" saved as .RData
/data_with_confs/ - a list of dataframes corresponding to each trio (list length = Number of trios) saved as .RData
