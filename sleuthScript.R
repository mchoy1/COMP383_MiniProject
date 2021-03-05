library(sleuth)
library(dplyr)
library(data.table)

#sleuth sample tables
stab <- read.table('sleuthInput.txt',header=TRUE, stringsAsFactors = FALSE,
                    sep='\t')

#initialize sleuth object
sleuth_object <- sleuth_prep(stab)

#test both conditions for differential expression
#model comparison of two conditions
sleuth_object <- sleuth_fit(~condition, 'full')

#reduced model
sleuth_object <- sleuth_fit(sleuth_object,~1,'reduced')

#likelihood ratio test for differential expression b/n conditions
sleuth_object <- sleuth_lrt(sleuth_object,'reduced','full')

#sleuth_object results
sleuth_table <- sleuth_results(sleuthObject, 'reduced:full','lrt',show_all=FALSE)

#define significance (FDR < 0.05) & fiter results
sleuth_significance <- dpplyr::filter(sleuth_table, qval <= 0.05) %>% dpylr:arrage(pval)

#select specified values
significant_sleuth <- sleught_significance %>% select(target_id, test_stat, pval, qval)

#write results to table with header and tab-delimit
write.table(significant_sleuth, file = 'sleuth_outfile.txt',quote=FALSE,row.names=FALSE)
