coef_result<-GetFoundercoefs('dataset.islet.rnaseq', 'ENSMUSG00000006732', '3')
mediate_results<-GetMediate('dataset.islet.rnaseq', 'ENSMUSG00000006732', '1_4533435')
expr_results<-GetExpression('dataset.islet.rnaseq', 'ENSMUSG00000006732')
snp_result<-GetSnpAssocMapping('dataset.islet.rnaseq', 'ENSMUSG00000006732', '10', 10000000)

coef_result<-GetFoundercoefs('dataset.clinical.phenotypes', 'partial_inflation', '3')
expr_results<-GetExpression('dataset.clinical.phenotypes', 'partial_inflation')
snp_result<-GetSnpAssocMapping('dataset.clinical.phenotypes', 'Glu_iAUC', '10', 10000000)
