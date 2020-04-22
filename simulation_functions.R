# enrichment function
perform_enrichment <- function(ewas_vars, gwas_vars, total_vars) 
{
  q <- sum(ewas_vars %in% gwas_vars) # overlap between ewas and gwas pathways
  m <- length(gwas_vars) - q # number of pathways identified by GWAS but not EWAS
  k <- length(ewas_vars) - q # number of pathways identified by EWAS but not by GWAS
  n <- length(total_vars) - q - m - k # all pathways minus pathways identified by EWAS and GWAS
  tab <- matrix(c(q, m, k, n), 2, 2)
  out <- fisher.test(tab, alternative = "greater")
  return(out)
}

sim <- function(prop_causal_genes, n_ewas_genes, n_gwas_genes, percent_ewas_causal, prop_consequent_genes)
{
	# 1. sample causal and consequential genes from total number of genes
	# 2. sample gwas and ewas genes and do gene overlap tests
	# 3. link these to pathways and do pathway overlap tests
	# 4. write out results

	# 1.
	ca_genes <- sample_n(all_genes, prop_causal_genes * nrow(all_genes))
	con_genes <- sample_n(all_genes, prop_consequent_genes * nrow(all_genes))

	# 2. 
	gwasg <- sample(ca_genes$ensembl_gene_id, n_gwas_genes)
	ewas_cag <- sample(ca_genes$ensembl_gene_id, n_ewas_genes * percent_ewas_causal)
	ewasg <- unique(c(ewas_cag, sample(con_genes$ensembl_gene_id, n_ewas_genes - length(ewas_cag))))

	outg <- perform_enrichment(ewasg, gwasg, all_genes$ensembl_gene_id)

	# 3.
	database <- c("go", "kegg") 
	outp <- lapply(1:2, function(x) {
		pathway_dat <- get(paste0(database[x], "_terms"))
		gwasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% gwasg) %>%
			pull(pathway_id) %>%
			unique
		ewasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% ewasg) %>%
			pull(pathway_id) %>%
			unique
		out_res <- perform_enrichment(gwasp, ewasp, unique(pathway_dat$pathway_id))
		overlap <- sum(ewasp %in% gwasp)
		return(list(res = out_res, overlap = overlap))
	})
	names(outp) <- database
	
	# 4.
	out <- list(gene = outg, 
				gene_overlap = sum(ewasg %in% gwasg),
				go = outp$go$res,  
				go_overlap = outp$go$overlap, 
				kegg = outp$kegg$res, 
				kegg_overlap = outp$kegg$overlap)
	return(out)
}

apply_sim <- function(params) {
	out <- lapply(split(params, 1:nrow(params)), function(x) {
		out <- sim3(x$prop_causal_genes, x$n_ewas_genes, x$n_gwas_genes, 
					x$percent_ewas_causal, x$prop_consequent_genes)
		x$gene_overlap <- out$gene_overlap
		x$go_overlap <- out$gene_overlap
		x$kegg_overlap <- out$kegg_overlap
		x$or_g <- out$gene$estimate
		x$or_go <- out$go$estimate
		x$or_kegg <- out$kegg$estimate
		x$p_g <- out$gene$p.value
		x$p_go <- out$go$p.value 
		x$p_kegg <- out$kegg$p.value
		return(x)
	})
	return(out)
}

