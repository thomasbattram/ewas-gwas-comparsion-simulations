# -------------------------------------------------------
# Simulating EWAS and GWAS pathway overlap
# -------------------------------------------------------

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

source("R/simulation_functions.R")

args <- commandArgs(trailingOnly = TRUE)
split <- args[1]
split <- as.numeric(unlist(strsplit(split, ",")))
message("split = ", split)

# PATH TO FILES: /projects/MRC-IEU/research/projects/ieu2/p1/058/working/data/

# gene ontology terms - table of terms + genes linked to the terms
# ensembl_gene_id pathway_id description
go_terms <- read_tsv("gene_ontology_terms.txt") %>% 
	filter(!is.na(pathway_id))
# same as above but for kegg
kegg_terms <- read_tsv("kegg_terms.txt") 

# gene names + positions
# ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
all_genes <- read_tsv("ensembl_genes.txt", guess_max = 1e6)
all_genes$size <- all_genes$end_position - all_genes$start_position
summary(all_genes$size)

# empirical data
trait_res <- read_tsv("fishers_test_empirical_overlap_res.txt")

# gwas n genes data
gwas_g_dat <- read_tsv("n_gwas_genes.txt")

# unique pathways
u_go <- unique(go_terms$pathway_id)
u_kegg <- unique(kegg_terms$pathway_id)

# -------------------------------------------------------
# Setup parameters for simulations
# -------------------------------------------------------

# combine pathway datasets
all_path_terms <- bind_rows(list(go = go_terms, kegg = kegg_terms), .id = "database")

trait_gene_dat <- trait_res %>%
	left_join(gwas_g_dat) %>%
	dplyr::select(trait, n_ewas_genes, n_gwas_genes) %>%
	distinct() 

max(trait_gene_dat$n_ewas_genes) / nrow(all_genes)
max(trait_gene_dat$n_gwas_genes) / nrow(all_genes)
# minimum prop = 0.04

params <- expand.grid(
	trait = unique(trait_res$trait), 
	prop_causal_genes = c(0.05, 0.075, 0.1, 0.125, 0.15),
	prop_consequent_genes = c(0.05, 0.075, 0.1, 0.125, 0.15),
	percent_ewas_causal = c(0, 0.01, 0.05, 0.1, 0.5, 1), 
	sim = c(1:1000),
	gene_overlap = NA, 
	go_overlap = NA, 
	kegg_overlap = NA, 
	or_g = NA,
	or_go = NA, 
	or_kegg = NA, 
	p_g = NA, 
	p_go = NA, 
	p_kegg = NA
	) %>%
	left_join(trait_gene_dat)


# -------------------------------------------------------
# run the simulations
# -------------------------------------------------------
params <- params %>%
	dplyr::filter(between(sim, split[1], split[2]))

print("starting simulations")

out <- apply_sim(params)

out <- bind_rows(out)

make_dir <- function(path) {
	system(paste("mkdir", path))
}

path <- "sim_temp/"
make_dir(path)
nam <- paste0("sim_", split[2], ".txt")
write.table(out, file = paste0(path, nam),
            row.names = F, col.names = T, quote = F, sep = "\t")

if (split[2] == max(params$sim)) {
    sim_files <- list.files(path)
    sim_res <- map_dfr(sim_files, function(fil) {
            x <- read_tsv(paste0(path, fil))
            return(x)
    })
    nam <- paste0("sim_res.txt")
    write.table(sim_res, file = nam,
                            row.names = F, col.names = T, quote = F, sep = "\t")
}