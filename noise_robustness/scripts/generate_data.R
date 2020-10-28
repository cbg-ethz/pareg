library(tidyverse)


# parameters
fname_in <- snakemake@input$fname

fname_geneset <- snakemake@output$fname_geneset
fname_termdatabase <- snakemake@output$fname_termdatabase

# read data
df_all <- read_csv(fname_in)
df_all %>%
  head

# compute geneset sizes
df_groupsizes <- df_all %>%
  group_by(gs_name) %>%
  summarize(n = n())

df_groupsizes %>%
  ggplot(aes(x = n)) +
    geom_histogram() +
    xlab("Geneset size") +
    scale_y_log10() +
    theme_minimal()

# filter for genesets of reasonable size
min_size <- 200
max_size <- 300

df_sub <- df_groupsizes %>%
  filter(min_size <= n & n <= max_size) %>%
  inner_join(df_all, by = "gs_name") %>%
  select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol)

df_sub %>%
  head

# select enriched terms
enriched_terms <- df_sub %>%
  distinct(term) %>%
  sample_n(5) %>%
  pull(term)

# create study geneset
geneset <- df_sub %>%
  filter(term %in% enriched_terms) %>%
  distinct(gene) %>%
  pull(gene)

# identify background geneset
all_genes <- df_sub %>%
  distinct(gene) %>%
  pull(gene)

background_genes <- setdiff(all_genes, geneset)

# add noise
fp_rate <- 0.1 # alpha
fn_rate <- 0.4 # beta

selection_mask_fn <- sample(c(FALSE, TRUE), length(geneset), replace = TRUE, prob = c(1 - fn_rate, fn_rate)) # remove from selected geneset
selection_mask_fp <- sample(c(FALSE, TRUE), length(background_genes), replace = TRUE, prob = c(1 - fp_rate, fp_rate)) # add from background set

noisy_geneset <- c(geneset[selection_mask_fn], background_genes[selection_mask_fp])

# save results
tibble(
  gene = c(noisy_geneset, all_genes),
  type = c(rep("study", length(noisy_geneset)), rep("all", length(all_genes)))
) %>%
  write_csv(fname_geneset)

df_sub %>%
  mutate(term_is_enriched = term %in% enriched_terms) %>%
  write_csv(fname_termdatabase)
