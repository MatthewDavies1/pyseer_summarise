#NOTE: this is written for bespoke use, it is not ready out-of-the-box.
#script to take pyseer annotated kmers output and group kmers into overlapping regions and either extract their annotation or assign a new annotation
#also summarises the statistics by gene or by cluster (continuous overlapping region of multiple genes)


library(tidyr)
library(IRanges)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("")
raw <- read.table(file = "annotated_kmers_nobadchisq.txt", sep = "\t", header = F) #output annotated kmers from pyseer. Remove bad-chisquared kmers first.


raw.long <- separate_longer_delim(raw, cols = 8, delim = ",") %>%
  filter(V8 != "")
raw.wide <- separate_wider_delim(raw.long, cols = 8, delim = ":", 
                                 too_few = "align_start", names_sep = "_") 
raw.wide <- separate_wider_delim(raw.wide, cols = 9, delim = ";", 
                                 too_few = "align_start", names = c("range","upstream","gene","downstream")) 

#group kmers by overlapping regions per chromosome. Create a list of genes of interest with its total region to be input into 'script summarise kmer stats by gene.R'
df <- raw.wide %>%
  separate(range, into = c("start", "end"), sep = "-", convert = T)

df.clusters <- df %>%
  group_by(V8_1) %>%
  group_modify(~ {
    ir <- IRanges(start = .x$start, end = .x$end)
    reduced <- reduce(ir)
    hits <- findOverlaps(ir, reduced)
    .x$cluster <- subjectHits(hits)
    .x
  }) %>%
  ungroup()

df.clusters <- df.clusters %>%
  mutate(gene = ifelse(
    gene == "", paste0(V8_1, "_intergenic", cluster),
    gene
  ))

goi <- df.clusters %>%
  group_by(contig = V8_1, cluster) %>%
  summarise(cluster_start = min(start),
            cluster_end = max(end),
            genes = paste(unique(gene), collapse = ";"),
            .groups = "drop")

goi <- goi %>%
  mutate(genes_unique = make.unique(genes, sep = "_"))


# #temp extract emm gene kmers
# emm <-df.clusters[df.clusters$gene =='emm',]
# write.csv(emm, file = "emm_kmers.csv")

####script summarise kmer stats by gene
colnames(df.clusters)[1:8] <- c("contig","kmer","af","filter.pvalue","lrt.pvalue","beta","beta.std.err","variant_h2")
df.clusters$maf <- ifelse(df.clusters$af > 0.5, 1 - df.clusters$af, df.clusters$af)

library(GenomicRanges)
library(dplyr)

# query <- IRanges(goi$extract_start, goi$extract_end)
# subject <- IRanges(raw.wide$start, raw.wide$end)
# 
# overlaps <- findOverlaps(query, subject)

goi.gr <- GRanges(
  seqnames = goi$contig,
  ranges = IRanges(start = goi$cluster_start,
                   end = goi$cluster_end),
  ID = goi$genes,
  cluster = goi$cluster
)

# cluster.gr <- GRanges(
#   seqnames = goi$contig,
#   ranges = IRanges(start = goi$cluster_start,
#                    end = goi$cluster_end),
#   ID = goi$genes_unique,
#   cluster = goi$cluster
# )

kmers.gr <- GRanges(
  seqnames = df.clusters$contig,
  kmernames = df.clusters$kmer,
  ranges = IRanges(start = df.clusters$start,
                   end = df.clusters$end),
  af = df.clusters$af,
  lrt.pvalue = df.clusters$lrt.pvalue,
  beta = df.clusters$beta,
  maf = df.clusters$maf,
  cluster = df.clusters$cluster
)

overlaps <- findOverlaps(kmers.gr, goi.gr, type = "any")
kmer.idx <- queryHits(overlaps)
gene.idx <- subjectHits(overlaps)

results.df<- tibble(
  gene_id = goi.gr$ID[gene.idx],
  kmer_id = df.clusters$kmer[kmer.idx],
  lrt.pvalue = kmers.gr$lrt.pvalue[kmer.idx],
  af = kmers.gr$af[kmer.idx],
  beta = kmers.gr$beta[kmer.idx],
  maf = kmers.gr$maf[kmer.idx],
  cluster = df.clusters$cluster[kmer.idx],
  contig = df.clusters$contig[kmer.idx]
) 

results.df2 <- results.df %>%
  group_by(gene_id) %>%
  summarise(hits = n(),
            max_lrt_pvalue = -log10(max(lrt.pvalue)),
            mean_af = mean(af),
            mean_maf = mean(maf),
            mean_beta = mean(beta))

results.cluster <- results.df %>%
  group_by(contig, cluster) %>%
  summarise(genes = paste(unique(gene_id), collapse = ";"),
            hits = n(),
            max_lrt_pvalue = -log10(max(lrt.pvalue)),
            mean_af = mean(af),
            mean_maf = mean(maf),
            mean_beta = mean(beta))

#unique(goi$gene..pyseer...me.annotation.[! goi$gene..pyseer...me.annotation. %in% unique(results.df2$gene_id)])


###############by gene
prelim_plot_genes <- ggplot(data = subset(results.df2, hits > 15), aes(x=mean_beta, y=max_lrt_pvalue, colour=mean_maf, size=hits, label=gene_id)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=30), show.legend = FALSE, colour='black', max.overlaps = Inf) +
  scale_size("Number of k-mers", range=c(1,10)) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive k-mers") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)")

ggsave(filename = "geom_point_kmers_genes.svg", plot = prelim_plot_genes, device = "svg",
       width = 7.5, height = 4.5, units = "in", dpi = 500)
# ###flip x axis with MAF
# ggplot(results.df2, aes(x=mean_maf, y=max_lrt_pvalue, colour=mean_beta, size=hits, label=gene_id)) +
#   geom_point(alpha=0.5) +
#   geom_text_repel(aes(size=10), show.legend = FALSE, colour='black', max.overlaps = Inf) +
#   scale_size("Number of k-mers", range=c(1,10)) +
#   scale_colour_gradient('Average effect size') +
#   theme_bw(base_size=14) +
#   ggtitle("Invasive k-mers") +
#   xlab("Average MAF") +
#   ylab("Maximum -log10(p-value)")
# 
#########by cluster
prelim_plot_clusters <- ggplot(results.cluster, aes(x=mean_beta, y=max_lrt_pvalue, colour=mean_maf, size=hits, label=genes)) +
  geom_point(alpha=0.5) +
  geom_text_repel(data = subset(results.cluster, max_lrt_pvalue > 5),
                  aes(size=10), show.legend = FALSE, colour='black', max.overlaps = Inf) +
  scale_size("Number of k-mers", range=c(1,10)) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive k-mers") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)")

ggsave(filename = "geom_point_kmers_clusters.svg", plot = prelim_plot_clusters, device = "svg",
       width = 7.5, height = 4.5, units = "in", dpi = 500)
# ###flip x axis with MAF
# ggplot(results.cluster, aes(x=mean_maf, y=max_lrt_pvalue, colour=mean_beta, size=hits, label=genes)) +
#   geom_point(alpha=0.5) +
#   geom_text_repel(aes(size=10), show.legend = FALSE, colour='black', max.overlaps = Inf) +
#   scale_size("Number of k-mers", range=c(1,10)) +
#   scale_colour_gradient('Average effect size') +
#   theme_bw(base_size=14) +
#   ggtitle("Invasive k-mers") +
#   xlab("Average MAF") +
#   ylab("Maximum -log10(p-value)")

#############################################################################
########summarise which isolates the kmers are present in, to overlay onto phylogenetic tree
# isolates <- read.csv(file = "list_isolates.csv", header = FALSE)
# isolates <- unlist(isolates, use.names = F)

kmer.extracts <- read.csv(file = "annotated_extracted_kmers.csv", 
                          header = FALSE, sep = "|")
kmer.extracts[] <- lapply(kmer.extracts, function(x) gsub("Spy:.", "Spy", x))
kmer.extracts[1] <- lapply(kmer.extracts[1], function(x) gsub(" ", "", x))

kmer.long <- separate_longer_delim(kmer.extracts, cols = 2, delim = " ") %>%
  filter(V2 != "")
colnames(kmer.long) <- c("kmer_id", "onsnummer")

resultsdf_with_kmers <- kmer.long %>%
  left_join(results.df %>% select(kmer_id, gene_id), by = "kmer_id", relationship = "many-to-many") %>%
  filter(!is.na(gene_id))

per.sample.kmers <-resultsdf_with_kmers %>%
  distinct(onsnummer, kmer_id, gene_id) %>%
  group_by(onsnummer, gene_id) %>%
  summarise(sample_kmers = n(), .groups = "drop") %>%
  left_join(resultsdf_with_kmers %>%
              distinct(kmer_id, gene_id) %>%
              group_by(gene_id) %>%
              summarise(total_kmers = n(), .groups = "drop"),
            by = "gene_id") %>%
  mutate(prop = sample_kmers / total_kmers,
         label = paste0(gene_id, " (", total_kmers, ")"))

df.sample.x.gene <- per.sample.kmers %>%
  select(onsnummer, prop, label) %>%
  pivot_wider(names_from = label,
              values_from = prop,
              values_fill = 0)

colnames(df.sample.x.gene) <- gsub(";", "|", colnames(df.sample.x.gene))

write.csv(df.sample.x.gene, file = "kmer_proportions_per_sample.csv", row.names = F)

# gene_onsnummers <- resultsdf_with_kmers %>%
#   group_by(gene_id) %>%
#   summarise(samples = list(unique(onsnummer)))
# 
# results.final <- results.df2 %>%
#   left_join(gene_onsnummers, by = "gene_id")





#####filter out singletons

# 
# ggplot(results.df2, aes(x=mean_beta, y=max_lrt_pvalue, colour=mean_maf, size=hits, label = gene_id)) +
#   geom_point(alpha=0.5) +
#   geom_text_repel(data = subset(results.df2, max_lrt_pvalue > 1),
#                   aes(size=30), show.legend = FALSE, colour='black', max.overlaps = Inf) +
#   scale_size("Number of k-mers", range=c(1,10)) +
#   scale_colour_gradient('Average MAF') +
#   theme_bw(base_size=14) +
#   ggtitle("Invasive k-mers") +
#   xlab("Average effect size") +
#   ylab("Maximum -log10(p-value)")
# ggplot(results.cluster, aes(x=mean_beta, y=max_lrt_pvalue, colour=mean_maf, size=hits, label=genes)) +
#   geom_point(alpha=0.5) +
#   geom_text_repel(data = subset(results.cluster, max_lrt_pvalue > 6.5),
#     aes(size=30), show.legend = FALSE, colour='black', max.overlaps = Inf) +
#   scale_size("Number of k-mers", range=c(1,10)) +
#   scale_colour_gradient('Average MAF') +
#   theme_bw(base_size=14) +
#   ggtitle("Invasive k-mers") +
#   xlab("Average effect size") +
#   ylab("Maximum -log10(p-value)")
# 
# 
# ggplot(results.df2, aes(x=mean_maf, y=max_lrt_pvalue, colour=mean_beta, size=hits, label=gene_id)) +
#   geom_point(alpha=0.5) +
#   geom_text_repel(data = subset(results.df2, max_lrt_pvalue > 6.5),
#                   aes(size=50), show.legend = FALSE, colour='black', max.overlaps = Inf) +
#   scale_size("Number of k-mers", range=c(1,10)) +
#   scale_colour_gradient('Average effect size') +
#   theme_bw(base_size=14) +
#   ggtitle("Invasive k-mers") +
#   xlab("Average MAF") +
#   ylab("Maximum -log10(p-value)")
# ggplot(results.cluster, aes(x=mean_maf, y=max_lrt_pvalue, colour=mean_beta, size=hits, label=genes)) +
#   geom_point(alpha=0.5) +
#   geom_text_repel(data = subset(results.cluster, max_lrt_pvalue > 6.5),
#                   aes(size=50), show.legend = FALSE, colour='black', max.overlaps = Inf) +
#   scale_size("Number of k-mers", range=c(1,10)) +
#   scale_colour_gradient('Average effect size') +
#   theme_bw(base_size=14) +
#   ggtitle("Invasive k-mers") +
#   xlab("Average MAF") +
#   ylab("Maximum -log10(p-value)")
