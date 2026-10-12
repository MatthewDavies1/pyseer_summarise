# pyseer_summarise
Expands upon pyseer annotated kmer clustering and statistics.

Takes pyseer's annotated kmer output and groups kmers into overlapping regions that now include intergenic kmers (instead of clustering by annotation alone). Extracts annotation when possible or assigns new (intergenic) annotations to unannotated regions. Improves pyseer's analysis of overlapping genes present on different reading frames by treating each individual region as separate 'genes'.

Pyseer's built-in statistical approach is applied to each overlapping region of kmers (now including upstream and downstream intergenic regions) but without grouping overlapping genes into one.
