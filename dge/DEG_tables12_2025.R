## Tables for A. Platt (12/2025)

library(tidyverse)
(base_path <- path.expand("~/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/rna-seq-host-virus/data_and_res_gh_repo/"))

# 
# 1. HHQ_inf vs. DMSO_inf, T1 (genes can be taken from "deseq_DEGSets_T1_76samp_host", Set A)
# 2. HHQ_inf vs. DMSO_inf, T2 (genes can be taken from "deseq_DEGSets_T2_76samp_host", Set A)
# 3. HHQ_inf vs. DMSO_inf, T3 (genes can be taken from "deseq_DEGSets_T3_76samp_host", Set A)
# 4. HHQ_inf vs. DMSO_inf, T4 (genes can be taken from "deseq_DEGSets_T4_76samp_host", Set A)
# 5. HHQ_inf vs. DMSO_inf, all time points together (lists 1-4 combined)
# 
# 6. A. (HHQ_inf vs. DMSO_inf, T1 without genes overlapping with HHQ_cntl vs. DMSO_cntl a.k.a Set C T1) (genes can be taken from "deseq_DEGSets_T1_76samp_host", Set A with Set B genes taken out)
# 7. A. HHQ_inf vs. DMSO_inf, T2 without genes overlapping with HHQ_cntl vs. DMSO_cntl a.k.a Set C T2) (genes can be taken from "deseq_DEGSets_T2_76samp_host", Set A with Set B genes taken out))
# 8. HHQ_inf vs. DMSO_inf, T3 without genes overlapping with HHQ_cntl vs. DMSO_cntl a.k.a Set C T3) (genes can be taken from "deseq_DEGSets_T3_76samp_host", Set A with Set B gene taken out)
# 9. HHQ_inf vs. DMSO_inf, T4 without genes overlapping with HHQ_cntl vs. DMSO_cntl a.k.a Set C) (genes can be taken from "deseq_DEGSets_T4_76samp_host", Set A with Set B genes taken out)
# 10. HHQ_inf vs. DMSO_inf, all time points together (lists 6-9 combined)
# 
# 11. HHQ_cntl vs. DMSO_cntl, T1 (genes can be taken from "deseq_DEGSets_T1_76samp_host", Set C)
# 12. HHQ_cntl vs. DMSO_cntl, T2 (genes can be taken from "deseq_DEGSets_T2_76samp_host", Set C)
# 13. HHQ_cntl vs. DMSO_cntl, T3 (genes can be taken from "deseq_DEGSets_T3_76samp_host", Set C)
# 14. HHQ_cntl vs. DMSO_cntl, T4 (genes can be taken from "deseq_DEGSets_T4_76samp_host", Set C)
# 15. HHQ_cntl vs. DMSO_cntl, all time points (lists 11-14 combined)
# 
# 16. HHQ_inf vs. HHQ_cntl, all time points, (genes can be taken from deseq_DEGSets_T1-4_76samp_host set D at all time points)
# 17. DMSO_inf vs. DMSO_cntl, all time points, (genes can be taken from deseq_DEGSets_T1-4_76samp_host set F at all time points)

#So it will be the genes from above (1-17) with their respective Log2FCshrnk_ashr and padj IHW for that comparison and time point, which I have been taking from "deseq_counts_contrastsOfInterest_76samp_host.xlxs". From "annot-host-ext.xlxs" I think we only need information from the columns listed below, but it can all go in if that's easier. 
#"annot-host-ext.xlxs"
#tab "GO+KEGG"  columns N,O,P,Q, and R
#tab "Harriet" column I
#tab "Pollara_SOM" columns O,P, Q, and R
#tab "EggNogg-Annotations" column J


