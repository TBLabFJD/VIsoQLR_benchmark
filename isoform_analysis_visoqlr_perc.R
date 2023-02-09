
library(bedr)
library(ggplot2)
library(ggforce)
library(stringr)
library(proxy)

rm(list=ls()) 


#===========#
# FUNCTIONS #
#===========#

transcript_data_loading <- function(input_path,  inputformat = "GTF", gene_paste = F){
  if (inputformat == "GTF") {separation = " "}
  if (inputformat == "GFF3") {separation = "="}
  
  gtf = read.delim(input_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  
  transcripts = gtf[gtf$V3 == "transcript","V9"]
  if (length(transcripts) == 0){transcripts = gtf[,"V9"]}
  transcriptinfo = gsub("; ", ";", transcripts)
  transcriptinfo = do.call("rbind", strsplit(x = transcriptinfo, split = ";"))
  cnames = do.call("rbind", strsplit(x = transcriptinfo[1,], split = separation))[,1]
  transcriptinfo = gsub("^.* ", "", transcriptinfo, perl = T)
  colnames(transcriptinfo) = cnames
  transcriptinfo = data.frame(transcriptinfo, stringsAsFactors = F)
  if(gene_paste){transcriptinfo$transcript_id = paste(transcriptinfo$gene_id, transcriptinfo$transcript_id, sep = "_")}
  transcriptinfo = transcriptinfo[!duplicated(transcriptinfo[,"transcript_id"]),]
  
  
  gtf = gtf[gtf$V3 == "exon", c("V1", "V4","V5","V9", "V6")]
  colnames(gtf) = c("gene", "start", "end", "id", "score")
  gtf$score = as.numeric(gtf$score)
  gtf$id=gsub("^.*transcript_id ", "", gtf$id, perl = TRUE)
  gtf$id=gsub(";.*$", "", gtf$id, perl = TRUE)
  if(gene_paste){gtf$id = paste(gtf$gene, gtf$id, sep = "_")}
  
  gtf = merge(gtf, transcriptinfo, by.x = "id", by.y = "transcript_id")
  
  
  # Size annotation
  gtf_split = split(gtf, gtf$id)
  gtf_lengths = unlist(lapply(gtf_split, function(x) sum(x$end - x$start + 1)))
  gtf$size = paste0(gtf_lengths[gtf$id], "bp")
  
  gtf = type.convert(gtf, as.is = T)
  return(gtf)
}

bed_data_loading <- function(input_path){
  bed6 = read.delim(input_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  colnames(bed6) = c("gene", "start", "end", "id", "score", "orientation")
  bed6$score = as.numeric(bed6$score)
  bed6$start = bed6$start + 1
  
  
  
  # Size annotation
  bed6_split = split(bed6, bed6$id)
  bed6_lengths = unlist(lapply(bed6_split, function(x) sum(x$end - x$start + 1)))
  bed6$size = paste0(bed6_lengths[bed6$id], "bp")
  
  return(bed6)
}







match_intersection_calculation <- function(gold, test){
  
  gold_split = split(gold, gold$id)
  
  test<- bedr(
    input = list(i = test), 
    method = "sort", 
    params = "",
    verbose = F,
    check.chr = FALSE
  )
  
  
  intersected_bed = do.call("rbind", lapply(gold_split, function(x) bedr(
    input = list(a = x, b = test), 
    method = "intersect", 
    params = "-wo",
    verbose = F,
    check.chr = FALSE
  )))
  
  colnames(intersected_bed) = c("gene", "start", "end", "id", "size", "gene.b", "start.b", "end.b", "id.b",     
                                "size.b", "V11" )
  
  
  
  intersected_bed_split = split(intersected_bed, paste(intersected_bed$id, intersected_bed$id.b))
  
  
  intersection_perc = do.call("rbind", lapply(intersected_bed_split, function(x){
    c(id = x$id[1],
      id2 = x$id.b[1],
      intersect = sum(as.numeric(x$V11))/max(as.numeric(gsub("bp", "", c(x$size, x$size.b))))
    )
  }))
  intersection_perc = data.frame(intersection_perc, stringsAsFactors = F)
  
  return(intersection_perc)
}





#==============#
# DATA LOADING #
#==============#


# GMAP
setwd("/home/gonzalo/tblab/mnt/tblab/gonzalo/VIsoQLR_article/uhrr/isoforms_SIRV_GMAP/")
viso_path25 = "/home/gonzalo/tblab/mnt/tblab/gonzalo/VIsoQLR_article/uhrr/isoforms_SIRV_GMAP/VIsoQLR_SIRV_0.25perc/visoqlr_results.gtf"
viso_path05 = "/home/gonzalo/tblab/mnt/tblab/gonzalo/VIsoQLR_article/uhrr/isoforms_SIRV_GMAP/VIsoQLR_SIRV_0.5perc/visoqlr_results.gtf"
viso_path1 = "/home/gonzalo/tblab/mnt/tblab/gonzalo/VIsoQLR_article/uhrr/isoforms_SIRV_GMAP/VIsoQLR_SIRV_1perc/visoqlr_results.gtf"
viso_path2 = "/home/gonzalo/tblab/mnt/tblab/gonzalo/VIsoQLR_article/uhrr/isoforms_SIRV_GMAP/VIsoQLR_SIRV_2perc/visoqlr_results.gtf"
viso_path3 = "/home/gonzalo/tblab/mnt/tblab/gonzalo/VIsoQLR_article/uhrr/isoforms_SIRV_GMAP/VIsoQLR_SIRV_3perc/visoqlr_results.gtf"

# Minimap2
# st2_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/stringtie_results_minimap/flnc.gtf"
# flair_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/flair_results_minimap/flnc.flair.collapse.isoforms.bed6"
# viso_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/visoqlr_results_minimap/visoqlr_results.gtf"
# flair_count_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/flair_results_minimap/flnc.flair.collapse.firstpass.q.counts"





gs_path = "/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/documentation/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"
gold = transcript_data_loading(gs_path)[,c("gene", "start", "end", "id", "size")]
gold$start = gold$start - 1



iso_intersection_viso = function(viso_path, perc_tag){
  viso = transcript_data_loading(viso_path, gene_paste = T)[,c("gene", "start", "end", "id", "size")]
  viso$start = viso$start - 1
  
  
  
  #==================#
  # RUN INTERSECTION #
  #==================#
  
  
  
  viso_intersection_perc = match_intersection_calculation(gold, viso)
  
  
  length(unique(viso$id))
  
  
  intersect_threshold = 0.99
  # intersect_threshold = 0.999
  # intersect_threshold = 0.9995
  
  
  # table(viso_intersection_perc$intersect>intersect_threshold)
  # length(unique(viso_intersection_perc[viso_intersection_perc$intersect>intersect_threshold,]$id))
  
  
  
  #================#
  # COUNTS LOADING #
  #================#
  
  
  # VISOQLR
  viso_count = transcript_data_loading(viso_path, gene_paste = T)[,c("id", "n_of_reads", "gene")]
  viso_count = viso_count[!duplicated(viso_count),]
  viso_count = do.call("rbind", lapply(split(viso_count, viso_count$gene), function(x){
    x$perc = x$n_of_reads/sum(x$n_of_reads)*100
    return(x)
  }))
  viso_count = viso_count[,c("id", "n_of_reads", "perc", "gene")]
  row.names(viso_count) = viso_count$id
  colnames(viso_count) = c("id", "cov", "perc", "gene")
  
  
  
  # 
  # 
  # # Number of isoforms 
  # # Total number
  # number_of_iso_total = rbind(
  #   paste(length(unique(viso_count$id)), " (", length(unique(viso_count$id[viso_count$perc > 1])), ")", sep = ""),
  #   paste(length(unique(flair_count$id)), " (", length(unique(flair_count$id[flair_count$perc > 1])), ")", sep = ""),
  #   paste(length(unique(st2_count$id)), " (", length(unique(st2_count$id[st2_count$perc > 1])), ")", sep = "")
  # )
  # 
  # # Per SIRV
  # number_of_iso_per_SIRV = rbind(
  #   paste(table(viso_count$gene), " (", table(viso_count$gene[viso_count$perc>1]), ")", sep = ""),
  #   paste(table(flair_count$gene), " (", table(flair_count$gene[flair_count$perc>1]), ")", sep = ""),
  #   paste(table(st2_count$gene), " (", table(st2_count$gene[st2_count$perc>1]), ")", sep = "")
  # )
  # 
  # number_of_iso = cbind(c("VIsoQLR", "FLAIR", "StringTie2"), number_of_iso_per_SIRV, number_of_iso_total)
  # colnames(number_of_iso) = c("Program", names(table(viso_count$gene)), "Total")
  # 
  # write.table(number_of_iso, "number_of_isoforms_table.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
  # 
  
  
  
  # Join coverage information to intersected isoforms
  viso_intersection_perc_2 = viso_intersection_perc[viso_intersection_perc$intersect>0.99,]
  viso_intersection_perc_2$coverage = viso_count[viso_intersection_perc_2$id2, "cov"]
  viso_intersection_perc_2$abundance = viso_count[viso_intersection_perc_2$id2, "perc"]
  viso_intersection_perc_2$gene = viso_count[viso_intersection_perc_2$id2, "gene"]
  viso_intersection_perc_2[duplicated(viso_intersection_perc_2$id2),]
  
  
  
  
  viso_intersection_perc_3 = do.call("rbind", lapply(sirv_list, function(x){
    selectedrows = viso_intersection_perc_2[viso_intersection_perc_2$id %in% x,]
    abundance = selectedrows[!duplicated(selectedrows$id2), "abundance"]
    coverage = selectedrows[!duplicated(selectedrows$id2), "coverage"]
    c(paste(x, collapse = " "), sum(abundance), perc_tag, sum(coverage))
  }))
  
  return(viso_intersection_perc_3)

}

viso_intersection25 = iso_intersection_viso(viso_path25, "0.25%")
viso_intersection05 = iso_intersection_viso(viso_path05, "0.5%")
viso_intersection1 = iso_intersection_viso(viso_path1, "1%")
viso_intersection2 = iso_intersection_viso(viso_path2, "2%")
viso_intersection3 = iso_intersection_viso(viso_path3, "3%")





# Create a df containing the theorical abundances of each transcript
sirv_list=as.list(unique(gold$id))
sirv_list = sirv_list[!sirv_list %in% c("SIRV701", "SIRV705", "SIRV604", "SIRV612", "SIRV502")]
sirv_list[[length(sirv_list)+1]] = c("SIRV701", "SIRV705")
sirv_list[[length(sirv_list)+1]] = c("SIRV604", "SIRV612")

sirv_df = data.frame(id = do.call("c", lapply(sirv_list, function(x) {paste(x, collapse = " ")})))
sirv_df$sirv_count = 1
sirv_df$sirv_count[sirv_df$id %in% c("SIRV701 SIRV705", "SIRV604 SIRV612")] = 2
sirv_df$sirv_count[sirv_df$id %in% c("SIRV502")] = 0

sirv_df$gene = gsub("..$", "", sub(" .*$", "", sirv_df$id))
sirv_df_perc = do.call("rbind", lapply(split(sirv_df, sirv_df$gene), function(x) {
  x$perc = x$sirv_count/sum(x$sirv_count)*100
  x$program = "Gold\nStandard"
  return(x[,c("id", "perc", "program")])
}))
sirv_df_perc <- unname(as.matrix(sirv_df_perc))




df_heatmap = data.frame(rbind(sirv_df_perc, viso_intersection25[,1:3], viso_intersection05[,1:3], viso_intersection1[,1:3], viso_intersection2[,1:3], viso_intersection3[,1:3]))
df_heatmap$X2 = as.numeric(as.character(df_heatmap$X2))
df_heatmap$X3 = factor(df_heatmap$X3, levels = c("Gold\nStandard", "0.25%", "0.5%", "1%", "2%", "3%"))
df_heatmap$X4 = gsub("..$", "", sub(" .*$", "", df_heatmap$X1))

ggplot(df_heatmap, aes(x = X3, y = X1, fill = X2)) + 
  geom_tile(colour = "gray") +
  scale_fill_gradient(low = "white", high = "turquoise4") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0,0)) +
  labs(fill = "Abundance (%)") +
  facet_col(~X4, scales = 'free_y', shrink=F, space = "free") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(colour="black", fill="white", linetype="solid"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# setwd("/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/")
ggsave(filename = "Supp.FigureX.pdf", device = "pdf", width = 4, height = 12)
ggsave(filename = "Supp.FigureX.jpeg", device = "jpeg", width = 4, height = 12)
ggsave(filename = "Supp.FigureX.png", device = "png", width = 4, height = 12)









# Coverage df creation
df_coverage = data.frame(cbind(viso_intersection25[,c(1,4)], viso_intersection05[,4], viso_intersection1[,4], viso_intersection2[,4], viso_intersection3[,4]))
colnames(df_coverage) = c("SIRV", "0.25%", "0.5%", "1%", "2%", "3%")

write.table(df_coverage, "coverage_table_visoqlr_percentages.tsv", col.names = T, row.names = F, quote = F, sep = "\t")






viso_count_fun = function(viso_path){
  
  viso_count = transcript_data_loading(viso_path, gene_paste = T)[,c("id", "n_of_reads", "gene")]
  viso_count = viso_count[!duplicated(viso_count),]
  viso_count = do.call("rbind", lapply(split(viso_count, viso_count$gene), function(x){
    x$perc = x$n_of_reads/sum(x$n_of_reads)*100
    return(x)
  }))
  viso_count = viso_count[,c("id", "n_of_reads", "perc", "gene")]
  row.names(viso_count) = viso_count$id
  colnames(viso_count) = c("id", "cov", "perc", "gene")
  
  return(viso_count)
}

viso_count25 = viso_count_fun(viso_path25)
viso_count05 = viso_count_fun(viso_path05)
viso_count1 = viso_count_fun(viso_path1)
viso_count2 = viso_count_fun(viso_path2)
viso_count3 = viso_count_fun(viso_path3)



# Number of isoforms
# Total number
number_of_iso_total = rbind(
  paste(length(unique(viso_count25$id)), " (", length(unique(viso_count25$id[viso_count25$perc > 1])), ")", sep = ""),
  paste(length(unique(viso_count05$id)), " (", length(unique(viso_count05$id[viso_count05$perc > 1])), ")", sep = ""),
  paste(length(unique(viso_count1$id)), " (", length(unique(viso_count1$id[viso_count1$perc > 1])), ")", sep = ""),
  paste(length(unique(viso_count2$id)), " (", length(unique(viso_count2$id[viso_count2$perc > 1])), ")", sep = ""),
  paste(length(unique(viso_count3$id)), " (", length(unique(viso_count3$id[viso_count3$perc > 1])), ")", sep = "")
)

# Per SIRV
number_of_iso_per_SIRV = rbind(
  paste(table(viso_count25$gene), " (", table(viso_count25$gene[viso_count25$perc>1]), ")", sep = ""),
  paste(table(viso_count05$gene), " (", table(viso_count05$gene[viso_count05$perc>1]), ")", sep = ""),
  paste(table(viso_count1$gene), " (", table(viso_count1$gene[viso_count1$perc>1]), ")", sep = ""),
  paste(table(viso_count2$gene), " (", table(viso_count2$gene[viso_count2$perc>1]), ")", sep = ""),
  paste(table(viso_count3$gene), " (", table(viso_count3$gene[viso_count3$perc>1]), ")", sep = "")
)

number_of_iso = cbind(c("0.25%", "0.5%", "1%", "2%", "3%"), number_of_iso_per_SIRV, number_of_iso_total)
colnames(number_of_iso) = c("Percentage", names(table(viso_count$gene)), "Total")

write.table(number_of_iso, "number_of_isoforms_table_visoqlr_percentages.tsv", col.names = T, row.names = F, quote = F, sep = "\t")









# Correlation calculation


abundance_split = split(df_heatmap, df_heatmap$X3)
abundance_split = lapply(abundance_split, function(x) x[order(x$X1),])

# cor(abundance_split$`Gold\nStandard`$X2, abundance_split$VIsoQLR$X2)
# cor(abundance_split$`Gold\nStandard`$X2, abundance_split$FLAIR$X2)
# cor(abundance_split$`Gold\nStandard`$X2, abundance_split$StringTie2$X2)
# 
# cor(abundance_split$FLAIR$X2, abundance_split$VIsoQLR$X2)

simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`0.25%`$X2), method = "cosine")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`0.5%`$X2), method = "cosine")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`1%`$X2), method = "cosine")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`2%`$X2), method = "cosine")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`3%`$X2), method = "cosine")

simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`0.25%`$X2), method = "euclidean")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`0.5%`$X2), method = "euclidean")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`1%`$X2), method = "euclidean")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`2%`$X2), method = "euclidean")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`3%`$X2), method = "euclidean")


df_cor_sirv = data.frame(row.names = c("0.25%", "0.5%", "1%", "2%", "3%"))
for(sirv in unique(df_heatmap$X4)){
  
  df_heatmap_subset = df_heatmap[df_heatmap$X4 == sirv,]
  abundance_split = split(df_heatmap_subset, df_heatmap_subset$X3)
  abundance_split = lapply(abundance_split, function(x) x[order(x$X1),])
  
  
  df_cor_sirv[,sirv]=c(
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`0.25%`$X2), method = "cosine"),
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`0.5%`$X2), method = "cosine"),
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`1%`$X2), method = "cosine"),
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`2%`$X2), method = "cosine"),
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$`3%`$X2), method = "cosine")
  )
  
}

write.table(df_cor_sirv, "cosine_dist_visoqlr_percentages.tsv", col.names = T, row.names = T, sep = "\t", quote = F)




