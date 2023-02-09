
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
st2_path = "/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/stringtie_results/flnc.gtf"
flair_path = "/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/flair_results/flnc.flair.collapse.isoforms.bed6"
viso_path = "/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/visoqlr_results/visoqlr_results.gtf"
flair_count_path = "/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/flair_results/flnc.flair.collapse.firstpass.q.counts"

# Minimap2
# st2_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/stringtie_results_minimap/flnc.gtf"
# flair_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/flair_results_minimap/flnc.flair.collapse.isoforms.bed6"
# viso_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/visoqlr_results_minimap/visoqlr_results.gtf"
# flair_count_path = "/home/gonzalo/Documents/VIsoQLR_article/revision/uhrr/flair_results_minimap/flnc.flair.collapse.firstpass.q.counts"





gs_path = "/home/gonzalo/Documents/old_computer/VIsoQLR_article/revision/uhrr/documentation/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"
gold = transcript_data_loading(gs_path)[,c("gene", "start", "end", "id", "size")]
gold$start = gold$start - 1


st2 = transcript_data_loading(st2_path)[,c("gene", "start", "end", "id", "size")]
st2$start = st2$start - 1


flair = bed_data_loading(flair_path)[,c("gene", "start", "end", "id", "size")]
flair$start = flair$start - 1
flair$id = gsub(";0", "", flair$id)

viso = transcript_data_loading(viso_path, gene_paste = T)[,c("gene", "start", "end", "id", "size")]
viso$start = viso$start - 1



#==================#
# RUN INTERSECTION #
#==================#

gold_intersection_perc = match_intersection_calculation(gold, gold)

st2_intersection_perc = match_intersection_calculation(gold, st2)

flair_intersection_perc = match_intersection_calculation(gold, flair)

viso_intersection_perc = match_intersection_calculation(gold, viso)



length(unique(gold$id))
length(unique(st2$id))
length(unique(flair$id))
length(unique(viso$id))



intersect_threshold = 0.99
# intersect_threshold = 0.999
# intersect_threshold = 0.9995

table(gold_intersection_perc$intersect>intersect_threshold)
length(unique(gold_intersection_perc[gold_intersection_perc$intersect>intersect_threshold,]$id))

table(st2_intersection_perc$intersect>intersect_threshold)
length(unique(st2_intersection_perc[st2_intersection_perc$intersect>intersect_threshold,]$id))

table(flair_intersection_perc$intersect>intersect_threshold)
length(unique(flair_intersection_perc[flair_intersection_perc$intersect>intersect_threshold,]$id))

table(viso_intersection_perc$intersect>intersect_threshold)
length(unique(viso_intersection_perc[viso_intersection_perc$intersect>intersect_threshold,]$id))







#================#
# COUNTS LOADING #
#================#



# STRINGTIE2
st2_count = transcript_data_loading(st2_path)[,c("id", "cov", "gene")]
st2_count = st2_count[!duplicated(st2_count),]
st2_count = do.call("rbind", lapply(split(st2_count, st2_count$gene), function(x){
  x$perc = x$cov/sum(x$cov)*100
  return(x)
}))
st2_count = st2_count[,c("id", "cov", "perc", "gene")]
row.names(st2_count) = st2_count$id
colnames(st2_count) = c("id", "cov", "perc", "gene")


# FLAIR
flair_count = read.delim(flair_count_path, stringsAsFactors = F, header = F, quote = "")
flair_tmp = flair[,c("gene", "id")]
flair_tmp = flair_tmp[!duplicated(flair_tmp),]
flair_count_tmp = merge(flair_count, flair_tmp, by.x = "V1", by.y = "id")
flair_count = do.call("rbind", lapply(split(flair_count_tmp, flair_count_tmp$gene), function(x){
  x$perc = x$V2/sum(x$V2)*100
  return(x)
}))
flair_count = flair_count[,c("V1", "V2", "perc", "gene")]
row.names(flair_count) = flair_count$V1
colnames(flair_count) = c("id", "cov", "perc", "gene")


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





# Number of isoforms 
# Total number
number_of_iso_total = rbind(
  paste(length(unique(viso_count$id)), " (", length(unique(viso_count$id[viso_count$perc > 1])), ")", sep = ""),
  paste(length(unique(flair_count$id)), " (", length(unique(flair_count$id[flair_count$perc > 1])), ")", sep = ""),
  paste(length(unique(st2_count$id)), " (", length(unique(st2_count$id[st2_count$perc > 1])), ")", sep = "")
)

# Per SIRV
number_of_iso_per_SIRV = rbind(
  paste(table(viso_count$gene), " (", table(viso_count$gene[viso_count$perc>1]), ")", sep = ""),
  paste(table(flair_count$gene), " (", table(flair_count$gene[flair_count$perc>1]), ")", sep = ""),
  paste(table(st2_count$gene), " (", table(st2_count$gene[st2_count$perc>1]), ")", sep = "")
)

number_of_iso = cbind(c("VIsoQLR", "FLAIR", "StringTie2"), number_of_iso_per_SIRV, number_of_iso_total)
colnames(number_of_iso) = c("Program", names(table(viso_count$gene)), "Total")

write.table(number_of_iso, "number_of_isoforms_table.tsv", col.names = T, row.names = F, quote = F, sep = "\t")




# Join coverage information to intersected isoforms
st2_intersection_perc_2 = st2_intersection_perc[st2_intersection_perc$intersect>0.99,]
st2_intersection_perc_2$coverage = st2_count[st2_intersection_perc_2$id2, "cov"]
st2_intersection_perc_2$abundance = st2_count[st2_intersection_perc_2$id2, "perc"]
st2_intersection_perc_2$gene = st2_count[st2_intersection_perc_2$id2, "gene"]
st2_intersection_perc_2[duplicated(st2_intersection_perc_2$id2),]

flair_intersection_perc_2 = flair_intersection_perc[flair_intersection_perc$intersect>0.99,]
flair_intersection_perc_2$coverage = flair_count[flair_intersection_perc_2$id2, "cov"]
flair_intersection_perc_2$abundance = flair_count[flair_intersection_perc_2$id2, "perc"]
flair_intersection_perc_2$gene = flair_count[flair_intersection_perc_2$id2, "gene"]
flair_intersection_perc_2[duplicated(flair_intersection_perc_2$id2),]

viso_intersection_perc_2 = viso_intersection_perc[viso_intersection_perc$intersect>0.99,]
viso_intersection_perc_2$coverage = viso_count[viso_intersection_perc_2$id2, "cov"]
viso_intersection_perc_2$abundance = viso_count[viso_intersection_perc_2$id2, "perc"]
viso_intersection_perc_2$gene = viso_count[viso_intersection_perc_2$id2, "gene"]
viso_intersection_perc_2[duplicated(viso_intersection_perc_2$id2),]




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



st2_intersection_perc_3 = do.call("rbind", lapply(sirv_list, function(x){
  selectedrows = st2_intersection_perc_2[st2_intersection_perc_2$id %in% x,]
  abundance = selectedrows[!duplicated(selectedrows$id2), "abundance"]
  coverage = selectedrows[!duplicated(selectedrows$id2), "coverage"]
  c(paste(x, collapse = " "), sum(abundance), "StringTie2", sum(coverage))
}))


flair_intersection_perc_3 = do.call("rbind", lapply(sirv_list, function(x){
  selectedrows = flair_intersection_perc_2[flair_intersection_perc_2$id %in% x,]
  abundance = selectedrows[!duplicated(selectedrows$id2), "abundance"]
  coverage = selectedrows[!duplicated(selectedrows$id2), "coverage"]
  c(paste(x, collapse = " "), sum(abundance), "FLAIR", sum(coverage))
}))


viso_intersection_perc_3 = do.call("rbind", lapply(sirv_list, function(x){
  selectedrows = viso_intersection_perc_2[viso_intersection_perc_2$id %in% x,]
  abundance = selectedrows[!duplicated(selectedrows$id2), "abundance"]
  coverage = selectedrows[!duplicated(selectedrows$id2), "coverage"]
  c(paste(x, collapse = " "), sum(abundance), "VIsoQLR", sum(coverage))
}))


df_heatmap = data.frame(rbind(sirv_df_perc, st2_intersection_perc_3[,1:3], flair_intersection_perc_3[,1:3], viso_intersection_perc_3[,1:3]))
df_heatmap$X2 = as.numeric(as.character(df_heatmap$X2))
df_heatmap$X3 = factor(df_heatmap$X3, levels = c("Gold\nStandard", "VIsoQLR", "FLAIR", "StringTie2"))
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
ggsave(filename = "Figure5.pdf", device = "pdf", width = 4, height = 12)
ggsave(filename = "Figure5.jpeg", device = "jpeg", width = 4, height = 12)
ggsave(filename = "Figure5.png", device = "png", width = 4, height = 12)

# 
# ggsave(filename = "Figure6.pdf", device = "pdf", width = 4, height = 12)
# ggsave(filename = "Figure6.jpeg", device = "jpeg", width = 4, height = 12)
# ggsave(filename = "Figure6.png", device = "png", width = 4, height = 12)







# Coverage df creation
df_coverage = data.frame(cbind(viso_intersection_perc_3[,c(1,4)], flair_intersection_perc_3[,4], st2_intersection_perc_3[,4]))
colnames(df_coverage) = c("SIRV", "VIsoQLR", "Flair", "StringTie2")
df_coverage$StringTie2 = round(as.numeric(as.character(df_coverage$StringTie2)), 2)

write.table(df_coverage, "coverage_table.tsv", col.names = T, row.names = F, quote = F, sep = "\t")








# Correlation calculation


abundance_split = split(df_heatmap, df_heatmap$X3)
abundance_split = lapply(abundance_split, function(x) x[order(x$X1),])

# cor(abundance_split$`Gold\nStandard`$X2, abundance_split$VIsoQLR$X2)
# cor(abundance_split$`Gold\nStandard`$X2, abundance_split$FLAIR$X2)
# cor(abundance_split$`Gold\nStandard`$X2, abundance_split$StringTie2$X2)
# 
# cor(abundance_split$FLAIR$X2, abundance_split$VIsoQLR$X2)

simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$VIsoQLR$X2), method = "cosine")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$FLAIR$X2), method = "cosine")
simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$StringTie2$X2), method = "cosine")

dist(rbind(abundance_split$`Gold\nStandard`$X2, abundance_split$VIsoQLR$X2), method = "euclidean")
dist(rbind(abundance_split$`Gold\nStandard`$X2, abundance_split$FLAIR$X2), method = "euclidean")
dist(rbind(abundance_split$`Gold\nStandard`$X2, abundance_split$StringTie2$X2), method = "euclidean")

df_cor_sirv = data.frame(row.names = c("VIsoQLR", "FLAIR", "StringTie2"))
for(sirv in unique(df_heatmap$X4)){

  df_heatmap_subset = df_heatmap[df_heatmap$X4 == sirv,]
  abundance_split = split(df_heatmap_subset, df_heatmap_subset$X3)
  abundance_split = lapply(abundance_split, function(x) x[order(x$X1),])
  
  # df_cor_sirv[,sirv]=c(
  #   cor(abundance_split$`Gold\nStandard`$X2, abundance_split$VIsoQLR$X2, method = "spearman"),
  #   cor(abundance_split$`Gold\nStandard`$X2, abundance_split$FLAIR$X2),
  #   cor(abundance_split$`Gold\nStandard`$X2, abundance_split$StringTie2$X2)
  # )
  
  df_cor_sirv[,sirv]=c(
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$VIsoQLR$X2), method = "cosine"),
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$FLAIR$X2), method = "cosine"),
    simil(list(abundance_split$`Gold\nStandard`$X2, abundance_split$StringTie2$X2), method = "cosine")
  )

}

write.table(df_cor_sirv, "cosine_dist.tsv", col.names = T, row.names = T, sep = "\t", quote = F)




