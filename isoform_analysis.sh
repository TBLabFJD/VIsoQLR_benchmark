#!/bin/bash




#================#
# Download  Data #
#================#

# PacBio RNA-seq Data
# wget https://downloads.pacbcloud.com/public/dataset/UHR_IsoSeq/PolishedMappedTranscripts/out.abundance.txt
# wget https://downloads.pacbcloud.com/public/dataset/UHR_IsoSeq/PolishedMappedTranscripts/out.fastq
# wget https://downloads.pacbcloud.com/public/dataset/UHR_IsoSeq/PolishedMappedTranscripts/out.gff
wget https://downloads.pacbcloud.com/public/dataset/UHR_IsoSeq/FullLengthReads/flnc.bam
bamToFastq -i flnc.bam -fq flnc.fastq

# Download SIRV information (sequence and metadata)
wget https://www.lexogen.com/wp-content/uploads/2021/06/SIRV_Set1_Norm_Sequences_20210507.zip
unzip SIRV_Set1_Norm_Sequences_20210507.zip


# Download all gene sequences (human)
# Downloaded with BioMart. Ensembl Genes 108 -> Human genes (GRCh38.p13) -> 
# -> Filters -> Gene type: protein_coding -> MANE Select: Only ->
# -> Attributes -> Sequences -> SEQUENCES -> Unspliced (Gene) ->
# -> HEADER INFORMATION -> Gene stable ID, Gene stable ID version, Gene name ->
# -> Results Export  all results to -> Compressed file (.gz) -> FASTA -> Unique results only -> Go


module load miniconda/3.7

#===========#
# Variables #
#===========#
working_dir="/mnt/tblab/gonzalo/VIsoQLR_article/uhrr"
treads="20"
fastq="${working_dir}/raw_data/flnc.fastq"
sample="$(basename ${fastq} .fastq)" 

# genomePrefix="SIRV"
# genomePrefix="PAX6"
# genomePrefix="TP53"
genomePrefix="human"

aligner="GMAP"
# aligner="Minimap2"



#============================#
# Reference Index Generation #
#============================#

fasta="${working_dir}/references/${genomePrefix}.fasta"
reference="${working_dir}/references/gmap_index_${genomePrefix}/"
mkdir ${reference}
gmap_build -D ${reference} -d ${genomePrefix} ${fasta}



#=========#
# Mapping #
#=========#

mapping_dir="${working_dir}/mapped_${genomePrefix}_${aligner}"
mkdir ${mapping_dir}

sam="${mapping_dir}/${sample}.unsorted.sam"
bam="${mapping_dir}/${sample}.unsorted.bam"
error="${mapping_dir}/log_${sample}.err"
sorted="${mapping_dir}/${sample}.sorted.bam"
bed="${mapping_dir}/${sample}.bed"



if [[aligner="GMAP"]] ; then



# BAM output

gmap -n1 -t ${treads} --cross-species -f samse -D ${reference} -d ${genomePrefix} ${fastq} > ${sam} 2> ${error}
samtools view -F 4 -Su ${sam} > ${bam}
samtools sort --threads ${treads} ${bam} > ${sorted}
samtools index ${sorted}

bamToBed -split -i ${bam} > ${bed}


elif [[aligner="Minimap2"]] ; then 


minimap2 -ax splice:hq -uf --secondary=no -t ${treads} ${fasta} ${fastq} | samtools view -Su | samtools sort > ${sorted}
samtools index ${sorted}


fi









#============#
# StringTie2 #
#============#

stringtie="/home/gonzalo/software/stringtie-2.2.1.Linux_x86_64/stringtie"
${stringtie} ${sorted} -f 0 -L -p ${treads} -o ${working_dir}/stringtie_results/${sample}.gtf






#=======#
# FLAIR #
#=======#

bamToBed -i ${sorted} -bed12 > ${mapping_dir}/${sample}.sorted.bed12

annot_file="${working_dir}/documentation/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"
flair correct -q ${mapping_dir}/${sample}.sorted.bed12 -f ${annot_file} -g ${fasta} -o ${working_dir}/flair_results/${sample}.flair.concat
flair collapse -g ${fasta} -q ${working_dir}/flair_results/${sample}.flair.concat_all_corrected.bed -r ${fastq} -o ${working_dir}/flair_results/${sample}.flair.collapse --keep_intermediate

bed12ToBed6 -i ${working_dir}/flair_results/${sample}.flair.collapse.isoforms.bed > ${working_dir}/flair_results/${sample}.flair.collapse.isoforms.bed6




# minimap2 -ax splice --secondary=no -t 4 ${working_dir}/references/gencode.v42.pc_transcripts.fa ${working_dir}/flnc.fastq

# flair align --threads 8 -g ${fasta} -r ${fastq} -o ${working_dir}/flair_results_minimap/${sample}.flair.align
# flair correct -q ${working_dir}/flair_results_minimap/${sample}.flair.align.bed -f ${annot_file} -g ${fasta} -o ${working_dir}/flair_results_minimap/${sample}.flair.concat
# flair collapse -g ${fasta} -q ${working_dir}/flair_results_minimap/${sample}.flair.concat_all_corrected.bed -r ${fastq} -o ${working_dir}/flair_results_minimap/${sample}.flair.collapse --keep_intermediate
# bed12ToBed6 -i ${working_dir}/flair_results_minimap/${sample}.flair.collapse.isoforms.bed > ${working_dir}/flair_results_minimap/${sample}.flair.collapse.isoforms.bed6




