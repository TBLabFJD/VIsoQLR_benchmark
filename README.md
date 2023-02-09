# VIsoQLR_benchmark

This repostory contains the scripts used for data preprocessing, isoform calling and isoform comparison of VIsoQLR, FLAIR and StringTie2.

SIRV transcripts were analyced on a public PacBio RNA-seq Data (https://github.com/PacificBiosciences/DevNet/wiki/Sequel-II-System-Data-Release:-Universal-Human-Reference-(UHR)-Iso-Seq).


Scripts:
 - `ìsoform_analysis.sh` contains the data preprocessing (data acquisition and instructions, GMAP intex generation, FLAIR and StringTie2 anlalysis with dafault parameters). 
 - `isoform_comparison.R` contains the isoform comparison between VIsoQLR, FLAIR, StringTie2. 
 - `isoform_analysis_visoqlr_perc.R` contains the isoform comparison between VIsoQLR with different read threshold values.
