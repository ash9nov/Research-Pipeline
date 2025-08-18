#!/bin/bash

#Date: 05.01.2025

alias fancy="column -s $'\t' -t|less -S"
alias grp="grep -v '##'"


#config:
data_path="/home/ashishks/3X1TB_storage/1TB_1st/PHD/WES_PROJECT/proj_DATA/ensembl-vep-data"
GRCh_v37_decoy="/home/ashishks/3X1TB_storage/1TB_2nd/NGS_MetaData/GRCh_v37_decoy/human_g1k_v37_decoy.fasta"
NOC=34 #number of computing cores

#Files
I_P="WGStest2_100022779910.filtered.vcf.gz"
O_P="~/3X1TB_storage/1TB_2nd/WGStest2_quadro_analysis_11.08.2023/vcf/Pro-100022779910/WGStest2_100022779910.filtered_annotated.vcf"

#Annotation
alias VEP="~/3X1TB_storage/1TB_2nd/my_tools/VEP_26.09.2023/ensembl-vep/vep"
VEP \
--dir /home/ashishks/.vep/ \
--fork $NOC --offline --show_ref_allele --variant_class  --everything --nearest symbol --force_overwrite --buffer_size 20000 --assembly GRCh37 --tab \
--fasta $GRCh_v37_decoy \
-i $I_P \
-o $O_P \
--plugin SpliceAI,snv=$data_path/spliceai_scores/spliceai_scores.raw.snv.hg19.vcf.gz,indel=$data_path/spliceai_scores/spliceai_scores.raw.indel.hg19.vcf.gz \
--plugin dbNSFP,$data_path/dbNSFP4.3/dbNSFP4.3a/dbNSFP4.3a_grch37.gz,ALL \
--plugin AlphaMissense,file=$data_path/Alphamissense/AlphaMissense_hg19.tsv.gz,cols=all

################################################################################################
#FILTERATION
################################################################################################
alias FILTER_VEP="~/3X1TB_storage/1TB_2nd/my_tools/VEP_26.09.2023/ensembl-vep/filter_vep"

I_P="WGStest2_100022779910.filtered.annotated.vcf" 
O_P="WGStest2_100022779910.VEP-filterd.vcf"
FILTER_VEP \
--input_file $I_P \
--output_file $O_P \
--format tab \
--filter \
"(AF < 0.05 or not AF) \
and (gnomADg_AF < 0.05 or not gnomADg_AF) \
and (gnomAD_genomes_AF < 0.05 or not gnomAD_genomes_AF) \
and (CANONICAL is YES) \
and ((BIOTYPE is protein_coding) or (BIOTYPE is nonsense_mediated_decay)) \
and (Consequence != downstream_gene_variant) \
and (Consequence != intron_variant) \
and (Consequence != intron_variant,non_coding_transcript_variant) \
and (Consequence != mature_miRNA_variant) \
and (Consequence != non_coding_transcript_exon_variant) \
and (Consequence != synonymous_variant) \
and (Consequence != upstream_gene_variant) \
and (Consequence != intron_variant,NMD_transcript_variant) \
and (Consequence != splice_polypyrimidine_tract_variant,intron_variant)"


