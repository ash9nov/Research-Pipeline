#! /bin/bash

#Loading the PATHs configurations:
        . /data/AKS_data/Research_pipleine/code/PATHs.conf

#Binder declartation (for path variables defined in PATHs.configs)
# BINDERS=(
        #         -B $data:/data \
        #         -B $HG38:/HG38 \
        #         -B $GRCh_v37_decoy:/GRCh_v37_decoy \
        #         -B $VEP_cache:/VEP_cache \
        #         -B $VEP_PLUGINS_DATA:/VEP_PLUGINS_DATA \
        #         -B $Functional_effect:/Functional_effect \
        #         -B $MaveDB:/MaveDB \
        #         -B $Gene_tolerance_to_change:/Gene_tolerance_to_change \
        #         -B $DosageSensitivity:/DosageSensitivity \
        #         -B $LOEUF:/LOEUF \
        #         -B $pLI:/pLI \
        #         -B $Motif:/Motif \
        #         -B $FunMotifs:/FunMotifs \
        #         -B $Nearby_features:/Nearby_features \
        #         -B $Downstream:/Downstream \
        #         -B $TSSDistance:/TSSDistance \
        #         -B $Pathogenicity_predictions:/Pathogenicity_predictions \
        #         -B $dbNSFP:/dbNSFP \
        #         -B $LOFTEE:/LOFTEE \
        #         -B $Phenotype_data_and_citations:/Phenotype_data_and_citations \
        #         -B $DisGeNET:/DisGeNET \
        #         -B $Mastermind:/Mastermind \
        #         -B $SatMutMPRA:/SatMutMPRA \
        #         -B $Regulatory_impact:/Regulatory_impact \
        #         -B $Enformer:/Enformer \
        #         -B $Splicing_predictions:/Splicing_predictions \
        #         -B $SpliceAI:/SpliceAI \
        #         -B $Structural_variant_data:/Structural_variant_data \
        #         -B $CNV_annotation_CUSTOM:/CNV_annotation_CUSTOM \
        #         -B $StructuralVariantOverlap:/StructuralVariantOverlap \
        #         -B $Transcript_annotation:/Transcript_annotation \
        #         -B $NMD:/NMD \
        #         -B $RiboseqORFs:/RiboseqORFs \
        #         -B $UTRAnnotator:/UTRAnnotator \
        #         -B $Variant_data:/Variant_data \
        #         -B $LOVD:/LOVD \
        #         -B $IP_DIR:/IP_DIR \
        #         -B $OP_DIR:/OP_DIR
        # )
        
#Loading the VEP Singularity
        singularity run "${BINDERS[@]}" $VEPsif
#Running the VEP Singularity
        singularity exec --cleanenv --no-home "${BINDERS[@]}" $VEPsif vep \
                --dir /VEP_cache/ --refseq --offline \
                --assembly GRCh38 --fasta /HG38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
                --fork $NOC \
                --buffer_size 400000 \
                --force_overwrite --tab \
                --show_ref_allele --variant_class  --everything --nearest symbol \
                --plugin MaveDB,file=/MaveDB/MaveDB_variants.tsv.gz,cols=all \
                --plugin Downstream \
                --plugin TSSDistance \
                --plugin dbNSFP,/dbNSFP/dbNSFP4.7a_grch38.gz,ALL \
                --plugin satMutMPRA,file=/SatMutMPRA/satMutMPRA_GRCh38_ALL.gz,cols=ALL \
                --plugin Enformer,file=/Enformer/enformer_grch38.vcf.gz \
                --plugin SpliceAI,snv=/SpliceAI/spliceai_scores/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/SpliceAI/spliceai_scores/spliceai_scores.raw.indel.hg38.vcf.gz \
                --plugin NMD \
                --plugin UTRAnnotator,file=/UTRAnnotator/uORF_5UTR_GRCh38_PUBLIC.txt \
                --plugin LOVD \
                -i /IP_DIR/WGStest_100016253159.mini-test-data.vcf \
                -o /OP_DIR/WGStest_100016253159.mini-test-OP.csv

                #--plugin dbNSFP,/dbNSFP/dbNSFP5.1a_grch38.gz,ALL \
                #--plugin DosageSensitivity,file=/DosageSensitivity/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz \
                #--plugin RiboseqORFs.file=/RiboseqORFs/Ribo-seq_ORFs.bed.gz \






#SNP/INDEL
#Copied from Annotation-Filteration_code_V5.sh
#Annotation
    singularity exec --cleanenv --no-home -B /data/NGS_AMG_misc/CNV-MT-data:/data,${SAMPLE_dir}:/SAMPLE_dir $VEPsif \
            vep \
            --dir /home/ashishks/.vep/ \
            --fork $NOC --offline --show_ref_allele --variant_class  --everything --nearest symbol --force_overwrite --buffer_size 20000 --assembly GRCh37 --tab \
            --fasta $GRCh_v37_decoy \
            -i $I_P \
            -o $O_P \
            --plugin SpliceAI,snv=$data_path/spliceai_scores/spliceai_scores.raw.snv.hg19.vcf.gz,indel=$data_path/spliceai_scores/spliceai_scores.raw.indel.hg19.vcf.gz \
            --plugin dbNSFP,$data_path/dbNSFP4.3/dbNSFP4.3a/dbNSFP4.3a_grch37.gz,ALL \
            --plugin AlphaMissense,file=$data_path/Alphamissense/AlphaMissense_hg19.tsv.gz,cols=all

#CNV/SV
#Copied from DRAGEN_CNV_Analysis_code_V2_Prod_sequencial_final.sh
#F4: VepAnnotating 
        VepAnnotating(){
            #singularity exec --cleanenv --no-home -B /data/NGS_AMG_stage/NGS_AMG_misc:/data,${SAMPLE_dir}:/SAMPLE_dir $VEPsif \
            singularity exec --cleanenv --no-home -B /data/NGS_AMG_misc/CNV-MT-data:/data,${SAMPLE_dir}:/SAMPLE_dir $VEPsif \
            vep \
                --dir /data/ENSEMBL-VEP/VEP112/vep_cache \
                --fork $NOC --offline --show_ref_allele --variant_class --everything --force_overwrite --buffer_size 100 --assembly GRCh37 --max_sv_size 1000000 --tab \
                --fasta /data/GRCh_resources/GRCh_v37_decoy/human_g1k_v37_decoy.fasta \
                --plugin StructuralVariantOverlap,file=/data/ENSEMBL-VEP/ensembl-vep-plugin-data/StructuralVariantOverlap/gnomad_v2.1_sv.sites.vcf.gz \
                -i /SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped.vcf \
                -o /SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf \
                --custom file=/SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped.vcf.gz,short_name=Sample,format=vcf,type=exact,coords=0,fields=FILTER    
        }

#MT
#Copied from MT_analysis_V3.sh
#F2: VepAnnotatingMT
        VepAnnotatingMT(){
            singularity exec --cleanenv --no-home -B ${data}:/data,${dataMT}:/dataMT,${SAMPLE_dirMT}:/SAMPLE_dirMT,${bed_MT}:/bed_MT $VEPsif \
            vep \
                --dir /data/ENSEMBL-VEP/VEP112/vep_cache \
                --refseq \
                --fork $NOC --offline --force_overwrite --assembly GRCh38 --tab \
                --fasta /data/GRCh_resources/HG38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
                -i /SAMPLE_dirMT/${sample}.${VType}.${Gpanel}-chopped.vcf \
                -o /SAMPLE_dirMT/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf \
                --custom file=/bed_MT/${Gpanel}_columns.bed.gz,short_name=GeneSymbol,format=bed \
                --custom file=/data/MT_MetaData/ClinVar_data/clinvar_20240611.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=GENEINFO%CLNSIG%CLNREVSTAT%CLNDN \
                --custom file=/data/MT_MetaData/gnomAD_data/gnomad.genomes.v3.1.sites.chrM.vcf.bgz,short_name=gnomAD,format=vcf,type=exact,coords=0,fields=AN%AC_hom%AC_het%AF_hom%AF_het%max_hl \
                --custom file=/dataMT/MT_MetaData/HelixMTdb/HelixMTdb_20200327.vcf.gz,short_name=HelixMTdb,format=vcf,type=exact,coords=0,fields=feature%gene%counts_hom%AF_hom%counts_het%AF_het%mean_ARF%max_ARF%haplogroups_for_homoplasmic_variants%haplogroups_for_heteroplasmic_variants \
                --custom file=/dataMT/MT_MetaData/MitoMapMTdb/MitoMap_disease.vcf.gz,short_name=MitoMap_disease,format=vcf,type=exact,coords=0,fields=AC%AF%aachange%homoplasmy%heteroplasmy%PubmedIDs%Disease%DiseaseStatus%HGFL \
                --custom file=/dataMT/MT_MetaData/MitoMapMTdb/MitoMap_polymorphisms.vcf.gz,short_name=MitoMap_polymorphisms,format=vcf,type=exact,coords=0,fields=AC%AF%HGFL \
                --custom file=/dataMT/MT_MetaData/MitoMapMTdb/MitoMap_mitotip-scores.vcf.gz,short_name=MitoTip,format=vcf,type=exact,coords=0,fields=MitoTIP_Score%Quartile%Count%Percentage%Mitomap_Status
        }