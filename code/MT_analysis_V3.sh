#Custom CNV_ANALYSIS for control samples|at-PROD-server|Date: 02.10.2024
#_____________________________________________________
#1.CONFIGs
data="/data/NGS_AMG_misc/CNV-MT-data"
dataMT="/data/AKS_data/my_tools/"

Gene_panels="${data}/Dragen_CNV_Analysis/Gene_panels/bed"
Gene_panels_MT="${data}/MT_MetaData/BioMart_genes/MT_gene_panels"
GENEPANEL_MT="${data}/MT_MetaData/BioMart_genes/Mitochondria_gene_list.bed"

RUNS_dir="/data/NGS_AMG/GATK_WGS/RUNS"
GRCh_v37_decoy="${data}/GRCh_resources/GRCh_v37_decoy/human_g1k_v37_decoy.fasta"
HG38="${data}/GRCh_resources/HG38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

#[SVDB]
SVDBsif="${data}/SVDB/SVDB_latest.sif"
#[VEP]
VEPsif="${data}/ENSEMBL-VEP/VEP112/vep.sif"
VEPcherry="${data}/ENSEMBL-VEP/VEP112/vep-cherrypicked-1115-issue"

NOC=32 #NumberOfComputingcores

#Config for VariantTypes
CNV="cnv_sv"
SV="sv"
MERGED="merged"
MT="mt"

#WORKINGDIR="/data/AKS_data/staging/ControlSample_SV_annotation"
WORKINGDIR="/data/AKS_data/staging/24GENOM5-241008"

#bed="${WORKINGDIR}/bed"
bed_MT="${WORKINGDIR}/bed_MT"
#______________________________________________________________________________
#1. SampleSheet
    # Do this work manually
#2.Functions
    #F0: BedProcessingMT
        BedProcessingMT(){
            cp ${Gene_panels_MT}/${Gpanel}.bed ${bed_MT}
            cat ${bed_MT}/${Gpanel}.bed |cut -f 1,2,3,5 > ${bed_MT}/${Gpanel}_columns.bed
            bgzip -c ${bed_MT}/${Gpanel}_columns.bed > ${bed_MT}/${Gpanel}_columns.bed.gz
            tabix -f -p bed ${bed_MT}/${Gpanel}_columns.bed.gz
        }
    #F1: PanelChoppingMT
        PanelChoppingMT(){
            bedtools intersect -header \
                -a ${SAMPLE_dirMT}/${sample}.${VType}.vcf \
                -b ${bed_MT}/${Gpanel}_columns.bed | \
            bcftools norm -d none \
                > ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped.vcf
        }
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
    #F3 PostVepAnnotatingMT
        PostVepAnnotatingMT(){
            join -a 1 -j 1 \
                <(cat \
                        <(echo -e "Location\tREF\tALT\tDP_full\tMQ\tGT\tAD(alt,ref)\tAF\tDP_vc\tFILTER\tSQ\tF1R2\tF2R1") \
                        <(paste \
                            <(grep -v '#' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped.vcf|cut -f 1,2|awk -v OFS='\t' '{print $1":"$2}') \
                            <(grep -v '#' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped.vcf|cut -f 4,5) \
                            <(grep -v '#' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped.vcf|cut -f 8|sed 's/;/\t/g'|sed 's/DP=//g'|sed 's/MQ=//g'|cut -f 1-2) \
                            <(grep -v '#' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped.vcf|cut -f 10|cut -d ':' -f 1-8|tr ':' '\t') \
                        ) \
                        |sort -k1 \
                ) \
                <(join -a 1 -j1 \
                    <(grep -v '##' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf|cut -f 2,3,8,24-36,38,40-45,49-53,55-56,59-60,63-67|sort -k1|uniq) \
                    <(cat \
                        <(echo -e "Location\tgDNA-for-ALAMUT\tcDNA-for-ALAMUT") \
                        <(join -a 1 \
                            <(grep -v '##' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf|cut -f 2,3,8,20,24|grep -v 'Location'|sed 's/:/\t/g'| awk -v OFS='\t' '{print$1":"$2,"m."$2$5">"$3}'|sort -k1|uniq) \
                            <(grep -v '##' ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf|cut -f 2,3,8,20,24|grep -v 'Location'|awk -v OFS='\t' '$3!="-" {print$1, "NC_012920.1("$5"):n."$3$4">"$2}'|sort -k1) \
                            |sed 's/ /\t/g'|sort -k1 \
                        )
                    )|sed 's/ /\t/g' \
                )|sed 's/ /\t/g' \
                |awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$16,$50,$51,$17,$19,$20,$21,$41,$42,$45,$46,$22,$23,$24,$25,$26,$27,$30,$31,$32,$33}'|uniq \
                > ${SAMPLE_dirMT}/${sample}.${VType}.${Gpanel}-prefinal.tsv
        }
        
    #F4: F0+F1+F2+F3 Combinedfun
        Combinedfun(){
            #MT
                VType="${MT}"
                    #BedProcessingMT
                    #PanelChoppingMT
                    #VepAnnotatingMT
                    PostVepAnnotatingMT
        }
#__________________________________________________________________________________________________________________________

#3.Chopping/Merging/Indexing/VepAnnotating/Filtering/AddingMissingTranscripts/ColumnChopping of CNV,SV &Merged for FocusGene/GenePanel
    #While looping for all the samples in run
        while IFS=$'\t' read -r -a SampleSheet; do
            #Variables:
                sample=${SampleSheet[0]}
                Gpanel=${SampleSheet[1]}
            #Sample config
                SAMPLE_dirMT="${WORKINGDIR}/${sample}/MT/RAW"
                echo "-----"$sample"_&_"$Gpanel"-----"
            #running all the samples
                Combinedfun
        done < ${WORKINGDIR}/SampleSheetMT
#_____________________________________________________
#4. Copy to FINAL
    while IFS=$'\t' read -r -a SampleSheet; do
        #Variables:
            sample=${SampleSheet[0]}
            Gpanel=${SampleSheet[1]}
        #Copy
            VType="${MT}"
        #    cp ${sample}/MT/RAW/${sample}.${VType}.${Gpanel}-final.tsv ${sample}/MT/FINAL
            
    done < SampleSheetMT