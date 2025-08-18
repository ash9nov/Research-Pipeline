#Custom CNV_ANALYSIS for control samples|at-PROD-server|Date: 02.10.2024
#_____________________________________________________
#0.A.CONFIGs
data="/data/NGS_AMG_misc/CNV-MT-data"

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

WORKINGDIR="/data/AKS_data/staging/ControlSample_SV_annotation"
bed="${WORKINGDIR}/bed"

#_____________________________________________________
#0.B.Functions
    #F1: PanelChopping
        PanelChopping(){
            bedtools intersect -header \
                -a ${SAMPLE_dir}/${sample}.${VType}.vcf \
                -b ${bed}/${Gpanel}.bed | \
            bcftools norm -d none \
                > ${SAMPLE_dir}${sample}.${VType}.${Gpanel}-chopped.vcf
        }
    #F2: VariantMerging
        VariantMerging(){
            singularity exec --cleanenv --no-home -B ${SAMPLE_dir}:/WDIR $SVDBsif \
            svdb \
                --merge --vcf \
                /WDIR/${sample}.${CNV}.${Gpanel}-chopped.vcf \
                /WDIR/${sample}.${SV}.${Gpanel}-chopped.vcf \
                --overlap 0.7 --bnd_distance 2500 \
                > ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped.vcf
        }
    #F3: Compressing
        Compressing(){
            bgzip -c \
                ${sample}.${VType}.${Gpanel}-chopped.vcf \
                > ${sample}.${VType}.${Gpanel}-chopped.vcf.gz
            tabix ${sample}.${VType}.${Gpanel}-chopped.vcf.gz
        }
    #F4: VepAnnotating
        VepAnnotating(){
            #singularity exec --cleanenv --no-home -B /data/NGS_AMG_stage/NGS_AMG_misc:/data,${SAMPLE_dir}:/SAMPLE_dir $VEPsif \
            singularity exec --cleanenv --no-home -B /data/NGS_AMG_misc/CNV-MT-data:/data,${SAMPLE_dir}:/SAMPLE_dir $VEPcherry \
            vep \
                --dir /data/ENSEMBL-VEP/VEP112/vep_cache \
                --fork $NOC --offline --show_ref_allele --variant_class --everything --force_overwrite --buffer_size 100 --assembly GRCh37 --max_sv_size 1000000 --tab \
                --fasta /data/GRCh_resources/GRCh_v37_decoy/human_g1k_v37_decoy.fasta \
                --plugin StructuralVariantOverlap,file=/data/ENSEMBL-VEP/ensembl-vep-plugin-data/StructuralVariantOverlap/gnomad_v2.1_sv.sites.vcf.gz \
                -i /SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped.vcf \
                -o /SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf \
                --custom file=/SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped.vcf.gz,short_name=Sample,format=vcf,type=exact,coords=0,fields=FILTER    
        }
    #F5: VepFiltering
        VepFiltering(){
            singularity exec --cleanenv --no-home -B /data/NGS_AMG_misc/CNV-MT-data:/data,${SAMPLE_dir}:/SAMPLE_dir,${WORKINGDIR}:/WORKINGDIR $VEPsif \
            filter_vep \
                -i /SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf \
                -o /SAMPLE_dir/${sample}.${VType}.${Gpanel}-chopped-annotated-filtered.vcf \
                --format tab \
                --force_overwrite \
                --filter $filter"(Feature in /WORKINGDIR/genepanelcolumn5.bed) or (not Feature)"
        }
    #F6: MissedTranscriptAddition 
        MissedTranscriptAddition(){
            cat \
                ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated-filtered.vcf \
                <(cat \
                        <(grep "#Uploaded_variation" ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated-filtered.vcf) \
                        <(egrep \
                            -f <(diff \
                                    <(grep -v '##' ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated-filtered.vcf|cut -f 1|uniq) \
                                    <(grep -v '##' ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf|cut -f 1|uniq)\
                                    |sed s'/> //g'\
                                )\
                            ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated.vcf\
                        )\
                    |singularity exec --cleanenv --no-home -B /data/NGS_AMG_misc/CNV-MT-data:/data,${SAMPLE_dir}:/SAMPLE_dir $VEPsif \
            filter_vep --format tab --filter "(Feature match NM)"|grep -v "#Uploaded_variation" \
                ) \
                > ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated-filtered-PlusMissedTranscriptVariants.vcf
        }
    #F7: ColumnChopping
        ColumnChopping(){
            grep -v '##' ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-chopped-annotated-filtered-PlusMissedTranscriptVariants.vcf|\
                        cut -f 1,2,5,7,13,19,20,24,42,43,46,47,48,49,77,86,87,88,90 \
                        > ${SAMPLE_dir}/${sample}.${VType}.${Gpanel}-final.tsv
        }
    #F8: F1+F2+F3+F4+F5+F6+F7 Combinedfun
        Combinedfun(){
            #CNV
                VType="${CNV}"
                PanelChopping
                Compressing
                VepAnnotating
                VepFiltering
                MissedTranscriptAddition
                ColumnChopping
            #SV
                VType="${SV}"
                PanelChopping
                Compressing
                VepAnnotating
                VepFiltering
                MissedTranscriptAddition
                ColumnChopping
            #MERGED
                VType="${MERGED}"
                VariantMerging
                Compressing
                VepAnnotating
                VepFiltering
                MissedTranscriptAddition
                ColumnChopping
        }
#_____________________________________________________    
#1.
    #SampleSheet DONE
#_____________________________________________________
#2.Making Sample Directories and Copying respective files form Q:/
    while IFS=$'\t' read -r -a SampleSheet; do
        #Variables:
            sample=${SampleSheet[0]}
        #Creating Sample directories
                mkdir -p ${WORKINGDIR}/${sample}/CNV/RAW
                mkdir -p ${WORKINGDIR}/${sample}/CNV/FINAL
            #Sample config
                SAMPLE_dir="${WORKINGDIR}/${sample}/CNV/RAW"
            #COPY    
                cp ${WORKINGDIR}/Samples/${sample}* ${SAMPLE_dir}
    done < ${WORKINGDIR}/SampleSheet
#_____________________________________________________
#3.Chopping/Merging/Indexing/VepAnnotating/Filtering/AddingMissingTranscripts/ColumnChopping of CNV,SV &Merged for FocusGene/GenePanel
    #While looping for all the samples in run
        while IFS=$'\t' read -r SampleSheet; do
            #Variables:
                sample=${SampleSheet[0]}
                Gpanel=${SampleSheet[1]}
            #Sample config
                SAMPLE_dir="${WORKINGDIR}/${sample}/CNV/RAW"
                echo $sample
            #running all the samples
            VType=${CNV}
            echo $VType
                PanelChopping
                #Combinedfun
            VType=${SV}
            echo $VType
                PanelChopping
                #Combinedfun
            #VType=${MERGED}
                #Combinedfun
        done < ${WORKINGDIR}/SampleSheet
#_____________________________________________________