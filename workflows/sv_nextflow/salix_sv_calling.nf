#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reference      = '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Reference/Salix_viminalis.fasta'
params.bam_dir        = '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Mapping_output/bam_salix_all/'
params.chr_list       = '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/chr_list.txt'
params.threads        = 8
params.pyenv_v        = ''
params.overlap_size   = 0.8
params.ref_bed_awk    = '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/make_reference_bed.awk'

params.conda_init     = '/cfs/klemming/home/s/saidash/miniconda3/etc/profile.d/conda.sh'
params.survivor_env   = 'survivor'


process makeExcludeBed {

    module 'bcftools'

    output:
    path 'exclude.bed', emit: exclude_bed
    path 'reference.bed', emit: reference_bed

    shell:
    '''
    awk -f !{params.ref_bed_awk} !{params.reference} > reference.bed
    awk 'NR==FNR {keep[$1]=1; next} !keep[$1]' !{params.chr_list} reference.bed > exclude.bed
    '''
}


process dellyCallRaw {

    module 'delly'
    module 'bcftools'
    cpus params.threads

    tag { "${sample}:${sv_type}" }

    input:
    tuple val(sample), path(bam), path(bai), val(sv_type)
    path exclude_bed

    output:
    tuple val(sv_type), path("${sample}_${sv_type}_filter.bcf"), emit: filtered_bcf_by_type

    shell:
    '''
    set -euo pipefail
    export OMP_NUM_THREADS=!{task.cpus}

    delly call -t !{sv_type} \
      -g !{params.reference} \
      -x !{exclude_bed} \
      -o "!{sample}_!{sv_type}.bcf" \
      !{bam}

    bcftools view -O b -f PASS "!{sample}_!{sv_type}.bcf" > "!{sample}_!{sv_type}_filter.bcf"
    '''
}


process dellyMergeSites {

    module 'delly'
    cpus params.threads

    tag { sv_type }
    stageInMode 'copy'

    input:
    tuple val(sv_type), path(bcfs, stageAs: 'bcf_*')

    output:
    tuple val(sv_type), path("sites_${sv_type}.bcf"), emit: sites_bcf

    shell:
    '''
    set -euo pipefail
    ls -1 bcf_* | sort > inputs.list
    delly merge -o "sites_!{sv_type}.bcf" $(cat inputs.list)
    '''
}


process dellyGenotype {

    module 'delly'
    module 'bcftools'
    cpus params.threads

    tag { "${sample}:${sv_type}" }

    input:
    tuple val(sv_type), path(sites_bcf), val(sample), path(bam), path(bai)
    path exclude_bed

    output:
    tuple val(sv_type), path("${sample}_${sv_type}.geno.bcf"), emit: geno_bcf_by_type
    tuple val(sv_type), path("${sample}_${sv_type}.geno.bcf.csi"), emit: geno_bcf_by_type_index

    shell:
    '''
    set -euo pipefail
    export OMP_NUM_THREADS=!{task.cpus}

    delly call -t !{sv_type} \
      -g !{params.reference} \
      -x !{exclude_bed} \
      -v !{sites_bcf} \
      -o "!{sample}_!{sv_type}.geno.bcf" \
      !{bam}

    bcftools index -f -c "!{sample}_!{sv_type}.geno.bcf"
    '''
}


process dellyMergeGenos {

    module 'bcftools'
    cpus params.threads

    tag { sv_type }
    stageInMode 'copy'

    input:
    tuple val(sv_type), path(geno_bcfs)

    output:
    tuple val(sv_type), path("merged_${sv_type}.bcf"), emit: merged_bcf
    path("merged_${sv_type}.bcf.csi"), emit: merged_bcf_csi, optional: true

    shell:
    '''
    set -euo pipefail

    printf "%s\n" !{geno_bcfs} | sort > inputs.list

    while read -r f; do
        bcftools index -f -c "$f"
    done < inputs.list

    bcftools merge -m id -O b -o "merged_!{sv_type}.bcf" $(cat inputs.list)
    bcftools index -f -c "merged_!{sv_type}.bcf"
    '''
}


process dellyConcat {

    module 'bcftools'
    cpus params.threads

    publishDir 'vcf_output', mode: 'copy'

    input:
    path merged_list

    output:
    path 'delly.vcf', emit: delly_vcf

    shell:
    '''
    set -euo pipefail

    echo "=== merged.list ===" >&2
    cat !{merged_list} >&2

    while read -r f; do
        [ -f "$f.csi" ] || bcftools index -f -c "$f"
    done < !{merged_list}

    bcftools concat -a -f !{merged_list} -O v -o delly.vcf
    '''
}


process dysguRunAndFilter {

    module 'cray-python/3.11.7'
    module 'bcftools'
    cpus params.threads

    tag { sample }

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}_filter.vcf", emit: filtered_vcf

    shell:
    '''
    set -euo pipefail

    dysgu run -p!{params.threads} !{params.reference} tmpdir !{bam} > "!{sample}.vcf"
    rm -rf tmpdir
    bcftools view -i 'FILTER="PASS"' "!{sample}.vcf" > "!{sample}_filter.vcf"
    '''
}


process dysguMerge {
    module 'cray-python/3.11.7'
    cpus params.threads

    publishDir 'vcf_output', mode: 'copy'

    input:
    path filtered_vcfs

    output:
    path 'dysgu.vcf', emit: dysgu_vcf

    shell:
    '''
    set -euo pipefail

    printf "%s\n" !{filtered_vcfs} > inputs.list
    mkdir -p wd

    dysgu merge -p8 --input-list inputs.list --wd wd --svs-out dysgu.vcf
    test -s dysgu.vcf
    '''
}


process mantaRun {
    module 'manta/1.6.0'
    module 'samtools/1.20'
    module 'bcftools'
    cpus params.threads

    publishDir 'vcf_output/manta_per_sample', mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.manta.vcf.gz", emit: manta_vcfs
    path "${sample}.manta.vcf.gz.tbi", emit: manta_vcfs_tbi

    shell:
    """
    set -euo pipefail

    echo "configManta.py: \$(which configManta.py)" >&2
    echo "convertInversion.py: \$(which convertInversion.py)" >&2
    echo "samtools: \$(which samtools)" >&2

    CHR_ARGS=\$(awk '{print "--region "\$0}' ${params.chr_list} | tr '\\n' ' ')

    configManta.py --referenceFasta ${params.reference} --bam ${bam} \$CHR_ARGS --runDir MantaWorkflow

    cd MantaWorkflow
    python2 runWorkflow.py -m local -j ${params.threads}

    SAMTOOLS=\$(which samtools)
    convertInversion.py "\$SAMTOOLS" ${params.reference} results/variants/diploidSV.vcf.gz \
      | bcftools sort -O z -o ../${sample}.manta.vcf.gz

    bcftools index -f -t ../${sample}.manta.vcf.gz
    """
}


process mantaMergeSurvivor {
    module 'bcftools'
    cpus params.threads
    publishDir 'vcf_output', mode: 'copy'

    input:
    path vcfs

    output:
    path 'manta.vcf', emit: manta_vcf

    shell:
    '''
    set -euo pipefail

    source !{params.conda_init}
    conda activate !{params.survivor_env}

    : > vcfs.list
    for f in !{vcfs}; do
        echo "$f" >> vcfs.list
    done

    SURVIVOR merge vcfs.list 500 1 1 1 0 50 manta.raw.vcf
    bcftools sort -O v -o manta.vcf manta.raw.vcf
    '''
}


process overlap_vcf {

    module 'bcftools'
    module 'bedtools'

    publishDir 'vcf_output', mode: 'copy'

    input:
    path delly_vcf
    path dysgu_vcf
    path manta_vcf

    output:
    path 'manta_overlap.vcf'

    shell:
    '''
    [[ ! -e delly.vcf ]] && cp !{delly_vcf} delly.vcf || true
    [[ ! -e dysgu.vcf ]] && cp !{dysgu_vcf} dysgu.vcf || true
    [[ ! -e manta.vcf ]] && cp !{manta_vcf} manta.vcf || true


    for file in *.vcf; do
        for sv_type in {DEL,DUP,INV}; do
            filename=$(basename $file .vcf)
            bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\n' $file \
                | awk '($5 == "'"$sv_type"'" && $6 == "PASS") || ($5 == "'"$sv_type"'" && $6 == ".")' \
                | awk -v OFS='\t' '{ if ($3 != "." && $3 < $2) { t=$2; $2=$3; $3=t } print }' \
                | sort -k1,1 -k2,2n \
                > "$filename"_"$sv_type".bed
        done
    done

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\t%INFO/INSLEN\n' delly.vcf \
        | awk '($5 == "INS" && $6 == "PASS") || ($5 == "INS" && $6 == ".")' \
        | awk '{print $0, $2+$7}' \
        | awk -v OFS='\t' '{print $1,$2,$8,$4,$5,$6,$7}' \
        | awk -v OFS='\t' '{ if ($3 != "." && $3 < $2) { t=$2; $2=$3; $3=t } print }' \
        | sort -k1,1 -k2,2n \
        > delly_INS.bed

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\t%INFO/SVLEN\n' dysgu.vcf \
        | awk '($5 == "INS" && $6 == "PASS") || ($5 == "INS" && $6 == ".")' \
        | awk '{print $0, $2+$7}' \
        | awk -v OFS='\t' '{print $1,$2,$8,$4,$5,$6,$7}' \
        | awk -v OFS='\t' '{ if ($3 != "." && $3 < $2) { t=$2; $2=$3; $3=t } print }' \
        | sort -k1,1 -k2,2n \
        > dysgu_INS.bed

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\t%INFO/SVLEN\n' manta.vcf \
        | awk '($5 == "INS" && $6 == "PASS") || ($5 == "INS" && $6 == ".")' \
        | awk '{print $0, $2+$7}' \
        | awk -v OFS='\t' '{print $1,$2,$8,$4,$5,$6,$7}' \
        | awk '$7 != "."' \
        | awk -v OFS='\t' '{ if ($3 != "." && $3 < $2) { t=$2; $2=$3; $3=t } print }' \
        | sort -k1,1 -k2,2n \
        > manta_INS.bed

    for sv_type in {DEL,DUP,INS,INV}; do
        bedtools multiinter -i dysgu_"$sv_type".bed delly_"$sv_type".bed manta_"$sv_type".bed > common_overlap_"$sv_type".bed
        bedtools intersect -a common_overlap_"$sv_type".bed -b dysgu_"$sv_type".bed delly_"$sv_type".bed manta_"$sv_type".bed \
            -f !{params.overlap_size} -r -wa -wb > common_overlap_with_ID_"$sv_type".bed
        bedtools groupby -g 1-8 -c 9 -o count_distinct -i common_overlap_with_ID_"$sv_type".bed > distinct_overlap_"$sv_type".bed
        awk '$9==3' distinct_overlap_"$sv_type".bed | awk -v OFS='\t' '{print $1,$2,$3}' > for_matching_"$sv_type".txt
        grep -F -f for_matching_"$sv_type".txt common_overlap_with_ID_"$sv_type".bed | awk '$9=="3"' | awk '{print $13}' > ID_list_"$sv_type".txt
    done

    cat ID_list_*.txt > ID_list.txt
    bcftools view -i'ID=@ID_list.txt' manta.vcf > manta_overlap.vcf
    '''
}


workflow {

    Channel
        .fromPath("${params.bam_dir}/*.bam")
        .ifEmpty { error "No BAMs found in: ${params.bam_dir}" }
        .map { bam ->
            def bai = file("${bam}.bai")
            if( !bai.exists() ) error "Missing BAM index for ${bam} (expected: ${bam}.bai)"
            tuple(bam.baseName, bam, bai)
        }
        .set { bam_ch }

    Channel
        .of('DEL','DUP','INS','INV')
        .set { svtype_ch }

    makeExcludeBed()

    def bam_x_sv = bam_ch.combine(svtype_ch)

    dellyCallRaw(bam_x_sv, makeExcludeBed.out.exclude_bed)

    dellyMergeSites(
        dellyCallRaw.out.filtered_bcf_by_type.groupTuple(by: 0)
    )

    def sites_x_bam = dellyMergeSites.out.sites_bcf.combine(bam_ch)

    dellyGenotype(sites_x_bam, makeExcludeBed.out.exclude_bed)

    dellyMergeGenos(
        dellyGenotype.out.geno_bcf_by_type
            .groupTuple(by: 0)
    )

    def order = ['DEL','DUP','INS','INV']

    def merged_list_ch = dellyMergeGenos.out.merged_bcf
      .map { svt, bcf -> [ svt: svt.toString(), bcf: bcf ] }
      .collect()
      .map { recs ->

        recs = recs.sort { a, b -> order.indexOf(a.svt) <=> order.indexOf(b.svt) }

        def text = recs.collect { it.bcf.toString() }.join('\n') + '\n'

        def f = file("merged_${workflow.runName}.list")
        f.text = text

        println "=== merged.list content ===\n${text}"
        return f
     }

    dellyConcat( merged_list_ch )

    dysguRunAndFilter(bam_ch)
    dysguMerge(dysguRunAndFilter.out.filtered_vcf.collect())

    mantaRun( bam_ch )
    mantaMergeSurvivor( mantaRun.out.manta_vcfs.collect() )

    overlap_vcf(
        dellyConcat.out.delly_vcf,
        dysguMerge.out.dysgu_vcf,
        mantaMergeSurvivor.out.manta_vcf
    )
}
