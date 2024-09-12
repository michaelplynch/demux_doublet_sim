id_ch = Channel.from(1,2,3,4,5,6)
seed_ch = Channel.from(1,2,3,4,5)
ngenes_ch = Channel.from(50,75,100,500,1000,5000,10000)

read_ch = Channel.fromPath(params.bam_path)
barcodes_ch = Channel.fromPath(params.barcodes_path)

params.dir="${projectDir}/data/output"

//params.tenx="/data/projects/yufei/210827_10X_KW9275_bcl/cellranger-6.1.1/GRCh38/BRI-1348/outs/filtered_feature_bc_matrix"
//params.ref="${projectDir}/data/input/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
//params.common_variants="${projectDir}/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf"


params.sim='./params_ccrcc.csv'
sims_ch = Channel.fromPath(params.sim)
    .splitCsv(header: true)
    .map { row ->
            tuple( row.key, row.pcdoub, row.n1, row.n2, row.n3, row.n4, row.n5, row.n6 )
        }
    .view()


log.info """\
    S N P   S I M   P I P E L I N E
    ===================================

    sims: ${params.sim}
    """
    .stripIndent(true)

// Merge to one bam file
process PARSE {
    publishDir "${params.dir}"
    input:
    path read_ch
    val id_ch

    output:
    path "parsed/${read_ch}_parsed.bam"
    
    shell:
    '''
    !{projectDir}/templates/parse_BAM_files.sh "!{read_ch}" "!{read_ch}_parsed.bam" K!{id_ch} parsed
    '''
}

process RENAME_BCS {
    publishDir "${params.dir}"

    input:
    path barcodes_ch
    val id_ch

    output:
    tuple val(id_ch), path("${barcodes_ch}_rn.tsv")

    shell:
    '''
    !{projectDir}/templates/parse_barcodes.sh "!{barcodes_ch}" !{id_ch}
    '''
}

process SUBSET {
    publishDir "${params.dir}"

    input:
    tuple  val(id1), path(h1), val(id2), path(h2), val(id3), path(h3), val(id4), path(h4), val(id5), path(h5), val(id6), path(h6)
    tuple val(key), val(doublets), val(n1), val(n2), val(n3), val(n4), val(n5), val(n6)
    
    output:
    tuple val(key), path("*.tsv")

    shell:
    '''
    Rscript !{projectDir}/templates/subset_barcodes.R !{h1} !{key} !{n1} !{id1} 
    Rscript !{projectDir}/templates/subset_barcodes.R !{h2} !{key} !{n2} !{id2}
    Rscript !{projectDir}/templates/subset_barcodes.R !{h3} !{key} !{n3} !{id3}
    Rscript !{projectDir}/templates/subset_barcodes.R !{h4} !{key} !{n4} !{id4}
    Rscript !{projectDir}/templates/subset_barcodes.R !{h5} !{key} !{n5} !{id5}
    Rscript !{projectDir}/templates/subset_barcodes.R !{h6} !{key} !{n6} !{id6}
    '''
}

process MERGE {
    cpus 4
    publishDir "${params.dir}"
    
    input:
    path parse_ch

    output:
    path "merged/merged.bam"

    shell:
    '''
    !{projectDir}/templates/merge_and_index_BAM.sh "!{parse_ch}" "merged.bam" "merged"
    '''
}

process MERGE_BCS {
    publishDir "${params.dir}"

    input:
    tuple val(key), path(barcodes)

    output:
    tuple val(key), path("barcodes_merged_${key}.tsv")

    shell:
    '''
    cat !{barcodes} > barcodes_merged_!{key}.tsv
    '''
}

process LOOKUP {
    publishDir "${params.dir}"

    input:
    tuple val(key), val(pcdoub), path(barcodes)

    output:
    tuple val(key), path("barcodes_merged_ccrcc_${key}_${pcdoub}pc.tsv"), emit: bcs_merged
    tuple val(key), val(pcdoub), path("lookup_table_doublets_ccrcc_${key}_${pcdoub}pc.tsv"), emit: lookup

    shell:
    '''
    Rscript !{projectDir}/templates/generate_awk_lookup_tables_doublets.R "!{projectDir}/data" "!{projectDir}/data" "!{projectDir}/data" "barcodes_merged_!{key}.tsv" !{pcdoub} !{key}
    '''
}

process SIMDOUB {
    publishDir "${params.dir}"

    input:
    val merge_ch
    tuple val(key), val(doublets), path(lookup_file)
    

    output:
    tuple val(key), val(doublets), path("doublets_ccrcc_${key}_${doublets}pc.bam"), path("doublets_ccrcc_${key}_${doublets}pc.bam.bai")
    
    shell:
    '''
    !{projectDir}/templates/parse_and_index_BAM_doublets.sh "!{projectDir}/data" "!{projectDir}/data" !{doublets} !{key} !{lookup_file} !{merge_ch}
    '''
}

// Run algorithms on sim data

process SOUP {
    cache true
    cpus 10
    publishDir "${params.dir}"
    module 'apptainer'
    container "${params.container__souporcell}"
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(lookup)

    output:
    tuple val(key), path("soup_${key}_${doublets}_1")

    shell:
    '''
    start=`date +%s`
     souporcell_pipeline.py \
        -i !{bam} \
        -b !{lookup} \
        -f !{params.ref} \
        -t 10 \
        -o "soup_!{key}_!{doublets}_1" \
        -k 6 \
        --common_variants !{params.common_variants} \
        --skip_remap SKIP_REMAP
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > soup_!{key}_!{doublets}_1/runtime.txt
    '''
}

process SOUPNGENES {
    cache true
    cpus 10
    publishDir "${params.dir}"
    module 'apptainer'
    container "${params.container__souporcell}"
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(lookup), val(ngenes)
    val(vcf)    
    val(workflow)

    output:
    tuple val(key), path("soup_${key}_${doublets}_1_${workflow}_${ngenes}")

    shell:
    '''
    start=`date +%s`
     souporcell_pipeline.py \
        -i !{bam} \
        -b !{lookup} \
        -f !{params.ref} \
        -t 10 \
        -o "soup_!{key}_!{doublets}_1" \
        -k 6 \
        --common_variants !{params.common_variants} \
        --skip_remap SKIP_REMAP
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > soup_!{key}_!{doublets}_1/runtime.txt
    cp -r soup_!{key}_!{doublets}_1 soup_!{key}_!{doublets}_1_!{workflow}_!{ngenes}
    '''
}

process DEMUXSNP {
    cache false
    cpus 10
    publishDir "${params.dir}"
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(barcodes), val(pcdoub), path(lookup), path(merged_bcs), val(seed), path(souporcell)
    val(ngenes)
    output:
    path "demuxSNP_${key}_${doublets}_${seed}"

    shell:
    '''
    start=`date +%s`
    Rscript !{projectDir}/templates/sce_sim2.R !{params.tenx} !{merged_bcs} !{lookup} !{key} 
    Rscript !{projectDir}/templates/demuxSNP_nf_pre.R !{souporcell}/common_variants_covered.vcf "sce.rdata" !{key} !{doublets} !{seed} !{ngenes}
    vartrix_linux --bam !{bam} \
        --cell-barcodes !{barcodes} \
        --fasta !{params.ref} \
        --threads 10 \
        --out-matrix demuxSNP_!{key}_!{doublets}_!{seed}/out_matrix.mtx \
        --vcf "!{souporcell}/common_variants_covered.vcf"
    Rscript !{projectDir}/templates/demuxSNP_nf_post.R sce.rdata demuxSNP_!{key}_!{doublets}_!{seed}/out_matrix.mtx !{barcodes} !{key} !{doublets} !{seed} !{souporcell}
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > demuxSNP_!{key}_!{doublets}_!{seed}/runtime.txt
    '''
}

process DEMUXSNPNGENES {
    cache true
    cpus 10
    publishDir "${params.dir}"
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(barcodes), val(pcdoub), path(lookup), path(merged_bcs), val(seed), path(souporcell), val(ngenes)
    val(workflow)
    output:
    path "demuxSNP_${key}_${doublets}_${seed}_${workflow}_${ngenes}"
    path "demuxSNP_${key}_${doublets}_${seed}_${workflow}_${ngenes}/vcf_sub.vcf", emit: vcf
    
    shell:
    '''
    start=`date +%s`
    Rscript !{projectDir}/templates/sce_sim2.R !{params.tenx} !{merged_bcs} !{lookup} !{key} 
    Rscript !{projectDir}/templates/demuxSNP_nf_pre.R !{souporcell}/common_variants_covered.vcf "sce.rdata" !{key} !{doublets} !{seed} !{ngenes}
    vartrix_linux --bam !{bam} \
        --cell-barcodes !{barcodes} \
        --fasta !{params.ref} \
        --threads 10 \
        --out-matrix demuxSNP_!{key}_!{doublets}_!{seed}/out_matrix.mtx \
        --vcf "demuxSNP_!{key}_!{doublets}_!{seed}/vcf_sub.vcf"
    Rscript !{projectDir}/templates/demuxSNP_nf_post.R sce.rdata demuxSNP_!{key}_!{doublets}_!{seed}/out_matrix.mtx !{barcodes} !{key} !{doublets} !{seed} !{souporcell}
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > demuxSNP_!{key}_!{doublets}_!{seed}/runtime.txt
    cp -r demuxSNP_!{key}_!{doublets}_!{seed} demuxSNP_!{key}_!{doublets}_!{seed}_!{workflow}_!{ngenes}
    '''
}

process VIREO {
    cache true
    cpus 5
    publishDir "${params.dir}"
    module 'conda'
    conda '/cm/shared/apps/conda/202308v2/envs/cellsnp'
    
    input:
    tuple val(key), val(doublets), path(bam), path(bai), path(lookup), val(seed)

    output:
    path "vireo_${key}_${doublets}_${seed}"

    shell:
    '''
    start=`date +%s`
    cellsnp-lite -s !{bam} \
        -O cellsnp_out \
        -R "!{projectDir}/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.vcf" \
        -b !{lookup} \
        -p 5 \
        --gzip
    vireo -c cellsnp_out \
        -N 6 \
        -o "vireo_!{key}_!{doublets}_!{seed}" \
        -p 5 \
        --randSeed !{seed}
    end=`date +%s`
    runtime=$((end-start))
    echo "$runtime" > vireo_!{key}_!{doublets}_!{seed}/runtime.txt
    '''
}

workflow {
    // parse and merge bams
    parse_ch = PARSE(read_ch, id_ch)
    merge_ch = MERGE(parse_ch.collect())

    // process and subset barcodes + create lookup
    rename_bcs_ch = RENAME_BCS(barcodes_ch,id_ch)
    rename_bcs_ch.toSortedList({ a, b -> a[0] <=> b[0] }).flatten().collect().view()
    rename_bcs_ch.toSortedList({ a, b -> a[0] <=> b[0] }).flatten().collect().set { barcodes_list_ch}
    
    //barcodes_list_ch.view()
    SUBSET(barcodes_list_ch,sims_ch)
    //SUBSET.out.view()
    merged_bcs_ch=MERGE_BCS(SUBSET.out)
    sims_ch.map {key,pcdoub,n1,n2,n3,n4,n5,n6 ->
            tuple(key,pcdoub)}
        .combine(merged_bcs_ch, by: 0)
        .set {doub_bcs_ch}
    LOOKUP(doub_bcs_ch)
    //LOOKUP.out.lookup.view()

    // simulate doublets
    SIMDOUB(merge_ch, LOOKUP.out.lookup)
    test_ch = SIMDOUB.out
        .combine(LOOKUP.out.bcs_merged,by:0)
    runseeds_ch=test_ch.combine(seed_ch)
    //test_ch.view()
    //runseeds_ch.view()

    // run souporcell and demuxSNP
    SOUP(test_ch)
    //VIREO(runseeds_ch)
    demuxsnp_ch = test_ch.combine(LOOKUP.out.lookup,by:0).combine(merged_bcs_ch,by:0).combine(seed_ch).combine(SOUP.out,by:0)
    demuxsnp_ch.view()
    //SOUP.out.view()

    // Run demuxSNP with different # genes
    DEMUXSNP(demuxsnp_ch,500)
    DEMUXSNPNGENES(demuxsnp_ch.filter{ it[0].matches(/^(e).*/) && it[8] == 1 }.combine(ngenes_ch),"ngenes")

    SOUPNGENES(test_ch.filter{ it[0].matches(/^(e).*/)}.combine(ngenes_ch),DEMUXSNPNGENES.out.vcf,'ngenes')
}
