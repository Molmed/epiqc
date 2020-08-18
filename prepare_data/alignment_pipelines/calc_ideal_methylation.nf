//this should tries to make an ideal but realistic methylation dataset
// steps:
// - proces unconverted fastq files (*.[12].fastq.gz), adapter trim, then align reads to 
//   the grch38 reference (downloaded by the tool) with bwa-mem
// - process input bed graphs to calculate an aggregate methylation frequency for each C
// - in-silico convert bwa-mem aligned reads 
// - re-map the converted reads to a converted reference using bwa-meth

// optionally --fastq_filter or bedgraph_filter can be used to specify a 
// regex matching samples to process (e.g. '(NoBS|WG).*)', default= all samples

params.outdir = 'output'
params.genome = 'genome/grch38_core+bs_controls.fa' //todo - move to zenodo?
params.unconverted_ = '.'
params.bedgraph_dir = '.'
params.aligner_cpus = 16
flowcell="find ${run_dir} -name '*.1.fastq*' | head -n 1 | xargs zcat | head -n 1 | cut -f 3 -d ':'".execute().text
params.fastq_filter = '.*'
params.bedgraph_filter = '^(?:(?!_ox.).)*$' //not ox libs from ox-bs libs

process cache_genome {
    publishDir File.dirname(params.genome)

    output: 
        val(file("*.fa")) into grch38_ref

    when:
        ! File.exists(params.genome)

    shell:
    '''
        if -f !{params.genome} ;

        else
            wget https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa 
            cp grch38_core+bs_controls.fa !{params.genome}
        fi
    '''
}

process bwameth_index {

    input: file(grch38_ref) //this is only used to cause this step to wait for download to complete

    genome_dir = File.dirname(params.genome)

    when:
        ! File.exists(params.genome + "/bwameth")

    conda 'bwameth=0.2.2'
    shell:
    '''
        mkdir -p !{genome_dir}/bwameth \
          && cd !{genome_dir}/bwameth \
          && ln -s !{params.genome} genome.fa \
          && bwameth.py --index genome.fa
    '''
}

process bwa_index {

    input: file(grch38_ref) //this is only used to cause this step to wait for download to complete

    genome_dir = File.dirname(params.genome)

    conda 'bwa=0.7.17'

    when:
        ! File.exists(params.genome + "/bwa-mem")

    shell:
    '''
        mkdir -p !{genome_dir}/bwa-mem \
          && cd !{genome_dir}/bwa-mem \
          && ln -s !{params.genome} genome.fa \
          && bwa index genome.fa
    '''
}

process prepare_metadata_for_bedGraphs {
    output: 'bedgraph_lib_info.csv' into bedgraph_files
    shell:
    '''
    for bg in `find !{params.bedgraph_dir} -maxdepth 1 -name '*.bedGraph.gz'`; do
        if ! [[ $bg =~ !{bedgraph_filter} ]]; then
            continue;
        fi
        dirname=`dirname $bg`
        library=`basename $bg sorted.merged.bedGraph.gz`
        method=`echo $library | cut -f 1 -d _`
        cell_line=`echo $library |cut -f 2 -d _`
        lab=`echo $library |cut -f 3 -d _`
        printf "${library},${method},${lab},${cell_line},${bg}\n">>'bedgraph_lib_info.csv'
    done
    '''
}

bedgraph_files
    .splitCsv(header: ['library', 'method', 'lab','cell_line', 'bedgraph'])
    .into {merge_bedgraphs}

//merge bedgraphs by cell line and recalculate % meth for all libs
process merge_bedgraphs {
    tag 
    input: 
    set library, method, lab, cell_line, flowcell, file(bedgraph) from merge_bedgraphs
}


process prepare_metadata_for_fastqs {
    output: 'fastq_lib_info.csv' into fastq_files

    shell:
    '''
    for read1_fq in `find !{params.fastq_dir} -maxdepth 1 -name '*1.fastq*'`; do

        if ! [[ $read1_fq =~ !{filter} ]]; then
            continue;
        fi
        dirname=`dirname $read1_fq`
        library=`basename $(basename $read1_fq .1.fastq) .1.fastq.gz`
        method=`echo $library | cut -f 1 -d _`
        cell_line=`echo $library |cut -f 2 -d _`
        lab=`echo $library |cut -f 3 -d _`
        tile=all
        header=`zcat -f $read1_fq | head -n 1`
        flowcell=`echo $header | cut -f 3 -d ':'`
        lane=`echo $header | head -n 1 | cut -f 4 -d ':'`
        read2_fq=$(ls ${dirname}/*${library}.2.*);

        if [ -f $read1_fq ] && [ -f $read2_fq ]; then
            printf "${library},${method},${lab},${cell_line},${read1_fq},${read2_fq}\n">>'fastq_lib_info.csv'
        else
            echo "Unable to find $read1_fq or $read2_fq"
            exit 255
        fi
    done
    '''

}

fastq_files
    .splitCsv(header: ['library', 'method', 'lab','cell_line', 'flowcell', 'read1', 'read2'])
    .into {fastq_for_bwamem}

process bwamem_md {
    cpus params.aligner_cpus
    errorStrategy 'retry'

    tag { [method, lab, cell_line] }

    conda 'bwa=0.7.17 fastp=0.20.1 seqtk=1.3 samblaster=0.1.24'

    input:
        set library, method, lab, cell_line, flowcell, read1, read2 from fastq_for_bwamem

    output:
        set val(library), method, lab, cell_line, flowcell, file("*.md.bam") into bwa-mem_aligned_files

    shell:
    '''
    inst_name=$(head -n 1 < <(zcat -f !{read1}) | cut -f 1 -d ':' | sed 's/^@//')
    fastq_barcode=$(head -n 1 < <(zcat -f !{read1}) | sed -r 's/.*://')
    if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || [[ "${inst_name:0:2}" == 'NB' ]]; then
        #these 2-color instruments should have poly-g stretches trimmed since they are likely indicative of inactive clusters
        trim_polyg='--trim_poly_g'
    else
        trim_polyg=''
    fi

    seqtk mergepe <(zcat -f !{read1}) <(zcat -f !{read2}) \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j "!{library}_combined_fastp.json" \
    | bwa-mem -p -t !{task.cpus} -R \"@RG\\\\tID:${fastq_barcode}\\\\tSM:!{library}\\\\tBC:${fastq_barcode}\\\\tCN:Multiple\\\\tPL:ILLUMINA\\\\tPU:!{flowcell}-${fastq_barcode}" --reference !{genome} /dev/stdin 2>  "!{library}_${fastq_barcode}!{flowcell}.log.bwamem" \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -S -f bam -o "!{library}_${fastq_barcode}.bwa.md.bam" /dev/stdin;
    '''

}

process insilico_convert {

    conda 'pysam'

    input: 
        tuple library, method, lab, cell_line, flowcell, file(bam) from bwa-mem_aligned_files
    
    output:
        tuple library, method, lab, cell_line, flowcell, file("*.1.fastq.gz"), file("*.2.fastq.gz") into fastq_for_bwameth

    shell:
    '''
        echo "magic tool here..."
    '''
}

process bwameth_md {
    cpus params.aligner_cpus
    errorStrategy 'retry'
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.{md.bam}*'

    tag { [method, lab, cell_line] }

    conda 'bwameth=0.2.2 fastp=0.20.1 seqtk=1.3 samblaster=0.1.24 '

    input:
         set library, method, lab, cell_line, flowcell, read1, read2 from fastq_for_bwameth

    output:
        set val(library), file("*.aln.bam") into bwameth_aligned_files

    shell:
    '''
    inst_name=$(head -n 1 < <(zcat -f !{fq_set.insert_read1}) | cut -f 1 -d ':' | sed 's/^@//')
    fastq_barcode=$(head -n 1 < <(zcat -f !{fq_set.insert_read1}) | sed -r 's/.*://')
    if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || [[ "${inst_name:0:2}" == 'NB' ]]; then
        #these 2-color instruments should have poly-g stretches trimmed since they are likely indicative of inactive clusters
        trim_polyg='--trim_poly_g'
    else
        trim_polyg=''
    fi

    seqtk mergepe <(zcat -f !{read1}) <(zcat -f !{read2}) \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j "!{library}_combined_fastp.json" \
    | bwameth.py -p -t !{task.cpus} \
       --read-group  \"@RG\\\\tID:${fastq_barcode}\\\\tSM:!{library}\\\\tBC:${fastq_barcode}\\\\tCN:Multiple\\\\tPL:ILLUMINA\\\\tPU:!{flowcell}-${fastq_barcode}" \
       --reference !{genome} /dev/stdin 2>  "!{library}_${fastq_barcode}!{flowcell}.log.bwamem" \
    | samblaster 2> !{library}.log.samblaster
    | sambamba view -t 2 -S -f bam -o "!{library}_${fastq_barcode}.bwameth.md.bam" /dev/stdin;
    '''
}

//todo: methyldackel

