//this should consume fastq files, adapter trim, then align reads to the grch38 reference using multiple aligment tools
//it assumes that all input files end with .[12].fastq.gz

// optionally --filter can be used to specify a POSIX regex matching samples to process (e.g. '(NoBS|WG).*)', default= all samples

params.outdir = 'output'
params.genome = 'genome/grch38_core+bs_controls.fa'
params.fastq_dir = '.'
params.aligner_cpus = 16
flowcell="find ${run_dir} -name '*.1.fastq*' | head -n 1 | xargs zcat | head -n 1 | cut -f 3 -d ':'".execute().text
params.filter = '.*'

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

    genome_dir = File.dirname(params.genome)/

    shell:
    '''
        mkdir -p !{genome_dir}/bwameth \
          && cd !{genome_dir}/bwameth \
          && ln -s !{params.genome} genome.fa \
          && bwameth.py --index genome.fa
    '''
}

process prepare_metadata_for_fastqs {

    shell:
    '''
    for read1_fq in `find !{params.fastq_dir} -maxdepth 1 -name '*1.fastq*'`; do

        if ! [[ $read1_fq =~ !{filter} ]]; then
            continue;
        fi
        dirname=`dirname $read1_fq`
        library=`basename $(basename $read1_fq .1.fastq) .1.fastq.gz`
        tile=all
        header=`zcat -f $read1_fq | head -n 1`
        flowcell=`echo $header | cut -f 3 -d ':'`
        lane=`echo $header | head -n 1 | cut -f 4 -d ':'`
        read2_fq=$(ls ${dirname}/*${library}.2.*);

        if [ -f $read1_fq ] && [ -f $read2_fq ]; then
            printf "${library},${flowcell},${lane},${tile},${read1_fq},${read2_fq}\n">>'tile_info.csv'
        else
            echo "Unable to find $read1_fq or $read2_fq"
            exit 255
        fi
    done
    '''

}

fastq_files
    .splitCsv(header: ['library', 'flowcell', 'lane','tile', 'read1', 'read2'])
    .into {bwameth;  bitmapperbs; bismark; gembs;}

process bwameth {
    cpus params.aligner_cpus
    errorStrategy 'retry'

    tag { [flowcell, library, tile] }

    conda 'bwameth=0.2.2 fastp=0.20.1 seqtk=1.3'

    input:
        set library, flowcell, lane, tile, read1, read2  from bwameth

    output:
        set val(library), file("*.aln.bam") into aligned_files

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

    seqtk mergepe <(zcat -f !{fq_set.insert_read1}) <(zcat -f !{fq_set.insert_read2}) \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j "!{fq_set.library}_combined_fastp.json" \
    | bwameth.py -p -t !{task.cpus} --read-group "@RG\\tID:!{fq_set.barcode}\\tSM:!{fq_set.library}\\tBC:!{fq_set.barcode}\\tCN:NEB\\tPL:ILLUMINA\\tPU:!{fq_set.flowcell}-!{fq_set.barcode}.!{fq_set.lane}" --reference !{genome} /dev/stdin 2>  "!{fq_set.library}_!{fq_set.barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.log.bwamem" \
    | sambamba view -t 2 -S -f bam -o "!{fq_set.library}_!{fq_set.barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.aln.bam" /dev/stdin;

    '''

}

process bwameth_markdups {
    cpus 8
    errorStrategy 'retry'
    tag { library }
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.{md.bam}*'
    conda "samtools=1.9 samblaster=0.1.24 sambamba=0.7.0"

    input:
        set val(library), file(libraryBam) from aligned_bams.groupTuple()

    output:
        set val(library), file('*.md.bam'), file('*.md.bam.bai') into bwameth_md_bams

    shell:
    '''
    samtools cat  -b <( find . -name '*.aln.bam' ) \
    | samtools view -h /dev/stdin \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=!{tmp_dir} -t !{task.cpus} -m 20GB -o !{library}.md.bam /dev/stdin

    '''
}

process gembs {

}

process gembs_markdups {

}

process bitmapperbs {

}

process bitmapperbs_markdups {

}

process bismark {

}

process bismark_markdups {

}


process summarize_results {

}


