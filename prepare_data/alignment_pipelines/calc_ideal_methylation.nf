// this should construct an idealize but realistic methylation dataset 
// from a standard DNA library and a set of methylation data (bedgraphs)
// steps:
// - proces unconverted fastq files (*.[12].fastq.gz), adapter trim, then align reads to 
//   the grch38 reference (downloaded by the tool) with bwa-mem
// - process input bed graphs to calculate an aggregate methylation frequency for each C
// - in-silico convert bwa-mem aligned reads 
// - re-map the converted reads to a converted reference using bwa-meth


params.outdir = 'output'
params.bwameth_genome = file('genome/grch38_core+bs_controls.fa').toAbsolutePath() //todo - move to zenodo?
params.bwamem_genome = file('/mnt/galaxy/data/genome/grch38_full/bwa/GRCh38_full_analysis_set_plus_decoy_hla.fa').toAbsolutePath() //todo - move to zenodo?
params.unconverted_fastqs_path = './*.R{1,2}.fastq.*'
params.bedgraphs_path = './*.bedGraph.gz'
params.aligner_cpus = 16

bedgraph_files = Channel.fromPath(params.bedgraphs_path)
unconverted_fastq_files = Channel.fromFilePairs(params.unconverted_fastqs_path)

process cache_genome {
    publishDir file(params.bwameth_genome).parent

    output: 
        val(file("*.fa")) into grch38_ref

    when:
        ! params.bwameth_genome.exists()

    shell:
    '''
        if [ -f !{params.bwameth_genome} ]; then
            echo "already present"
        else
            wget https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa 
            mkdir -p $(dirname !{params.bwameth_genome})
            cp grch38_core+bs_controls.fa !{params.bwameth_genome}
        fi
    '''
}

process bwameth_index {

    input: file(grch38_ref) //this is only used to cause this step to wait for the genome download to complete

    genome_dir = file(params.bwameth_genome).parent

    when:
        ! (params.bwameth_genome + "/bwameth").exists()

    conda 'bwameth=0.2.2'

    shell:
    '''
        mkdir -p !{genome_dir}/bwameth \
          && cd !{genome_dir}/bwameth \
          && ln -s !{params.bwameth_genome} genome.fa \
          && bwameth.py --index genome.fa
    '''
}

process bwa_index {

    input: file(grch38_ref) //this is only used to cause this step to wait for download to complete

    genome_dir = file(params.bwameth_genome).parent

    conda 'bwa=0.7.17'

    when:
        ! file(params.bwameth_genome + "/bwa-mem").exists

    shell:
    '''
        mkdir -p !{genome_dir}/bwa-mem \
          && cd !{genome_dir}/bwa-mem \
          && ln -s !{params.bwameth_genome} genome.fa \
          && bwa index genome.fa
    '''
}

process prepare_metadata_for_bedGraphs {
    input: file(bedgraph) from bedgraph_files.toList()
    output: file('bedgraph_lib_info.csv') into bedgraph_info

    shell:
    '''
    for bg in `find "$(pwd)" -maxdepth 1 -name '*.bedGraph.gz'`; do
        dirname=`dirname $bg`
        cell_line=`echo $(basename $bg) | cut -f 2 -d _`
        printf "${cell_line},${bg}\n" >> 'bedgraph_lib_info.csv'
    done
    '''
}

bedgraph_info.splitCsv()
             .set{library_bedgraphs}


// merge bedgraphs by cell line and recalculate % meth for all libs
process merge_bedgraphs {
    cpus 1
    tag {cell_line}

    input: 
        tuple cell_line, file(bed) from library_bedgraphs.groupTuple()
    output: 
    file("combined_bedgraph") into combined_bedgraphs

    conda 'bedtools=2.29.2 pigz=2.3.4 parallel=20200722 sed=4.8'

    shell:
    '''
        for file in input.*; do ln -s $(cat $file) ./; done
        ungz_bedgs=$(
            for f in *.bedGraph.gz; do 
                echo -n "<(pigz -d < $f | sed 's/\\t/:/g4') "
            done
        ) && echo " bedtools unionbedg -filler '0:0:0' -i ${ungz_bedgs}" | bash \
        | tr ':' "\\t" | \
        awk -v OFS='\\t' '{  
                offset=3; meth=0; unmeth=0
                for(i=0; i<=NF-offset+1; i++) {
                    if( i%3 == 1) { meth += $(i+offset+1) }
                    else if( i%3 == 2) { unmeth += $(i+offset+1) }
                }; 
                avg_meth=meth/(meth+unmeth);
                print $1,$2,$3,avg_meth,meth,unmeth
            }' > combined_bedgraph

    '''
}

// process prepare_metadata_for_fastqs {
//     input: 
//         tuple base, file(r1), file(r2) from unconverted_fastq_files.flatten().toList()
//     output: 
//         file('fastq_lib_info.csv') into fastq_files

//     shell:
//     '''
//     for read1_fq in `find . -maxdepth 1 -name '*1.fastq*'`; do

//         dirname=`dirname $read1_fq`
//         library=!{base} #`basename $(basename $read1_fq .1.fastq) .1.fastq.gz`
//         method=`echo $library | cut -f 1 -d _`
//         cell_line=`echo $library |cut -f 2 -d _`
//         lab=`echo $library |cut -f 3 -d _`
//         tile=all
//         header=`zcat -f $read1_fq | head -n 1`
//         flowcell=`echo $header | cut -f 3 -d ':'`
//         lane=`echo $header | head -n 1 | cut -f 4 -d ':'`
//         read2_fq=$(ls ${dirname}/*${library}*2.*);

//         if [ -f $read1_fq ] && [ -f $read2_fq ]; then
//             printf "${library},${method},${lab},${cell_line},${read1_fq},${read2_fq}\n" >> 'fastq_lib_info.csv'
//         else
//             echo "Unable to find $read1_fq or $read2_fq"
//             exit 255
//         fi
//     done
//     '''

// }

// fastq_files
//     .splitCsv(header: ['library', 'method', 'lab','cell_line', 'flowcell', 'read1', 'read2'])
//     .set {fastq_for_bwamem}

process bwamem_md {
    cpus params.aligner_cpus
    errorStrategy 'retry'
    publishDir "output", mode: 'copy', pattern: '*.{sorted.bam}*'

    conda 'bwa=0.7.17 fastp=0.20.1 seqtk=1.3 samblaster=0.1.24 sambamba'

    input:
        tuple val(library), file(reads) from unconverted_fastq_files

    output:
        tuple file('*.sorted.bam'), file('*.sorted.bam.bai') into bwa_mem_aligned_files

    shell:
    '''
    inst_name=$(head -n 1 < <(zcat -f *1.fastq.gz) | cut -f 1 -d ':' | sed 's/^@//')
    fastq_barcode=$(head -n 1 < <(zcat -f *1.fastq.gz) | sed -r 's/.*://')
    flowcell=$(head -n 1 < <(zcat -f *1.fastq.gz) | cut -f 3 -d ':')
    if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || [[ "${inst_name:0:2}" == 'NB' ]]; then
        #these 2-color instruments should have poly-g stretches trimmed since they are likely indicative of inactive clusters
        trim_polyg='--trim_poly_g'
    else
        trim_polyg=''
    fi

    seqtk mergepe <(zcat -f *1.fastq.gz) <(zcat -f *2.fastq.gz) \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j "!{library}_combined_fastp.json" \
    | bwa mem -p -t !{task.cpus} -R \"@RG\\\\tID:${fastq_barcode}\\\\tSM:!{library}\\\\tBC:${fastq_barcode}\\\\tCN:Multiple\\\\tPL:ILLUMINA\\\\tPU:${flowcell}-${fastq_barcode}" !{params.bwamem_genome} /dev/stdin 2>  "!{library}_${fastq_barcode}${flowcell}.log.bwamem" \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -S -f bam -o "!{library}_${fastq_barcode}.bwa.md.bam" /dev/stdin;

    sambamba sort "!{library}_${fastq_barcode}.bwa.md.bam"
    sambamba index *.sorted.bam

    '''

}

bwa_mem_aligned_files.into{bwa_mem_contigs; bwa_mem_split}

process capture_contigs {

    conda 'samtools'

    input:
        tuple file(bam), file(bai) from bwa_mem_contigs
    
    output:
        file('contigs.txt') into contigs

    shell:
    '''
    samtools view -H !{bam} | grep SN | grep -v HLA | grep -v chrUn \
    | grep -v alt | grep -v random | cut -f 2 | cut -f 2 -d ":" > contigs.txt
    '''

}

contigs.splitText().set{chrom_list}
// chrom_list = ['chr1', 'chr2', 'chr3']

process split_bams {
    cpus 2
    conda 'sambamba'

    input:
        tuple file(bam), file(bai) from bwa_mem_split
        each chrom from chrom_list

    output:
        tuple chrom, file('*.split.bam'), file('*.split.bam.bai') into split_bams

    shell:
    '''
    my_chrom=$(echo "!{chrom}" | tr -d \\n)
    sambamba view -h -f bam -t 2 -o "${my_chrom}.bwa.md.split.bam" *.sorted.bam ${my_chrom}
    sambamba index "${my_chrom}.bwa.md.split.bam"
    '''
}

process insilico_convert {

    conda 'pysam'

    input: 
        tuple chrom, file(bam), file(bai) from split_bams
        each file(bedgraph) from combined_bedgraphs
    
    output:
        file('*.converted.bam') into converted_bams

    shell:
    '''
    my_chrom=$(echo "!{chrom}" | tr -d \\n)
    grep -w ${my_chrom} !{bedgraph} > ${my_chrom}_bedgraph    
    /mnt/home/mcampbell/20201001_epiqc/trial/convert_reads.py --bam !{bam} --bed ${my_chrom}_bedgraph --out !{bam}.converted.bam
    '''
}

process combine_sort_bams {
    conda 'samtools'

    input:
        file(bam) from converted_bams.collect()

    output:
        file('merged_converted.sorted.bam') into combined_bam

    shell:
    '''
    samtools cat -o merged_converted.bam *.bam
    samtools sort -n -T ./ -o merged_converted.sorted.bam merged_converted.bam
    '''

}

process bam_to_fastq {
    conda 'samtools'

    input:
        file(bam) from combined_bam

    output:
        tuple file('read1.fastq.gz'), file('read2.fastq.gz') into fastq_for_bwameth

    shell:
    '''
    samtools fastq -n -1 read1.fastq.gz -2 read2.fastq.gz !{bam}
    '''
}

process bwameth_md {
    cpus params.aligner_cpus
    // errorStrategy 'retry'
    publishDir "output", mode: 'copy', pattern: '*.{md.bam}*'

    conda 'bwameth=0.2.2 fastp=0.20.1 seqtk=1.3 samblaster=0.1.24 sambamba'

    input:
        tuple file(read1), file(read2) from fastq_for_bwameth

    output:
        tuple file("*.bwameth.md.bam"), file("*.bwameth.md.bam.bai") into bwameth_aligned_files

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
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j "combined_fastp.json" \
    | bwameth.py -p -t !{task.cpus} \
       --read-group  \"@RG\\\\tID:${fastq_barcode}\\\\tBC:${fastq_barcode}\\\\tCN:Multiple\\\\tPL:ILLUMINA\\\\tPU:${fastq_barcode}" \
       --reference !{params.bwameth_genome} /dev/stdin 2>  "${fastq_barcode}.log.bwamem" \
    | samblaster 2> log.samblaster \
    | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=./ -t !{task.cpus} -m 20GB -o "${fastq_barcode}.bwameth.md.bam" /dev/stdin
    '''
}

bwameth_aligned_files.into{aligned_for_mbias; aligned_for_extract}

process methylDackel_mbias {
    cpus 8
    errorStrategy 'retry'
    conda "methyldackel=0.4.0 samtools=1.9"

    input:
        tuple file(md_file), file(md_bai) from aligned_for_mbias

    output:
        file('*.svg') into mbias_output_svg
        file('*.tsv') into mbias_output_tsv

    shell:
    '''
    echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > UII_combined_mbias.tsv
    chrs=(`samtools view -H !{md_file} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

    for chr in ${chrs[*]}; do
        for context in CHH CHG CpG; do
            arg=''
            if [ $context = 'CHH' ]; then
            arg='--CHH --noCpG'
            elif [ $context = 'CHG' ]; then
            arg='--CHG --noCpG'
            fi
            # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
            # not sure why we need both --keepDupes and -F, probably a bug in mbias
            join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <( \
                MethylDackel mbias --noSVG $arg -@ !{task.cpus} -r $chr !{params.bwameth_genome} !{md_file} | \
                tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
            ) \
            <( \
                MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ !{task.cpus} -r $chr !{params.bwameth_genome} !{md_file} | \
                tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
            ) \
            | sed "s/^/${chr}\t${context}\t/" \
            >> UII_combined_mbias.tsv
        done
    done
    # makes the svg files for trimming checks
    MethylDackel mbias -@ !{task.cpus} --noCpG --CHH --CHG -r ${chrs[0]} !{params.bwameth_genome} !{md_file} UII_chn
    for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

    MethylDackel mbias -@ !{task.cpus} -r ${chrs[0]} !{params.bwameth_genome} !{md_file} UII_cpg
    for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;

    '''

}

process methylDackel_extract {
    cpus 8
    // publishDir "${default_dest_path}/${email}/${flowcell}${dest_modifier}", mode: 'copy'
    conda "methyldackel=0.4.0 pigz=2.4"

    input:
        tuple file(md_file), file(md_bai) from aligned_for_extract

    output:
        file('*.methylKit.gz') into extract_output

    shell:
    '''
    MethylDackel extract --methylKit --OT 0,0,0,95 --OB 0,0,5,0 -@ !{task.cpus} --CHH --CHG -o UII !{params.bwameth_genome} !{md_file}
    pigz -p !{task.cpus} *.methylKit
    '''

}
