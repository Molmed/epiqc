
params.bams_glob = '*.bam'
params.samstats_all_glob = 'samstats/*.samstats'
params.samstats_filtered_glob = 'samstats_noDovetail_filter3844/*.samstats'
params.samstats_q10_filtered_glob = 'samstats_noDovetail_filter3884_minq10/*.samstats'
params.outdir = 'output'

samstats_all_reads = Channel
    .fromPath(params.samstats_all_glob)
    .map { file -> tuple(file.baseName, file) }
samstats_F3844 = Channel
    .fromPath(params.samstats_filtered_glob)
    .map { file -> tuple(file.baseName, file) }
samstats_F3844_q10 = Channel
    .fromPath(params.samstats_q10_filtered_glob)
    .map { file -> tuple(file.baseName, file) }

// process samstats {
//     conda "samtools=1.10 sed"

//     input:
//         tuple library, file(bam) from bams
//     output: 
//         tuple library, file('*.samstats') into samstats_all
//         tuple library, file('*.samstats_F3844') into samstats_F3844
//         tuple library, file('*.samstats_F3844_q10') into samstats_F3844_q10

//     shell:
//     '''
//     for filter in '' '-F 3844' '-F 3844 -q 10'; do 
//         ext=$(echo ${filter} | sed 's/-(.) /\1/g' | tr ' ' '_') 
//         #running all the samstats steps in one process to take advantage of caching
//         samtools view -bS ${flags} !{bam} \
//         | samtools stats -p - \
//         > !{library}.samstats_${ext}
//     done
//     '''
// }

//extract relevant values from all reads stats files
//could not think of a way to make this reusable across stats files
process extract_from_all_reads_stats {
    conda 'grep'
    tag { library }
    input: 
        tuple library, file(samstat_file) from samstats_all_reads
    output: 
        tuple library, file('*raw_total_sequences') into all_reads_total
        tuple library, file('*reads_mapped') into all_reads_mapped

    shell:
    '''
    for entry in 'raw total sequences:' 'reads mapped:' 'bases mapped:'; do 
        ext=$(echo "${entry}" | tr -d ':' | tr ' ' '_')
        grep -E "^SN\\s+${entry}" '!{samstat_file}' \
        | cut -f 3 > "!{library}.all_reads_${ext}"
    done
    '''
}
//extract relevant values from samstats_F3844 filtered files
process extract_from_3844_filtered_stats {
    conda 'grep'
    tag { library }
    input: 
        tuple library, file(samstat_file) from samstats_F3844
    output: 
        tuple library, file('*raw_total_sequences') into F3844_filtered_total
        tuple library, file('*reads_properly_paired') into F3844_filtered_proper
        tuple library, file('*bases_mapped') into F3844_filtered_bases
        tuple library, file('*bases_mapped_cigar') into F3844_filtered_bases_cigar
    shell:
    '''
    for entry in 'raw total sequences:' 'reads properly paired:' 'bases mapped:' 'bases mapped (cigar)'; do 
        ext=$(echo "${entry}" | sed 's/[\\(\\):]//g' | tr ' ' '_')
        grep -E "^SN\\s+${entry}" '!{samstat_file}' \
        | cut -f 3 > "!{library}.F3844_${ext}"
    done
    '''
}

//extract relevant values from samstats_F3844_q10 filtered files
process extract_from_F3844_q10_filtered_stats {
    conda 'grep'
    tag { library }
    input: 
        tuple library, file(samstat_file) from samstats_F3844_q10
    output: 
        tuple library, file('*raw_total_sequences') into F3844_q10_filtered_total
        tuple library, file('*reads_properly_paired') into F3844_q10_filtered_proper
        tuple library, file('*bases_mapped') into F3844_q10_filtered_bases
        tuple library, file('*bases_mapped_cigar') into F3844_q10_filtered_bases_cigar

    shell:
    '''
    for entry in 'raw total sequences:' 'reads properly paired:' 'bases mapped:' 'bases mapped (cigar)'; do 
        ext=$(echo "${entry}" | sed 's/[\\(\\):]//g' | tr ' ' '_')
        grep -E "^SN\\s+${entry}" '!{samstat_file}' \
        | cut -f 3 > "!{library}.F3844_q10_${ext}"
    done
    '''
}
combined_metrics = all_reads_total.mix(all_reads_mapped,
    F3844_filtered_total,
    F3844_filtered_proper,
    F3844_filtered_bases,
    F3844_filtered_bases_cigar,
    F3844_q10_filtered_total,
    F3844_q10_filtered_proper,
    F3844_q10_filtered_bases,
    F3844_q10_filtered_bases_cigar).groupTuple()
    //.subscribe{println "File: ${it.name} => ${it.text}" }

//pastes all the individual entries together
process combined_metrics_by_library {
    input:
        tuple library, file('*') from combined_metrics
    output:
        file('*.lib_stats') into library_stats

    shell:
    '''
        paste <(echo !{library}) \
              *.all_reads_raw_total_sequences \
              *.all_reads_reads_mapped \
              *.F3844_raw_total_sequences \
              *.F3844_reads_properly_paired \
              *.F3844_bases_mapped \
              *.F3844_bases_mapped_cigar \
              *.F3844_q10_raw_total_sequences \
              *.F3844_q10_reads_properly_paired \
              *.F3844_q10_bases_mapped \
              *.F3844_q10_bases_mapped_cigar \
          > !{library}.lib_stats
    '''
}

//cats all the lines and calc losses
process concatenate_all {
    publishDir params.outdir, mode:'copy'
    
    input: 
        file('*') from library_stats.collect();
    output:
        file('alignment_summary_stats.csv') into alignment_summary_stats

    shell:
    '''
        echo -e "library\tall_reads_raw_total_sequences\tall_reads_reads_mapped\tF3844_raw_total_sequences\tF3844_reads_properly_paired\tF3844_bases_mapped\tF3844_bases_mapped_cigar\tF3844_q10_raw_total_sequences\tF3844_q10_reads_properly_paired\tF3844_q10_bases_mapped\tF3844_q10_bases_mapped_cigar\t"\
            > alignment_summary_stats.csv
        cat *.lib_stats >> alignment_summary_stats.csv
    '''
}



