// Script parameters
fwd_primers = file( params.fwd_primers )
rev_primers = file( params.rev_primers )

params.dev = false
params.number_of_inputs = 1

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/
Channel
    .fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .take( params.dev ? params.number_of_inputs : -1 )
	.set { read_pairs }

// cutadapt based quality filtering and QC
// Filter to only reads with primer detected and trim poly-g
process cutadapt {
    
	publishDir "${params.outdir}",
		saveAs: {filename ->
			if (filename.indexOf(".log") > 0) "${pair_id}/logs/${filename}"
			else if (filename.indexOf(".html") > 0) "${pair_id}/htmls/${filename}"
			else "${pair_id}/${filename}"
		}

	input:
	set pair_id, file(reads) from read_pairs
	file fwd_primers
	file rev_primers

	output:
        file("trimmed_noprimers")
        file("trimmed_primers")
	file("*.cutadapt.log")
        file("*.html")

 	conda 'bioconda::cutadapt=3.4 bioconda::fastp=0.20.1'

        time '120m'
	cpus 4
	penv 'smp'
	memory '8 GB'
	
	script:
	"""
	#!/usr/bin/env bash
	set -e

	mkdir untrimmed_primers

	cutadapt \
	    --action=retain \
	--discard-untrimmed \
        -g file:${fwd_primers} \
		-G file:${rev_primers} \
		--pair-adapters \
		-j \$NSLOTS \
        -o untrimmed_primers/{name}_${pair_id}_R1_primers.fastq.gz \
        -p untrimmed_primers/{name}_${pair_id}_R2_primers.fastq.gz \
        ${reads[0]} \
        ${reads[1]} \
        > ${pair_id}.primers.cutadapt.log

	mkdir trimmed_primers

	for i in `ls untrimmed_primers/*_R1_primers.fastq.gz | grep -oP "(?<=untrimmed_primers/).*(?=_R1_primers)"`;
	do
	fastp \
		-i untrimmed_primers/\${i}_R1_primers.fastq.gz \
		-I untrimmed_primers/\${i}_R2_primers.fastq.gz \
		-o trimmed_primers/\${i}_R1_trimmed_primers.fastq.gz \
		-O trimmed_primers/\${i}_R2_trimmed_primers.fastq.gz \
		-w \$NSLOTS \
		-h \${i}_primers.html
	done

	mkdir untrimmed_noprimers

	cutadapt \
	    --action=trim \
        --discard-untrimmed \
        -g file:${fwd_primers} \
		-G file:${rev_primers} \
		--pair-adapters \
		-j \$NSLOTS \
        -o untrimmed_noprimers/{name}_${pair_id}_R1_noprimers.fastq.gz \
        -p untrimmed_noprimers/{name}_${pair_id}_R2_noprimers.fastq.gz \
        ${reads[0]} \
        ${reads[1]} \
        > ${pair_id}.noprimers.cutadapt.log

        mkdir trimmed_noprimers
	
	for i in `ls untrimmed_noprimers/*_R1_noprimers.fastq.gz | grep -oP "(?<=untrimmed_noprimers/).*(?=_R1_noprimers)"`;
	do
	fastp \
		-i untrimmed_noprimers/\${i}_R1_noprimers.fastq.gz \
		-I untrimmed_noprimers/\${i}_R2_noprimers.fastq.gz \
		-o trimmed_noprimers/\${i}_R1_trimmed_noprimers.fastq.gz \
		-O trimmed_noprimers/\${i}_R2_trimmed_noprimers.fastq.gz \
		-w \$NSLOTS \
		-h \${i}_noprimers.html
	done
	
	"""
}

