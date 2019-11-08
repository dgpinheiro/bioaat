#!/bin/bash

indir=$1


# SE ${indir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi

# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

outdir=$2

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

mkdir -p ${outdir}/processed/fastqc/pre
mkdir -p ${outdir}/processed/atropos
mkdir -p ${outdir}/processed/prinseq
mkdir -p ${outdir}/processed/fastqc/pos

for r1 in `ls ${indir}/*_R1.fastq`; do
	
	r2=`echo ${r1} | sed 's/_R1.fastq/_R2.fastq/'`
	
	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi

	name=`basename ${r1} | sed 's/_R1.fastq//'`
	echo -e "FastQC pre-evaluation using sample ${name}: ${r1} & ${r2} ...\n"

	fastqc -t 2 \
   		${r1} \
   		-o ${outdir}/processed/fastqc/pre/ \
		 > ${outdir}/processed/fastqc/pre/${name}_R1.log.out.txt \
		2> ${outdir}/processed/fastqc/pre/${name}_R1.log.err.txt

	fastqc -t 2 \
   		${r2} \
		-o ${outdir}/processed/fastqc/pre/ \
		 > ${outdir}/processed/fastqc/pre/${name}_R2.log.out.txt \
		2> ${outdir}/processed/fastqc/pre/${name}_R2.log.err.txt
	
	echo -e "Running atropos (insert) for adapter trimming using sample ${name}: ${r1} & ${r2} ...\n"
	
	# necessário inicializar o ambiente pyenv dentro desta sessão
	eval "$(pyenv init -)"
	# ativando ambiente Python para execução do atropos
	pyenv activate atropos

	atropos trim --aligner insert \
             -e 0.1 \
             -n 2 \
             -m 1 \
             --op-order GAWCQ \
             --match-read-wildcards \
             -O 20 \
             -q 25 \
             -T 2 \
             --correct-mismatches conservative \
             --pair-filter any \
             -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT  \
             -o ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
             -p ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
             -pe1 ${r1} \
             -pe2 ${r2} \
             --untrimmed-output        ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
             --untrimmed-paired-output ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \
               > ${outdir}/processed/atropos/${name}.atropos.log.out.txt \
              2> ${outdir}/processed/atropos/${name}.atropos.log.err.txt
	
	echo -e "Running atropos (adapter) for adapter trimming using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq ...\n"
	atropos trim    --aligner adapter \
                -e 0.1 \
                -n 2 \
                -m 1 \
                --match-read-wildcards \
                -O 3 \
                -q 20 \
                --pair-filter both \
                -pe1 ${outdir}/processed/atropos/${name}_R1.atropos_untrimmed.fastq \
                -pe2 ${outdir}/processed/atropos/${name}_R2.atropos_untrimmed.fastq \
                -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
                -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT  \
                -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                -G CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
                -T 2 \
                -o  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq  \
                -p  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
                 >  ${outdir}/processed/atropos/${name}.atropos_adapter.log.out.txt \
                2>  ${outdir}/processed/atropos/${name}.atropos_adapter.log.err.txt
	
	echo -e "Merging atropos adapter trimming results using sample ${name}: ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq + ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq and ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq ...\n"
	
	cat       ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
        	  ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
   		> ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq

	cat       ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
        	  ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq \
   		> ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq

	echo -e "Removing useless atropos results ...\n"

	rm -f ${outdir}/processed/atropos/${name}_R1.atropos_insert.fastq \
	      ${outdir}/processed/atropos/${name}_R1.atropos_adapter.fastq \
	      ${outdir}/processed/atropos/${name}_R2.atropos_insert.fastq \
	      ${outdir}/processed/atropos/${name}_R2.atropos_adapter.fastq

	echo -e "PrinSeq processing: ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq & ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq ...\n"

	prinseq-lite.pl -fastq  ${outdir}/processed/atropos/${name}_R1.atropos_final.fastq \
			-fastq2 ${outdir}/processed/atropos/${name}_R2.atropos_final.fastq \
			-out_format 3 \
			-trim_qual_window 3 \
			-trim_qual_step 1 \
			-trim_qual_right 30 \
			-trim_qual_type mean \
			-trim_qual_rule lt \
			-out_good ${outdir}/processed/prinseq/${name}.atropos_final.prinseq \
			-out_bad  null \
			-lc_method dust \
			-lc_threshold 30 \
			-min_len 20 \
			-trim_tail_right 5 \
			-trim_tail_left 5 \
			-ns_max_p 80 \
			-noniupac \
			 > ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.out.log \
			2> ${outdir}/processed/prinseq/${name}.atropos_final.prinseq.err.log

	echo -e "FastQC pos-evaluation using sample ${name}: ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq & ${outdir}/processed/prinseq/${name}_2.atropos_final.prinseq.fastq ...\n"

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1.log.err.txt

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}_2.atropos_final.prinseq.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2.log.err.txt

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_1_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_1_singletons.fastq \
	   	   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_1_singletons.log.err.txt
	fi

	# SE EXISTIR <SAMPLE_NAME>.atropos_final.prinseq_2_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.atropos_final.prinseq_2_singletons.fastq \
		   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.atropos_final.prinseq_2_singletons.log.err.txt
	fi


done
