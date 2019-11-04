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

# SE ${outdir} EXISTE, SE FOI PASSADO ARGUMENTO NA LINHA DE COMANDO
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
mkdir -p ${outdir}/processed/prinseq
mkdir -p ${outdir}/processed/fastqc/pos

for r1 in `ls ${indir}/*_R1.fastq`; do
	r2=`echo ${r1} | sed 's/_R1.fastq/_R2.fastq/'`

	name=`basename ${r1} | sed 's/_R1.fastq//'`
	echo "FastQC pre-evaluation using sample ${name}: ${r1} & ${r2} ..."

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

	echo "PrinSeq processing: ${r1} & ${r2} ..."

	prinseq-lite.pl -fastq  ${r1} \
			-fastq2 ${r2} \
			-out_format 3 \
			-trim_qual_window 3 \
			-trim_qual_step 1 \
			-trim_qual_right 30 \
			-trim_qual_type mean \
			-trim_qual_rule lt \
			-out_good ${outdir}/processed/prinseq/${name}.prinseq \
			-out_bad  null \
			-lc_method dust \
			-lc_threshold 30 \
			-min_len 20 \
			-trim_tail_right 5 \
			-trim_tail_left 5 \
			-ns_max_p 80 \
			-noniupac \
			 > ${outdir}/processed/prinseq/${name}.prinseq.out.log \
			2> ${outdir}/processed/prinseq/${name}.prinseq.err.log

	echo "FastQC pos-evaluation using sample ${name}: ${outdir}/processed/prinseq/${name}.prinseq_1.fastq & ${outdir}/processed/prinseq/${name}_2.prinseq.fastq ..."

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}.prinseq_1.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.prinseq_1.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.prinseq_1.log.err.txt

	fastqc -t 2 \
	   ${outdir}/processed/prinseq/${name}_2.prinseq.fastq \
	   -o ${outdir}/processed/fastqc/pos/ \
	    > ${outdir}/processed/fastqc/pos/${name}.prinseq_2.log.out.txt \
	   2> ${outdir}/processed/fastqc/pos/${name}.prinseq_2.log.err.txt

	# SE EXISTIR <SAMPLE_NAME>.prinseq_1_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.prinseq_1_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.prinseq_1_singletons.fastq \
	   	   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.prinseq_1_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.prinseq_1_singletons.log.err.txt
	fi

	# SE EXISTIR <SAMPLE_NAME>.prinseq_2_singletons.fastq
	if [ -e "${outdir}/processed/prinseq/${name}.prinseq_2_singletons.fastq" ]; then
		fastqc -t 2 \
		   ${outdir}/processed/prinseq/${name}.prinseq_2_singletons.fastq \
		   -o ${outdir}/processed/fastqc/pos/ \
	            > ${outdir}/processed/fastqc/pos/${name}.prinseq_2_singletons.log.out.txt \
	           2> ${outdir}/processed/fastqc/pos/${name}.prinseq_2_singletons.log.err.txt
	fi
done
