#!/bin/bash

mkdir -p ./output/processed/fastqc/pre
mkdir -p ./output/processed/prinseq
mkdir -p ./output/processed/fastqc/pos

for r1 in `ls ./raw/*_R1.fastq`; do
	r2=`echo ${r1} | sed 's/_R1.fastq/_R2.fastq/'`

	name=`basename ${r1} | sed 's/_R1.fastq//'`
	echo "FastQC pre-evaluation using sample ${name}: ${r1} & ${r2} ..."

	fastqc -t 2 \
   		${r1} \
   		-o ./output/processed/fastqc/pre/ \
		 > ./output/processed/fastqc/pre/${name}_R1.log.out.txt \
		2> ./output/processed/fastqc/pre/${name}_R1.log.err.txt

	fastqc -t 2 \
   		${r2} \
		-o ./output/processed/fastqc/pre/ \
		 > ./output/processed/fastqc/pre/${name}_R2.log.out.txt \
		2> ./output/processed/fastqc/pre/${name}_R2.log.err.txt

	echo "PrinSeq processing: ${r1} & ${r2} ..."

	prinseq-lite.pl -fastq  ${r1} \
			-fastq2 ${r2} \
			-out_format 3 \
			-trim_qual_window 3 \
			-trim_qual_step 1 \
			-trim_qual_right 30 \
			-trim_qual_type mean \
			-trim_qual_rule lt \
			-out_good ./output/processed/prinseq/${name}.prinseq \
			-out_bad  ./output/processed/prinseq/${name}.prinseq.bad \
			-lc_method dust \
			-lc_threshold 30 \
			-min_len 20 \
			-trim_tail_right 5 \
			-trim_tail_left 5\
			-ns_max_p 80 \
			-noniupac \
			 > ./output/processed/prinseq/${name}.prinseq.out.log \
			2> ./output/processed/prinseq/${name}.prinseq.err.log

	echo "FastQC pos-evaluation using sample ${name}: ./output/processed/prinseq/${name}.prinseq_1.fastq & ./output/processed/prinseq/${name}_2.prinseq.fastq ..."

	fastqc -t 2 \
	   ./output/processed/prinseq/${name}.prinseq_1.fastq \
	   -o ./output/processed/fastqc/pos/ \
	    > ./output/processed/fastqc/pos/${name}.prinseq_1.log.out.txt \
	   2> ./output/processed/fastqc/pos/${name}.prinseq_1.log.err.txt

	fastqc -t 2 \
	   ./output/processed/prinseq/${name}_2.prinseq.fastq \
	   -o ./output/processed/fastqc/pos/ \
	    > ./output/processed/fastqc/pos/${name}.prinseq_2.log.out.txt \
	   2> ./output/processed/fastqc/pos/${name}.prinseq_2.log.err.txt

	# SE EXISTIR <SAMPLE_NAME>.prinseq_1_singletons.fastq
	if [ -e "./output/processed/prinseq/${name}.prinseq_1_singletons.fastq" ]; then
		fastqc -t 2 \
		   ./output/processed/prinseq/${name}.prinseq_1_singletons.fastq \
	   	   -o ./output/processed/fastqc/pos/ \
	            > ./output/processed/fastqc/pos/${name}.prinseq_1_singletons.log.out.txt \
	           2> ./output/processed/fastqc/pos/${name}.prinseq_1_singletons.log.err.txt
	fi

	# SE EXISTIR <SAMPLE_NAME>.prinseq_2_singletons.fastq
	if [ -e "./output/processed/prinseq/${name}.prinseq_2_singletons.fastq" ]; then
		fastqc -t 2 \
		   ./output/processed/prinseq/${name}.prinseq_2_singletons.fastq \
		   -o ./output/processed/fastqc/pos/ \
	            > ./output/processed/fastqc/pos/${name}.prinseq_2_singletons.log.out.txt \
	           2> ./output/processed/fastqc/pos/${name}.prinseq_2_singletons.log.err.txt
	fi
done
