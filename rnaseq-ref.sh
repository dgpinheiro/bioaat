#!/bin/bash

num_threads=2

indir=$1


# SE ${indir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 1 NA LINHA DE COMANDO
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

# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 2 NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi

refgtf=$3
# SE ${refgtf} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO
if [ ! ${refgtf} ]; then
	echo "Missing GTF file."
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "Not found GTF file (${refgtf})."
	exit
fi

refseq=$4
# SE ${refseq} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO
if [ ! ${refseq} ]; then
	echo "Missing GENOME fasta file."
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})."
	exit
fi

./preprocess2.sh "${indir}" "${outdir}"

mkdir -p ${outdir}/star_index
mkdir -p ${outdir}/star_out_pe
mkdir -p ${outdir}/star_out_se
mkdir -p ${outdir}/star_out_final

for r1 in `ls ${outdir}/processed/prinseq/*.atropos_final.prinseq_1.fastq`; do
	r1_singletons=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_1_singletons.fastq/'`
	if [ ! -e "${r1_singletons}" ]; then
		touch ${r1_singletons}
	fi
	
	r2=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_2.fastq/'`
	
	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi
	
	r2_singletons=`echo ${r2} | sed 's/prinseq_2.fastq/prinseq_2_singletons.fastq/'`
	if [ ! -e "${r2_singletons}" ]; then
		touch ${r2_singletons}
	fi

	name=`basename ${r1} | sed 's/.atropos_final.prinseq_1.fastq//'`
	
	if [ ! -e "${outdir}/star_index/SAindex" ]; then
		echo "Indexing genome (${refseq}) ..."
		
		STAR 	--runThreadN        ${num_threads} \
     			--runMode           genomeGenerate \
     			--genomeFastaFiles  ${refseq} \
     			--genomeDir         ./${outdir}/star_index \
    	 		--sjdbGTFfile       ./genome.gtf \
			--genomeSAindexNbases 12 \
     			--sjdbOverhang      149 \
		 > ./${outdir}/star_index/STAR.index.log.out.txt \
		2> ./${outdir}/star_index/STAR.index.log.err.txt

	fi

	echo "STAR alignment PE with sample ${name}: ${r1} & ${r2} ..."
	
	mkdir -p ${outdir}/star_out_pe/${name}
	
	STAR 	--runThreadN  	    ${num_threads} \
         	--genomeDir	    ${outdir}/star_index \
         	--readFilesIn       ${r1} ${r2} \
         	--sjdbGTFfile	    ${refgtf} \
         	--outFileNamePrefix ./$outdir/star_out_pe/${name} \
		 > ./${outdir}/star_index/STAR.alignment_pe.log.out.txt \
		2> ./${outdir}/star_index/STAR.alignment_pe.log.err.txt

	echo "STAR alignment SE with sample ${name}: ${r1_singletons} & ${r2_singletons} ..."
	
	STAR 	--runThreadN  	    ${num_threads} \
         	--genomeDir	    ${outdir}/star_index \
         	--readFilesIn       ${r1_singletons},${r2_singletons} \
         	--sjdbGTFfile	    ${refgtf} \
         	--outFileNamePrefix ./$outdir/star_out_se/${name} \
		 > ./${outdir}/star_index/STAR.alignment_se.log.out.txt \
		2> ./${outdir}/star_index/STAR.alignment_se.log.err.txt

	echo "Merging STAR alignment PE & SE ..."

        # Combinar resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons)       
        samtools merge -@ ${num_threads} -f -n  ${outdir}/star_out_final/${name}/Aligned.out.bam \
                                                ${outdir}/star_out_pe/${name}/Aligned.out.bam \
                                                ${outdir}/star_out_se/${name}/Aligned.out.bam

	echo "Sorting STAR alignment final ..."
        # Ordenando o resultado do alinhamento por coordenadas genômicas
        # - exigência para executar o cufflinks
        samtools sort -@ ${num_threads} -o      ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
                                                ${outdir}/star_out_final/${name}/Aligned.out.bam

	echo "Collecting alignment statistics ..."
	
	SAM_nameSorted_to_uniq_count_stats.pl ${outdir}/star_out_final/${name}/Aligned.out.bam > ${outdir}/star_out_final/${name}/Aligned.stats.txt
	
done
