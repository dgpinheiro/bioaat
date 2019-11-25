#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

# as linhas que iniciam com cerquilha são comentários

# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi

# output - diretório para armazenar o resultado do processo de montagem
output=$2

# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi


num_threads=23

###
# Arquivos e diretórios de entrada (input)
#

###
# Arquivos e diretórios de saída (output) 
#

basedir_out="${output}"

trinity_out="${basedir_out}/gg-assembled"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir


mkdir -p ${trinity_out}

bamfiles=()

bamfiles=( $( find ${input} -name 'Aligned.out.sorted.bam' ) )

samtools merge -f ${basedir_out}/All.sorted.bam ${bamfiles[*]}

if [ ! -d ${trinity_out}/Trinity.timing ]; then
	
	echo -e "Assembling step (Trinity) ..."

	Trinity --KMER_SIZE 27 \
		--output ${trinity_out} \
		--seqType fq \
		--max_memory 160G \
		--CPU ${num_threads} \
		--min_per_id_same_path 95 \
		--max_diffs_same_path  5 \
		--path_reinforcement_distance 5 \
		--group_pairs_distance 500 \
		--min_glue 5 \
		--min_contig_length 600 \
               	--min_kmer_cov 3 \
		--genome_guided_bam ${basedir_out}/All.sorted.bam \
		--genome_guided_max_intron 10000 \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
fi
