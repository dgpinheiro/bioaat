#!/bin/bash

# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

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

num_threads="8"
mem_gb="10G"

###
# Arquivos e diretórios de saída (output) 
#

basedir_out="${output}/"

renamed_out="${basedir_out}/renamed"

trinity_out="${basedir_out}/trinity_assembled"

mkdir -p ${renamed_out}


left=()
left_singleton=()

right=()
right_singleton=()

echo "Renaming step ..."

mkdir -p ${trinity_out}

for fastq in `ls ${input}/*.fastq`; do
	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq}`;
	if [[ ! $fastqbn =~ \.bad_ ]]; then
		renamed_fastq="${renamed_out}/${fastqbn}"
		if [ ! -e ${renamed_fastq} ]; then
			echo -e "\tRenaming ${fastqbn} ..."
			if [[ ${fastqbn} =~ _1[\._] ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/1$/) { print $1"/1" } else { print $0 } } else if (NR%4==3) { print "+" } else { print $0 } }' ${fastq} > ${renamed_fastq}
			elif [[ ${fastqbn} =~ _2[\._]  ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/2$/) { print $1"/2" } else { print $0 } } else if (NR%4==3) { print "+" } else { print $0 } }' ${fastq} > ${renamed_fastq}
			else 
				echo "Warning: ${fastqbn} discarded!"
			fi
		fi
		
		if [[ ${fastqbn} =~ _1[\._] ]]; then
			if [[ ${fastqbn} =~ singletons ]]; then
				left_singleton=($(printf "%s\n" ${left_singleton[@]} ${renamed_fastq} | sort -u ))
			else
				left=($(printf "%s\n" ${left[@]} ${renamed_fastq}  | sort -u ))
			fi
		elif [[ ${fastqbn} =~ _2[\._] ]]; then
			if [[ ${fastqbn} =~ singleton ]]; then
				right_singleton=($(printf "%s\n" ${right_singleton[@]} ${renamed_fastq}  | sort -u ))
			else
				right=($(printf "%s\n" ${right[@]} ${renamed_fastq}  | sort -u ))
			fi
		else
			echo "Warning: ${fastqbn} discarded!"
		fi
	fi
done


if [ ! -d ${trinity_out}/Trinity.timing ]; then
	
	echo -e "Assembling step (Trinity) ..."

	Trinity --KMER_SIZE 27 \
		--output ${trinity_out} \
		--seqType fq \
		--max_memory ${mem_gb} \
		--CPU ${num_threads} \
		--min_per_id_same_path 95 \
		--max_diffs_same_path  5 \
		--path_reinforcement_distance 5 \
		--group_pairs_distance 500 \
		--min_glue 5 \
		--min_contig_length 600 \
		--min_kmer_cov 3 \
		--left $(IFS=, ; echo "${left[*]},${left_singleton[*]}") \
		--right $(IFS=, ; echo "${right[*]},${right_singleton[*]}") \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
fi
