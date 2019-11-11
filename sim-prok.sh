#!/bin/bash
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2012  Universidade de São Paulo
#
#  Universidade de São Paulo
#  Laboratório de Biologia do Desenvolvimento de Abelhas
#  Núcleo de Bioinformática (LBDA-BioInfo)
#
#  Daniel Guariz Pinheiro
#  dgpinheiro@gmail.com
#  http://zulu.fmrp.usp.br/bioinfo 
#

rm -f transcriptoma.fa

rm -f transcriptoma.fa

IFS=$'\n'
for accline in $(cat ./ACCS.txt); do
        acc=`echo ${accline} | cut -f 1`
        seqref=`echo ${accline} | cut -f 2`
        chr_start=`echo ${accline} | cut -f 3`
        chr_stop=`echo ${accline} | cut -f 4`
        strand=`echo ${accline} | cut -f 5`

        echo "Pegando FASTA para ${acc}  [${seqref}:${chr_start}-${chr_stop}(${strand})] ..."

        efetch -db nucleotide -id ${seqref}  -format fasta \
        -chr_start ${chr_start} \
        -chr_stop ${chr_stop} \
        -strand ${strand} | \
        sed "s/^>.*/>${acc}/" \
        >> transcriptoma.fa
done

for biogroup in A B; do
	for rep in 1 2; do
		echo "Gerando reads para amostra ${biogroup} réplica ${rep} ..."

		generate_fragments.py -r transcriptoma.fa \
		   -a ./abundance_${biogroup}.txt \
		   -o ./tmp.frags_${biogroup}_${rep} \
		   -t 25000 \
		   -i 300 \
		   -s 30

		cat ./tmp.frags_${biogroup}_${rep}.1.fasta | renameSeqs.pl \
		   -if FASTA \
		   -of FASTA \
		   -p SAMPLE${biogroup}${rep} \
		   -w 1000 | \
		   sed 's/^>\(\S\+\).*/>\1/' \
		   > ./frags_${biogroup}${rep}.fa

		cat ./frags_${biogroup}${rep}.fa | simNGS -a \
		AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-p paired \
		/usr/local/bioinfo/simNGS/data/s_4_0099.runfile \
		-n 151 > ./SAMPLE${biogroup}${rep}.fastq 2> SAMPLE${biogroup}${rep}.err.txt

		mkdir -p ./raw

		deinterleave_pairs SAMPLE${biogroup}${rep}.fastq \
		   -o ./raw/SAMPLE${biogroup}${rep}_R1.fastq \
		      ./raw/SAMPLE${biogroup}${rep}_R2.fastq

		rm -f ./tmp.frags_${biogroup}_${rep}.1.fasta ./frags_${biogroup}${rep}.fa ./SAMPLE${biogroup}${rep}.fastq ./SAMPLE${biogroup}${rep}.err.txt

		echo "Número de reads ${biogroup}${rep} R1:" $(echo "$(cat raw/SAMPLE${biogroup}${rep}_R1.fastq | wc -l)/4" | bc)
		echo "Número de reads ${biogroup}${rep} R2:" $(echo "$(cat raw/SAMPLE${biogroup}${rep}_R2.fastq | wc -l)/4" | bc)

	done
done
