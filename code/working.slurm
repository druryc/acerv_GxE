#!/bin/bash
#SBATCH --job-name=angsd
#SBATCH --partition=shared
#SBATCH --time=0-24:00:00 ## time format is DD-HH:MM:SS
#SBATCH --cpus-per-task=12
#SBATCH --mem=10000 ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node
#SBATCH --error=hello-%A.err ## %A - filled with jobid
#SBATCH --output=hello-%A.out ## %A - filled with jobid
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=druryc@hawaii.edu
#srun -I30 -p sandbox -N 1 -c 1 --mem=6G -t 0-01:00:00 --pty /bin/bash

source ~/.bash_profile
#module load lang/Java
#module load lang/Anaconda3/5.1.0
#module load lang/Python/3.7.2-intel-2018.5.274   
#module load  bio/Biopython/1.73-intel-2018.5.274-Python-2.7.15
#module load bio/BLAST+/2.2.31
#module load lang/R
#source activate tools_seq

#ALIGN HOST DATA
	#bwa index ~/reference/amillepora.fasta
	#samtools faidx ~/reference/amillepora.fasta
	#gunzip *.gz
	#ls *.fastq >list.txt; sed -i 's/.fastq//g' list.txt
	
	#while read f; do
	#	bwa mem -t 8 ~/reference/amillepora.fasta "$f".fastq > "$f".sam;
	#	samtools view -bh "$f".sam | samtools sort -o "$f".bam;
	#	samtools index "$f".bam;
	#done<list.txt

#ANALYZE HOST DATA
	#ls *.bam >bamlist.txt
	#angsd -bam bamlist.txt -GL 1 -nThreads 12 -doPost 1 -doGeno 8 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-4  -minInd 10 -minmapQ 30 -minQ 20 -doCounts 1 -dumpCounts 2 -geno_minDepth 7 -doGlf 2 -doIBS 1 -makeMatrix 1 -out genolike_depth
	#python ~/pcangsd/pcangsd.py -beagle genolike.beagle.gz -o pcangsd_out -threads 12 -indf_save -sites_save 

#NGSADMIX
	#~/NGSadmix -likes genolike.beagle.gz -K 1 -P 1 -o output1_1 -minMaf 0.05 -printInfo 1
	#~/NGSadmix -likes genolike.beagle.gz -K 2 -P 1 -o output2_1 -minMaf 0.05 -printInfo 1
	#~/NGSadmix -likes genolike.beagle.gz -K 3 -P 1 -o output3_1 -minMaf 0.05 -printInfo 1
	#~/NGSadmix -likes genolike.beagle.gz -K 4 -P 1 -o output4_1 -minMaf 0.05 -printInfo 1
	#grep 'best like' *.log > comparison.txt
	#nano comparison.txt

#ALIGN AND SUMMARISE SYMBIONT DATA
	#cat *.fasta >symbiont_combined.fasta
	#bwa index ~/reference/sym/symbiont_combined.fasta
	#samtools faidx ~/reference/sym/symbiont_combined.fasta
	#while read f; do
	#	 bwa mem -t 8 ~/reference/sym/symbiont_combined.fasta "$f".fastq > "$f"_symb.sam;
	#     samtools view -F 260 -q 30 "$f"_symb.sam -o "$f"_filt.sam;
	#     grep 'kb8' "$f"_filt.sam |wc -l >> "$f"_summary.txt; 
	#     grep 'isogroup' "$f"_filt.sam |wc -l >> "$f"_summary.txt; 
	#     grep 'comp' "$f"_filt.sam |wc -l >> "$f"_summary.txt; 
	#     grep 'TRINITY' "$f"_filt.sam|wc -l >> "$f"_summary.txt;
	#done<list.txt
	#cat *summary.txt >symbiont_summary.txt

#ALIGN AND ANALYZE ADDITIONAL HOST DATA 
	#'short_list_additional.bed' produced in R script from 58 loci of interest to limit angsd sites for computational/filtering efficiency
	
	#ls *.fastq >list.txt; sed -i 's/.fastq.gz//g' list.txt
	#while read f; do
	#	gunzip "$f".fastq.gz;
	#	bwa mem -t 12 ~/reference/amillepora.fasta "$f".fastq > "$f".sam;
	#	samtools view -bh -@ 11 "$f".sam -L short_list_additional.bed -@ 11 | samtools sort -@ 11 -o "$f".bam;
	#	samtools index -@ 11"$f".bam;
	#done<list.txt
	#ls *.bam >bamlist.txt
	#angsd -bam bamlist.txt -GL 1 -nThreads 12 -doPost 1 -doGeno 8 -doMajorMinor 1 -doMaf 2 -minInd 5 -minmapQ 30 -minQ 20 -doCounts 1 -dumpCounts 2 -doGlf 2 -out acerv_extra

#PREP BLAST DB
	#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
	#gunzip uniprot_sprot.fasta.gz
	#makeblastdb -in uniprot_sprot.fasta -dbtype prot

#BLAST ALL, LOOK UP GO, REFORMAT
	#blastx -query list-aa  -db ./uniprot_sprot.fasta -evalue 0.0001 -num_threads 4 -num_alignments 5 -out blast_output-aa.txt -outfmt 6
	#grep -o "sp|[A-Z0-9]*" blast_output-aa.txt | sed 's/sp|//' >gene_names-aa.txt

	#mkdir GOfromweb
	#while read f; do
	#	wget "https://www.uniprot.org/uniprot/?query=id:"$f"&columns=id,go&format=tab" --output-document=./GOfromweb/"$f".txt;
	#done<gene_names-aa.txt

	#cd GOfromweb
	#scp to laptop

#FINAL GO FORMATTING FOR GO_MWU
	#cat *.txt > master_GO.txt
	#sed '/^Entry/d' master_GO.txt > master_GO_reformat.txt
	#sed '/\t$/d' master_GO_reformat.txt > int.txt
	#awk -F '\t' '{print $1}' int.txt >gene_list.txt
	#sed '/\t$/d' int.txt | tr '\n' 'x' |grep -oE 'GO:.{9}'|tr "x" "\n" |sed 's/\]//g' |perl -p -e 's/;\n/;/'|sed '/^$/d' >GOterms.txt
	#paste gene_list.txt GOterms.txt |awk '{ print length, $0 }' | sort -nr -s |sed -E 's/[0-9]+\s//' >GO_output.txt
