all: multiqc_report.html

bt_b_anthracis*: ../data/genomes/Bacillus_anthracis.fa
	bowtie2-build ../data/genomes/Bacillus_anthracis.fa bt_b_anthracis

metagenome.1.fastq metagenome.2.fastq: ../data/b_anthracis.csv
	python ../adrsm -p 0.8 -M 0.04 -m 0.01  ../data/short_genome_list.csv

metagenome_collapse.fq *_fastqc.html metagenome.settings: metagenome.1.fastq metagenome.2.fastq
	fastqc metagenome.*.fastq
	AdapterRemoval --basename metagenome --file1 metagenome.1.fastq --collapse --file2 metagenome.2.fastq --outputcollapsed metagenome_collapse.fq

aligned.bam: metagenome_collapse.fq bt_b_anthracis*
	bowtie2 -x bt_b_anthracis -U metagenome_collapse.fq | samtools view -bS -F 4 - | samtools sort - > aligned.bam

tmp/aligned/DamagePlot.pdf: aligned.bam data/genomes/Bacillus_anthracis.fa
	damageprofiler -i aligned.bam -r data/genomes/Bacillus_anthracis.fa -o tmp

multiqc_report.html: metagenome.settings tmp/aligned/DamagePlot.pdf
	multiqc .

clean:
	rm -rf metagenome* aligned.bam tmp stats.csv bt_b_anthracis* multiqc* 
