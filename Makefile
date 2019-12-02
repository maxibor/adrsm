all: tmp/aligned/DamagePlot.pdf

bt_b_anthracis*: data/genomes/Bacillus_anthracis.fa
	bowtie2-build data/genomes/Bacillus_anthracis.fa bt_b_anthracis

metagenome.1.fastq metagenome.2.fastq: data/b_anthracis.csv
	python adrsm -p 0.8 -M 0.6 -m 0.1  ./data/b_anthracis.csv

metagenome.{1,2}.trimmed.fastq: metagenome.1.fastq metagenome.2.fastq
	AdapterRemoval --basename metagenome --file1 metagenome.1.fastq --file2 metagenome.2.fastq --output1 metagenome.1.trimmed.fastq --output2 metagenome.2.trimmed.fastq

aligned.bam: metagenome.{1,2}.trimmed.fastq bt_b_anthracis*
	bowtie2 -x bt_b_anthracis -1 metagenome.1.trimmed.fastq -2 metagenome.2.trimmed.fastq | samtools view -bS - | samtools sort - > aligned.bam

tmp/aligned/DamagePlot.pdf: aligned.bam data/genomes/Bacillus_anthracis.fa
	damageprofiler -i aligned.bam -r data/genomes/Bacillus_anthracis.fa -o tmp

clean:
	rm -rf metagenome* aligned.bam tmp stats.csv bt_b_anthracis*