source activate TSDfinder
conda env export > doc/environment.yml

if [ ! -d fastq ]; then
	mkdir -p bam clean fastq fastqc/raw fastqc/clean output stat TE/index track
	cp example/example.TE_coordinates.bed TE/TE_coordinates.bed
fi
