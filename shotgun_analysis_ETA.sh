#!/bin/bash

eval "$(conda shell.bash hook)"

#Concatenate all fastq files from NextSeq run

#cd 00_fastq_uncat
#for d in *L001_R1*; \
#	do cat $d ${d%%_L001*}_L002_R1_001.fastq.gz ${d%%_L001*}_L003_R1_001.fastq.gz ${d%%_L001*}_L004_R1_001.fastq.gz > ../01_fastq_raw/${d}; \
#done
#for d in *L001_R2*; \
#	do cat $d ${d%%_L001*}_L002_R2_001.fastq.gz ${d%%_L001*}_L003_R2_001.fastq.gz ${d%%_L001*}_L004_R2_001.fastq.gz > ../01_fastq_raw/${d}; \
#done

mkdir 02_fastqc_raw 03_fastq_trimmed 04_fastqc_trimmed 05_fastq_unmapped 06_fastqc_unmapped 07_kraken 08_bracken 09_ARGs-OAP 10_HUMAnN 11_mumame 12_metaSPAdes 13_metaQUAST 14_prokka 15_BacPipe

#Assess quality of the reads with FastQC

cd 01_fastq_raw
fastqc *.fastq.gz -o ../02_fastqc_raw/ -t 16
cd ../02_fastqc_raw/
mkdir zips
for d in *.zip; \
	do unzip $d -d zips/; \
done
cd zips/
mkdir txt
for d in *; \
	do cd $d; \
	mv "fastqc_data.txt" ../txt/"$d.txt"; \
	cd ../; \
done
cd txt/
sed -s -n '4p;7p' *_R1_001_fastqc.txt > summary.txt
cd ../../../01_fastq_raw

#Adapter trimming with TrimGalore

conda activate trimgalore 
for d in *_R1_001.fastq.gz; \
    do trim_galore --output_dir ../03_fastq_trimmed/ --quality 30 --cores 8 --paired $d ${d%%_R1*}_R2_001.fastq.gz; \
done 
conda deactivate

#Assess quality of trimmed reads and calculate number of reads removed

cd ../03_fastq_trimmed
fastqc *.fq.gz -o ../04_fastqc_trimmed/ -t 16
cd ../04_fastqc_trimmed/
mkdir zips
for d in *.zip; \
	do unzip $d -d zips/; \
done
cd zips/
mkdir txt
for d in *; \
	do cd $d; \
	mv "fastqc_data.txt" ../txt/"$d.txt"; \
	cd ../; \
done
cd txt/
sed -s -n '4p;7p' *_val_1_fastqc.txt > summary.txt
cd ../../../03_fastq_trimmed

#Host reads removal with bowtie2

conda activate mapping
for d in *1.fq.gz; \
    do bowtie2 -p 16 -x ~/reference_database/GRCh38_noalt_as/GRCh38_noalt_as -1 $d -2 ${d%%_R1*}_R2_001_val_2.fq.gz --un-conc-gz ../05_fastq_unmapped/${d%%_S*}_unmapped; \
done
cd ../05_fastq_unmapped
for d in *.1; \
    do mv $d ${d%%.1}_R1.fastq.gz; \
done
for d in *.2; \
    do mv $d ${d%%.2}_R2.fastq.gz; \
done

#Additional removal with SNAP - needs uncompressed files

gunzip *.gz
for d in *_R1.fastq; \
    do snap-aligner paired ~/reference_database/snap_index/ $d ${d%%_R1*}_R2.fastq -o ./${d%%_u*}_SNAP.bam -t 16 -F u; \
    samtools sort -n ${d%%_u*}_SNAP.bam -o ${d%%_u*}_SNAP.sorted.bam; \
    bedtools bamtofastq -i ${d%%_u*}_SNAP.sorted.bam -fq ${d%%_u*}_SNAP_1.fq -fq2 ${d%%_u*}_SNAP_2.fq; \
    rm *.bam; \
done

#Merge paired-end reads for HUMAnN analysis before compressing

#for d in *_1.fq; \
#	do cat $d ${d%%_1*}_2.fq > ../10_HUMAnN/${d%%_1*}.fastq; \
#done

#Compress the files to save space (and remove bowtie2 output)

rm *.fastq
gzip *.f*
#mkdir bowtie_out
#mv *.fastq.gz bowtie_out

conda deactivate

#FastQC to calculate number of reads removed in host removal step

fastqc *.f*.gz -o ../06_fastqc_unmapped
cd ../06_fastqc_unmapped/
mkdir zips
for d in *.zip; \
	do unzip $d -d zips/; \
done
cd zips/
mkdir txt
for d in *; \
	do cd $d; \
	mv "fastqc_data.txt" ../txt/"$d.txt"; \
	cd ../; \
done
cd txt/
sed -s -n '4p;7p' *_SNAP_1_fastqc.txt > summary.txt
cd ../../../05_fastq_unmapped

#Taxonomical classification of unmapped files

conda activate kraken
for d in *_1.fq.gz; \
    do kraken2 --paired --threads 16 --db ~/reference_database/refseq_kraken/ --output ../07_kraken/${d%%_S*}_refseq --report ../07_kraken/${d%%_S*}_refseq.report $d ${d%%_1.*}_2.fq.gz; \
	kraken2 --paired --threads 16 --db ~/reference_database/uhgg_database/ --output ../07_kraken/${d%%_S*}_uhgg --report ../07_kraken/${d%%_S*}_uhgg.report $d ${d%%_1.*}_2.fq.gz; \
done
cd ../07_kraken

#Combine kraken results into one file

combine_kreports.py -r *_refseq.report -o kraken_refseq_summary.txt
combine_kreports.py -r *_uhgg.report -o kraken_uhgg_summary.txt

#Run bracken for microbial abundance estimation

for d in *_refseq.report; \
	do est_abundance.py -i $d -k ~/reference_database/refseq_kraken/database150mers.kmer_distrib -o ../08_bracken/${d%%.report}.bracken; \
done
for d in *_uhgg.report; \
	do est_abundance.py -i $d -k ~/reference_database/uhgg_database/database150mers.kmer_distrib -o ../08_bracken/${d%%.report}.bracken; \
done

#Combine bracken output

cd ../08_bracken
python2 ~/tools/Bracken-2.5/analysis_scripts/combine_bracken_outputs.py --files *_refseq.bracken -o bracken_refseq_summary.txt
python2 ~/tools/Bracken-2.5/analysis_scripts/combine_bracken_outputs.py --files *_uhgg.bracken -o bracken_uhgg_summary.txt
conda deactivate

#Run ARGs-OAP for acquired resistance genes detection on the reads

cd ../05_fastq_unmapped
argoap_pipeline_stageone_version2 -i ./ -m ../09_ARGs-OAP/metadata.txt -o ../09_ARGs-OAP -n 16 -z -s 
gzip *.fq
cd ../09_ARGs-OAP 
argoap_pipeline_stagetwo_version2 -i extracted.fa -m meta_data_online.txt -n 16 -o ARGs_OAP_output -b
cd ../05_fastq_unmapped
rm *.fa

#Run HUMAnN for metabolic pathway prediction

#conda activate humann
#cd ../10_HUMAnN
#mkdir output
#for d in *.fastq; \
#	do humann -i $d -o output/${d%%_SNAP*} --threads 16 --remove-temp-output; \
#done

#Retrieve, normalize and merge all HUMAnN ouptut

#cd output
#mkdir genefamilies pathcoverage pathabundance
#for d in *; \
#	do humann_renorm_table --input $d/*_genefamilies.tsv --output genefamilies/${d}_genefamilies_relab.tsv --units relab; \
#	cp $d/*_pathcoverage.tsv pathcoverage/; \
#	humann_renorm_table --input $d/*_pathabundance.tsv --output pathabundance/${d}_pathabundance_relab.tsv --units relab; \
#done
#humann_join_tables --input genefamilies --output humann_genefamilies.tsv --file_name genefamilies_relab
#humann_join_tables --input pathcoverage --output humann_pathcoverage.tsv --file_name pathcoverage
#humann_join_tables --input pathabundance --output humann_pathabundance.tsv --file_name pathabundance_relab
#conda deactivate

#Run mumame to detect mutations leading to AMR

mumame -d ~/reference_database/mumame/mutation_database -o ../11_mumame/mumame_output *.fq.gz

#Assembly of unmapped reads with metaSPAdes

conda activate spades
for d in *_1.fq.gz; \
    do spades.py -1 $d -2 ${d%%_1.*}_2.fq.gz -o ../12_metaSPAdes/${d%%_S*} --meta -t 16; \
done
cd ../12_metaSPAdes
mkdir scaffolds gfa
for d in *; \
    do cd $d; \
    mv "scaffolds.fasta" ../scaffolds/"$d.fasta"; \
    mv "assembly_graph_with_scaffolds.gfa" ../gfa/"$d.gfa"; \
    cd ../; \
done
conda deactivate

#Quality assessment of assemblies with metaQUAST

#cd scaffolds
#conda activate quast
#metaquast.py *.fasta -o ../../13_metaQUAST -t 16
#conda deactivate

#Metagenome annotation with prokka

cd scaffolds
conda activate prokka
for d in *.fasta; \
	do prokka --outdir ../../14_prokka/${d%%.fasta} --prefix ${d%%.fasta} $d --cpus 16 --metagenome --centre X --compliant; \
done
conda deactivate

#Move fasta and gff to BacPipe folder

cp *.fasta ../../15_BacPipe
cd ../../14_prokka
for d in *; \
	do cp $d/*.gff ../15_BacPipe; \
done


