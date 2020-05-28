##Codes are from [Taj Azarian Gist Page](https://gist.githubusercontent.com/EpiDemos82/ad030f26f1df2d5e1a7c799f0987ca02/) 

##GENERAL TEXT OR FILE MANIPULATION

#Find lines in a list (e.g. file names) that are not present in another list
#This is good for checking whether downstream files are present (i.e. pipeline ran susscessfully)

`comm -23 <(sort All.txt) <(sort Finished.txt)`

#Looping over anything
`for f in $(cat names.txt); do whatever to ${f}; done`

#renaming file extensions using bash code
`for i in *.old; do mv -- "$i" "${i%.old}.new"; done`

#Rename any text string in a file using a tab delimitted lookup table where col1=old name and col2=new name
#This is useful for renaming taxa in trees
`awk 'NR==FNR{a[$1]=$2;next}{ for (i in a) gsub(i ":",a[i] ":")}1' names.txt RAxML.tre > RAxML_renamed.tree`


#You can chance ":" to ">" if you are working with FASTA files, but here is a MUCH faster method using Python (Thanks Karel Brinda)

`python3 -c 'import csv;a=dict(csv.reader(open("names.txt"),delimiter="\t"));lines=open("input.fasta").read().split("\n");lines=map(lambda x: ">"+a[x[1:]] if x and x[0]==">" else x,lines);print("\n".join(lines));' > output.fasta`


#Rename file names based on a lookup table (see above)  - useful for renaming fastq/fasta/gff files
`perl -lane 'rename("$F[0].fasta", "$F[1].fasta")' rename.txt`
`perl -lane 'rename("$F[0].fastq.gz", "$F[1].fastq.gz")' rename.txt`

#Creating a file with a list of stem names (for example fasta files or fastq files)
`ls *.fasta | awk -F.fasta '{print $1}' > names.txt`
`ls *_1.fastq.gz | awk -F_1.fastq.gz '{print $1}' > names.txt`
`ls *.fasta | sed 's/.fasta//' > names.txt # another alternative`

#Radomly subsample (e.g.200) files from a directory and copy them to another directory 
`shuf -zen200 /source/*.ext | xargs -0 cp -t /destination`

##WORKING WITH FASTA FILES

#counting number of sequences in a fasta file
`grep -c "^>" file.fa`

#add something to end of all header lines
`sed 's/>.*/&WHATEVERYOUWANT/' file.fa > outfile.fa`

#Renaming contig names in a de novo assembly using a simple count
#Can also use fastx_renamer from the fastX toolkit
`awk '/^>/{print ">" ++i; next}{print}' < in.fasta > renamed_contigs.fasta`

#Extract taxa ids
`grep -o -E "^>\w+" file.fasta | tr -d ">"`

#linearize your sequences (i.e. remove the sequence wrapping)
`sed -e 's/\(^>.*$\)/#\1#/' file.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d'`

#Remove duplicated sequences
`sed -e '/^>/s/$/@/' -e 's/^>/#/' file.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort -u -t $'\t' -f -k 2,2  | sed -e 's/^/>/' -e 's/\t/\n/'`

#extracting a taxa group from a fasta file using samtools (must have samtools installed)
`xargs samtools faidx all_taxa.fasta < list_of_taxa_to_sub.txt > sub_taxa.fasta`

#Calculate the proptions of Ns in your multi-fasta alignment using Seqtk
`seqtk comp input.fa | awk '{x=$3+$4+$5+$6;y=$2;print $1,y-x,y,(y-x)/y}'`

#Extract fasta sequence from gff file (requires gff to have sequence at end)
`for i in $(ls *.gff); do cat ${i} | sed '1,/##FASTA/d' > ${i}.fasta; done`

#Counting CDSs (annotated coding sequences) in gff file
`for i in $(ls *.gff); do echo ${i}; cat ${i} | grep -c CDS; done > CDS.txt`



##SPECIAL TASKS FOR FASTQ OR FASTA FILES

#Renaming the ridiculous fastq file names that Illumina machines spit out
#When you confirm the chnage is correct, remove "-n" to make the change permanent - note: learn how to use "rename"!
`rename -n 's/_S.*R1_001/_1/' *.fastq.gz`
`rename -n 's/_S.*R2_001/_2/' *.fastq.gz`

#Evaluate assembly statistics - you can also use quast if you want more extensive output
`conda install assembly-stats`
`assembly-stats -t assembly.fasta` 

#Print the number of sequences in a gzipped fastq file (.fastq.gz)
`zcat your.fastq.gz | echo $(('wc -l'/4))`

#Print the number of bases in a gzipped fastq file (.fastq.gz)
`zcat file.fq.gz | paste - - - - | cut -f2 | wc -c`

#Print the number of bases in a fasta file
`awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' fasta_file.fasta`

#Subsample a FASTQ (change K for the number of reads to sample)
`cat ${1}_1.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' | awk -v k=300000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > sub.fastq}'`

#Convert FASTQ to FASTA
`cat file.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > file.fasta`

#Using Sed
`sed '/^@/!d;s//>/;N' sample1.fq > sample1.fa`

#remove the description information from a fasta file and just keep the identifier
`perl -p -i -e 's/>(.+?) .+/>$1/g' sample1.fa`

#Generate random reads from a draft genomes using bbmap random reads
`bbmap/randomreads.sh ref=draft_genome.fasta reads=2000000 length=150 illuminanames=t paired=t interleaved=f maxq=36 midq=32 minq=28 out=reads`

#Subsampline two taxa from the pan-genome-sequences ROARY output folder (I use this to blast against if looking for specific COGS)
`ls *.fa.aln > names.txt` #Creates a name of all COG alignments 
`for f in $(cat names.txt); do cat ${f} | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | head -n 2 | awk '{printf("%s\n%s\n",$1,$2)}' > ${f}.fasta; done #Subsampling 2 taxa from each COG file`
`find . -path './*.fa.aln.fasta' -prune -type f -exec cat {} + > combined.fasta` #Concatenating those sequences

#extracting a genomic region from a de novo assembly using a query - for example if you want to pull out the CPS loci in pneumo
`mkdir region`
`blastn -subject denovo_assembly.fasta -query region_of_interest.fasta -outfmt 6 -perc_identity 90 -max_target_seqs 3 -evalue 0.1e-200 > blastout.txt; sort -k1,1 -k12,12nr -k11,11n  blastout.txt | awk '{print $2}' > contigs.txt; awk '!seen[$1]++' contigs.txt > contigs.txt; xargs samtools faidx denovo_assembly.fasta < contigs.txt > ./region/region_of_interest_from_denovo.fasta; done`

#This can also be looped
`or f in $(cat names.txt); do blastn -subject ${f}.fasta -query region_of_interest.fasta -outfmt 6 -perc_identity 90 -max_target_seqs 3 -evalue 0.1e-200 > blastout.txt; sort -k1,1 -k12,12nr -k11,11n  blastout.txt | awk '{print $2}' > contigs.txt; awk '!seen[$1]++' contigs.txt > contigs.txt; xargs samtools faidx ${f}.fasta < contigs.txt > ./region/${f}_region_of_interest.fasta; done`
  
  
  
#PHYLOGENETICS
#Using SNP-sites to extract only polymorphic sites from a Roary output file
`snp-sites -o SNPs.fasta -c core_gene_alignment.aln`

#RAxML with ascertainment bias correction and 100 bs replicates (for SNP Alignments)
`raxmlHPC-AVX -f a -p 12345 -x 12345 -# 100 -n SNPs -s SNPs.fasta -m ASC_GTRGAMMA --asc-corr=lewis` 

#SNP/Ascertainment bias correction with IQ-TREE and automatic model selection (this method is almost as accurate as RAxML as less computationally taxing) - bootstrapping can be added as well
`iqtree -s SNPs.fasta -m MFP+ASC -nt AUTO`




`#Example pipeline
 #!/bin/bash

echo "Starting FastQC analysis on sample "${1}
 fastqc ${1}_1.fastq.gz ${1}_2.fastq.gz
echo "Finished FastQC analysis on sample "${1}

echo "Starting Kraken with sample "${1}
 kraken2 --db minikraken2_v1_8GB --output ${1}.output --report ${1}.report --paired --use-names --gzip-compressed ${1}_1.fastq.gz ${1}_2.fastq.gz
echo "Finished kraken on sample "${1}

echo "Simplifying kraken output for sample "${1}
 awk '$1 >= 5 {print}' sample1.report | awk '$4 == "S" {print}' > ${1}.short.report
 awk 'FNR == 1 {print $6}' ${1}.short.report > genus
 awk 'FNR == 1 {print $7}' ${1}.short.report > species
echo "Finished simplifying output for sample "${1}

echo "Cleaning up with trim"
 trimmomatic PE -threads 1 -phred33 -trimlog ${1}.log ${1}_1.fastq.gz ${1}_2.fastq.gz ${1}_1.fq.gz ${1}_unpaired.fq.gz ${1}_2.fq.gz ${1}_unpaired.fq.gz SLIDINGWINDOW:10:20 MINLEN:31 TRAILING:20
echo "Finished cleaning up with trim"

echo "Assemble with Spades"
 spades.py -1 ${1}_1.fq.gz -2 ${1}_2.fq.gz -o ${1} -k 55 --only-assembler -t 1
echo "Finished assembly with Spades"`

  
