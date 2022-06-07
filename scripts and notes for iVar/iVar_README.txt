iVar looping

Summary: 
1.	Initial prep
	o	 prep all needed files and packages
2.	Run first iVar script: ivar.sh 
	o	loops through the first half of iVar commands, taking the paired R1 and R2 fastq files, trimming and aligning them
3.	Run second iVar script: ivar2.sh (only if you ran technical replicates for each sample)
	o	loop through the second half of ivar commands, merging the technical replicates (_R) for each sample and calling variants
4.	Final important step: the iVar output files don’t have the sample ID inside the file, only in the file name. I used an awk command to insert a column in each variant table with the file name

Instructions:
1.	Initial prep:
	•	install iVar and dependencies: might want to use conda and these instructions: https://andersen-lab.github.io/ivar/html/installpage.html
	•	get viral reference sequence, .gff3 annotation file, .bed primer file, and all your sample fastq files (compress them if needed: gzip *.fastq)
	•	index reference sequence (get bwa if you need it):

bwa index MN985325.fasta

	•	prepare file lists
		o	InputFileList.txt lists all fastq file names you want to process
		o	InputFileList2.txt lists just the sample IDs without listing the technical replicates (ie just sample1 instead of sample1 and sample1_R)
	•	you will need to modify script to be specific to your directory and data
		o	reference files, output path, etc.
2.	Run ivar.sh
	•	to ensure script is recognized as a bash script, may need to run chmod u+x ivar.sh
	•	command to run: ./ivar.sh
3.	Run ivar2.sh
	•	make sure ivar2.sh is in the output folder you made, and cd into it
	•	command to run: ./ivar2.sh
4.	Add sample ID to all .tsv variant tables you have created:

awk -i inplace -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' *.tsv
