#!/bin/sh

umask 002;
export READSDIR=$1;
export THREADS=$2;
CURRWORKDIR=`pwd`;

#If no thread count was passed, used all available threads on the cluster
if [ -z "${THREADS}" ];
	then export THREADS=`grep -c ^processor /proc/cpuinfo`;
fi

#If no reads directory was passed, sets it to the current directory
if [ -z "${READSDIR}" ];
	then export READSDIR=".";
fi

#Verify the existence of the temporary directory
if [ -d $(readlink -e $TMPDIR) ];
	then echo "Temporary Directory: ${TMPDIR}";
	else echo "Temporary Directory does not exist";
fi

#go into Reads dir
export READSDIR=`readlink -e $READSDIR`;
cd ${READSDIR};
export POOL=`echo ${READSDIR} | cut -f1-6 -d "/"`;
export PROJECTID=`basename $(readlink -e ${READSDIR}/..) | sed 's/WorkDir//g'`;
echo ${PROJECTID};

## prepping the fastqs for processing ##
#decompress fastqs to temporary directory
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.1.fq.bz2 | sed 's/1:N:0:.*/1:N:0:/g' > ${TMPDIR}/{}.1.fq";
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.2.fq.bz2 | sed 's/2:N:0:.*/3:N:0:/g' > ${TMPDIR}/{}.2.fq";

#merge both reads
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Merge:'; usearch70 -fastq_mergepairs ${TMPDIR}/{}.1.fq -reverse ${TMPDIR}/{}.2.fq -fastq_allowmergestagger -fastq_minovlen 50 -fastqout ${TMPDIR}/{}.Merged.fq; echo";
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo {}; perl ~mcwong/mergeReads.pl ${TMPDIR}/{}.1.fq ${TMPDIR}/{}.2.fq ${TMPDIR}/{}.Merged.fq > ${TMPDIR}/{}.Temp.fq";

#seperate into raw and standard fastq files for each read
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} "echo {}; usearch70 -fastq_filter ${TMPDIR}/{}.Temp.fq -fastqout ${TMPDIR}/{}.FilteredRaw.fq -eeout -fastq_minlen 200";
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} "echo {}; usearch70 -fastq_filter ${TMPDIR}/{}.Temp.fq -relabel "{}_" -fastqout ${TMPDIR}/{}.Filtered.fq -eeout -fastq_minlen 200";
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} "echo {}; perl ~mcwong/filterSeqs.pl ${TMPDIR}/{}.FilteredRaw.fq > ${TMPDIR}/{}.FilteredRaw2.fq";
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} "echo {}; perl ~mcwong/filterSeqs.pl ${TMPDIR}/{}.Filtered.fq > ${TMPDIR}/{}.Filtered2.fq";

#create the seqs.fq files
cat ${TMPDIR}/*.FilteredRaw2.fq > ${TMPDIR}/seqs.raw.fq &
cat ${TMPDIR}/*.Filtered2.fq > ${TMPDIR}/seqs.fq &
wait;

#run bowtie to strip out PhiX
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.raw.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.raw.filtered.fq -S /dev/null 2>${READSDIR}/../../Logs/phix.raw.bleed.txt;
bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.fq --end-to-end --very-sensitive --reorder -p ${THREADS} --un ${TMPDIR}/seqs.filtered.fq -S /dev/null 2>${READSDIR}/../../Logs/phix.bleed.txt;

#Make archived copies of the fastq files
cat ${READSDIR}/../SampleList | xargs -I {} mkdir ${READSDIR}/{};
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} 'cat ${TMPDIR}/{}.1.fq | pbzip2 -p1 -c > ${READSDIR}/{}/{}.1.fq.bz2; cat ${TMPDIR}/{}.2.fq | pbzip2 -p1 -c > ${READSDIR}/{}/{}.2.fq.bz2';
mkdir ${READSDIR}/../split_libraries;
fq2fa ${TMPDIR}/seqs.filtered.fq ${READSDIR}/../split_libraries/seqs.fna;

mkdir ${READSDIR}/../../Deliverables;

#Uparse for both 97% and 99% similarity
for j in {97,99}; 
do echo $j;
done | parallel -I {} '
if [ {} -eq 97 ];
then list="0.4 0.7 1.0 1.3 1.6 1.9 2.2 2.5 2.8 3.0";
else list="0.4 0.7 1.0";
fi
mkdir ${TMPDIR}/uparse{};
cd ${TMPDIR}/uparse{};
usearch70 -derep_fulllength ${READSDIR}/../split_libraries/seqs.fna -output ${TMPDIR}/uparse{}/derep.fna -sizeout -uc ${TMPDIR}/uparse{}/derep.uc 2>&1;
usearch70 -sortbysize ${TMPDIR}/uparse{}/derep.fna -output ${TMPDIR}/uparse{}/sorted.fa -minsize 2;
cp ${TMPDIR}/uparse{}/sorted.fa ${TMPDIR}/uparse{}/temp.fa;
for i in $list;
do usearch70 -cluster_otus ${TMPDIR}/uparse{}/temp.fa -otus ${TMPDIR}/uparse{}/temp1.fa -otu_radius_pct $i -uc ${TMPDIR}/uparse{}/cluster_$i.uc -fastaout ${TMPDIR}/uparse{}/clustering.$i.fasta.out;
cat ${TMPDIR}/uparse{}/clustering.$i.fasta.out | grep "^>" | grep chimera | sed "s/^>//g" | sed -re "s/;n=.*up=/\t/g" | sed "s/;$//g" | tee -a ${TMPDIR}/uparse{}/chimeras.txt > ${TMPDIR}/uparse{}/chimeras.$i.txt;
cat ${TMPDIR}/uparse{}/clustering.$i.fasta.out | grep "^>" > ${TMPDIR}/uparse{}/uparse{}ref.decisions.$i.txt;
rm ${TMPDIR}/uparse{}/clustering.$i.fasta.out;
mv ${TMPDIR}/uparse{}/temp1.fa ${TMPDIR}/uparse{}/temp.fa;
done;
mv ${TMPDIR}/uparse{}/temp.fa ${TMPDIR}/uparse{}/otus1.fa;
usearch80 -uchime_ref ${TMPDIR}/uparse{}/otus1.fa -db /gpfs1/db/GenbankEukaryotes/ITSDbV8.udb -strand plus -nonchimeras ${TMPDIR}/uparse{}/otus.fa -uchimeout ${TMPDIR}/uparse{}/uchimeref.uc;
cat ${TMPDIR}/uparse{}/uchimeref.uc | cut -f2,17 | grep -v "Y$" | cut -f1 | /users/mcwong/getSeq ${TMPDIR}/uparse{}/otus1.fa > ${TMPDIR}/uparse{}/otus.fa;
usearch80 -usearch_global ${TMPDIR}/uparse{}/otus.fa -db /gpfs1/db/GenbankEukaryotes/ITSDbV8.udb -maxaccepts 0 -maxrejects 0 -strand plus -id .97 -query_cov .95 -threads $THREADS -uc ${TMPDIR}/uparse{}/{}.centroids.uc -gapopen 5.0I/0.0E -gapext 1.0I/0.0E -top_hits_only;
perl ${WONGGITREPO}/ITS_workflows/cleanHitsTableITS.pl ${TMPDIR}/uparse{}/{}.centroids.uc /gpfs1/db/GenbankEukaryotes/TaxaITSDb.txt > ${TMPDIR}/uparse{}/{}.clean.uc;
mv ${TMPDIR}/uparse{}/NewTaxa.txt ${TMPDIR}/uparse{}/NewTaxa.{}.txt;
cat ${TMPDIR}/uparse{}/{}.centroids.uc | grep "*$" | cut -f9 | ~mcwong/getSeq ${TMPDIR}/uparse{}/otus.fa > ${TMPDIR}/uparse{}/missed.fa;
usearch80 -usearch_global ${TMPDIR}/uparse{}/missed.fa -db /gpfs1/db/GenbankEukaryotes/ITSDbV8.udb -maxaccepts 0 -maxrejects 0 -strand both -id .80 -query_cov .95 -threads $THREADS -uc ${TMPDIR}/uparse{}/missed.centroids.uc -gapopen 5.0I/0.0E -gapext 1.0I/0.0E -top_hits_only;
cat ${TMPDIR}/uparse{}/{}.clean.uc | grep -v "*$" | sed -re "s/\t1$//g" | cut -f10 | sed "s/.*/*/g" > ${TMPDIR}/uparse{}/second;
cat ${TMPDIR}/uparse{}/{}.clean.uc | grep -v "*$" | sed -re "s/\t1$//g" | cut -f1-9 > ${TMPDIR}/uparse{}/first;
paste first second >> ${TMPDIR}/uparse{}/missed.centroids.uc;
perl ${WONGGITREPO}/ITS_workflows/cleanHitsTableITS.pl ${TMPDIR}/uparse{}/missed.centroids.uc /gpfs1/db/GenbankEukaryotes/TaxaITSDb.txt > ${TMPDIR}/uparse{}/missed.clean.uc;
mv ${TMPDIR}/uparse{}/NewTaxa.txt ${TMPDIR}/uparse{}/NewTaxa.missed.{}.txt;
cat ${TMPDIR}/uparse{}/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | ~mcwong/getSeq ${TMPDIR}/uparse{}/derep.fna > ${TMPDIR}/uparse{}/singletons.fna;
usearch70 -usearch_global ${TMPDIR}/uparse{}/singletons.fna -db ${TMPDIR}/uparse{}/sorted.fa -id .99 -uc ${TMPDIR}/uparse{}/singletons2otus.uc -strand plus -threads $THREADS -maxaccepts 32 -maxrejects 128 -query_cov .85 -wordlength 12;
perl ${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_* ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/{}.clean.uc --derep ${TMPDIR}/uparse{}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --taxonomy ${TMPDIR}/uparse{}/NewTaxa.{}.txt --uchime ${TMPDIR}/uparse{}/uchimeref.uc;
mv ${TMPDIR}/uparse{}/otu_table.biom ${TMPDIR}/uparse{}/otu_table.{}.biom;
mv ${TMPDIR}/uparse{}/reads2otus.txt ${TMPDIR}/uparse{}/reads2otus.{}.txt;
perl ${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_* ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/missed.clean.uc --perc --otus ${TMPDIR}/uparse{}/missed.clean.uc --derep ${TMPDIR}/uparse{}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --taxonomy ${TMPDIR}/uparse{}/NewTaxa.missed.{}.txt; 
mv ${TMPDIR}/uparse{}/otu_table.biom ${TMPDIR}/uparse{}/otu_table_missed.{}.biom;
mv ${TMPDIR}/uparse{}/reads2otus.txt ${TMPDIR}/uparse{}/reads2otus.missed.{}.txt;
perl ${WONGGITREPO}/ITS_workflows/make_uc_for_unmapped.pl ${TMPDIR}/uparse{}/missed.clean.uc > ${TMPDIR}/uparse{}/missed.unmapped.uc;
mv ${TMPDIR}/uparse{}/NewTaxa.unmapped.txt ${TMPDIR}/uparse{}/NewTaxa.unmapped.{}.txt;
perl ${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_* ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/missed.unmapped.uc --derep ${TMPDIR}/uparse{}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --taxonomy ${TMPDIR}/uparse{}/NewTaxa.unmapped.{}.txt;
mv reads2otus.txt ${TMPDIR}/uparse{}/reads2otus.unmapped.{}.txt;
mv ${TMPDIR}/uparse{}/otu_table.biom ${TMPDIR}/uparse{}/otu_table_unmapped.{}.biom;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table.{}.biom -o ${TMPDIR}/uparse{}/stats.otu_table.{}.txt;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table_missed.{}.biom -o ${TMPDIR}/uparse{}/stats.otu_table_missed.{}.txt;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table_unmapped.{}.biom -o ${TMPDIR}/uparse{}/stats.otu_table_unmapped.{}.txt;
cat ${READSDIR}/../split_libraries/seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c > ${TMPDIR}/uparse{}/Stats.{}.MergedReads.txt;
cat ${TMPDIR}/uparse{}/stats.otu_table.{}.txt | tail -n +17 | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${TMPDIR}/uparse{}/Stats.{}.MappedReads.txt;
cat ${TMPDIR}/uparse{}/stats.otu_table_missed.{}.txt | tail -n +17 | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${TMPDIR}/uparse{}/Stats.{}.MappedReads.missed.txt;
cat ${TMPDIR}/uparse{}/stats.otu_table_unmapped.{}.txt | tail -n +17 | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${TMPDIR}/uparse{}/Stats.{}.MappedReads.unmapped.txt;
perl ${WONGGITREPO}/ITS_workflows/StatsComparisonMergedVsMapped.pl ${TMPDIR}/uparse{}/Stats.{}.MergedReads.txt ${TMPDIR}/uparse{}/Stats.{}.MappedReads.{}.txt ${TMPDIR}/uparse{}/Stats.{}.MappedReads.missed.txt ${TMPDIR}/uparse{}/Stats.{}.MappedReads.unmapped.txt > ${READSDIR}/../Stats.{}.Combined.txt;
cat ${TMPDIR}/uparse{}/missed.centroids.uc | grep "^N" | cut -f9 | ~mcwong/getSeq ${TMPDIR}/uparse{}/missed.fa > ${TMPDIR}/uparse{}/OTUsFailedToMap2ndTime.fa;
for i in `ls ${TMPDIR}/uparse{}/ | grep .biom`;
	do cp ${TMPDIR}/uparse{}/${i} ${READSDIR}/../../Deliverables/${PROJECTID}.${i};
done;
tar -cvf ${TMPDIR}/uparse{}.tar.bz2 -C ${TMPDIR}/uparse{} --use-compress-program=pbzip2;
mv ${TMPDIR}/uparse{}.tar.bz2 ${READSDIR}/..;
cd ${READSDIR};
chmod -R 777 ${READSDIR}/../uparse{}.tar.bz2;
' 1>> ${READSDIR}/../../Logs/Parallel.log 2>> ${READSDIR}/../../Logs/Parallel.err &
bigJob=`jobs -p`;

## construct the deliverables ##

#zip up Read and Merge Filetar.bz2s
cat ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.1.fq" > ${TMPDIR}/Read1.fq;
cat ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.2.fq" > ${TMPDIR}/Read2.fq;
wait $pidlist;
mv ${TMPDIR}/seqs.raw.filtered.fq ${TMPDIR}/MergedRaw.fq;
mv ${TMPDIR}/seqs.filtered.fq ${TMPDIR}/MergedStandard.fq;
pbzip2 -p${THREADS} ${TMPDIR}/Read1.fq;
pbzip2 -p${THREADS} ${TMPDIR}/Read2.fq;
pbzip2 -p${THREADS} ${TMPDIR}/MergedRaw.fq;
pbzip2 -p${THREADS} ${TMPDIR}/MergedStandard.fq;

#recover barcodes for deliverables
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/Read1.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADS} -c > ${TMPDIR}/RawReadsBarcodes.fq.bz2;
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/MergedRaw.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADS} -c > ${TMPDIR}/MergedRawBarcodes.fq.bz2;
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/MergedStandard.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADS} -c > ${TMPDIR}/MergedStandardBarcodes.fq.bz2;

#make the deliverables directory and move files into it

for i in `ls ${TMPDIR}/*.fq.bz2`; 
	do name=`basename $i`; 
	cp $i ${READSDIR}/../../Deliverables/${PROJECTID}.${name};
done;

cp ${TMPDIR}/uparse99/otu_table.biom ${READSDIR}/../../Deliverables/${PROJECTID}.99.otu_table.biom;
cp ${READSDIR}/../SampleList ${READSDIR}/../../Deliverables/${PROJECTID}.SampleList;
cat ${READSDIR}/../../samplesheet.${PROJECTID}.csv | grep -f ${READSDIR}/../SampleList | cut -f3,5 -d "," | tr "," "\t" > ${READSDIR}/../../Deliverables/${PROJECTID}.SampleSheet.txt;
head -1 ${GITREPO}/Miscellaneous/IlluminaHeaderExample > ${READSDIR}/../../Deliverables/${PROJECTID}.ExampleQiimeMappingFile.txt;
tail -n+1 ${READSDIR}/../../Deliverables/${PROJECTID}.SampleSheet.txt | sed -re 's/(.*)\t(.*)/\1\t\2\tGGACTACHVGGGTWTCTAAT\tGTGCCAGCMGCCGCGGTAA\t\1/g' >> ${READSDIR}/../../Deliverables/${PROJECTID}.ExampleQiimeMappingFile.txt;
cat ${GITREPO}/Miscellaneous/Versions.txt > ${READSDIR}/../../Deliverables/${PROJECTID}.SoftwareVersionInformation.txt;

#return to working directory when script was launched
cd $CURRWORKDIR;

#exit without error status once completed
exit 0;
