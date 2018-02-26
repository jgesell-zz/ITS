#/bin/sh

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
export THREADSPLIT=$[THREADS / 2];
echo ${PROJECTID};

#Create a readcount file if there i not one already
if [ ! -f ${READSDIR}/../../${PROJECTID}.barcodeCounts.txt ];
then for i in `cat ${READSDIR}/../SampleList`;
do count=`cat ${TMPDIR}/${i}.1.fq | wc -l`;
count=$[count / 4];
echo -e "${i}\t${count}";
done > ${READSDIR}/../../${PROJECTID}.barcodeCounts.txt;
fi;

## prepping the fastqs for processing ##
#decompress fastqs to temporary directory
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.1.fq.bz2 | sed 's/1:N:0:.*/1:N:0:/g' > ${TMPDIR}/{}.1.fq";
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "bzcat {}.2.fq.bz2 | sed 's/2:N:0:.*/3:N:0:/g' > ${TMPDIR}/{}.2.fq";

#merge both reads
cat ${READSDIR}/../SampleList | parallel -j${THREADS} -I {} "echo '{} Merge:'; usearch81 -fastq_mergepairs ${TMPDIR}/{}.1.fq -reverse ${TMPDIR}/{}.2.fq -threads 1 -fastq_minovlen 50 -fastqout ${TMPDIR}/{}.Merged.fq; echo";

#seperate into raw and standard fastq files for each read
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} "echo {}; usearch81 -fastq_filter ${TMPDIR}/{}.Merged.fq -fastqout ${TMPDIR}/{}.FilteredRaw.fq -fastq_minlen 200 -fastq_maxee_rate .005 -threads 1";
cat ${READSDIR}/../SampleList | parallel -j$THREADS -I {} "echo {}; usearch81 -fastq_filter ${TMPDIR}/{}.Merged.fq -relabel \"{}_\" -fastqout ${TMPDIR}/{}.Filtered.fq -fastq_minlen 200 -fastq_maxee_rate .005 -threads 1 && rm ${TMPDIR}/{}.Merged.fq";

#create the seqs.fq files
(cat ${TMPDIR}/*.FilteredRaw.fq > ${TMPDIR}/seqs.raw.fq && rm ${TMPDIR}/*.FilteredRaw.fq) &
(cat ${TMPDIR}/*.Filtered.fq > ${TMPDIR}/seqs.fq && rm ${TMPDIR}/*.Filtered.fq) &
wait;

#run bowtie to strip out PhiX
(bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.raw.fq --end-to-end --very-sensitive --reorder -p ${THREADSPLIT} --un ${TMPDIR}/Merged_Reads_Raw.fq -S /dev/null 2>${READSDIR}/../../Logs/phix.raw.bleed.txt && rm ${TMPDIR}/seqs.raw.fq) &
(bowtie2 -x ${PHIXDB} -U ${TMPDIR}/seqs.fq --end-to-end --very-sensitive --reorder -p ${THREADSPLIT} --un ${TMPDIR}/Merged_Reads.fq -S /dev/null 2>${READSDIR}/../../Logs/phix.bleed.txt && rm ${TMPDIR}/seqs.fq ) &

#Setup for iterative uparse
mkdir -p ${READSDIR}/../split_libraries &
mkdir -p ${READSDIR}/../../Deliverables &
wait;

fq2fa ${TMPDIR}/Merged_Reads.fq ${READSDIR}/../split_libraries/seqs.fna;
usearch70 -derep_fulllength ${READSDIR}/../split_libraries/seqs.fna -output ${TMPDIR}/derep.fna -sizeout -uc ${TMPDIR}/derep.uc 2>&1;
usearch70 -sortbysize ${TMPDIR}/derep.fna -output ${TMPDIR}/sorted.fa -minsize 2;
#Uparse for both 97% and 99% similarity

for j in {97,99};
do echo $j;
done | parallel -I {} '
mkdir ${TMPDIR}/uparse{};
cd ${TMPDIR}/uparse{};
cp ${TMPDIR}/sorted.fa ${TMPDIR}/uparse{}/temp.fa;
if [ {} -eq 97 ];
then for i in {0.4,0.7,1.0,1.3,1.6,1.9,2.2,2.5,2.8,3.0};
do usearch70 -cluster_otus ${TMPDIR}/uparse{}/temp.fa -otus ${TMPDIR}/uparse{}/temp1.fa -otu_radius_pct ${i} -uc ${TMPDIR}/uparse{}/cluster_${i}.uc -fastaout ${TMPDIR}/uparse{}/clustering.${i}.fasta.out;
cat ${TMPDIR}/uparse{}/clustering.${i}.fasta.out | grep "^>" | grep chimera | sed "s/^>//g" | sed -re "s/;n=.*up=/\t/g" | sed "s/;$//g" | tee -a ${TMPDIR}/uparse{}/chimeras.txt > ${TMPDIR}/uparse{}/chimeras.${i}.txt;
cat ${TMPDIR}/uparse{}/clustering.${i}.fasta.out | grep "^>" > ${TMPDIR}/uparse{}/uparse{}ref.decisions.${i}.txt;
rm ${TMPDIR}/uparse{}/clustering.${i}.fasta.out;
mv ${TMPDIR}/uparse{}/temp1.fa ${TMPDIR}/uparse{}/temp.fa;
done;
else for i in {0.4,0.7,1.0};
do usearch70 -cluster_otus ${TMPDIR}/uparse{}/temp.fa -otus ${TMPDIR}/uparse{}/temp1.fa -otu_radius_pct ${i} -uc ${TMPDIR}/uparse{}/cluster_${i}.uc -fastaout ${TMPDIR}/uparse{}/clustering.${i}.fasta.out;
cat ${TMPDIR}/uparse{}/clustering.${i}.fasta.out | grep "^>" | grep chimera | sed "s/^>//g" | sed -re "s/;n=.*up=/\t/g" | sed "s/;$//g" | tee -a ${TMPDIR}/uparse{}/chimeras.txt > ${TMPDIR}/uparse{}/chimeras.${i}.txt;
cat ${TMPDIR}/uparse{}/clustering.${i}.fasta.out | grep "^>" > ${TMPDIR}/uparse{}/uparse{}ref.decisions.${i}.txt;
rm ${TMPDIR}/uparse{}/clustering.${i}.fasta.out;
mv ${TMPDIR}/uparse{}/temp1.fa ${TMPDIR}/uparse{}/temp.fa;
done;
fi;
mv ${TMPDIR}/uparse{}/temp.fa ${TMPDIR}/uparse{}/otus1.fa;
usearch80 -uchime_ref ${TMPDIR}/uparse{}/otus1.fa -db /gpfs1/db/GenbankEukaryotes/ITSDbV8.udb -strand plus -nonchimeras ${TMPDIR}/uparse{}/otus.fa -uchimeout ${TMPDIR}/uparse{}/uchimeref.uc;
cat ${TMPDIR}/uparse{}/uchimeref.uc | cut -f2,17 | grep -v "Y$" | cut -f1 | /users/mcwong/getSeq ${TMPDIR}/uparse{}/otus1.fa > ${TMPDIR}/uparse{}/otus.fa;
usearch80 -usearch_global ${TMPDIR}/uparse{}/otus.fa -db /gpfs1/db/GenbankEukaryotes/ITSDbV8.udb -maxaccepts 0 -maxrejects 0 -strand plus -id .97 -query_cov .95 -threads ${THREADSPLIT} -uc ${TMPDIR}/uparse{}/{}.centroids.uc -gapopen 5.0I/0.0E -gapext 1.0I/0.0E -top_hits_only;
perl ${WONGGITREPO}/ITS_workflows/cleanHitsTableITS.pl ${TMPDIR}/uparse{}/{}.centroids.uc /gpfs1/db/GenbankEukaryotes/TaxaITSDb.txt > ${TMPDIR}/uparse{}/{}.clean.uc;
mv ${TMPDIR}/uparse{}/NewTaxa.txt ${TMPDIR}/uparse{}/NewTaxa.{}.txt;
cat ${TMPDIR}/uparse{}/{}.centroids.uc | grep "*$" | cut -f9 | ~mcwong/getSeq ${TMPDIR}/uparse{}/otus.fa > ${TMPDIR}/uparse{}/missed.fa;
usearch80 -usearch_global ${TMPDIR}/uparse{}/missed.fa -db /gpfs1/db/GenbankEukaryotes/ITSDbV8.udb -maxaccepts 0 -maxrejects 0 -strand both -id .80 -query_cov .95 -threads ${THREADSPLIT} -uc ${TMPDIR}/uparse{}/missed.centroids.uc -gapopen 5.0I/0.0E -gapext 1.0I/0.0E -top_hits_only;
cat ${TMPDIR}/uparse{}/{}.clean.uc | grep -v "*$" | sed -re "s/\t1$//g" | cut -f10 | sed "s/.*/*/g" > ${TMPDIR}/uparse{}/second;
cat ${TMPDIR}/uparse{}/{}.clean.uc | grep -v "*$" | sed -re "s/\t1$//g" | cut -f1-9 > ${TMPDIR}/uparse{}/first;
paste first second >> ${TMPDIR}/uparse{}/missed.centroids.uc;
perl ${WONGGITREPO}/ITS_workflows/cleanHitsTableITS.pl ${TMPDIR}/uparse{}/missed.centroids.uc /gpfs1/db/GenbankEukaryotes/TaxaITSDb.txt > ${TMPDIR}/uparse{}/missed.clean.uc;
mv ${TMPDIR}/uparse{}/NewTaxa.txt ${TMPDIR}/uparse{}/NewTaxa.missed.{}.txt;
cat ${TMPDIR}/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | ~mcwong/getSeq ${TMPDIR}/derep.fna > ${TMPDIR}/uparse{}/singletons.fna;
usearch70 -usearch_global ${TMPDIR}/uparse{}/singletons.fna -db ${TMPDIR}/sorted.fa -id .99 -uc ${TMPDIR}/uparse{}/singletons2otus.uc -strand plus -threads ${THREADSPLIT} -maxaccepts 32 -maxrejects 128 -query_cov .85 -wordlength 12;
perl ${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_* ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/{}.clean.uc --derep ${TMPDIR}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --taxonomy ${TMPDIR}/uparse{}/NewTaxa.{}.txt --uchime ${TMPDIR}/uparse{}/uchimeref.uc;
mv ${TMPDIR}/uparse{}/otu_table.biom ${TMPDIR}/uparse{}/otu_table.{}.biom;
mv ${TMPDIR}/uparse{}/reads2otus.txt ${TMPDIR}/uparse{}/reads2otus.{}.txt;
perl ${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_* ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/missed.clean.uc --perc --otus ${TMPDIR}/uparse{}/missed.clean.uc --derep ${TMPDIR}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --taxonomy ${TMPDIR}/uparse{}/NewTaxa.missed.{}.txt; 
mv ${TMPDIR}/uparse{}/otu_table.biom ${TMPDIR}/uparse{}/otu_table_missed.{}.biom;
mv ${TMPDIR}/uparse{}/reads2otus.txt ${TMPDIR}/uparse{}/reads2otus.missed.{}.txt;
perl ${WONGGITREPO}/ITS_workflows/make_uc_for_unmapped.pl ${TMPDIR}/uparse{}/missed.clean.uc > ${TMPDIR}/uparse{}/missed.unmapped.uc;
mv ${TMPDIR}/uparse{}/NewTaxa.unmapped.txt ${TMPDIR}/uparse{}/NewTaxa.unmapped.{}.txt;
perl ${WONGGITREPO}/16S_workflows/resolveIterativeUparse.pl ${TMPDIR}/uparse{}/cluster_* ${TMPDIR}/uparse{}/singletons2otus.uc ${TMPDIR}/uparse{}/missed.unmapped.uc --derep ${TMPDIR}/derep.uc --chimeras ${TMPDIR}/uparse{}/chimeras.txt --taxonomy ${TMPDIR}/uparse{}/NewTaxa.unmapped.{}.txt;
mv reads2otus.txt ${TMPDIR}/uparse{}/reads2otus.unmapped.{}.txt;
mv ${TMPDIR}/uparse{}/otu_table.biom ${TMPDIR}/uparse{}/otu_table_unmapped.{}.biom;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table.{}.biom -o ${TMPDIR}/uparse{}/stats.otu_table.{}.txt;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table_missed.{}.biom -o ${TMPDIR}/uparse{}/stats.otu_table_missed.{}.txt;
biom summarize-table -i ${TMPDIR}/uparse{}/otu_table_unmapped.{}.biom -o ${TMPDIR}/uparse{}/stats.otu_table_unmapped.{}.txt;
cat ${READSDIR}/../split_libraries/seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c > ${TMPDIR}/uparse{}/Stats.{}.MergedReads.txt;
cat ${TMPDIR}/uparse{}/stats.otu_table.{}.txt | tail -n +16 | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${TMPDIR}/uparse{}/Stats.{}.MappedReads.txt;
cat ${TMPDIR}/uparse{}/stats.otu_table_missed.{}.txt | tail -n +16 | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${TMPDIR}/uparse{}/Stats.{}.MappedReads.missed.txt;
cat ${TMPDIR}/uparse{}/stats.otu_table_unmapped.{}.txt | tail -n +16 | sed "s/^ //g" | sed -re "s/: /\t/g" | sed "s/\.0$//g" > ${TMPDIR}/uparse{}/Stats.{}.MappedReads.unmapped.txt;
perl ${WONGGITREPO}/ITS_workflows/StatsComparisonMergedVsMapped.pl ${TMPDIR}/uparse{}/Stats.{}.MergedReads.txt ${TMPDIR}/uparse{}/Stats.{}.MappedReads.txt ${TMPDIR}/uparse{}/Stats.{}.MappedReads.missed.txt ${TMPDIR}/uparse{}/Stats.{}.MappedReads.unmapped.txt > ${READSDIR}/../../Deliverables/Read_QC.{}.txt;
cat ${TMPDIR}/uparse{}/missed.centroids.uc | grep "^N" | cut -f9 | ~mcwong/getSeq ${TMPDIR}/uparse{}/missed.fa > ${TMPDIR}/uparse{}/OTUsFailedToMap2ndTime.fa;
merge_otu_tables.py -i ${TMPDIR}/uparse{}/otu_table.{}.biom,${TMPDIR}/uparse{}/otu_table_missed.{}.biom,${TMPDIR}/uparse{}/otu_table_unmapped.{}.biom -o ${TMPDIR}/uparse{}/otu_table_merged.{}.biom;
biom convert -i ${TMPDIR}/uparse{}/otu_table.{}.biom -o ${READSDIR}/../../Deliverables/OTU_Table.{}.txt --header-key taxonomy --to-tsv &
biom convert -i ${TMPDIR}/uparse{}/otu_table_merged.{}.biom -o ${READSDIR}/../../Deliverables/OTU_Table.Merged.{}.txt --header-key taxonomy --to-tsv &
biom convert -i ${TMPDIR}/uparse{}/otu_table_missed.{}.biom -o ${READSDIR}/../../Deliverables/OTU_Table.Missed.{}.txt --header-key taxonomy --to-tsv &
biom convert -i ${TMPDIR}/uparse{}/otu_table_unmapped.{}.biom -o ${READSDIR}/../../Deliverables/OTU_Table.Unmapped.{}.txt --header-key taxonomy --to-tsv &
cp otu_table.{}.biom ${READSDIR}/../../Deliverables/OTU_Table.{}.biom &
cp otu_table_merged.{}.biom ${READSDIR}/../../Deliverables/OTU_Table.Merged.{}.biom &
cp otu_table_missed.{}.biom ${READSDIR}/../../Deliverables/OTU_Table.Missed.{}.biom &
cp otu_table_unmapped.{}.biom ${READSDIR}/../../Deliverables/OTU_Table.Unmapped.{}.biom &
cp ${TMPDIR}/derep.fna ${TMPDIR}/uparse{}/ &
cp ${TMPDIR}/derep.uc ${TMPDIR}/uparse{}/ &
cp ${TMPDIR}/sorted.fa ${TMPDIR}/uparse{}/ &
wait;
tar -cvf ${READSDIR}/../uparse{}.tar.bz2 -C ${TMPDIR} uparse{} -I pbzip2 && cd ${READSDIR} && rm -rf ${TMPDIR}/uparse{};
' 1>> ${READSDIR}/../../Logs/IterativeUparse.log 2>> ${READSDIR}/../../Logs/IterativeUparse.err &
bigJob=`jobs -p`;

## construct the deliverables ##

#zip up Read and Merge Filetar.bz2s
cat ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.1.fq && rm ${TMPDIR}/{}.1.fq && rm ${TMPDIR}/{}.1.fq" > ${TMPDIR}/Raw_Read1.fq &
cat ${READSDIR}/../SampleList | parallel -j1 -I {} "cat ${TMPDIR}/{}.2.fq && rm ${TMPDIR}/{}.2.fq && rm ${TMPDIR}/{}.2.fq" > ${TMPDIR}/Raw_Read3.fq &
smallJobs=`jobs -p | grep -v "${bigJob}"`;
wait ${smallJobs};
pbzip2 -p${THREADS} ${TMPDIR}/Raw_Read1.fq;
pbzip2 -p${THREADS} ${TMPDIR}/Raw_Read3.fq;
pbzip2 -p${THREADS} ${TMPDIR}/Merged_Reads_Raw.fq;
pbzip2 -p${THREADS} ${TMPDIR}/Merged_Reads.fq;

#recover barcodes for deliverables
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/Raw_Read1.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADSPLIT} -c > ${TMPDIR}/Raw_Read2_Barcodes.fq.bz2 &
${WONGGITREPO}/16S_workflows/recoverBarcodesForRaw.pl ${TMPDIR}/Merged_Reads_Raw.fq.bz2 ${READSDIR}/../../${PROJECTID}Barcodes/Project_${PROJECTID}/Sample_${PROJECTID}/${PROJECTID}_NoIndex_L001_R2_001.fastq.gz | pbzip2 -p${THREADSPLIT} -c > ${TMPDIR}/Merged_Barcodes_Raw.fq.bz2 &
smallJobs=`jobs -p | grep -v "${bigJob}"`;
wait ${smallJobs};

#make the deliverables directory and move files into it
(cp ${TMPDIR}/*.fq.bz2 ${READSDIR}/../../Deliverables/ && rm ${TMPDIR}/*.fq.bz2) &
(head -1 ${GITREPO}/Miscellaneous/IlluminaHeaderExample > ${READSDIR}/../../Deliverables/Demultiplex_Sheet.txt && cat ${READSDIR}/../../samplesheet.${PROJECTID}.csv | grep -f ${READSDIR}/../SampleList | cut -f3,5 -d "," | tr "," "\t" | tail -n+1 | sed -re 's/(.*)\t(.*)/\1\t\2\tGGACTACHVGGGTWTCTAAT\tGTGCCAGCMGCCGCGGTAA\t\1/g' >> ${READSDIR}/../../Deliverables/Demultiplex_Sheet.txt) &
cat ${GITREPO}/ITS/Versions.txt > ${READSDIR}/../../Deliverables/${PROJECTID}.SoftwareVersionInformation.txt &
wait;
chmod -R 755 ${READSDIR}/../../Deliverables;
if [ -r "${READSDIR}/../../Deliverables/OTU_Table.99.biom" -a -r "${READSDIR}/../../Deliverables/OTU_Table.97.biom" ];
then collab=`readlink -e ${READSDIR} | cut -f5 -d "/"`;
pool=`readlink -e ${READSDIR} | cut -f6 -d "/"`;
if [ "${collab}" != "StatsProject" ];
then echo -e "${collab} ${pool} has completed running thru the 16S V4 pipeline.  Attached are the read statistics for this run.\nAll other deliverables can be found on the CMMR cluster at the following location:\t`readlink -e ${READSDIR}/../../Deliverables`" | mail -a ${READSDIR}/../../Deliverables/Read_QC.txt -s "${collab} ${pool} has completed" gesell@bcm.edu,dls1@bcm.edu,mcross@bcm.edu,Jacqueline.O\'Brien@bcm.edu,Nadim.Ajami@bcm.edu, carmical@bcm.edu, jcope@diversigen.com;
elif [ "${collab}" = "StatsProject" ];
then echo -e "${collab} ${pool} has completed running thru the 16S V4 pipeline.  Attached are the read statistics for this run.\nAll other deliverables can be found on the CMMR cluster at the following location:\t`readlink -e ${READSDIR}/../../Deliverables`" | mail -a ${READSDIR}/../../Deliverables/Read_QC.txt -s "${collab} ${pool} has completed" gesell@bcm.edu, carmical@bcm.edu;
fi;
else echo -e "${collab} ${pool} run failed, please check reason" | mail -a ${READSDIR}/../../Deliverables/Read_QC.txt -s "${collab} ${pool} has failed" ${USER}@bcm.edu;
fi;

#return to working directory when script was launched
cd ${CURRWORKDIR};

#exit without error status once completed
exit 0;
