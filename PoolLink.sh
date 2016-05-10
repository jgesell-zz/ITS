#!/bin/sh

#This program takes the first name, last name and pool number for a MiSeq run and creates the file structure to run our current pipeline.  It assumes that you have already demultiplexed the whole pool, but not the individual collaborator.

firstName=$1;
lastName=$2;
pool=$3;
prefix=$4

#make the necessary directory structures
mkdir ${firstName}${lastName};
mkdir ${firstName}${lastName}/Pool${pool};
cd ${firstName}${lastName}/Pool${pool};
mkdir ${lastName}Pool${pool}WorkDir;
mkdir ${lastName}Pool${pool}WorkDir/Reads;
mkdir ${lastName}Pool${pool}Reads;
mkdir ${lastName}Pool${pool}Reads/Project_${lastName}Pool${pool}
mkdir ${lastName}Pool${pool}Barcodes;
mkdir ${lastName}Pool${pool}Barcodes/Project_${lastName}Pool${pool};
mkdir ${lastName}Pool${pool}Barcodes/Project_${lastName}Pool${pool}/Sample_${lastName}Pool${pool};
mkdir Logs;

if [ -z "$prefix" ];
	then prefix=`echo $lastName`;
fi;

#cat in the information for files that cannot be softlinked
cat ../../StatsProject/ITS/Pool${pool}/Pool${pool}WorkDir/SampleList | grep "${prefix}" > ${lastName}Pool${pool}WorkDir/SampleList;
cat ../../StatsProject/ITS/Pool${pool}/samplesheet.*${pool}.csv | grep "${prefix}" > samplesheet.${lastName}Pool${pool}.csv

#softlink the required demultiplexed reads into the Reads/Project_* directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/ | grep -f ${lastName}Pool${pool}WorkDir/SampleList`; do ln -s ../../../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/$i ${lastName}Pool${pool}Reads/Project_${lastName}Pool${pool}/$i; done;

#softlink the other files in the master pool Reads directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/`; do name=`echo $i | sed "s:Pool${pool}:${lastName}Pool${pool}:g" | sed "s:Overall::g" | sed "s:ReagentTest::g" `; ln -s ../../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/$i ${lastName}Pool${pool}Reads/$name; done;

#softlink the items in the individual reads into the WorkDir/Reads directory
for i in `find ${lastName}Pool${pool}Reads/Project_${lastName}Pool${pool}/Sample_*/*.bz2`; do name=`echo $i | cut -f4 -d "/" | cut -f1 -d "_"`; num=`echo $i | cut -f6 -d "_" | cut -c2`; ln -s ../../${i} ${lastName}Pool${pool}WorkDir/Reads/${name}.${num}.fq.bz2; done

#softlink the un-demultiplexed reads
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../../../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/$i ${lastName}Pool${pool}Barcodes/Project_${lastName}Pool${pool}/Sample_${lastName}Pool${pool}/$name; done;

#softlink the files in the barcodes directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/ | grep -v "Project_Pool${pool}"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" | sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/$i ${lastName}Pool${pool}Barcodes/$name; done;

#softlink the log files into the Logs directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Logs/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../../StatsProject/ITS/Pool${pool}/Logs/$i Logs/$name; done;

#softlink any remaining files into the base directory, substituting the PoolID where appropriate
for i in `ls ../../StatsProject/ITS/Pool${pool}/ | grep -v "Logs" | grep -v " Pool${pool}Barcodes" | grep -v "Pool${pool}Reads" | grep -v "Pool${pool}WorkDir" | grep -v "Deliverables" | grep -v "samplesheet.${pool}.csv"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}Pool${pool}:g"`; ln -s ../../StatsProject/ITS/Pool${pool}/$i $name; done;

#echo the number of samples found
numSamples=`cat ${lastName}Pool${pool}WorkDir/SampleList | wc -l`;
echo "Total samples found: ${numSamples}"; 

#automatically launch the processing job
link=`readlink -e ${lastName}Pool${pool}WorkDir/Reads/;`
echo "${GITREPO}/ITS/fullPipelineSplit.sh $link 40" | qsub -l ncpus=20 -q batch -N ${lastName}Pool${pool}.Process -d `pwd -P` -V;
