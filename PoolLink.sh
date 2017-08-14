#!/bin/sh

#This program takes the first name, last name and pool number for a MiSeq run and creates the file structure to run our current pipeline.  It assumes that you have already demultiplexed the whole pool, but not the individual collaborator.

collabName=$1;
pool=$2;
prefix=$3;
poolName=$4;
sampleList=$5;

if [ ! -e "${sampleList}" ];
then sampleList="";
else sampleList=`readlink -e ${sampleList}`;
fi;

if [ -z "$poolName" ];
then poolName=Pool`echo $pool`;
fi;

firstName=`echo ${collabName} | cut -f1 -d "_"`;
lastName=`echo ${collabName} | cut -f2 -d "_"`;

if [ "${firstName}" = "${lastName}" ];
then firstName="";
fi;

cd /gpfs1/projects/jgesell;

if [ -d "${firstName}${lastName}/${poolName}" ];
then poolName=`echo ${poolName}.Pool${Pool}`;
fi;

#make the necessary directory structures
mkdir -p ${firstName}${lastName}/${poolName};
cd ${firstName}${lastName}/${poolName};
mkdir -p ${lastName}${poolName}WorkDir/Reads;
mkdir -p ${lastName}${poolName}Reads/Project_${lastName}${poolName}
mkdir -p ${lastName}${poolName}Barcodes/Project_${lastName}${poolName}/Sample_${lastName}${poolName};
mkdir Logs;

if [ -z "$prefix" ];
	then prefix=`echo $lastName`;
fi;

#cat in the information for files that cannot be softlinked
if [ -e "${sampleList}" ];
then cat ../../StatsProject/ITS/Pool${pool}/Pool${pool}WorkDir/SampleList | grep -f "${sampleList}" > ${lastName}${poolName}WorkDir/SampleList;
cat ../../StatsProject/ITS/Pool${pool}/samplesheet.*${pool}.csv | grep -f "${sampleList}" > samplesheet.${lastName}${poolName}.csv
cat ../../StatsProject/ITS/Pool${pool}/Pool${pool}.barcodeCounts.txt | grep -f "${sampleList}" >  ${lastName}${poolName}.barcodeCounts.txt
else cat ../../StatsProject/ITS/Pool${pool}/Pool${pool}WorkDir/SampleList | grep "${prefix}" > ${lastName}${poolName}WorkDir/SampleList;
cat ../../StatsProject/ITS/Pool${pool}/samplesheet.*${pool}.csv | grep "${prefix}" > samplesheet.${lastName}${poolName}.csv
cat ../../StatsProject/ITS/Pool${pool}/Pool${pool}.barcodeCounts.txt | grep "${prefix}" >  ${lastName}${poolName}.barcodeCounts.txt
fi;

#softlink the required demultiplexed reads into the Reads/Project_* directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/ | grep -f ${lastName}${poolName}WorkDir/SampleList`; do ln -s ../../../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/Project_Pool${pool}/$i ${lastName}${poolName}Reads/Project_${lastName}${poolName}/$i; done;

#softlink the other files in the master pool Reads directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/`; do name=`echo $i | sed "s:Pool${pool}:${lastName}${poolName}:g" | sed "s:Overall::g" | sed "s:ReagentTest::g" `; ln -s ../../../StatsProject/ITS/Pool${pool}/Pool${pool}Reads/$i ${lastName}${poolName}Reads/$name; done;

#softlink the items in the individual reads into the WorkDir/Reads directory
for i in `find ${lastName}${poolName}Reads/Project_${lastName}${poolName}/Sample_*/*.bz2`; do name=`echo $i | cut -f4 -d "/" | cut -f1 -d "_"`; num=`echo $i | cut -f6 -d "_" | cut -c2`; ln -s ../../${i} ${lastName}${poolName}WorkDir/Reads/${name}.${num}.fq.bz2; done

#softlink the un-demultiplexed reads
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -s ../../../../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/Project_Pool${pool}/Sample_Pool${pool}/$i ${lastName}${poolName}Barcodes/Project_${lastName}${poolName}/Sample_${lastName}${poolName}/$name; done;

#softlink the files in the barcodes directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/ | grep -v "Project_Pool${pool}"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" | sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -s ../../../StatsProject/ITS/Pool${pool}/Pool${pool}Barcodes/$i ${lastName}${poolName}Barcodes/$name; done;

#softlink the log files into the Logs directory
for i in `ls ../../StatsProject/ITS/Pool${pool}/Logs/`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -s ../../../StatsProject/ITS/Pool${pool}/Logs/$i Logs/$name; done;

#softlink any remaining files into the base directory, substituting the PoolID where appropriate
for i in `ls ../../StatsProject/ITS/Pool${pool}/ | grep -v "Logs" | grep -v " Pool${pool}Barcodes" | grep -v "Pool${pool}Reads" | grep -v "Pool${pool}WorkDir" | grep -v "Deliverables" | grep -v "samplesheet.${pool}.csv" | grep -v "barcodeCounts.txt"`; do name=`echo $i | sed "s:Overall::g" | sed "s:ReagentTest::g" |  sed "s:Pool${pool}:${lastName}${poolName}:g"`; ln -s ../../StatsProject/ITS/Pool${pool}/$i $name; done;

#echo the number of samples found
numSamples=`cat ${lastName}${poolName}WorkDir/SampleList | wc -l`;
echo "Total samples found: ${numSamples}"; 

#Allocate the threads needed for the pipeline
if [ ${numSamples} -lt 40 ];
then PROCS=$[ $[numSamples / 2] + $[numSamples % 2]];
THREADS=${numSamples};
else PROCS=20;
THREADS=40;
fi;

#automatically launch the processing job
link=`readlink -e ${lastName}${poolName}WorkDir/Reads/;`
echo "${GITREPO}/ITS/fullPipelineSplit.sh ${link} ${THREADS} " | qsub -l ncpus=${PROCS} -q batch -N ${lastName}${poolName}.Process -d `pwd -P` -V;
