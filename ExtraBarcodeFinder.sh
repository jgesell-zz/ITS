#!/bin/sh

if [ -z ${Pool} ];
then Pool=`basename \`pwd -P\``;
fi;

cat Barcodes.txt | cut -f2 > ExpectedBarcodes;
~mcwong/reversecomp.pl ExpectedBarcodes > ExpectedBarcodesReversed;
zcat ${Pool}Barcodes/Project_${Pool}/Sample_${Pool}/${Pool}_NoIndex_L001_R2_001.fastq.gz | grep -A1 "^@M70287" | grep -v "^@M70287" | grep -v "\-\-" | cut -c1-12 | sort | uniq -c | sort -n | perl ~gesell/Programs/checknonbarcodes.pl ~gesell/Programs/revGolay > ${Pool}NonGolay.txt &
zcat ${Pool}Barcodes/Project_${Pool}/Sample_${Pool}/${Pool}_NoIndex_L001_R2_001.fastq.gz | grep -A1 "^@M70287" | grep -v "^@M70287" | grep -v "\-\-" | cut -c1-12 | sort | uniq -c | sort -n | perl ~gesell/Programs/checkbarcodes.pl ~gesell/Programs/revGolay > ${Pool}Golay.txt &
wait;
cat ${Pool}NonGolay.txt | grep -v -f ExpectedBarcodesReversed | grep -P "\d\d\d\d"> NonGolayFoundBarcodes.txt && rm  ${Pool}NonGolay.txt &
cat ${Pool}Golay.txt | grep -v -f  ExpectedBarcodesReversed | grep -P "\d\d\d\d"> AdditionalGolayFoundBarcodes.txt && rm ${Pool}Golay.txt &
wait;
rm ExpectedBarcodes ExpectedBarcodesReversed;
echo -e "Non Golay Barcodes found ( > 1000 Reads): `cat NonGolayFoundBarcodes.txt | wc -l`\nUnexpected Barcodes Found ( > 1000 Reads): `cat AdditionalGolayFoundBarcodes.txt | wc -l`" > HiddenBarcodeStats.txt;
