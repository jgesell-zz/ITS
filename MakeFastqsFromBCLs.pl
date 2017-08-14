#!/usr/bin/env perl




#CMMR Pipeline for converting BCLs to a .biom file
#Author: Matt Wong
#Date Updated: 2-12-14
#Version: 0.1
#
#
#
#Comments need throughout
#Take in a BCL directory path and an HGSC 16S samplesheet, and generates a .biom file with the samples in it. 
#



use warnings;
use strict;
use DataBrowser qw(browse);
use Spreadsheet::XLSX;
use Getopt::Long;
use File::Basename;
use Spreadsheet::ParseExcel;
use Cwd qw(abs_path);
use IPC::System::Simple qw(capture system);
use threads;
use XML::Hash;

my $bclDir;
my $sampleSheet;
my $help;
my $barcodes;
my $projectName;
my $threads;
my $forceAscii;
my $noRC;
my $only;
my $barcodeMaskCL;

GetOptions ("sampleSheet=s" => \$sampleSheet,
	    "bclDir=s"      => \$bclDir,
	    "barcodes=s"    => \$barcodes,
	    "help"          => \$help,
	    "projectName=s" => \$projectName,
	    "threads=i"     => \$threads,
	    "forceAscii"    => \$forceAscii,
	    "noRC"          => \$noRC,
	    "onlyBarcodes"  => \$only,
	    "baseMask=s" => \$barcodeMaskCL);
	    
sub usage {
	my $usage = <<END;
usage: 16SIlluminaPipeline.pl --BclDir /path/to/Data/Intensities/BaseCalls --sampleSheet samplesheet.xls --projectName MyProjectName --barcodes /path/to/barcodes.fastq --threads <n_threads>
END
	
	print STDERR "$usage";
	exit()
}

sub findBclDir {
	
	#Find the BCL (BaseCalls) directory given the input from the command line
	#Automatically search down if you gave the top level BCL path
	#Accepts windows formatted paths

	my $nothing = 1 + 1;	
	my $dir = $_[0];

	$dir =~ s/\\/\//g;
	$dir =~ s/.*://g;
	$dir =~ s/\/$//g;
	if (-f $dir) {
		$dir = readlink($dir);
	}
	unless (-d $dir) {
		die "$dir does not exist!\nDid you put it in quotes?\n";
	}
	$dir = abs_path($dir);
	my $BclDir = "$dir/Data/Intensities/BaseCalls";
	if (-d $BclDir) {
		return($BclDir);
	} elsif ((-d $dir) and $dir =~ m/Data\/Intensities\/BaseCalls/) {
		return($dir);
	} else {
		die "directory is not basecall directory";
	}
}

sub cleanSampleNames {
	
	#Cleans samples names that might have characters that interfere with filenames.

	my @old = @{$_[0]};
	my @new;
	foreach my $name (@old) {
		$name =~ s/\s+$//g;
		$name =~ s/\.$//g;
		$name =~ s/^\s+//g;
		$name =~ s/\s+\t/\t/g;
		$name =~ s/\'/x/g;
		$name =~ s/\-/./g;
		$name =~ s/[^\w\-\+\=\ \.]//g;
		$name =~ s/\s+/./g;
		$name =~ s/\.\././g;
		$name =~ s/_/./g;
		push @new, $name;
	}
	return(\@new);
}

sub reverseComp {

	#Quick reverse complement
	#This should really be in a separate package

	my @sequences = @{$_[0]};
	my $noRC = $_[1];
	if ($noRC) {
		return(\@sequences);
	}
	my @reverse;
	foreach my $seq (@sequences) {
		$seq =~ tr/ACGT/TGCA/;
		$seq = reverse($seq);
		push @reverse, $seq;
	}
	return (\@reverse);
}

sub parseSampleSheetXls {
	
	#This parses an HGSC sample sheet in xls format.	
	#Grabs all sample and barcodes that occur after line 10

	my %samplesToBarcodes;
	my %barcodesToSamples;
	my $sampleIDPos;
	my $barcodePos;
	my $sampleSheet = $_[0];
	my $noRC = $_[1];
	my $parser = Spreadsheet::ParseExcel->new();
	my $workbook = $parser->parse($sampleSheet);
	if ( !defined $workbook ) {
		die $parser->error(), ".\n";
	}
	my $worksheet = $workbook->worksheet(0);
	my ($row_min, $row_max) = $worksheet->row_range();
	my ($col_min, $col_max) = $worksheet->col_range();
	for my $col ( $col_min .. $col_max ) {
		my $cell = $worksheet->get_cell(3, $col);
		my $val;
		if ($cell) {
			$val = $cell->value();
		} else {
			next;
		}
		if ($val =~ m/Sample.*External.*ID/) {
			$sampleIDPos = $col;
		}
		if ($val =~ m/Molecular.*Barcode.*ID/g) {
			$barcodePos = $col;
		}
	}
	my @samples;
	my @barcodes;
	foreach my $row (9 .. $row_max) {
		my $sampleCell = $worksheet->get_cell($row, $sampleIDPos);
		my $barcodeCell = $worksheet->get_cell($row, $barcodePos);
		unless ($sampleCell && $barcodeCell) {
			next
		}

		my $sampleId = $worksheet->get_cell($row, $sampleIDPos)->value();
		my $barcode = $worksheet->get_cell($row, $barcodePos)->value();
		#print "$sampleId = $barcode\n"; 
		push @samples, $sampleId;
		push @barcodes, $barcode;
	}
	my @sampleNames = @{cleanSampleNames(\@samples)};
	my @revBarcodes = @{reverseComp(\@barcodes, $noRC)};
	if (scalar @sampleNames != scalar @revBarcodes) {
		die "sample names and barcodes don't match\n";
	}
	for (my $i = 0; $i < scalar(@sampleNames); $i++) {
		unless($sampleNames[$i]) {next}
		if (exists($samplesToBarcodes{$sampleNames[$i]})) {
			die "sample $sampleNames[$i] exists more than once\n";
		}
		if (exists($barcodesToSamples{$revBarcodes[$i]})) {
			die "barcode $revBarcodes[$i] exists more tha once\n";
		}
		$samplesToBarcodes{$sampleNames[$i]} = $revBarcodes[$i];
		$barcodesToSamples{$revBarcodes[$i]} = $sampleNames[$i];
	}
	return (\%samplesToBarcodes);
}

sub parseSampleSheetXlsx {

	#This parses an HGSC sample sheet in xlsx format.
	#Grab all samples and barcodes that start after line 10.

	my %samplesToBarcodes;
	my %barcodesToSamples;
	my $sampleIDPos;
	my $barcodePos;
	my $sampleSheet = $_[0];
	my $noRC = $_[1];
	my $excel = Spreadsheet::XLSX -> new ($sampleSheet);
	my $sheet = ${$excel -> {Worksheet}}[0];
	foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
		my $cell = $sheet -> {Cells} [3] [$col];
		if ($cell) {
			my $val = $cell -> {Val};
			if ($val =~ m/Sample.*External.*ID/) {
				$sampleIDPos = $col;
			}
			if ($val =~ m/Molecular.*Barcode.*ID/g) {
				$barcodePos = $col;
			}
		}	
	}
	my @samples;
	my @barcodes;
	foreach my $row (9 .. $sheet -> {MaxRow}) {
		my $sampleId = $sheet -> {Cells} [$row] [$sampleIDPos] -> {Val};
		my $barcode = $sheet -> {Cells} [$row] [$barcodePos] -> {Val};
		push @samples, $sampleId;
		push @barcodes, $barcode;

	}
	my @sampleNames = @{cleanSampleNames(\@samples)};
	my @revBarcodes = @{reverseComp(\@barcodes, $noRC)};
	if (scalar @sampleNames != scalar @revBarcodes) {
		die "sample names and barcodes don't match\n";
	}
	for (my $i = 0; $i < scalar(@sampleNames); $i++) {
		unless($sampleNames[$i]) {next}
			if (exists($samplesToBarcodes{$sampleNames[$i]})) {
			die "sample $sampleNames[$i] exists more than once\n";
		}
		if (exists($barcodesToSamples{$revBarcodes[$i]})) {
			die "barcode $revBarcodes[$i] exists more tha once\n";
		}
		$samplesToBarcodes{$sampleNames[$i]} = $revBarcodes[$i];
		$barcodesToSamples{$revBarcodes[$i]} = $sampleNames[$i];
	}
	return (\%samplesToBarcodes);
}

sub parseSampleSheetText {
	
	#This reads a file with <Sample>\t<barcode> as each line
	#Reverse complements the barcodes.

	my %samplesToBarcodes;
	my %barcodesToSamples;
	my @samples;
	my @barcodes;
	my $sampleSheet = $_[0];
	my $noRC = $_[1];
	open IN, "$sampleSheet";
	while (my $line = <IN>) {
		chomp $line;
		my @parts = split /\t/, $line;
		unless (scalar(@parts) == 2) {
			@parts = split /,/, $line;
		}
		unless (scalar(@parts) == 2) {
			die "Text samplesheet not delimited by tab or comma\n";
		}
		if ($parts[1] =~ m/[^ACGT]/) {
			die "barcode $parts[1] contains non-ACGT characters. Line: $line\n";
		}
		push @samples, $parts[0];
		push @barcodes, $parts[1];
	}
	my @sampleNames = @{cleanSampleNames(\@samples)};
	my @revBarcodes = @{reverseComp(\@barcodes, $noRC)};
	if (scalar @sampleNames != scalar @revBarcodes) {
		die "sample names and barcodes don't match\n";
	}
	for (my $i = 0; $i < scalar(@sampleNames); $i++) {
		unless($sampleNames[$i]) {next}
			if (exists($samplesToBarcodes{$sampleNames[$i]})) {
			die "sample $sampleNames[$i] exists more than once\n";
		}
		if (exists($barcodesToSamples{$revBarcodes[$i]})) {
			die "barcode $revBarcodes[$i] exists more tha once\n";
		}
		$samplesToBarcodes{$sampleNames[$i]} = $revBarcodes[$i];
		$barcodesToSamples{$revBarcodes[$i]} = $sampleNames[$i];
	}
	return (\%samplesToBarcodes);
}

sub getFlowcellId {

	#Grabs the flowcell ID form the BCL directory for easy samplesheet generation

	my $dir = $_[0];
	my $configFile = "$dir/config.xml";
	my $xml_converter = XML::Hash->new();
	
	unless (-f $configFile) {
		die "FATAL: Can't find config.xml in Basecall directory: $dir\n";
	}
	open IN, "$configFile";
	my $xmlLines = "";
	while (my $line = <IN>) {
		chomp $line;
		$xmlLines .= $line;
	}
	my $xmlHash = $xml_converter->fromXMLStringtoHash($xmlLines);
	my $flowcellId = $xmlHash->{'BaseCallAnalysis'}->{'Run'}->{'RunParameters'}->{'RunFlowcellId'}->{'text'};
	my $read1Length = ($xmlHash->{'BaseCallAnalysis'}->{'Run'}->{'RunParameters'}->{'Reads'}->[0]->{'LastCycle'}->{'text'} + 1) - ($xmlHash->{'BaseCallAnalysis'}->{'Run'}->{'RunParameters'}->{'Reads'}->[0]->{'FirstCycle'}->{'text'});
	my $read2Length = ($xmlHash->{'BaseCallAnalysis'}->{'Run'}->{'RunParameters'}->{'Reads'}->[-1]->{'LastCycle'}->{'text'} + 1) - ($xmlHash->{'BaseCallAnalysis'}->{'Run'}->{'RunParameters'}->{'Reads'}->[-1]->{'FirstCycle'}->{'text'});
	if ($read1Length != $read2Length) {
		warn ("Read1 and Read2 are not the right length; Read1 $read1Length, Read2 $read2Length\n");
	}
	my $barcodeMask = $read1Length;
	if ($read2Length < $read1Length) {
		$barcodeMask = $read2Length;
	}
	my $laneId = $xmlHash->{'BaseCallAnalysis'}->{'Run'}->{'TileSelection'}->{'Lane'}->{'Index'};
	return ($flowcellId, $laneId, $barcodeMask);
}

sub callBcl {

	#This calls the BCL to fastq function in CASAVA and generates separate fastqs per sample
	

	my $projectName = $_[0];
	my $bclDir = $_[1];
	my $samplesheet = $_[2];
	my $threads = $_[3];
	my $barcodeMask = $_[4];
	my $projDir = "$projectName" . "Reads";
	

	#Should this be enabled? If the project directory exists, configureBclToFastq.pl will throw an error.
	#Perhaps an automatic way to move the old directory out of the way?
	
	#if (-d $projDir) {
	#	die "$projDir already exists, this will fail\n";
	#}
	my $cmd = "configureBclToFastq.pl --input-dir $bclDir --output-dir $projDir --use-bases-mask $barcodeMask,i12n,$barcodeMask --sample-sheet $samplesheet --fastq-cluster-count 0 --mismatches 1 2>Logs/$projectName.configureBcl.err";
	system($cmd);
	$projDir = `readlink -e $projDir`;
	chomp $projDir;
	#print STDERR "projDir = $projDir\n";
	#die;
	$cmd = "make -j $threads -C $projDir 1>Logs/make.$projectName.out 2>Logs/make.$projectName.err";
	system($cmd);
	$cmd = "rm -rf $projDir/Undetermined_indices $projDir/Temp";
	system($cmd);
	return($projDir);
}

sub callBclRaw {
	my $projectName = $_[0];
	my $bclDir = $_[1];
	my $samplesheet = $_[2];
	my $laneId = $_[3];
	my $threads = $_[4];
	my $barcodeMask = $_[5];
	my $projDir = "$projectName" . "Barcodes";
	
	#Same issue as in callBcl
	#if (-d $projDir) {
	#	die "$projDir already exists, this will fail\n";
	#}
	my $cmd = "configureBclToFastq.pl --input-dir $bclDir --output-dir $projDir --use-bases-mask $barcodeMask,y12n,$barcodeMask --sample-sheet $samplesheet --fastq-cluster-count 0 2>Logs/$projectName.configureBcl.barcodes.err";
	capture($cmd);
	$cmd = "make -j $threads -C $projDir 1>Logs/make.$projectName.barcodes.out 2>Logs/make.$projectName.barcodes.err";
	capture($cmd);
	return("$projDir/Project_$projectName/Sample_$projectName/${projectName}_NoIndex_L00${laneId}_R1_001.fastq.gz","$projDir/Project_$projectName/Sample_$projectName/${projectName}_NoIndex_L00${laneId}_R2_001.fastq.gz","$projDir/Project_$projectName/Sample_$projectName/${projectName}_NoIndex_L00${laneId}_R3_001.fastq.gz");
}

sub makeBclSampleSheet {
	
	#Makes the samplesheet for individual sample fastqs

	my %samplesToBarcodes = %{$_[0]};
	my $flowcellId = $_[1];
	my $laneId = $_[2];
	my $projectName = $_[3];
	my $header = "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject";
	my @sampleStrs;
	push @sampleStrs, $header;
	foreach my $sample (keys %samplesToBarcodes) {
		my $sampleStr = "$flowcellId,$laneId,$sample,,$samplesToBarcodes{$sample},,N,$laneId,MW,$projectName";
		push @sampleStrs, $sampleStr;
	}
	my $printStr = join "\n", @sampleStrs;
	open OUT, ">samplesheet.$projectName.csv";
	print OUT "$printStr\n";
	close OUT;
	return("samplesheet.$projectName.csv")
}

sub concatenateReads {

	#Unzips reads into a working directory
	#bzip2's the original reads for better compression
	
	my $projDir = $_[0];
	my $projName = $_[1];
	my $lane = $_[2];
	my $threads = $_[3];
	my %samplesToBarcodes = %{$_[4]};
	my %links;
	my $workDir = "$projName" . "WorkDir";
	unless (-d $workDir) {
		system("mkdir $workDir")
	}
	my $readDir = "$workDir/Reads";
	unless (-d $readDir) {
		system("mkdir $readDir");
	}
	open COUNT, ">$projName.barcodeCounts.txt" or die "can't open file to write barcode counts\n";
	my @compressFiles;
	foreach my $sample (keys %samplesToBarcodes) {
		system("echo $sample >> $workDir/SampleList");
		my $cmd = "zcat $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq.gz | tee $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq | wc -l";		
#my $cmd = "zcat $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq.gz | tee $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq | sed 's/1:N:0:.*/1:N:0:/g' | grep \"^\@HWI\" | wc -l";
		my $cmd2 = "zcat $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R2_001.fastq.gz > $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R2_001.fastq";
		my $output = capture($cmd);
		$output = $output / 4;
		system($cmd2);
		$cmd = "rm $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq.gz $projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R2_001.fastq.gz";
		system($cmd);
		push @compressFiles, "$projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq";
		push @compressFiles, "$projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R2_001.fastq";
		print COUNT "$sample\t$output\n";
		$links{"${readDir}/${sample}.1.fq.bz2"} = "$projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R1_001.fastq.bz2";
		$links{"${readDir}/${sample}.2.fq.bz2"} = "$projDir/Project_$projName/Sample_$sample/${sample}_$samplesToBarcodes{$sample}_L00${lane}_R2_001.fastq.bz2";
	}
	my $compressStr = join " ", @compressFiles;
	my $cmd = "pbzip2 -p$threads $compressStr";
	system($cmd);

	foreach my $link (keys %links){
		my $target = $links{$link};
		$cmd = "ln -s $target $link";
		system($cmd)
	}
	return($workDir);

}

sub getPhiXPercRaw {
	
	#Calculate PhiX percentage of the raw machine output. 

	my $read1 = $_[0];
	my $read2 = $_[1];
	my $projectName = $_[2];
	my $threads = $_[3];
	
	#database shouldn't be hardcoded.
	my $cmd = "bowtie2 -x \${PHIXDB} -U $read1,$read2 --end-to-end --very-sensitive --reorder -p $threads -S /dev/null 2>Logs/$projectName.phixOverall.stats.txt";
	system($cmd);
}

sub getPhiXPercBleed {

	#Calculate PhiX percentage that bled over into the demultiplexed samples.	

	my $workDir = $_[0];
	my $threads = $_[1];
	my $projectName = $_[2];
	#database shouldn't be hardcoded.
	my $cmd = "bowtie2 -x \${PHIXDB} -1 $workDir/Read1.fq -2 $workDir/Read2.fq --end-to-end --very-sensitive --reorder -p $threads -S /dev/null --un-conc $workDir/Reads.filtered.fq 2>Logs/$projectName.phix.bleed.stats";
	system($cmd);
	$cmd = "rm $workDir/Read1.fq $workDir/Read2.fq";
	system($cmd);
	$cmd = "pbzip2 -p$threads $workDir/Reads.filtered*";
	system($cmd);

}

sub makeBarcodes {

	#Runs bcl2fastq in non-demultiplexing mode
	#Basically, every read that passes filter is dumped into a monolithic fastq


	my $bclDir = $_[0];
	my $flowcellId = $_[1];
	my $laneId = $_[2];
	my $projName = $_[3];
	my $threads = $_[4];
	my $barcodeMask = $_[5];
	my $samplesheet = "sampleSheet.notDemultiplexed.$projName.csv";
	open OUT, ">$samplesheet";
	print OUT "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n";
	print OUT "$flowcellId,$laneId,$projectName,not_demultiplexed,,$projectName,N,$flowcellId,MCW,$projectName\n";
	close OUT;
	return(callBclRaw($projName, $bclDir, $samplesheet, $laneId, $threads, $barcodeMask));
}

if ($help) {
	usage()
}
unless ($bclDir and $sampleSheet and $projectName) {
	usage()
}

unless ($threads) {
	$threads = 1;
}
unless (-d "Logs") {
	system("mkdir Logs");
}
sub barcodeRun {
	
	#mini workflow to generate files needed to reconstitute barcodes files.

 	my $bclDir = $_[0];
	my $flowcellId = $_[1];
	my $laneId = $_[2];
	my $projectName = $_[3];
	my $threads = $_[4];
	my $barcodeMask = $_[5];
	my $barcodes;
	my $read1Raw;
	my $read2Raw;	
	($read1Raw, $barcodes, $read2Raw) = makeBarcodes($bclDir, $flowcellId, $laneId, $projectName, $threads, $barcodeMask);
	getPhiXPercRaw($read1Raw, $read2Raw, $projectName, $threads);
	return 0;
}


######## MAIN ########


my $filetype = `file "$sampleSheet"`;
my %samplesToBarcodes;
if ($filetype =~ m/cannot open/) {
	die "$sampleSheet doesn't exist.\n";
} elsif ($filetype =~ m/CDF V2 Document/) {
	%samplesToBarcodes = %{parseSampleSheetXls($sampleSheet, $noRC)};
} elsif ($filetype =~ m/Zip archive data/) {
	%samplesToBarcodes = %{parseSampleSheetXlsx($sampleSheet, $noRC)};
} elsif ($filetype =~ m/ASCII text/ || $forceAscii) {
	%samplesToBarcodes = %{parseSampleSheetText($sampleSheet, $noRC)};
} else {
	die "$sampleSheet not in recognizable format, must be xls or xlsx.\n";
}

$bclDir = findBclDir($bclDir);

my ($flowcellId, $laneId, $barcodeMask) = getFlowcellId($bclDir);

if ($barcodeMaskCL and $barcodeMaskCL < $barcodeMask) {
	$barcodeMask = "y" . $barcodeMaskCL . "n*";
} elsif ($barcodeMaskCL and $barcodeMaskCL > $barcodeMask) {
	die "--barcodeMask cannot be greater than the amount of cycles run. Detected Cycles = $barcodeMask\n";
} else {
	$barcodeMask = "y" . $barcodeMask;
}



my $read1Raw;
my $read2Raw;
my $barcodeThread;
unless ($barcodes) {
	if ($threads > 1) {
		$barcodeThread = threads->create(\&barcodeRun, $bclDir, $flowcellId, $laneId, $projectName, $threads, $barcodeMask);
	} else {
		#($read1Raw, $barcodes, $read2Raw) = makeBarcodes($bclDir, $flowcellId, $laneId, $projectName, $threads, $barcodeMask);
		#getPhiXPercRaw($read1Raw, $read2Raw, $projectName, $threads);
		barcodeRun($bclDir, $flowcellId, $laneId, $projectName, $threads, $barcodeMask);
	}
	if ($only) {
		exit(0);
	}
}

my $samplesheet = makeBclSampleSheet(\%samplesToBarcodes, $flowcellId, $laneId, $projectName);
my $projDir = callBcl($projectName, $bclDir, $samplesheet, $threads, $barcodeMask);
my $workDir = concatenateReads($projDir, $projectName, $laneId, $threads, \%samplesToBarcodes);
unless($barcodes) {
	if ($threads > 1) {
		$barcodeThread->join();
	}
}
