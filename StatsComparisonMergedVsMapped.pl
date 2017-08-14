#!/usr/bin/env perl
#
use warnings;
use strict;


my $raw = $ARGV[0];
my $merged = $ARGV[1];
my $picked = $ARGV[2];
my $low = $ARGV[3];
my $missed = $ARGV[4];

my $reads = {};

#Get Raw Reads
open IN, "$raw";
open STATSRRAW, ">tempRaw";
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/^\s+//g;
	my @parts = split /\s+/, $line;
	$reads->{$parts[0]}->{'raw'} = $parts[1];
	print STATSRAW "$parts[1]\n";
}
close STATSMERGED;
close IN;

#Get Merged Reads
open IN, "$merged";
open STATSMERGED, ">tempMerged";
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/^\s+//g;
	my @parts = split /\s+/, $line;
	$reads->{$parts[1]}->{'merged'} = $parts[0];
	print STATSMERGED "$parts[0]\n";
}
close STATSMERGED;
close IN;

#Get Mapped Reads
open IN, "$picked";
open STATSPICK, ">temppick";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
#	print "sample = $parts[0]\n";
	$reads->{$parts[0]}->{'picked'} = $parts[1];
	print STATSPICK "$parts[1]\n";
}
close STATSPICK;
close IN;

#Get Low-Identity Mapped Reads
open IN, "$low";
open STATSLOW, ">templow";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	$reads->{$parts[0]}->{'low'} = $parts[1];
	print STATSLOW "$parts[1]\n";
}
close STATSLOW;
close IN;

#Get Reads Missed in Initial Run Thru UParse
open IN, "$missed";
open STATSMISSED, ">tempmissed";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	$reads->{$parts[0]}->{'missed'} = $parts[1];
	print STATSMISSED "$parts[1]\n";
}
close IN;
close STATSMISSED;

#Correct for any reads that do not show up in one of the categories
foreach my $sample (keys %{$reads}) {
	unless (exists $reads->{$sample}->{'raw'}) {
		$reads->{$sample}->{'raw'} = 0
	}
	unless (exists $reads->{$sample}->{'picked'}) {
		$reads->{$sample}->{'picked'} = 0
	}
	unless (exists $reads->{$sample}->{'merged'}) {
		$reads->{$sample}->{'merged'} = 0;
	}
	unless(exists $reads->{$sample}->{'low'}) {
		$reads->{$sample}->{'low'} = 0;
	}
	unless(exists $reads->{$sample}->{'missed'}) {
		$reads->{$sample}->{'missed'} = 0;
	}

}

#Print out the stats on a per-sample basis
print "SampleID\tRaw Reads\tMapped\tLow Identity Mapped\tUnmapped\tRemainder=Chimeras+Singletons\n";
foreach my $sample (sort {$reads->{$a}->{'raw'} <=> $reads->{$b}->{'raw'}} keys %{$reads}) {
	my $raw = $reads->{$sample}->{'raw'};
	my $picked = $reads->{$sample}->{'picked'};
	my $low = $reads->{$sample}->{'low'};
	my $unmapped = $reads->{$sample}->{'missed'};
	my $diff = $reads->{$sample}->{'merged'} - $picked - $low - $unmapped;
	print "$sample\t$raw\t$picked\t$low\t$unmapped\t$diff\n";
}

#Print out the stats by read type
print "\n\n";
print "Raw Stats:\n";
my $capture = `cat tempRaw | ~mcwong/listStats.pl`;
print "$capture";
`rm tempRaw`;
print "\n\n";
print "Merged Stats:\n";
my $capture = `cat tempMerged | ~mcwong/listStats.pl`;
print "$capture";
`rm tempMerged`;
print "\n\n";
print "Mapped Stats:\n";
$capture = `cat temppick | ~mcwong/listStats.pl`;
print "$capture";
`rm temppick`;
print "\n\n";
print "Low Identity Hits:\n";
$capture = `cat templow | ~mcwong/listStats.pl`;
print "$capture";
`rm templow`;
print "\n\n";
print "Unmapped Reads:\n";
$capture = `cat tempmissed | ~mcwong/listStats.pl`;
print "$capture";
`rm tempmissed`;
