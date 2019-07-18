#!/usr/bin/perl -w
use strict;
use IO::Zlib;
use IPC::Run 'run';
use IPC::Pipeline;

# bowtie does crazy things for read lenght 10; lots of mappings with good quality to wrong positions.

# config:
my $minLen = 20;
my $maxLen = 150;
my $nThreads = 56;
my $assembly = "danRer10";
my $faF = "/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa.gz";

#-rw-r--r-- 1 piotr lenhard 81021094371 Jul 18 12:36 alreadyMappable.danRer10.txt

# input args:
my $alreadyMappableF = shift @ARGV;
my %needToCheck = ();

## read chromosomes in:
# just first scaffold:
my %genome = (); #keys {chr}, value: seq
my $seq = "";
my $chr = "";
my $faH = IO::Zlib->new($faF, 'rb') or die;
while(<$faH>)
{
	chomp;
	if ($_ =~ /^>(\S+)/)
	{
		$genome{$chr} = $seq;
		$seq = "";
		$chr = $1;
	}
	else
	{
		$seq .= uc $_;
	}
}
$faH->close() ;
$genome{$chr} = $seq;

### Optimisation: solve really bad regions by finding unmapped positions at
### maxReadlength+1. Save mappable positions to %needToCheck

mapWithBowtie($maxLen+1, spitOut($maxLen+1, 1), 1, "/dev/null");

### the main loop over read length:
### - chop the genome
### - map it
### - update noNeedToCheck, and the file
for(my $rLen = $minLen; $rLen <= $maxLen; $rLen++)
{
	print STDERR "Main Loop iter $rLen\n";
	my $fastqF = spitOut($rLen, 0);
	print STDERR "Mapping\n";
	mapWithBowtie($rLen, $fastqF, 0, $alreadyMappableF);
	# unlink($fastqF);
}

sub mapWithBowtie
{
	my($rLen, $infile, $reverse, $outF) = @_;
	my @pids = pipeline(my $bowtieIn, my $out, my $err,
	    #["srun", "-c", $nThreads, "--mem", "20G", "-J", "mappab.$rLen", "-w", "kraken",
		["bowtie2", "--phred33-quals", "--threads", $nThreads, "--no-head",
		"-x", "/mnt/biggles/data/alignment_references/bowtie2/$assembly",
		"--no-unal", "-U", $infile],
		["samtools", "view", "-S", "-q", 10, "-"],
		["./bam2mappablePositions2"]);

	close($bowtieIn);
	open(my $oH, ">>", $outF) or die;
	while(<$out>)
	{
		chomp;
		my($chr, $pos) = split(/\t/);
		# update the the hash:
		if($reverse)
		{
			$needToCheck{$chr}{$pos} = 1;
		}
		else
		{
			# we just mapped it successfully. No need anymore:
			delete($needToCheck{$chr}{$pos});
			# update the big file:
			print $oH $chr, "\t", $pos, "\t", $rLen, "\n";
		}
	}
	close($out);
	close($oH);
	close($err);

	#print STDERR "waiting for pids @pids\n";
	foreach my $pid (@pids)
	{
		waitpid($pid, 1);
	}
}

sub spitOut
{
	# the output uses Granges, sam coordinates :(
	my($readLen, $doWholeGenomeBool) = @_;
	my $outF = "chop".$readLen.".fastq";
	open(my $oH, ">", $outF) or die;
	foreach my $chr (sort keys %genome)
	{
		my $seq = $genome{$chr};
		for(my $i=1; $i<=length($seq) - $readLen + 1; $i++)
		{
			if($doWholeGenomeBool || defined $needToCheck{$chr}{$i})
			{
				print $oH '@'.$chr.":".($i)."-".($i+$readLen-1), "\n";
				print $oH substr($seq, $i-1, $readLen), "\n";
				print $oH "+\n";
				for(my $j=$readLen; $j; $j--)
				{
					print $oH "I";
				}
				print $oH "\n"
			}
		}
	}
	close($oH);
	return($outF);
}
