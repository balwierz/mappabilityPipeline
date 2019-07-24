#!/usr/bin/perl -w
use strict;
use IO::Zlib;
use File::Temp qw/ tempfile tempdir /;

# uses 85GB for danRer19
#14:16 start mapping 31
#14:28 start chopping 32	// 12 min of mapping
#14:43 start mapping 32		// 15 min of chopping

# bowtie2 does crazy things for read length 10; lots of mappings with good quality to wrong positions.
# note that botwie2 by default uses seed length 21 or 22 bp

# input args:
my $assembly = shift @ARGV;  # eg "danRer10"
my $alreadyMappableF = shift @ARGV;  # output file

# config:
my $minLen = 20;
my $maxLen = 150;
my $nThreads = 8;	# bowtie threads. There will be 3 more processes: samtools, bam2mappablePositions and this script.
my $faF = (-e "/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa.gz") ?
			"/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa.gz" :
			"/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa";
my $scratchDir = "/mnt/scratch/" . $ENV{USER};

# params to programs:
my $bowtieParams   = 	"--mm --phred33-quals --threads $nThreads" .
						" -x /mnt/biggles/data/alignment_references/bowtie2/$assembly --no-unal";
my $samtoolsParams = "-S -q 10";	# you can increase quality to 20 or even 30 if you wish

### read chromosomes in:
my %genome = %{readGenome($faF)}; #keys {chr}, value: seq

### init:
die "Cannot fine $faF" if ! -e $faF;
my %needToCheck = ();
my $tmpDir = tempdir( TEMPLATE => 'mappability.XXXXX', DIR => $scratchDir, CLEANUP => 1 );

### Optimisation: solve really bad regions by finding unmapped positions at
### maxReadlength+1. Save mappable positions to %needToCheck
my $wholeGenomeLongChopF = spitOut($maxLen+1, 1);
mapWithBowtie($maxLen+1, $wholeGenomeLongChopF, 1, "/dev/null");

### the main loop over read length:
### - chop the genome
### - マップをしながら　ファイルとneedToCheckを アップデート する
for(my $rLen = $minLen; $rLen <= $maxLen; $rLen++)
{
	print STDERR "Main Loop iter $rLen\n";
	my $fastqF = spitOut($rLen, 0);
	print STDERR "Mapping\n";
	mapWithBowtie($rLen, $fastqF, 0, $alreadyMappableF);
	unlink($fastqF);
}
### end of main()

### SUBs:
sub mapWithBowtie
{
	my($rLen, $infile, $reverse, $outF) = @_;
	open(BOWTIE,	"bowtie2 $bowtieParams -U $infile | " .
					"samtools view $samtoolsParams | " .
					"./bam2mappablePositions2 |") or die "cannot start mapping";

	open(my $oH, ">>", $outF) or die;
	while(<BOWTIE>)
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
	close(BOWTIE);
	close($oH);
}

sub spitOut
{
	# the output uses Granges, sam coordinates :(
	my($readLen, $doWholeGenomeBool) = @_;
	print STDERR "chopping $readLen\n";
	my $outF = "$tmpDir/" . join('.', "chop", $assembly, $readLen, "fastq");
	open(my $oH, ">", $outF) or die;
	if($doWholeGenomeBool)
	{
		foreach my $chr (sort keys %genome)
		{
			my $seq = $genome{$chr};
			for(my $i=1; $i<=length($seq) - $readLen + 1; $i++)
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
	else
	{
		foreach my $chr (keys %needToCheck)
		{
			my $seq = $genome{$chr};
			foreach my $i (keys %{$needToCheck{$chr}})
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

sub readGenome
{
	my($faF) = @_;	
	my $seq = "";
	my $chr = "";
	my %genome;
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
	return \%genome;
}
