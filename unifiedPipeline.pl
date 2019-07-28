#!/usr/bin/perl -w
use strict;
use IO::Zlib;
use File::Temp qw/ tempfile tempdir /;

# bowtie2 does crazy things for read length 10; lots of mappings with good quality to wrong positions.
# note that botwie2 by default uses seed length 21 or 22 bp

# input args:
my $assembly = shift @ARGV;  # eg "danRer10"
my $alreadyMappableF = shift @ARGV;  # output file
my $nThreads = shift @ARGV or die;	# bowtie threads. 
					# There will be 3 more processes: samtools, bam2mappablePositions and this script.

# config:
my $minLen = 20;
my $maxLen = 150;
my $faF = (-e "/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa.gz") ?
			"/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa.gz" :
			"/mnt/biggles/data/UCSC/goldenpath/$assembly/bigZips/$assembly.fa";
my $scratchDir = "/mnt/scratch/" . $ENV{USER};

# params to programs:
my $bowtieParams   = 	" --phred33-quals --threads $nThreads" .
						" -x /mnt/biggles/data/alignment_references/bowtie2/$assembly --no-unal";
my $samtoolsParams = "-S -q 10";	# you can increase quality to 20 or even 30 if you wish

### init:
die "Cannot fine $faF" if ! -e $faF;
my %needToCheck = (); # {chr}[pos]
	#meaning: -1=unmappable even at maxLen; 0=need to check; non-zero=the result
my $tmpDir = tempdir( TEMPLATE => 'mappability.XXXXX', DIR => $scratchDir, CLEANUP => 0 );

### read chromosomes in:
my %genome = %{readGenome($faF, \%needToCheck)}; #keys {chr}, value: seq

### Optimisation: solve really bad regions by finding unmapped positions at
### maxReadlength. Save mappable positions to %needToCheck
my $wholeGenomeLongChopF = spitOut($maxLen, 1);
mapWithBowtie($maxLen, $wholeGenomeLongChopF, 1, "/dev/null");
unlink($wholeGenomeLongChopF);

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

# post-process:
produceBigWig();

### end of main()

### SUBs:
sub produceBigWig
{
	my $sortOutF = "$tmpDir/$alreadyMappableF.bedGraph";
	my $chrLenF  = "$tmpDir/chromInfo.txt";
	open(my $sortOutH, ">", $sortOutF) or die;
	open(my $chrLenH, ">", $chrLenF) or die;
	foreach my $chr (sort keys %genome)
	{
		print $chrLenH $chr, "\t", length($genome{$chr}), "\n";
		for(my $i=1; $i<=length($genome{$chr}); $i++)
		{
			my $thisVal = (-1 == $needToCheck{$chr}[$i]) ? $needToCheck{$chr}[$i] : $maxLen;
			print $sortOutH join("\t", $chr, $i-1, $i, $thisVal), "\n";
		}
	}
	close $sortOutH;
	close $chrLenH;

	system("bedGraphToBigWig", $sortOutF, $chrLenF, "$alreadyMappableF.bw");
}

sub mapWithBowtie
{
	my($rLen, $infile, $reverse, $outF) = @_;
	$rLen = int($rLen);  #maybe it saves us some memory instead of strings???
	open(BOWTIE,	"bowtie2 $bowtieParams -U $infile | " .
					"samtools view $samtoolsParams | " .
					"./bam2mappablePositions2 |") or die "cannot start mapping: $!";

	#open(my $oH, ">>", $outF) or die;
	while(<BOWTIE>)
	{
		chomp;
		my($chr, $pos) = split(/\t/);
		# update the the hash:
		if($reverse)
		{
			$needToCheck{$chr}[$pos] = 0;
		}
		else
		{
			# we just mapped it successfully. Let's save the result
			$needToCheck{$chr}[$pos] = $rLen;
			# update the big file (obsolete right now!):
			# print $oH $chr, "\t", $pos, "\t", $rLen, "\n";
		}
	}
	close(BOWTIE);
	#close($oH);
}

sub spitOut
{
	# the output uses Granges, sam coordinates :(
	my($readLen, $doWholeGenomeBool) = @_;
	print STDERR "chopping $readLen\n";
	my $outF = "$tmpDir/" . join('.', "chop", $assembly, $readLen, "fastq");
	open(my $oH, ">", $outF) or die;
	foreach my $chr (sort keys %genome)
	{
		my $seq = $genome{$chr};
		for(my $i=1; $i<=length($seq) - $readLen + 1; $i++)
		{
			if($doWholeGenomeBool || $needToCheck{$chr}[$i] == 0)
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
	my($faF, $h) = @_;
	my $seq = "";
	my $chr = "";
	my %genome;
	my $faH = IO::Zlib->new($faF, 'rb') or die;
	while(<$faH>)
	{
		chomp;
		if ($_ =~ /^>(\S+)/)
		{
			if($chr)
			{
				$genome{$chr} = $seq;
				@{$h->{$chr}} = (-1) x (1+length($seq));
			}
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
	@{$h->{$chr}} = (-1) x (1 + length($seq));
	return \%genome;
}
