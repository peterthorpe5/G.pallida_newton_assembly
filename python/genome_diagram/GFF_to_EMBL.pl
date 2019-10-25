#!/usr/local/bin/perl -w
# Convert GFF files (v3) to EMBL
#
use strict;
use lib "/nfs/users/nfs_a/ar11/scripts";
use Seq;

$| = 1;

my($fasta, $gff_file) = @ARGV;

if(!defined $fasta)
{
	print STDERR "Usage: gff2embl_jigsaw.pl <fasta> <gff_file>\n";
	exit;
}

open(FA, "<$fasta") or die "$!";
my @seqs = <FA>;
my $seqs = Seq::parse_fasta_id(\@seqs);
close FA;


my %cds = ();
my %cds_strand = ();
my %gene_desc = ();

if(defined $gff_file && -e $gff_file)
{
	open(GFF, "<$gff_file") || die "$!";

	while(<GFF>)
	{
		chomp;
	
		next if /^#/;
	
		my($seq, $type, $feature, $start, $end, $score, $strand, $frame, $gene_name) = split /\t/;
	
		next unless $feature eq 'CDS';

		#$gene_name =~ /gene_id \"(\S+)\"/;
		$gene_name =~ /Parent=(\S+)/;

		my $id = $1;

		#print "$gene_name\t$id\t$start\t$end\t$frame\t$strand\n";
	
		push @{$cds{$seq}->{$id}}, [$start, $end, $frame+1];
		$cds_strand{$id} = $strand;
	}
	close GFF;
}

foreach my $contig (sort keys %$seqs)
{
	my $filename = $contig;

	print  "ID   $contig\n";
	print  "FH   Key             Location/Qualifiers\n";
	print  "FH\n";

	foreach my $cds_id (sort keys %{$cds{$contig}})
	{
		my @exons = ();

		my $frame;

		my $first_frame;
		my $last_frame;

		foreach my $pos (sort {$a->[0]<=>$b->[0]} @{$cds{$contig}->{$cds_id}})
		{
			push @exons, "$pos->[0]\.\.$pos->[1]";

			$first_frame = $pos->[2] unless defined $first_frame;
			$last_frame = $pos->[2];
		}

		$frame = $cds_strand{$cds_id} eq '+' ? $first_frame : $last_frame; 

		my $exon_string = join ',', @exons;

		if(scalar @exons > 1)
		{
			$exon_string = 'join('.$exon_string.')';
		}

		if($cds_strand{$cds_id} eq '-')
		{
			$exon_string = 'complement('.$exon_string.')';
		}

		print "FT   CDS             $exon_string\n";
		print "FT                   \/locus_tag\=\"$cds_id\"\n";
		print "FT                   \/desc\=\"$gene_desc{$cds_id}\"\n" if exists $gene_desc{$cds_id};
		print "FT                   \/codon_start\=$frame\n";
	}

	my $length = length $seqs->{$contig};
	print "SQ   Sequence $length BP\;\n";

	my $seq_lines = roundup($length / 60);

	my $embl_format = '';
	my $end_flag = 0;

	for my $line (1..$seq_lines)
	{
		$embl_format .= '     ';

		for my $col(1..6)
		{
			next if $end_flag;

			my $start = (($line - 1) * 60 ) + (($col - 1) * 10);

			my $seq_length = 10;

			if(($start + $seq_length) > $length)
			{
				$seq_length = $length - $start;
				$end_flag = 1;
			}

			$embl_format .= substr($seqs->{$contig}, $start, $seq_length);
			$embl_format .= ' ';

		}
	
		my $bp_label;
		my $spacing;

		if($end_flag)
		{
			$bp_label = $length;
			my $label_length = length $bp_label;
			my $left_over = ($length % 60);
			my $extra_gaps = int($left_over / 10);
			$spacing = 79 - 5 - $extra_gaps - $left_over - $label_length;

		}
		else
		{
			$bp_label = $line * 60;
			my $label_length = length $bp_label;
			$spacing = 79 - 10 - 60 - $label_length;
		}

		$embl_format .= ' ' x $spacing ."$bp_label\n";

		my $embl_length = length $embl_format;
	}
	$embl_format .= "\/\/\n";

	print "$embl_format";
}

sub roundup {
    my $n = shift;
        return(($n == int($n)) ? $n : int($n + 1))
}

