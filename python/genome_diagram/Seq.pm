package Seq;

sub revcomp
{
	my ($seq) = @_;

	my $revcomp = reverse($seq);
  	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
  	return $revcomp;
}

sub print_multiline
{
	my($seqs, $length) = @_;

	$length = 100 unless defined $length;

	my $start = 0;
	my $end = $start + $length - 1;

	foreach my $id (sort keys %$seqs)
	{
		my $start = 0;
		my $end = $start + $length - 1;

		print ">$id\n";

		my @seq = split //, $seqs->{$id};

		while($start <= $end)
		{
			my $line = join "", @seq[$start..$end];

			print "$line\n";

			$start += $length;
			$end = $end + $length > @seq  - 1 ? scalar @seq - 1: $end + $length;
		}
	}
	
}

sub embl2seqs
{
	my($embl_file) = @_;

	my %seqs = ();

	my $stream = Bio::SeqIO->new(-file => $embl_file, -format => 'EMBL');

	while ( (my $seq = $stream->next_seq()) )
	{
		my @features = $seq->get_SeqFeatures();

		my $current_id = '';

		foreach my $feat (@features)
		{
			next unless $feat->primary_tag eq 'CDS';

			foreach my $tag ( $feat->get_all_tags() )
			{
				if ($tag eq 'systematic_id' || $tag eq 'locus_tag')
				{
					$current_id = join(' ',$feat->get_tag_values($tag));
				}
			}

				my $cds_obj = $feat->spliced_seq;
				my $cds_seq = $cds_obj->seq;

				#print STDERR "$current_id\t$cds_seq\n";

				my $seq_length = length $cds_seq;

				if($seq_length % 3 != 0)
				{
					print STDERR "The sequence for $current_id is $seq_length bp long and is not a whole number of codons!\n";
					next;
				}

				my $substr_length = $seq_length - 3;

				#print $substr_length."\n";

				$cds_seq = substr($cds_seq, 0, $substr_length);


				my $translation = $cds_obj->translate();

				my $aa = $translation->seq();

				my $aa_length = length $aa;

				$aa = substr($aa, 0, $aa_length - 1);

				#Skip those sequences containing STOP codons
				if ($aa =~ /\*/)
				{
					print STDERR "$current_id contains a STOP codon and is being skipped\n";
				}
				else
				{
					$seqs{$current_id}->{'aa'} = $aa;
					$seqs{$current_id}->{'nt'} = $cds_seq;
				}
		}
	}

	return \%seqs;
}

sub get_aa2codon
{
	my %aa2codon = (
		'PHE'	=> ['UUU', 'UUC'],
		'LEU'	=> ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
		'ILE'	=> ['AUU', 'AUC', 'AUA'],
		'MET'	=> ['AUG'],
		'VAL'	=> ['GUU', 'GUC', 'GUG', 'GUA'],
		'SER'	=> ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
		'PRO'	=> ['CCU', 'CCC', 'CCA', 'CCG'],
		'THR'	=> ['ACU', 'ACC', 'ACA', 'ACG'],
		'ALA'	=> ['GCU', 'GCC', 'GCA', 'GCG'],
		'TYR'	=> ['UAU', 'UAC'],
		'CYS'	=> ['UGU', 'UGC'],
		'TRP'	=> ['UGG'],
		'HIS'	=> ['CAU', 'CAC'],
		'GLN'	=> ['CAA', 'CAG'],
		'ASN'	=> ['AAU', 'AAC'],
		'LYS'	=> ['AAA', 'AAG'],
		'ASP'	=> ['GAU', 'GAC'],
		'GLU'	=> ['GAA', 'GAG'],
		'ARG'	=> ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
		'GLY'	=> ['GGU', 'GGC', 'GGA', 'GGG']
	);

	return \%aa2codon;
}

sub translate
{
	my $seq = shift;

	my %codon2aa = (
		'TTT' => 'F',
		'TTC' => 'F',
		'TTA' => 'L',
		'TTG' => 'L', 
		'CTT' => 'L',
		'CTC' => 'L',
		'CTA' => 'L', 
		'CTG' => 'L',
		'ATT' => 'I', 
		'ATC' => 'I', 
		'ATA' => 'I',
		'ATG' => 'M',
		'GTT' => 'V', 
		'GTC' => 'V',
		'GTG' => 'V', 
		'GTA' => 'V',
		'TCT' => 'S', 
		'TCC'  => 'S', 
		'TCA'  => 'S',
		'TCG'  => 'S',
		'AGT'  => 'S', 
		'AGC'  => 'S',
		'CCT' => 'P', 
		'CCC' => 'P',
		'CCA' => 'P', 
		'CCG' => 'P',
		'ACT' => 'T',
		'ACC' => 'T', 
		'ACA' => 'T', 
		'ACG' => 'T',
		'GCT' => 'A', 
		'GCC' => 'A', 
		'GCA' => 'A', 
		'GCG' => 'A',
		'TAT' => 'Y',
		'TAC' => 'Y',
		'TGT' => 'C', 
		'TGC' => 'C',
		'TGG' => 'W',
		'CAT' => 'H', 
		'CAC' => 'H',
		'CAA' => 'Q', 
		'CAG' => 'Q',
		'AAT' => 'N', 
		'AAC' => 'N',
		'AAA' => 'K', 
		'AAG' => 'K',
		'GAT' => 'D', 
		'GAC' => 'D',
		'GAA' => 'E', 
		'GAG' => 'E',
		'CGT' => 'R', 
		'CGC' => 'R', 
		'CGA' => 'R', 
		'CGG' => 'R', 
		'AGA' => 'R', 
		'AGG' => 'R',
		'GGT' => 'G', 
		'GGC' => 'G', 
		'GGA' => 'G', 
		'GGG' => 'G',
		'NNN' => 'X',
		'TGA' => '*',
		'TAG' => '+',
		'TAA' => '#'
	);

	my @seq = split //, $seq;

	my $len = length $seq;

	my $remainder = $len % 3;
	
	my $its = int($len / 3);

	my $aa_seq = '';

	for(0..$its)
	{
		my $codon = substr($seq, $_*3, 3);

		if(exists $codon2aa{$codon})
		{
			$aa_seq .= $codon2aa{$codon};
		}
		elsif($codon eq 'TGA' || $codon eq 'TAG' || $codon eq 'TAA' || $codon eq '')
		{
			# Stop codon or no codon
		}
		else
		{
			print STDERR "No aa for $codon\n";
		}

	}
	return $aa_seq;
}

sub get_codon2aa
{
	my %codon2aa = (
		'UUU' => 'PHE',
		'UUC' => 'PHE',
		'UUA' => 'LEU',
		'UUG' => 'LEU', 
		'CUU' => 'LEU',
		'CUC' => 'LEU',
		'CUA' => 'LEU', 
		'CUG' => 'LEU',
		'AUU' => 'ILE', 
		'AUC' => 'ILE', 
		'AUA' => 'ILE',
		'AUG' => 'MET',
		'GUU' => 'VAL', 
		'GUC' => 'VAL',
		'GUG' => 'VAL', 
		'GUA' => 'VAL',
		'UCU' => 'SER', 
		'UCC'  => 'SER', 
		'UCA'  => 'SER',
		'UCG'  => 'SER',
		'AGU'  => 'SER', 
		'AGC'  => 'SER',
		'CCU' => 'PRO', 
		'CCC' => 'PRO',
		'CCA' => 'PRO', 
		'CCG' => 'PRO',
		'ACU' => 'THR', 
		'ACC' => 'THR', 
		'ACA' => 'THR', 
		'ACG' => 'THR',
		'GCU' => 'ALA', 
		'GCC' => 'ALA', 
		'GCA' => 'ALA', 
		'GCG' => 'ALA',
		'UAU' => 'TYR',
		'UAC' => 'TYR',
		'UGU' => 'CYS', 
		'UGC' => 'CYS',
		'UGG' => 'TRP',
		'CAU' => 'HIS', 
		'CAC' => 'HIS',
		'CAA' => 'GLN', 
		'CAG' => 'GLN',
		'AAU' => 'ASN', 
		'AAC' => 'ASN',
		'AAA' => 'LYS', 
		'AAG' => 'LYS',
		'GAU' => 'ASP', 
		'GAC' => 'ASP',
		'GAA' => 'GLU', 
		'GAG' => 'GLU',
		'CGU' => 'ARG', 
		'CGC' => 'ARG', 
		'CGA' => 'ARG', 
		'CGG' => 'ARG', 
		'AGA' => 'ARG', 
		'AGG' => 'ARG',
		'GGU' => 'GLY', 
		'GGC' => 'GLY', 
		'GGA' => 'GLY', 
		'GGG' => 'GLY',
		'---' => 'XAA'
	);

	return \%codon2aa;
}

sub parse_fastq
{
	# Assume that sequences take up only one line because quality score lines can start with the same letters
	# as header lines!
	my $input = shift;
	my %seqs;
	my %quals;
	my $header;
	my $seq;

	my $c = 0;

	foreach my $line (@$input)
	{
		chomp $line;

		if($c == 0)
		{
			$line =~ /^\@(\S+)/;
			$header = $1;
		}
		elsif($c == 1)
		{
			$seqs{$header} = $line;
		}
		elsif($c == 3)
		{
			$quals{$header} = $line;
			$c = -1;
		}

		$c++;
	}

	return (\%seqs, \%quals);
}
sub parse_fasta
{
        $/ = "\n";

	my $input = shift;
	my %seqs;
        my $header;
	my $seq;

	foreach my $line (@$input)
        {
                chomp $line;

                if ($line =~ /^>(.*)/ && defined $seq)
                {
                        $seqs{$header} = $seq;
										                        $header = $1;
													$seq = '';
										                }
										                elsif ($line =~ /^>(.*)/)
		{
										                        $header = $1;
	        }
                else
		{
	                $seq .= $line;
	        }
	}

        $seqs{$header} = $seq;

        return \%seqs;
}

sub parse_fasta_id
{
        $/ = "\n";

	        my $input = shift;
		        my %seqs;
			        my $header;
				        my $seq;

					        foreach my $line (@$input)
						        {
							                chomp $line;

									                if ($line =~ /^>(\S+)/ && defined $seq)
											                {
													                        $seqs{$header} = $seq;
																                                                                                                        $header = $1;
																													                                                                                                        $seq = '';
																																										                                                                                                }
																																																						                                                                                                elsif ($line =~ /^>(\S+)/)
																																																																		                {
																																																																				                                                                                                        $header = $1;
																																																																																	                }
																																																																																			                else
																																																																																					                {
																																																																																							                        $seq .= $line;
																																																																																										                }
																																																																																												        }

																																																																																													        $seqs{$header} = $seq;

																																																																																														        return \%seqs;
																																																																																															}


sub three2one
{
	chomp(my $code = shift);

	$code = uc $code;

	my %code = (
			'ALA'	=> 'A',
			'CYS'	=> 'C',
			'ASP'	=> 'D',
			'GLU'	=> 'E',
			'PHE'	=> 'P',
			'GLY'	=> 'G',
			'HIS'	=> 'H',
			'ILE'	=> 'I',
			'LYS'	=> 'K',
			'LEU'	=> 'L',
			'MET'	=> 'M',
			'ASN'	=> 'N',
			'PRO'	=> 'P',
			'GLN'	=> 'Q',
			'ARG'	=> 'R',
			'SER'	=> 'S',
			'THR'	=> 'T',
			'VAL'	=> 'V',
			'TRP'	=> 'W',
			'TYR'	=> 'Y'
	);

	if(exists $code{$code})
	{
		return $code{$code};
	}
	else
	{
		return 0;
	}	
}

sub fasta2emblseq
{
	my $seq = shift;

	my $length = length($seq);

	my $embl_format = "SQ   Sequence $length BP\;\n";

	my $seq_lines = roundup($length / 60);

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


1;

