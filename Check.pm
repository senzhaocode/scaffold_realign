package Check;
use strict;
use warnings;

	sub read_length {
		my ($fastq_1, $fastq_2) = @_;
		my $tmp_1 = `head -n 2 $fastq_1`; chomp $tmp_1; my $seq_1 = (split /\n/, $tmp_1)[1];
		my $tmp_2 = `head -n 2 $fastq_2`; chomp $tmp_2; my $seq_2 = (split /\n/, $tmp_2)[1];
		my $len_1 = length($seq_1); my $len_2 = length($seq_2);

		if ( $len_1 == $len_2 ) {
			return($len_1);
		} else {
			print "read length of fastq file is error\n"; exit;
			return(0);
		}
	}

	sub extract_scaffold { # Load scaffold seq; at current stage, only one breakpoint sequence was allowed to involve in alignment 
		my ($read_l, $scaffold, $ref) = @_;
		$/ = ">";
		open (IN, "$scaffold") || die "cannot load scaffold sequence dataset:$!\n";
		while ( <IN> ) {
			chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
			my @array = (split /\n/, $_);
			my $header = shift @array; $header =~s/\>//g;
			my $seq = join '', @array;
			my $front; my $back;
			if ( $seq =~/\|/ ) { 
				($front, $back) = (split /\|/, $seq)[0, 1];
			} elsif ( $seq =~/\*/ ) {
				($front, $back) = (split /\*/, $seq)[0, 1];
			} else {
				print "Breakpoint sequence has wrong split symbol, please use '|' or '*' as breakpoint symbol\n"; exit;
			}
			
			# Compare the length of given breakpoint seq with read length
			my $front_frag; my $front_tag = 0;
			if ( length($front) >= $read_l ) { # if given length of front seq >= read_length; trim length of front seq identical to read_length
				$front_frag = substr($front, -$read_l);
			} else { # if given length of front seq < read_length; mark $front_tag as value of "$read_l - length($front)". 
				$front_frag = $front;
				$front_tag = $read_l - length($front);
			}
			my $back_frag; my $back_tag = 0;
			if ( length($back) >= $read_l ) { #if given length of back seq >= read_length; trim length of back seq identical to read_length
				$back_frag = substr($back, 0, $read_l);
			} else { #if given length of back seq < read_length; mark $back_tag as value of "$read_l - length($back)". 
				$back_frag = $back;
				$back_tag = $read_l - length($back);
			}
			if ( exists($ref->{$header}) ) { print "\nNote: $header have replicates in scaffold sequence list\n"; }
			$ref->{$header}[0] = [$front_tag, $front_frag]; # breakpoint seq of GeneA
			$ref->{$header}[1] = [$back_tag, $back_frag]; # breakpoint seq of GeneB

			#set reverse complementary sequence
			$front_frag =~tr/ATCG/TAGC/; $front_frag=reverse($front_frag);
			$back_frag =~tr/ATCG/TAGC/; $back_frag=reverse($back_frag);
			$ref->{$header}[2] = [$front_tag, $front_frag]; # reverse complementary breakpoint seq of GeneA
			$ref->{$header}[3] = [$back_tag, $back_frag]; # reverse complementary breakpoint seq of GeneB
			#print "original: $ref->{$header}[0][1] $ref->{$header}[1][1]; reverse_complementary: $ref->{$header}[2][1] $ref->{$header}[3][1]\n"; #--testing
		}
		close IN; 
	}

	sub judge_scaffold {
		my ($scaffold, $tran_seq, $gene, $geneA, $geneB) = @_;

		# Load gene symbol file -- download from ensembl Biomart, Dec 4th 2016
		my %name;
		$/ = "\n";
		open (IN, "cut -s -f1,5 $gene | sort | uniq |") || die "cannot open gene symbol / ensembl id file:$!\n";
		while ( <IN> ) {
			chomp $_; my ($ensembl, $symbol) = (split /\t/, $_)[0, 1];
			$name{$symbol} = $ensembl; # key => gene_symbol; value => ensembl_id
		}
		close IN;

		# Load transcript (cdna seq) -- download from ensembl Biomart, Feb 2017
		my %cdna;
		$/ = ">";
		open (IN, "$tran_seq") || die "cannot load transcript sequence dataset:$!\n";
		while ( <IN> ) {
			chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
			my @array = (split /\n/, $_);
			my $header = shift @array; $header =~s/\>//g;
			my ($gene, $tran) = (split /\|/, $header);
			if (! $gene =~/ENSG/ ) { print "gene name $gene is wrong\n"; exit; }
			if (! $tran =~/ENST/ ) { print "transcript name $tran is wrong\n"; exit; }
			my $seq = join '', @array; my $len = length($seq);
			$cdna{$gene}{$len} = [$tran, $seq]; #print "$gene $len $tran\n$seq\n"; #--testing
		}
		close IN;
		$/ = "\n";

		# Judge whether input gene names are validated in the database; and whether breakpoint sequences match to the cDNA sequences
		my %scaff_final; # %scaff_final data structure: 
				 # $scaff_final{"GeneA"} = [sequence, transcript_id, break_seq_geneB]
				 # $scaff_final{"GeneB"} = [sequence, transcript_id, break_seq_geneB]
		my $tag_scaff_A; # if == 0, scaffold sequence of geneA fails to match cDNA of geneA; if == 1, pass
		my $tag_scaff_B; # if == 0, scaffold sequence of geneA fails to match cDNA of geneB; if == 1, pass 
		if ( $geneA =~/ENSG/ && $geneB =~/ENSG/ )  { # gene names as ensembl id
			if ( exists($cdna{$geneA}) ) {
				print "The input of $geneA as Ensembl gene id\n";	my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
				$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
			} else {
				print "Ensembl gene id $geneA is not in the list, please be sure the valid gene id\n"; exit;
			}

			if ( exists($cdna{$geneB}) ) {
				print "The input of $geneB as Ensembl gene id\n";	my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
				$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
			} else {
				print "Ensembl gene id $geneB is not in the list, please be sure the valid gene id\n"; exit;
			}
		} else { # gene names as gene symbol
			if ( exists($name{$geneA}) ) {
				if ( exists($cdna{$name{$geneA}}) ) {
					print "The input of $geneA as Refseq gene name\n";	my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
					$tag_scaff_A = &output_geneA($cdna{$name{$geneA}}, $scaffold, $scaff_final{$geneA});
				} else {
					print "Refseq gene name $geneA is not in the list, please be sure the correct gene name\n"; exit;
				}
			} else {
				print "Refseq gene name $geneA is not in the list, please be sure the correct gene name\n"; exit;
			}

			if ( exists($name{$geneB}) ) {
				if ( exists($cdna{$name{$geneB}}) ) {
					print "The input of $geneB as Refseq gene name\n";	my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
					$tag_scaff_B = &output_geneB($cdna{$name{$geneB}}, $scaffold, $scaff_final{$geneB});
				} else {
					print "Refseq gene name $geneB is not in the list, please be sure the correct gene name\n"; exit;
				}
			} else {
				print "Refseq gene name $geneB is not in the list, please be sure the correct gene name\n"; exit;
			}
		}

		return($tag_scaff_A, $tag_scaff_B, \%scaff_final);
	}

	sub output_geneB {
		my ($cdna_seq, $scaff_ref, $scaff_seq) = @_;
		
		my $tag = 0; my %hash; # collect all cdna with matched breakpoint seq for geneB 
		foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

			my $split_tag = $scaff_ref->[1][1];
			my @cdna_geneB = split /$split_tag/, $cdna_seq->{$rank}[1]; # print ">$cdna_seq->{$rank}[0]\n$cdna_seq->{$rank}[1]\n";
			if ( scalar(@cdna_geneB) == 2 ) {
				$tag = 1;
				if ( $scaff_ref->[1][0] == 0 ) { # length of breakpoint == length of read
					my $len_geneB = length($scaff_ref->[1][1]);
					push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
				} else {
					my $add = substr($cdna_geneB[1], 0, $scaff_ref->[1][0]); $add = $scaff_ref->[1][1].$add;
					my $len_geneB = length($add);
					push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneB]
				}
				next;
			} else {
				my $len = length($split_tag) - 3;
				$split_tag = substr($split_tag, -$len); # trim the frist 3 base of breakpoint sequence
				@cdna_geneB = split /$split_tag/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneB) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[1][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneB = length($scaff_ref->[1][1]);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
					} else {
						my $add = substr($cdna_geneB[1], 0, $scaff_ref->[1][0]); $add = $scaff_ref->[1][1].$add;
						my $len_geneB = length($add);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneB]
					}
					next;
				}
			}
		}

		if (! %hash ) {
			$tag = 0; # re-set $tag = 0;
			foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

				my $split_tag_r = $scaff_ref->[3][1];
				my @cdna_geneB_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneB_r) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[3][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneB = length($scaff_ref->[1][1]);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
					} else {
						my $add = substr($cdna_geneB_r[0], -$scaff_ref->[3][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $scaff_ref->[1][1].$add;
						my $len_geneB = length($add);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; #[sequence, transcript_id, break_seq_geneB]
					}
					next;
				} else {
					my $len_r = length($split_tag_r) - 3;
					$split_tag_r = substr($split_tag_r, 0, $len_r); # trim the last 3 base of breakpoint sequence
					@cdna_geneB_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
					if ( scalar(@cdna_geneB_r) == 2 ) {
						$tag = 1;
						if ( $scaff_ref->[3][0] == 0 ) { # length of breakpoint == length of read
							my $len_geneB = length($scaff_ref->[1][1]);
							push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
						} else {
							my $add = substr($cdna_geneB_r[0], -$scaff_ref->[3][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $scaff_ref->[1][1].$add;
							my $len_geneB = length($add);
							push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneB]
						}
						next;
					}
				}
			}
		}

		if ( $tag == 1 ) {
			foreach my $len ( reverse sort {$a <=> $b} keys %hash ) {
				foreach my $ref ( @{$hash{$len}} ) { # select the first element [0] of %hash (the longest transcript match to breakpoint region of geneB)
					$scaff_seq->[0] = $ref->[0]; $scaff_seq->[1] = $ref->[1]; $scaff_seq->[2] = $ref->[2]; # [sequence, transcript_id, break_seq_geneB]
					return(1);
				}
			}
		} else {
			return(0);
		}
	}

	sub output_geneA {
		my ($cdna_seq, $scaff_ref, $scaff_seq) = @_;

		my $tag = 0; my %hash; # collect all cdna with matched breakpoint seq for geneA
		foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

			my $split_tag = $scaff_ref->[0][1];
			my @cdna_geneA = split /$split_tag/, $cdna_seq->{$rank}[1]; # print ">$cdna_seq->{$rank}[0]\n$cdna_seq->{$rank}[1]\n";
			if ( scalar(@cdna_geneA) == 2 ) {
				$tag = 1;
				if ( $scaff_ref->[0][0] == 0 ) { # length of breakpoint == length of read
					my $len_geneA = length($scaff_ref->[0][1]);
					push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
				} else { # length of breakpoint < length of read; need to get extension
					my $add = substr($cdna_geneA[0], -$scaff_ref->[0][0]); $add = $add.$scaff_ref->[0][1]; 
					my $len_geneA = length($add);
					push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
				}
				next;
			} else {
				my $len = length($split_tag) - 3; 
				$split_tag = substr($split_tag, 0, $len); # trim the last 3 base of breakpoint sequence
				@cdna_geneA = split /$split_tag/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneA) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[0][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneA = length($scaff_ref->[0][1]);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
					} else {
						my $add = substr($cdna_geneA[0], -$scaff_ref->[0][0]); $add = $add.$scaff_ref->[0][1];
						my $len_geneA = length($add);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
					}
					next;
				}
			}
		}

		if (! %hash) {
			$tag = 0;
			foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

				my $split_tag_r = $scaff_ref->[2][1];
				my @cdna_geneA_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneA_r) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[2][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneA = length($scaff_ref->[0][1]);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
					} else { # length of breakpoint < length of read; need to get extension
						my $add = substr($cdna_geneA_r[1], 0, $scaff_ref->[2][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $add.$scaff_ref->[0][1];
						my $len_geneA = length($add);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
					}
					next;
				} else {
					my $len_r = length($split_tag_r) - 3; 
					$split_tag_r = substr($split_tag_r, -$len_r); # trim the first 3 base of breakpoint sequence
					@cdna_geneA_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
					if ( scalar(@cdna_geneA_r) == 2 ) {
						$tag = 1;
						if ( $scaff_ref->[2][0] == 0 ) { # length of breakpoint == length of read
							my $len_geneA = length($scaff_ref->[0][1]);
							push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
						} else { # length of breakpoint < length of read; need to get extension
							my $add = substr($cdna_geneA_r[1], 0, $scaff_ref->[2][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $add.$scaff_ref->[0][1];
							my $len_geneA = length($add);
							push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
						}
						next;
					}
				}
			}
		}

		if ( $tag == 1 ) {
			foreach my $len ( reverse sort {$a <=> $b} keys %hash ) {
				foreach my $ref ( @{$hash{$len}} ) { # select the first element [0] of %hash (the longest transcript match to breakpoint region of geneA)
					$scaff_seq->[0] = $ref->[0]; $scaff_seq->[1] = $ref->[1]; $scaff_seq->[2] = $ref->[2]; # [sequence, transcript_id, break_seq_geneA]
					return(1);
				}
			}
		} else {
			return(0);
		}
	}
1;
