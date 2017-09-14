package Genome_align;
use strict;
use warnings;

	sub discordant_specif {
		my ($dis_ref, $multiple_ref, $path, $gene, $read_type) = @_; #read_type: spanning / discordant
		foreach my $id ( keys %{$dis_ref} ) {
			my $hit_1 = `grep "\@$id " -A 3 $path/${read_type}_1.txt`; chomp $hit_1;
                        my $hit_2 = `grep "\@$id " -A 3 $path/${read_type}_2.txt`; chomp $hit_2;
                        my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
                        my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
                        $dis_ref->{$id}[0][1] = $seq_1;
                        $dis_ref->{$id}[0][2] = $seq_1_tr; #reverse complementary
			$dis_ref->{$id}[0][3] = $hit_1;
                        $dis_ref->{$id}[1][1] = $seq_2;
                        $dis_ref->{$id}[1][2] = $seq_2_tr; #reverse complementary
			$dis_ref->{$id}[1][3] = $hit_2;
		}

		# collect mutiple hit reads or bad quality reads
		open (IN, "awk -F '\t' -v OFS='\t' '(\$8>=0 && \$5<40)' $path/tmp/${read_type}_sec.sam |") || die "cannot $read_type sam file:$!\n";
                while ( <IN> ) {
                        chomp $_; my ($id, $flag, $chr, $pos, $quality, $cigar, $seq) = (split /\t/, $_)[0,1,2,3,4,5,9];
                        push @{$multiple_ref->{$id}{$seq}}, [$chr, $pos, $quality, $cigar, $flag]; #
                }
                close IN;
		my %mapped; # mapped read with good quality (and remove the other end of reads with multiple hits)
                open (IN, "awk -F '\t' -v OFS='\t' '(\$8>=0 && \$5>=40)' $path/tmp/${read_type}_sec.sam |") || die "cannot $read_type sam file:$!\n";
                while ( <IN> ) {
                        chomp $_; my ($id, $flag, $chr, $pos, $quality, $cigar, $seq) = (split /\t/, $_)[0,1,2,3,4,5,9];
                        next if ( exists($multiple_ref->{$id}) );
                        push @{$mapped{$id}{$seq}}, [$chr, $pos, $quality, $cigar, $flag]; #
                }
                close IN;

		foreach my $id ( keys %{$dis_ref} ) {
			if ( exists($mapped{$id}) ) {
				my $tag = 1; # judgement: 1 => fail; 0 => accept
				foreach my $seq ( keys %{$mapped{$id}} ) {
					# read with multiple mapping positions (>=2)
					if ( scalar(@{$mapped{$id}{$seq}}) > 1 ) {
						if ( $mapped{$id}{$seq}[0][2] <= $mapped{$id}{$seq}[1][2] ) { # read has multiple hits with equal quality - filtered out
							print "${read_type}-filtered: $id show multiple hits with equal quality: $seq\n";
							delete($dis_ref->{$id}); last;
						}	
					}
					# read not primary sequence mapped - filtered out
					if ( $mapped{$id}{$seq}[0][4] >= 256 ) { 
						print "${read_type}-filtered: $id mapped not primary sequence: $seq\n";
						delete($dis_ref->{$id}); last;
					}
					# check whether read mapped to GeneA/GeneB if it has unique mapping position (tag == 1, filtered out)
					if ( $seq eq $dis_ref->{$id}[0][1] ) {
						&discordant_func($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); #print "$read_type: $id - $dis_ref->{$id}[0][0]\t$tag\n";
                                	} elsif ( $seq eq $dis_ref->{$id}[0][2] ) {
						&discordant_func($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); #print "$read_type: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
                                	}
                                	if ( $seq eq $dis_ref->{$id}[1][1] ) {
						&discordant_func($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); #print "$read_type: $id - $dis_ref->{$id}[1][0]\t$tag\n";
                                	} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
						&discordant_func($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); #print "$read_type: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
                                	}       
                                	if ( $tag == 1 ) {
                                        	delete($dis_ref->{$id}); last;
                                	}	
                                } 
                        } else {
                                if ( exists($multiple_ref->{$id}) ) {
                                        print "${read_type}-filtered: $id show unspecific mapping in genome alignment\n";
                                        delete($dis_ref->{$id});
                                } else {
                                        print "${read_type}: $id ($dis_ref->{$id}[0][0]-$dis_ref->{$id}[1][0]) show no mapping hits in genome alignment, but probably correct\n"; # no-splicing alignment is more sensitive than splicing
                                }       
                        }       
                }
	}

	sub singlton_specif {
		my ($dis_ref, $multiple_ref, $path, $gene, $read_type) = @_;
		foreach my $id ( keys %{$dis_ref} ) {
                        my $hit_1 = `grep "\@$id " -A 3 $path/${read_type}_1.txt`; chomp $hit_1;
                        my $hit_2 = `grep "\@$id " -A 3 $path/${read_type}_2.txt`; chomp $hit_2;
                        my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
                        my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
                        $dis_ref->{$id}[0][1] = $seq_1;
                        $dis_ref->{$id}[0][2] = $seq_1_tr; #reverse complementary
			$dis_ref->{$id}[0][3] = $hit_1;
                        $dis_ref->{$id}[1][1] = $seq_2;
                        $dis_ref->{$id}[1][2] = $seq_2_tr; #reverse complementary
			$dis_ref->{$id}[1][3] = $hit_2;
		}
		# collect all the singlton split reads
		open (IN, "awk  -F '\t' -v OFS='\t' '(\$8>=0)' $path/tmp/${read_type}_sec.sam |") || die "cannot ${read_type} sam file:$!\n";
                while ( <IN> ) {
                        chomp $_; my ($name, $flag, $chr, $pos, $quality, $cigar, $seq) = (split /\t/, $_)[0,1,2,3,4,5,9];
                        push @{$multiple_ref->{$name}{$seq}}, [$chr, $pos, $quality, $cigar, $flag];
                }
                close IN;

		foreach my $id ( keys %{$dis_ref} ) {
                        if ( exists($multiple_ref->{$id}) ) {
                                my $tag = 1; # judgememnt: 1 => fail; 0 => accept
                                foreach my $seq ( keys %{$multiple_ref->{$id}} ) {
					# read with multiple mapping positions (>=2)
					if ( scalar(@{$multiple_ref->{$id}{$seq}}) > 1 ) {
						if ( $multiple_ref->{$id}{$seq}[0][2] <= $multiple_ref->{$id}{$seq}[1][2] ) { # read has multiple hits with equal quality - filtered out
							print "${read_type}-filtered: $id show multiple hits with equal quality: $seq\n"; 
							delete($dis_ref->{$id}); last;
						}
					}
					# read not primary sequence mapped - filtered out
					if ( $multiple_ref->{$id}{$seq}[0][4] >= 256 ) { 
						print "${read_type}-filtered: $id mapped not primary sequence: $seq\n";
						delete($dis_ref->{$id}); last;
					}

                                        if ( $seq eq $dis_ref->{$id}[0][1] ) { 
                                                &single_func($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #print "singlton: $id - $dis_ref->{$id}[0][0]\t$tag\n";
                                        } elsif ( $seq eq $dis_ref->{$id}[0][2] ) { 
                                                my $name_1 = $dis_ref->{$id}[0][0];
                                                &single_func($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
						if ( $name_1 eq "NULL" ) {
							$dis_ref->{$id}[0][0] = $dis_ref->{$id}[0][0]."(reverse_complement)";
						} #print "singlton: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
                                        }
                                        if ( $seq eq $dis_ref->{$id}[1][1] ) { 
                                                &single_func($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #print "singlton: $id - $dis_ref->{$id}[1][0]\t$tag\n";
                                        } elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
                                                my $name_2 = $dis_ref->{$id}[1][0];
                                                &single_func($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
						if ( $name_2 eq "NULL" ) {
							$dis_ref->{$id}[1][0] = $dis_ref->{$id}[1][0]."(reverse_complement)";
						} #print "singlton: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
                                        }
                                        if ( $tag == 1 ) { delete($dis_ref->{$id}); last; }
                                }
                        } else {
                                print "singlton: $id - not mapped in genome\n";
                        }
                }
	}

#####################################
# substitune function within module #
#####################################
	sub discordant_func {
		my ($id_ref, $read_ref, $mapped_ref, $partner_ref, $tag_ref, $read_type) = @_; # Note: $tag_ref is a refrerence

		if ( $partner_ref->[0][0] eq $partner_ref->[1][0] ) { # geneA and geneB are in the same chromosome
			if ( $mapped_ref->[0][0] eq $partner_ref->[0][0] ) {
				if ( $partner_ref->[0][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					$$tag_ref = 0;
				} elsif ( $partner_ref->[1][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					$$tag_ref = 0;
				} else {
					print "${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				print "${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		} else { # geneA and geneB are in the different chromosome
			if ( $mapped_ref->[0][0] eq $partner_ref->[0][0] ) { # the chromosome of mapping read == that of geneA
				if ( $partner_ref->[0][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					$$tag_ref = 0;
				} else {
					print "${read_type}-filtered: $id_ref ($read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[0][3]\n";
					$$tag_ref = 1;
				}
			} elsif ( $mapped_ref->[0][0] eq $partner_ref->[1][0] ) { # the chromosome of mapping read == that of geneB
				if ( $partner_ref->[1][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					$$tag_ref = 0;
				} else {
					print "${read_type}-filtered: $id_ref ($read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				print "${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		}
	}

	sub single_func { # whether unmapped read in the frist round (no-splicing) was mapped to geneA or geneB in the second round (splicing)
                my ($id_ref, $read_ref, $multiple_ref, $partner_ref, $tag_ref, $read_type) = @_; #Note: $read_ref and $tag_ref are two references

		if ( $multiple_ref->[0][2] < 40 ) { # set mapping quality filtering
			print "$read_type-filtered: $id_ref ($$read_ref) with bad mapping quality\n";
			$$tag_ref = 1;
		}

		if ( $partner_ref->[0][0] eq $partner_ref->[1][0] ) { # geneA and geneB are in the same chromosome
			if ( $multiple_ref->[0][0] eq $partner_ref->[0][0] ) {
				if ( $partner_ref->[0][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[0][3]; } # replace geneA to NULL
					$$tag_ref = 0;
				} elsif ( $partner_ref->[1][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[1][3]; } # replace geneB to NULL
					$$tag_ref = 0;
				} else {
					print "${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				print "${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		} else {	
                	if ( $multiple_ref->[0][0] eq $partner_ref->[0][0] ) { # the chromosome of mapping read == that of geneA
                        	if ( $partner_ref->[0][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[0][3]; } # replace geneA to NULL
					$$tag_ref = 0;
                        	} else {
                                	print "${read_type}-filtered: $id_ref ($$read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[0][3]\n";
                                        $$tag_ref = 1;
                                }
                        } elsif ( $multiple_ref->[0][0] eq $partner_ref->[1][0] ) { #the chromosome of mapping read == that of geneB
                        	if ( $partner_ref->[1][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[1][3]; } # replace geneB to NULL
					$$tag_ref = 0;
                                } else {
                                	print "${read_type}-filtered: $id_ref ($$read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[1][3]\n";
                                        $$tag_ref = 1;
                                }
                        } else {
                        	print "${read_type}-filtered: $id_ref ($$read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
                                $$tag_ref = 1;
                        }
		}
        }
1;
