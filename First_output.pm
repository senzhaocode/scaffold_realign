package First_output;
use strict;
use warnings;

	sub grep_spanning {
		my ($dis_ref, $fastq_1, $fastq_2, $path) = @_;
		open (OUT, ">>$path/read_mapped_info") || die "cannot open summary output path:$!\n";
		open (OUT1, ">$path/spanning_1.txt") || die "cannot open first round spanning_1 output path:$!\n";
		open (OUT2, ">$path/spanning_2.txt") || die "cannot open second round spanning_2 output path:$!\n";
		foreach my $id ( keys %{$dis_ref} ) {
			my $header = "@".$id." ";
			my $hit_1 = `grep "$header" -A 3 $fastq_1`; chomp $hit_1;
			my $hit_2 = `grep "$header" -A 3 $fastq_2`; chomp $hit_2;
			my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
			my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
                	if ( scalar(@{$dis_ref->{$id}}) == 2 ) {
                        	if ( $dis_ref->{$id}[0][3] eq $seq_1 ) { #first_end
                                	if ( $dis_ref->{$id}[1][3] eq $seq_2 ) { #second_end
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2_tr ) { #second(reverse_complement)
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->$id}[1][2](mapping orient wrong)\n";
                                	} 
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_1_tr ) { #first_end(reverse_complement)
                                	if ( $dis_ref->{$id}[1][3] eq $seq_2 ) { #second_end
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2_tr ) {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](mapping orient wrong)\n";
                                	} 
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2 ) { #second_end
                                	if ( $dis_ref->{$id}[1][3] eq $seq_1 ) { #first_end
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_1_tr ) { #first_end(reverse_complement)
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](mapping orient wrong)\t$dis_ref->{$id}[0][2]\n";
                                	} 
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2_tr ) { #second_end(reverse_complement)
                                	if ( $dis_ref->{$id}[1][3] eq $seq_1 ) { #first_end
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_1_tr ) { #first_end(reverse_complement)
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](mapping orient wrong)\t$dis_ref->{$id}[0][2](reverse_complement)\n";
                                	}
                        	} else {
                                	if ( $dis_ref->{$id}[1][3] eq $seq_1 ) {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](mapping orient wrong)\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_1_tr ) {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](mapping orient wrong)\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2 ) {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](mapping orient wrong)\t$dis_ref->{$id}[1][2]\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2_tr ) {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](mapping orient wrong)\t$dis_ref->{$id}[1][2](reverse_complement)\n";
                                	} else {
                                        	print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](mapping orient wrong)\t$dis_ref->{$id}[1][2](mapping orient wrong)\n";
                                	}
                        	}
			}
		}
		close OUT1;
		close OUT2;
		close OUT;
	}

	sub grep_discordant {
		my ($dis_ref, $fastq_1, $fastq_2, $path) = @_;
		open (OUT, ">$path/read_mapped_info") || die "cannot open summary output path:$!\n";
        	open (OUT1, ">$path/discordant_split_1.txt") || die "cannot open first round discordant_1 output path:$!\n";
        	open (OUT2, ">$path/discordant_split_2.txt") || die "cannot open second round discordant_2 output path:$!\n";
		
		foreach my $id ( keys %{$dis_ref} ) {
                	my $header = "@".$id." ";
                	my $hit_1 = `grep "$header" -A 3 $fastq_1`; chomp $hit_1;
                	my $hit_2 = `grep "$header" -A 3 $fastq_2`; chomp $hit_2;
                	my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
                	my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
                	if ( scalar(@{$dis_ref->{$id}}) == 2 ) {
                        	if ( $dis_ref->{$id}[0][3] eq $seq_1 ) { #first_end
                                	if ( $dis_ref->{$id}[1][3] eq $seq_2 ) { #second_end
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2_tr ) { #second(reverse_complement)
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->$id}[1][2](mapping orient wrong)\n";
                                	} 
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_1_tr ) { #first_end(reverse_complement)
                                	if ( $dis_ref->{$id}[1][3] eq $seq_2 ) { #second_end
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2_tr ) {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](mapping orient wrong)\n";
                                	} 
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2 ) { #second_end
                                	if ( $dis_ref->{$id}[1][3] eq $seq_1 ) { #first_end
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_1_tr ) { #first_end(reverse_complement)
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](mapping orient wrong)\t$dis_ref->{$id}[0][2]\n";
                                	} 
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2_tr ) { #second_end(reverse_complement)
                                	if ( $dis_ref->{$id}[1][3] eq $seq_1 ) { #first_end
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_1_tr ) { #first_end(reverse_complement)
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                                	} else {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](mapping orient wrong)\t$dis_ref->{$id}[0][2](reverse_complement)\n";
                                	}
                        	} else {
                                	if ( $dis_ref->{$id}[1][3] eq $seq_1 ) {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](mapping orient wrong)\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_1_tr ) {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](mapping orient wrong)\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2 ) {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](mapping orient wrong)\t$dis_ref->{$id}[1][2]\n";
                                	} elsif ( $dis_ref->{$id}[1][3] eq $seq_2_tr ) {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](mapping orient wrong)\t$dis_ref->{$id}[1][2](reverse_complement)\n";
                                	} else {
                                        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](mapping orient wrong)\t$dis_ref->{$id}[1][2](mapping orient wrong)\n";
                                	}
                        	}
			} elsif ( scalar(@{$dis_ref->{$id}}) == 1 ) { # Not included 
                        	#if ( $dis_ref->{$id}[0][2] eq "scaffold" ) {
                                #	if ( $dis_ref->{$id}[0][3] eq $seq_1 ) {
                                #        	print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\tGeneA/GeneB\n";
                                #	} elsif ( $dis_ref->{$id}[0][3] eq $seq_1_tr ) {
                                #        	print OUT "discordant_split\t$id\t$dis_ref->$id}[0][2](reverse_complement)\tGeneA/GeneB\n";
                                #	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2 ) {
                                #        	print OUT "discordant_split\t$id\tGeneA/GeneB\t$dis_ref->{$id}[0][2]\n";
                                #	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2_tr ) {
                                #        	print OUT "ddiscordant_split\t$id\tGeneA/GeneB\t$dis_ref->{$id}[0][2](reverse_complement)\n";
                                #	} else {
                                #        	print OUT "discordant_split\t$id\tmapping orient wrong\n";
                                #	}
                                #	print OUT1 "$hit_1\n";
                                #	print OUT2 "$hit_2\n";
                        	#} else {
                                #	print OUT "$id pairend read mapping wrong discordant read\n";
                        	#}
                	}
        	}
        	close OUT1;
        	close OUT2;
		close OUT;
	}
	
	sub grep_singlton {
		my ($dis_ref, $fastq_1, $fastq_2, $path) = @_;
		open (OUT, ">>$path/read_mapped_info") || die "cannot open summary output path:$!\n";
		open (OUT1, ">$path/singlton_split_1.txt") || die "cannot open first round singlton_1 output path:$!\n";
		open (OUT2, ">$path/singlton_split_2.txt") || die "cannot open first round singlton_2 output path:$!\n";
	
		foreach my $id ( keys %{$dis_ref} ) {
                	my $header = "@".$id." ";
                	my $hit_1 = `grep "$header" -A 3 $fastq_1`;
                	my $hit_2 = `grep "$header" -A 3 $fastq_2`;
                	my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
                	my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
                	if ( scalar(@{$dis_ref->{$id}}) == 1 ) {
                        	if ( $dis_ref->{$id}[0][3] eq $seq_1 ) {
                                	print OUT "singlton_split\t$id\t$dis_ref->{$id}[0][2]\tNULL\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_1_tr ) {
                                	print OUT "singlton_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\tNULL\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2 ) {
                                	print OUT "singlton_split\t$id\tNULL\t$dis_ref->{$id}[0][2]\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                        	} elsif ( $dis_ref->{$id}[0][3] eq $seq_2_tr ) {
                                	print OUT "singlton_split\t$id\tNULL\t$dis_ref->{$id}[0][2](reverse_complement)\n"; print OUT1 "$hit_1\n"; print OUT2 "$hit_2\n";
                        	} else {
                                	print OUT "singlton_split\t$id\tmapping orient wrong\n";
                        	}
                	}
        	}
        	close OUT1;
        	close OUT2;
        	close OUT;
	}
1;
