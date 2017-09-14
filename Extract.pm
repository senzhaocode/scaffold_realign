package Extract; 
use strict;
use warnings;

	sub Spanning {
		my ($dis_ref, $geneA, $geneB, $sam_file) = @_;

		# select mapped spanning reads (one mapping to GeneA/GeneB; the other mapping to GeneB/GeneA)
		open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7!=\"=\" && \$5>=40 && \$8>0)' $sam_file/tmp/hisats_noclip.sam | ") || die "cannot run awk discordant script 1:$!\n";
		while ( <IN> ) {
			chomp $_; #
			my ($name, $one_match, $one_pos, $other_match, $other_pos, $seq) = (split /\t/, $_)[0,2,3,6,7,9];
			if ( $one_match eq $geneA && $other_match eq $geneB ) {
				push @{$dis_ref->{$name}}, [$one_pos, "spanning", $one_match, $seq];
			} elsif ( $one_match eq $geneB && $other_match eq $geneA ) {
				push @{$dis_ref->{$name}}, [$one_pos, "spanning", $one_match, $seq];
			}
		}
		close IN;
		
		foreach my $name ( keys %{$dis_ref} ) {
			if ( scalar(@{$dis_ref->{$name}}) != 2 ) {
				print "Step2-0: $name does not show pairend spanning read (2)\n";
				delete($dis_ref->{$name});
			}
		}
	}
		
	sub Discordant {
		my ($dis_ref, $archor, $breakpoint_scaffold, $read_length, $sam_file) = @_;

		# select mapped discordant split reads (one mapping to GeneA/GeneB; the other mapping to scaffold)
		#*****/ 1. mapping quality >=40; ($3!=$7; $7!="=" -- discordant: one mapped to scaffold; the other mapped to GeneA/GeneB) /
		open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7!=\"=\" && \$5>=40 && \$8>0)' $sam_file/tmp/hisats_noclip.sam | grep 'caffold' |") || die "cannot run awk discordant script 1:$!\n";
		while ( <IN> ) {
			chomp $_; #
			my ($name, $one_match, $one_pos, $other_match, $other_pos, $seq) = (split /\t/, $_)[0,2,3,6,7,9];
			my $string = ""; # record mismatching information
			if ( $_ =~/MD\:Z\:([\w\^]+)/ ) { 
				$string = $1;
			} else {
				print "mismatch information string contains unexpected symbol.\n";
			}
			$string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

			if ( $one_match eq "scaffold" ) {
				if ( ($one_pos+$archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($one_pos+$read_length-$archor) ) { # read cross breakpoint
					if ( $num_mis < 6 ) { # the default setting -- mismatching number
						push @{$dis_ref->{$name}}, [$one_pos, "discordant_split", $one_match, $seq]; # read across breakpoint
					} else {
						push @{$dis_ref->{$name}}, [$one_pos, "mismatch_fail", $one_match, $seq]; # read with many mismatch base
					}
				} else { # read does not cross breakpoint
					push @{$dis_ref->{$name}}, [$one_pos, "anchor_fail", $one_match, $seq]; # read outside breakpoint
				}
			} elsif ( $other_match eq "scaffold" ) {
				if ( ($other_pos+$archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($other_pos+$read_length-$archor) ) { # read cross breakpoint
					if ( $num_mis < 6 ) { # the default setting -- mismatching number
						push @{$dis_ref->{$name}}, [$other_pos, "discordant_split", $one_match, $seq]; # read across breakpoint
					} else {
						push @{$dis_ref->{$name}}, [$one_pos, "mismatch_fail", $one_match, $seq]; # read with many mismatch base
					}
				} else { # read does not cross breakpoint
					push @{$dis_ref->{$name}}, [$other_pos, "anchor_fail", $one_match, $seq]; #read outside breakpoint
				}
			} else {
				print "Step2-1: No hit match to scaffold for discoardant split reads mapping (one mapping to GeneA/GeneB; the other mapping to scaffold)\n";
			}
		}
		close IN;

		#*****/2. mapping quality >=40; ($7=="=" -- discordant: both read mapped to scaffold)
		open (IN, "awk  -F '\t' -v OFS='\t' '(\$7==\"=\" && \$8>0 && \$5>=40 && \$9!=0)' $sam_file/tmp/hisats_noclip.sam | grep 'caffold' |") || die "cannot run awk discordant script 2:$!\n";
		while ( <IN> ) {
			chomp $_; #
			my ($name, $one_match, $one_pos, $other_match, $other_pos, $seq) = (split /\t/, $_)[0,2,3,6,7,9];
			my $string = ""; # record mismatching information
			if ( $_ =~/MD\:Z\:([\w\^]+)/ ) {
				$string = $1;
			} else {
				print "mismatch information string contains unexpected symbol.\n";
			}
			$string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

			if ( $one_match eq "scaffold" ) {
				if ( ($one_pos+$archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($one_pos+$read_length-$archor) ) {
					if ( $num_mis < 6 ) { # the default setting -- mismatching number
						push @{$dis_ref->{$name}}, [$one_pos, "discordant_split", $one_match, $seq]; # read across breakpoint
					} else {
						push @{$dis_ref->{$name}}, [$one_pos, "mismatch_fail", $one_match, $seq]; # read with many mismatch base
					}
				} else { # read does not cross breakpoint
					push @{$dis_ref->{$name}}, [$one_pos, "anchor_fail", $one_match, $seq]; # read outside breakpoint
				}
			} else {
				print "Step2-2: No hit match to scaffold for discordant split reads mapping (both read mapped to scaffold)\n";
			}
		}
		close IN;

		foreach my $name ( keys %{$dis_ref} ) {
                	if ( scalar(@{$dis_ref->{$name}} == 2) ) {
                        	if ( $dis_ref->{$name}[0][1] eq "discordant_split" and $dis_ref->{$name}[1][1] eq "discordant_split" ) {
                                	#print "pairwise: $name\t$dis_ref->{$name}[0][0]\t$dis_ref->{$name}[0][1]\t$dis_ref->{$name}[1][0]\t$dis_ref->{$name}[1][1]\n";
                        	} else {
                                	print "Step2: $name does not meet criteria of discordant split read mapping (2)\n";
                                	delete($dis_ref->{$name});
                        	}
                	} else {
                        	if ( $dis_ref->{$name}[0][1] eq "discordant_split" ) {
                                	#print "pairwise: $name\t$dis_ref->{$name}[0][0]\t$dis_ref->{$name}[0][1]\n";
                        	} else {
                                	print "Step2: $name does not meet criteria of discordant split read mapping (1)\n";
                                	delete($dis_ref->{$name});
                       		}
                	}
        	}
	}

	sub Singlton {
		my ($dis_ref, $archor, $breakpoint_scaffold, $read_length, $sam_file) = @_;
		
		# select mapped singlton split read (one mapping to scaffold; the other no mapping hit)
		#*****/ mapping quality >= 40; ($9=0 -- singlton mapped; $2=73|89|137|153 -- singlton mapping feature filtering) /
		open (IN, "awk -F '\t' -v OFS='\t' '(\$7==\"=\" && \$3==\"scaffold\" && \$5>=40 && \$8>0 && \$9==0 && (\$2==73 || \$2==89 || \$2==137 || \$2==153))' $sam_file/tmp/hisats_noclip.sam |") || die "cannot run awk singlton script:$!\n";
        	while ( <IN> ) {
                	chomp $_; #
                	my ($name, $one_match, $one_pos, $other_match, $other_pos, $seq) = (split /\t/, $_)[0,2,3,6,7,9];
			my $string = ""; # record mismatching information
                        if ( $_ =~/MD\:Z\:([\w\^]+)/ ) {
                                $string = $1;
                        } else {
                                print "mismatch information string contains unexpected symbol.\n";
                        }
                        $string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

                	if ( $one_match eq "scaffold" ) {
                        	if ( ($one_pos+$archor+2) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($one_pos+$read_length-$archor-2) ) { # read cross breakpoint: get extension (2, -2)
					if ( $num_mis < 6 ) {
                                		push @{$dis_ref->{$name}}, [$one_pos, "singlton_split", $one_match, $seq];
                                		#print "singlton: $name\t$singlton{$name}[0][0]\t$singlton{$name}[0][1]\n";
                                	} else {
						print "Step2-3: $name does not meet criteria of singlton split read mapping (1)\n";
					}
                        	} else {
                                	print "Step2-3: $name does not meet criteria of singlton split read mapping (1)\n";
                        	}
                	} else {
                        	print "Step2-3: No hit match to scaffold for singlton split read mapping (one to scaffold; the other unmapped)\n";
                	}
        	}
        	close IN;
	}
1;
