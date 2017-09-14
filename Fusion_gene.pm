package Fusion_gene;

use strict;
use warnings;

	sub cDNA {
		my ($ref, $path, $a, $b) = @_;
		my $tag_a = 0; my $tag_b = 0;
		open (IN, "cut -s -f1,3,4,5,9,10 $path | uniq |") || die "cannot open this path for cDNA path:$!\n";
		while ( <IN> ) {
			chomp $_; my ($gene_en, $start, $end, $symbol, $type, $chr) = (split /\t/, $_)[0,1,2,3,4,5];

			if ( $a eq $gene_en ) {
				push @{$ref}, [$chr, $start, $end, $a]; $tag_a = 1;
			} elsif ( $a eq $symbol ) {
				push @{$ref}, [$chr, $start, $end, $a]; $tag_a = 1;
			}

			if ( $b eq $gene_en ) {
				push @{$ref}, [$chr, $start, $end, $b]; $tag_b = 1;
			} elsif ( $b eq $symbol ) {
				push @{$ref}, [$chr, $start, $end, $b]; $tag_b = 1;
			}
		}
		close IN;
		if ( $tag_a == 1 && $tag_b == 1 ) {
		} else {
			print "Gene name does not match to annotation\n"; exit;
		}
	}
1;
