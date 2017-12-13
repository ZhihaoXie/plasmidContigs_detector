#!/usr/bin/perl -w

# Author:    Zhihao Xie  \(#^o^)/
# Date:      2017.08.03
# Version:   v1.0.0
# CopyRight: Copyright ©Zhihao Xie, All rights reserved.

use strict;
use Bio::SeqIO;

if (@ARGV < 3) {
    print STDERR "usage: perl $0 <glimmer_predict> <genome> <out_prefix> <codon table[int]>\n";
    print STDERR "\t<codon table> default is 11, supports [1, 11, 4, 25]\n";
    exit 1;
}
my ($predict_file, $genome, $out_tag, $gcode) = @ARGV;
unless ($gcode) {
    $gcode = 11;
}
my %hash;
my $fasta = Bio::SeqIO->new(-file=>$genome,-format=>'fasta');
while (my $seq = $fasta->next_seq()) {
    my $id = $seq->id;
    $hash{$id} = $seq;
}
$fasta->close();

# difined protein and gene output
my $faa_obj = Bio::SeqIO->new(-file=>">$out_tag.protein.faa",-format=>'fasta');
my $fna_obj = Bio::SeqIO->new(-file=>">$out_tag.gene.fnn",-format=>'fasta');
my $seqid;
my $seq;
open FILE, $predict_file || die $!;
open OUT, "> $out_tag.gene.gff";
print OUT "##gff-version  3\n";
while (<FILE>) {
    chomp;
    if (/^>(\S+)/) {
        $seqid = $1;
        $seq = $hash{$seqid};
        next;
    }
    my @a = split;
    my $gene_id = "$seqid\_$a[0]";
    if ($a[1]>$a[2]) {
        print OUT "$seqid\tGlimmer3\tgene\t$a[2]\t$a[1]\t.\t-\t.\tlocus_tag=$gene_id\n";
        my $str = $seq->subseq($a[2],$a[1]);
        my $new_seq = Bio::Seq->new(-seq=>$str,-id=>"$gene_id");
        $new_seq = $new_seq->revcom;
        my $faa = $new_seq->translate(-codontable_id=>$gcode,-complete=>1);
        my $faa_seq = $faa->seq;
        print OUT "$seqid\tGlimmer3\tCDS\t$a[2]\t$a[1]\t.\t-\t.\tlocus_tag=$gene_id;translation=$faa_seq\n";
        $faa_obj->write_seq($faa);
        $fna_obj->write_seq($new_seq);
    }
    else {
        print OUT "$seqid\tGlimmer3\tgene\t$a[1]\t$a[2]\t.\t+\t.\tlocus_tag=$gene_id\n";
        my $str = $seq->subseq($a[1],$a[2]);
        my $new_seq = Bio::Seq->new(-seq=>$str,-id=>"$gene_id");
        my $faa = $new_seq->translate(-codontable_id=>$gcode,-complete=>1);
        my $faa_seq = $faa->seq;
        print OUT "$seqid\tGlimmer3\tCDS\t$a[1]\t$a[2]\t.\t+\t.\tlocus_tag=$gene_id;translation=$faa_seq\n";
        $faa_obj->write_seq($faa);
        $fna_obj->write_seq($new_seq);
    }
}
close FILE;
close OUT;

