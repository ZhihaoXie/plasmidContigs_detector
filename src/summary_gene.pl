#!/usr/bin/perl -w
use Bio::SeqIO;
use File::Basename;

if (@ARGV < 3) {
    die "Usage:\n perl $0 genome_seq gene_seq output_file\n";
}
my $genome_seq = shift;
my $gene_seq = shift;
my $output = shift;
my $basename = basename $genome_seq;
$basename =~ s/\.genome\.fasta$//;

# genome
my $all_length;
my $in_obj = Bio::SeqIO->new(-file=>"$genome_seq", -format=>'fasta');
while (my $seq = $in_obj->next_seq()) {
    my $id = $seq->id;
    $all_length += $seq->length;
}

# gene summary
my ($gene_num, $gene_all_length, $gene_avg_length);
my $in_obj2 = Bio::SeqIO->new(-file=>"$gene_seq", -format=>'fasta');
while (my $seq = $in_obj2->next_seq()) {
    $gene_num++;
    $gene_all_length += $seq->length;
}
$in_obj2->close();

if($gene_num > 0) {
    $gene_avg_length = sprintf("%.2f", $gene_all_length/$gene_num);
} else {
    $gene_num = 0;
    $gene_avg_length = 0;
    $gene_all_length = 0;
}
my $rate = sprintf("%.2f", $gene_all_length/$all_length*100);
my $density = sprintf("%.2f", $gene_num/$all_length*1000);

open FO,"> $output";
print FO "Sample\tGenome size (bp)\tGene number\tGene total length (bp)\tGene average length (bp)\tGene length / Genome (%)\tGene density(per Kb)\n";
print FO "$basename\t$all_length\t$gene_num\t$gene_all_length\t$gene_avg_length\t$rate\t$density\n";
close FO;
