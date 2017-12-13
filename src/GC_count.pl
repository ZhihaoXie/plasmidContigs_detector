#!/usr/bin/perl -w
use Bio::SeqIO;
my $Usage="\n\tUsage:\n\t\tperl $0 <fasta_file>\n\n";
if(@ARGV<1){
	print STDERR $Usage;
	exit;
}
my $input =shift;
my $obj = Bio::SeqIO->new(-file=>$input,-format=>'fasta');
print "$input\n";
while(my $seq = $obj->next_seq){
	my $id = $seq->id;
	my $str = $seq->seq;
	my $n=$str=~tr/GCgc/GCgc/;
	my $length = $seq->length;
	my $pr = $n/$length*100;
	$len_n += $n;
	$len_g += $length;

	print "$id:\t";
	printf "%0.2f%%\n",$pr;
}
$pr = $len_n/$len_g*100;
print "$pr\n";
