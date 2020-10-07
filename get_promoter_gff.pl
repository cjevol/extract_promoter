#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use constant UPL => 2000;

my $ref=shift;
my $gff=shift;
my $list=shift;
# my $output=shift;

my %list=();
my %genes=();
open(IN, $list) or die "Cannot open $list: $!";
# $_=<IN>;
my @content=<IN>;
chomp @content;
close(IN);

foreach (@content){
	my $gene_name=(split/\s+/)[0];
	# print STDERR "$gene_name\n";
	$genes{$gene_name}=1 if($gene_name!~/^\s+$/);
}
print STDERR scalar keys %genes, " genes have been selected\n";

open(GFF, $gff) or die "Cannot open $gff: $!";
my %complete=();
while(<GFF>){
	# chomp;
	next if(/^#/);
	my ($chr, $type, $start, $end, $strand, $feat) = (split/\s+/)[0,2,3,4,6,8];
	next unless($type eq 'gene');
	my ($gene_name)=$feat=~/Name=(\w+);/;
	# print STDERR "$gene_name, ";
	if(exists $genes{$gene_name}){
		# print STDERR "$gene_name\n";
		$list{$chr}->{$gene_name}=[$start, $end, $strand];
		$complete{$gene_name}=1;
	}
}
close(GFF);
print STDERR "coords of ", scalar keys %complete, " genes have been found\n";
foreach (keys %complete){
	delete $genes{$_};
}
print STDERR "Cannot find info for genes: ", join(", ", keys %genes),"\n" if(scalar keys %genes >0);


my $in=Bio::SeqIO->new(-file =>$ref, -format=>'fasta');
my $out=Bio::SeqIO->new(-fh => \*STDOUT, -format=>'fasta');

while(my $seq=$in->next_seq){
	my $chr=$seq->id;
	my $length=$seq->length;
	if(exists $list{$chr}){
		foreach my $gene_name (keys %{$list{$chr}}){
			my ($start, $end, $strand)=@{$list{$chr}->{$gene_name}};

			if($strand eq '+'){
				my $new_start=$start - UPL;
				$new_start=1 if($new_start<1);
				my $subseq=$seq->subseq($new_start, $end);
				my $subseqObj=Bio::Seq->new(-id=>$gene_name, -seq=>$subseq);
				$out->write_seq($subseqObj);
			}else{
				my $new_end=$end + UPL;
				$new_end=$length if($new_end>$length);
				my $subseq=$seq->subseq($start, $new_end);
				my $subseqObj=Bio::Seq->new(-id=>$gene_name, -seq=>$subseq);
				$subseqObj=$subseqObj->revcom;
				$out->write_seq($subseqObj);
			}
		}
	}
}








