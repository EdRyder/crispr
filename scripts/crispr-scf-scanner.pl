#!/usr/bin/perl
#
# This script takes a directory of .scf sequence trace files and pattern matches against a list of CRISPRs 
# Useful for off-target validation of sequenced PCRs where you don't want to manually look through them all
# Output file shows the original crispr, the off-target match, strand and scf file
#
# USAGE perl crispr-scf-scanner.pl INPUT_LIST OUTPUT_FILE SCF_DIRECTORY

use Bio::Seq;
use Bio::SeqIO;
use File::Slurp;

## file name for crisprs input file on command line as first argument
$crisprs_in = $ARGV[0];
open (OUT, ">$ARGV[1]");
print OUT "Crispr Sequence	Off-target sequence	Strand	File\n";
$dir = $ARGV[2];

## get input file of crisprs and read into an array
@crisprs = read_file($crisprs_in, chomp => 1);

### used for testing
#$input = "TATTATTTGTGTCAAGTTGACAAGGGG";
#$input_rev = reverse($input);
#$input_rev =~ tr/ACGTacgt/TGCAtgca/;
####


## get directory of scf files 
#$dir = ".";
my @files = grep ( -f ,<$dir/*.scf>);
## Process each file and look for matches in the sequence to the list of crisprs 
foreach $file (@files) {
    $r = 0;
    print "processing $file\n";
 
    $testseq = Bio::SeqIO->new ('-format' => 'SCF',
                            '-file' => $file
                           );

    while ($s = $testseq->next_seq() ){ ## add in loop to check all crisprs in input file
       $sequence = $s->seq();
       #print ( "sequence is ", $sequence, "\n");
       foreach $crispr (@crisprs){
           @crispr_split = split(/\t/, $crispr);
           $crispr_off = uc($crispr_split[1]);
           $crispr_off_lc = $crispr_split[1];
           $crispr_rev = reverse($crispr_off);
           $crispr_rev =~ tr/ACGTacgt/TGCAtgca/;
           $original = $crispr_split[0];
           #print "line = $crispr rev = $crispr_rev off = $crispr_off \n";
           if (index($sequence, $crispr_off) != -1) {
              print OUT "$original	$crispr_off_lc	plus	$file\n";
	      $r++;   
           }
           if (index($sequence, $crispr_rev) != -1) {
              print OUT "$original	$crispr_off_lc	minus	$file\n";
              $r++;
           }
        }
    }       
    if ($r == 0){
       print  OUT "no match	no match	no match	$file\n";;
    }
} # end of each file
close OUT;
