#!/usr/bin/perl

### Info
# This script takes a tab-delimited file of crispr sequences and primers and creates a genbank file of the WT region
# Used for exon deletions containing 4 crisprs and 3 genotyping primers but can be edited for other combinations
# Currently set up for mice but can be adapted for other species by editing the script
# Crisprs are currently all on the plus strand (as that's how they are in iMits but the strand for primers are auto-detected accordingly for the genbank file

### file format:
# gene_id	chr	region_loc_start	region_loc_end	crispr1_seq	crispr2_seq	crispr3_seq	crispr4_seq	primer1_seq	primer2_seq	primer3_seq

### Usage
# perl crispr_genbank_creator.pl infile_name

### Modules
use Bio::EnsEMBL::Registry; # I had to install the ensembl API to my home dir as gen1 doesn't seem to have it
use Bio::SeqIO;
use File::Slurp;

$infile = $ARGV[0];

my $registry = 'Bio::EnsEMBL::Registry';
#$registry->load_registry_from_db( -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org' -user => 'anonymous');
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);
my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' ); ## Get mouse sequence object


### get coordinates from file
my @array = read_file($infile);
###

foreach $name (@array) { # loop through each file line and create a genbank file
   #chomp($name);
   $name =~ s/\r|\n//g;
   @split = split(/\t/, $name);
   $gene = $split[0];
   $file = ">".$gene.".gb"; ## get gene name
   $chr =  $split[1];
   $start = $split[2] - 1000; ## to get 5' flanking sequence from ensembl
   $stop = $split[3] + 1000; ## to get 3' flanking sequence from ensembl
   #print OUT "$split[0] $split[3]       $split[5]       $chr    $start  $stop   ";
   
   $out = Bio::SeqIO->new(-file   => $file, -format => 'GenBank');

   # get handle type - in this case a genome region
   #my @crispr_array = ("9","88407119","88407661");
   $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr,  $start, $stop); # get sequence from ensembl.org

   my $sequence = $slice->seq();
   $seq_obj = Bio::Seq->new(-seq => "$sequence", -display_id => "$gene", -desc => "CRISPR allele design",-alphabet => "dna" );

   # create crispr features
   $c=1;
   for ($n = 4; $n<8; $n++){ 
      ###  get locations of sequences
       $loc_start = index ($sequence, $split[$n])+1;
       $loc_length = length ($split[$n]); 
       $loc_end = $loc_start + $loc_length -1;
       $feat = new Bio::SeqFeature::Generic(-start       => $loc_start,
                                        -end         => $loc_end,
                                        -primary_tag => 'primer_bind',
                                        -tag => {note     => $c } );
      $seq_obj->add_SeqFeature($feat);      
      $c++;
   }

   # create primer features
   $c=0;   
   $lock_start = 0;
   @primer_array =("DF1","DR1","ER1");
   @primer_strand_array = ("1","-1","-1");
   
   $loc_length = 0;
   for ($n = 8; $n<11; $n++){ 
      ###  get locations of sequences
      $loc_start = index ($sequence, $split[$n]) +1;
      print "loc_start =  $loc_start ";
      ## not found? Try revcomp instead
      if ($loc_start == 0){  
          $split[$n] = reverse($split[$n]);
          $split[$n] =~ tr/ACGTacgt/TGCAtgca/;
          $loc_start = index ($sequence, $split[$n]) +1;
      }
      $loc_length = length($split[$n]); 
      $loc_end = $loc_start + $loc_length -1;
      print " $split[$n], $loc_start, $loc_end, $loc_length \n";
      $feat = new Bio::SeqFeature::Generic(-start       => $loc_start,
                                        -end         => $loc_end,
                                        -strand      => $primer_strand_array[$c],
                                        -primary_tag => 'PCR primer',
                                        -tag => {note     => $primer_array[$c] } );  
      $seq_obj->add_SeqFeature($feat);
      $c++;
   }

   ### write out the final file
   $out->write_seq($seq_obj);
}
