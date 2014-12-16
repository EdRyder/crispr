#!/usr/bin/perl

### The following script takes an output file from Cas-off finder  http://www.rgenome.net/cas-offinder/
### Parses the file, uses the ensembl perl API to grab a slice of sequence flanking each target site
### Uses bioperl::Run::Primer3 to design the primers based on the sequence and input parameters	
### Grabs the output of the top primer pair 
### Outputs it all to a tab file cripsr.out

### USAGE perl crispr-offtarget-primers2.pl INPUT_FILE OUTPUT_FILE

### Modules
use Bio::EnsEMBL::Registry; 
use File::Slurp;
use Bio::Tools::Primer3;
use Bio::Tools::Run::Primer3; # This doesn't seem to be part of the basic bioperl package
### WARNING - the default Run::Primer3.pm has some variables that don't work with the primer3_core version on Sanger and so the module needs editing. 
### SEQUENCE is now SEQUENCE_TEMPLATE and TARGET is now SEQUENCE_TARGET for example.

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

# get handle type - in this case a genome region
my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );

#read file into an array

$filename = $ARGV[0];
my @array = read_file($filename);
open (OUT, ">$ARGV[1]");
print OUT "Crispr	Off-target	Mismatches	Chr	Flanking Start	Flanking Stop	Left Primer	Right Primer	Left Tm	Right Tm	Size(bp)	Sequence\n";

foreach $name (@array) {
   chomp($name);
   @split = split(/\t/, $name);
   $start = $split[2] - 400; ## to get 5' flanking sequence from ensembl
   $stop = $split[2] + 420; ## to get 3' flanking sequence from ensembl
   $chr = $split[1];
   $chr =~ s/chr//; ## remove chr text from off-target output chromosome identifier or it breaks ensembl lookup
   print OUT "$split[0]	$split[3]	$split[5]	$chr	$start	$stop	";
   $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $stop); # get sequence from ensembl.org
   my $sequence = $slice->seq();

   ### PRIMER 3 PROGRAM
   my $primer3 = Bio::Tools::Run::Primer3->new( -outfile => "./temp.out",-path => "/usr/bin/primer3_core");
   $primer3->add_targets('SEQUENCE_TEMPLATE'=>$sequence);
   $primer3->add_targets('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>65);
   $primer3->add_targets('SEQUENCE_TARGET'=>"350,120", 'PRIMER_PRODUCT_SIZE_RANGE'=>"300-500");
   $results = $primer3->run;
   my $result1 = $results->primer_results(0);
   @resultkeys = ("PRIMER_LEFT_SEQUENCE","PRIMER_RIGHT_SEQUENCE","PRIMER_LEFT_TM","PRIMER_RIGHT_TM","PRIMER_PAIR_PRODUCT_SIZE");
   foreach $key (@resultkeys){  
      print OUT "${$result1}{$key}\t";
   }
   print OUT "$sequence\n";

   #### NEXT SEQUENCE

 }
close OUT;
