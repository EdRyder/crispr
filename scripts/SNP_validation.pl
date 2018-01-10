#!/usr/bin/perl

### The following script takes an input file of SNP locations and designs PCR primers for validation (capillary or MiSeq).
### Parses the file, uses the ensembl perl API to grab a slice of sequence flanking each target site
### Uses bioperl::Run::Primer3 to design the primers based on the sequence and imput parameters	
### Outputs it all to a tab file

### USAGE perl SNP_validation.pl input_file output_file
### File format (not all are needed for this script, this was just output file I got): 
### Variant_no	chrom	pos	ref	alt	indel	variant_allele_fraction	alternat_allele_detph	sample	TDNV_genotype_quality	MSAT_or_HOMP	Comments	Type
### eg
### Variant_no	chrom	pos	ref	alt	indel	variant allele fraction	alternat allele detph	sample	TDNV genotype quality	MSAT or HOMP	Comments	Type
###	1	5	16828975	A	C	0	0.155555556	7	sample1	10.12	ok		SNP


### Modules - needs both the Ensembl API and Bioperl, which can be downloaded from  http://www.ensembl.org/info/docs/api/api_installation.html
use Bio::EnsEMBL::Registry; 
use File::Slurp;
use Bio::Tools::Primer3;
use Bio::Tools::Run::Primer3; # I had to install this into my home dir as gen1 doesn't seem to have it
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
print OUT "ID	Mouse_ID	Chr	Start	Ref	Var	Flanking Start	Flanking Stop	Left Primer Name	Right Primer Name	Left Primer Sequence	Right Primer Sequence	Left Tm	Right Tm	Size(bp)	Sequence\n";

foreach $name (@array) {
   chomp($name);
   @split = split(/\t/, $name);
   $start = $split[2] - 400; ## to get 5' flanking sequence from ensembl
   $stop = $split[2] + 420; ## to get 3' flanking sequence from ensembl
   $chr = $split[1];
   $type =  $split[0]."_".$split[8]."_".$split[12];
   $chr =~ s/chr//; ## remove chr text from off-target output chromosome identifier or it breaks ensembl lookup
   print OUT "$split[0]	$split[8]	$split[1]	$split[2]	$split[3]	$split[4]	$start	$stop	".$type."_F	".$type."_R	";
   $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $stop); # get sequence from ensembl.org
   #$slice = $slice->get_repeatmasked_seq(); ## get hard repeatmasked sequence
   $slice = $slice->get_repeatmasked_seq(undef, 1); ## get softmarked sequence 
   my $sequence = $slice->seq();

   ### PRIMER 3 PROGRAM
   my $primer3 = Bio::Tools::Run::Primer3->new( -outfile => "./temp.out",-path => "/usr/bin/primer3_core");
   $primer3->add_targets('SEQUENCE_TEMPLATE'=>$sequence);
   $primer3->add_targets('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>65);
   $primer3->add_targets('SEQUENCE_TARGET'=>"330,140", 'PRIMER_PRODUCT_SIZE_RANGE'=>"250-400");
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
