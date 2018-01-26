#!/usr/bin/perl

### Gene Structure notes
# Takes a list of genes and crispr locations and outputs a file with the consequences to the gene
# Works at the exon level purely on locations. Does not calculate whether a stop codon is introduced from an indel, for example, as that requires sequence and allele input. 
# Meant for gene-specific exon deletions.
# Usage: perl CRISPR_deletion_consequences.pl infile outfile 
# File format (tabbed) = gene_name, start location, end location, source, project
# eg Zfp712  67041242        67042371        Experimental      WTSI_project_A
# Chr is not included as that information will be extracted from the gene object
###

### Modules - needs both the Ensembl API and Bioperl, which can be downloaded from  http://www.ensembl.org/info/docs/api/api_installation.html
use Bio::EnsEMBL::Registry; 
use Bio::SeqIO;
use Bio::Seq;
use List::Util qw(min max);
#use DBI;
#use DBD::mysql;

use File::Slurp; ## this may need installing
use DDP; ##used for debugging. Comment this out if you don't need or have it

$infile = $ARGV[0]; # input file
$outfile = $ARGV[1]; # output file. Two files will be created: .out and .summary


### Connect to the Ensembl databases. Here we will need the gene for transcript info, and a chromosome slice for the sequence and location
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);
my $gene_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'gene' ); ## Get mouse sequence object
my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
###

### create tab delimeted output files 
open (OUT, ">$outfile.out");
open (OUT2, ">$outfile.summary");
print (OUT "Gene	Transcript	Exon_Ensembl	Chr	Exon start	Exon stop	Strand	start_frame	end_frame	exon_rank	crispr_min	crispr_max	State\n");
print (OUT2 "Project	Gene	Transcript_ID	Transcript_name	Transcript_length	Type	Strand	Annotation_source	Chr	Del_start	Del_End	CRISPR_Deletion_size	Exons_affected	Start_frame	End_frame	Transcript_deleted	Percent_transcript_deleted	CDS_size	CDS_deleted	Percent_CDS deleted	Exon_consequence	Flanking_exon_consequence	ATG_deleted	Source	Domains	Sequence	Peptide_coordinates	Peptide_sequence	Peptide_deleted	Transcript_hit	Summary	URL\n");
### get coordinates from file
my @array = read_file($infile);
###

foreach $name (@array) { ## for each gene in the input file: 
   $name =~ s/\r|\n//g;
   @split = split(/\t/, $name);
   $in_gene = $split[0];
   $crispr_min = $split[1];
   $crispr_max = $split[2];
   $crispr_source = $split[3];
   $project_name = $split[4];
   ### hard-coded for testing
   #$in_gene = "Tacc3";
   #$crispr_min = 33657586;
   #$crispr_max = 33659727;
   ###

   print "processing gene $in_gene...\n"; ## all of the print statements can be removed for speed, but they do show the end user that something is happening...
   my $gene = $gene_adaptor->fetch_by_display_label($in_gene); # Get the gene info based on name. If the input list is different (eg MGI IDs) then this will have to be changed
   my @transcripts = @{ $gene->get_all_Transcripts }; ## get all transcripts for the gene and stick them in an array
   #p(@transcripts);
   my @exons = @{ $gene->get_all_Exons };
   print "Gene ", $gene->display_id, "\n"; 
   print scalar(@transcripts), " Transcripts\n";
   print scalar(@exons), " Exons\n";
   print $gene->description, "\n";
   my $exon_count = 0;
   
  foreach my $tr (@transcripts) { # loop through each transcript for this gene 
     # reset a load of variables from last time, so start anew.
     my $url = "";
     my @domain_summary = ();
     my @pep_start_array =();
     my @pep_stop_array =();
     my @cds_start_array =();
     my @cds_stop_array =();
     my $cds_ATG = "";
     my $hit = 0;
     my $summary = "";
     my @crispr_stat_array = ();
     my @exon_id_array = ();
     my $last_exon = "";
     my $domain_summ_string = "";
     my $temp_start = $tr->start();;
     my $temp_stop = 0;
     my $has_domain = "0";
     my $temp_start_neg = $tr->end();
     my $crispr_stat;
     my $crispr_start_frame = "";
     my $crispr_end_frame = "";     
     my $tr_id = $tr->stable_id();
     my $tr_display_id = $tr->display_id();
     my $tr_ens = $tr->external_name();
     my $tr_chr = $tr->seq_region_name();
     my $tr_type = $tr->biotype();
     my $tr_source = $tr->source();
     my $tr_strand = $tr->strand(); 
     my $tr_length = $tr->length();
     my $tr_cds_start = $tr->cdna_coding_start();
     my $tr_cds_end = $tr->cdna_coding_end();
     my $tr_cds_region_start = $tr->coding_region_start();
     my $tr_cds_region_end = $tr->coding_region_end();
     my $tr_cds_length = $tr_cds_end - $tr_cds_start + 1;
     if ($tr_cds_length == 1){ #no CDS
        $tr_cds_length = ""; ## otherwise CDS will report as 1bp and look stupid
     }
     my $deleted_length = 0; # exon material deleted
     my $url_start = $temp_start - 500;
     my $url_end = $temp_start_neg + 500;
     ## create URL
     $url = "http://www.ensembl.org/Mus_musculus/Location/View?r=".$tr_chr.":".$url_start."-".$url_end.";mr=".$tr_chr.":".$crispr_min."-".$crispr_max."";
     #$url = "http://www.ensembl.org/Mus_musculus/Location/View?g=".$gene->display_id.";mr=".$tr_chr.":".$crispr_min."-".$crispr_max."";
     ##


	 ### show the user that something happened.
     print "Processing $tr_display_id $tr_ens\n";
     print "biotype = $tr_type\n";
     print "transcription length = $tr_length, cds length = $tr_cds_length\n";
     my $tr_cds_region_length = $tr_cds_region_end - $tr_cds_region_start +1 ; 

     ### get cds co-ordinates and deleted coding sequences
	 # This used to a shed-load of if/else loops in the exon section, but this is much simpler and less likely to break
     my $cds_coord_start = "";
     my $cds_coord_end = "";
     my $cds_sequence = "";
     my $cds_deleted = "";
     if ($tr->translation()) {
        my @cds_coords = $tr->get_TranscriptMapper->genomic2cds($crispr_min, $crispr_max, $tr_strand);
        foreach my $cds_coord (@cds_coords) { ## loop through the coordinates and put them in an array - needed for multiple exons
           if ($cds_coord->can(id) ){ ## filter out gap objects
               push @cds_start_array,$cds_coord->start();
               push @cds_stop_array,$cds_coord->end();
           }
        }
        #p@cds_coords;
        #p@cds_start_array;
        #p@cds_stop_array;
        $cds_coord_start = $cds_start_array[0];
        $cds_coord_end = $cds_stop_array[-1];
        print "cds coordinates = $cds_coord_start \- $cds_coord_end\n";
        $cds_deleted = $cds_coord_end - $cds_coord_start + 1;
        if ($cds_deleted == 1){
           $cds_deleted = "";
        }
        if ($cds_coord_start == 1){
           $cds_ATG = "yes";
        }
        print "cds deleted  = $cds_deleted\n";
     }
     ###


     ### get protein co-ordinates and deleted peptide sequences
     my $pep_coord_start = "";
     my $pep_coord_end = "";
     my $pep_sequence = "";
     my $pep_deleted = "";
     if ($tr->translation()) {
        my @pep_coords = $tr->genomic2pep($crispr_min, $crispr_max, $tr_strand);
        foreach my $pep_coord (@pep_coords) { ## loop through the coordinates and put them in an array - needed for multiple exons
           if ($pep_coord->can(id) ){ ## filter out gap objects
               push @pep_start_array,$pep_coord->start();
               push @pep_stop_array,$pep_coord->end();
           }
        }
        $pep_coord_start = $pep_start_array[0];
        $pep_coord_end = $pep_stop_array[-1];
        print "Pep coordinates = $pep_coord_start \- $pep_coord_end\n";
        #p@pep_coords;
        $pep_sequence = $tr->translate->seq;
        $pep_deleted = substr($pep_sequence, $pep_coord_start-1, $pep_coord_end - $pep_coord_start + 1);
        if (length $pep_deleted == 1){  
           $pep_deleted = "";
        }
        print "Pep sequence = $pep_sequence\n";
        print "Pep deleted  = $pep_deleted\n";

     }
     ###


     ##### protein domains

      if ($tr->translation()) {
        my $transl = $tr->translation();
        #print "  ", $tr->translation()->stable_id(), "\n";

        my @domain_feats = @{$transl->get_all_DomainFeatures()};
        #p(@domain_feats);
        foreach my $df(@domain_feats) { #cycle through features
        my @coord_name_array = ();
        my @coord_type_array =();


          my @gc_array = ();
          my $dbid         = $df->display_id();
          my $coord_type = "";
          my $domain_start = $df->start();
          my $domain_end   = $df->end();
          my $domain_desc  = $df->idesc();
          my $domain_desc2  = $df->hdescription();
          my $domain_source = $df->hseqname();
          #print "    $dbid $domain_start $domain_end $domain_desc $domain_desc2\n";

          my @genomic_coords = $tr->pep2genomic($domain_start,$domain_end);
          foreach my $gc (@genomic_coords){ 
              push @gc_array, $gc->start;
              push @gc_array, $gc->end;
          }

          #p(@genomic_coords);
          if ($domain_source =~ /^PF/ || $domain_source =~ /^SSF/){ ## only get pfam domains and SuperFamily  ##EDIT THIS TO INCLUDE OTHER DOMAIN TYPES, ALTHOUGH THINGS MIGHT START LOOKING A BIT CROWDED
             $coord_start = min(@gc_array);
             $coord_end = max(@gc_array);
             #p(@gc_array);
             print "domain_source = $domain_source, coord_start=$coord_start,coord_end =$coord_end \n";
             if ($crispr_min <  $coord_start  && $crispr_max > $coord_end){
                  $coord_type = "flanked";
                  push @coord_name_array, $domain_desc2;                   
                  push @coord_type_array, $coord_type;
                  $has_domain = 1;
              }
             elsif( $crispr_min > $coord_start && $crispr_max < $coord_end){
                   $coord_type = "internal_disrupted";
                   push @coord_name_array, $domain_desc2;
                   push @coord_type_array, $coord_type;
                   $has_domain = 1;
                                 }
             elsif( $crispr_min < $coord_start && $crispr_max > $coord_start && $crispr_max < $coord_end){ # 5' part of exon deleted
                   $coord_type = "5_disrupted";
                   push @coord_name_array, $domain_desc2;
                   push @coord_type_array, $coord_type;
                   $has_domain = 1;
              }
             elsif( $crispr_min > $coord_start && $crispr_max > $coord_end && $crispr_min < $coord_end){ # 3' part of exon deleted
                   $coord_type = "3_disrupted";
                   push @coord_name_array, $domain_desc2;
                   push @coord_type_array, $coord_type;
                   $has_domain = 1;
              }
            
              #p(@coord_name_array);
              for ($n =0; $n<@coord_name_array; $n++){ ## create a summary array for this as one line, so it can be added in later.
                 $domain_summ_string = "$domain_source, $domain_desc2, $domain_start - $domain_end, ".$coord_name_array[$n].", ".$coord_type_array[$n].";"; 
                 push  @domain_summary, $domain_summ_string;
                 print "$domain_source $domain_desc2 found in CRISPR deletion, $domain_start - $domain_end, ".$coord_name_array[$n].", ".$coord_type_array[$n]." CRISPR $crispr_min -  $crispr_max, DOMAIN $coord_start - $coord_end \n";
              }
               #p(@domain_summary);
            } #end of pfam domains
        } #end of features
      } # end of translation
     ##### END OF PROTEIN DOMAIN

     # my $tr_strand = $tr->strand();
     # my $ccds_xrefs = @{ $tr->get_all_DBEntries('CCDS') };
     # print "CCDS=".$ccds_xrefs."\n";
     #print "CCDS info for $tr_id ";

     #print "ccds = $ccds_name\n";
     #if($ccds_xrefs > 0) { ### only get information from genes which have a CCDS ref. Commented out in this build as it's better to get everything and filter later. But if you really only want curated transcripts, but it back in.
        ### initiate arrays
        my @start_frame_array = ();
        my @end_frame_array = ();
        my @rank_array = ();
        my $outofframe = "";
        my $deleted_region = 0; # see how much of cds and transcript is gone.
        ###

        ### Cycle through each exon and see what the consequences of the deletion are.
        my @t_exons = @{ $tr->get_all_Exons };
        my $last_exon = @t_exons[-1]->stable_id();
        print "last exon = $last_exon\n";
        foreach my $exon ( @t_exons ) {
           my $rank = $tr->exon_rank($exon);
           #print "exon rank $rank\n";
           #my $estring = feature2string($exon);
           my $revised_length = 0;   
           my $e_stable_id  = $exon->stable_id();
           my $e_seq_region = $exon->slice->seq_region_name();
           my $e_start      = $exon->start();
           my $e_end        = $exon->end();
           my $e_strand     = $exon->strand();
           my $e_phase      = $exon->phase();
           my $e_end_phase  = $exon->end_phase();
           my $ii = 0;

           if ($e_strand == 1){
              ### positive strand
              #print "temp_start = $temp_start, temp_stop = $temp_stop, e_start = $e_start, e_end = $e_end\n";
              ## Deletion flanks the exon
              if ($crispr_min <  $e_start  && $crispr_max > $e_end){ 
                 print "loop 1 ";
                 $crispr_stat = "Flanking";
                 push @start_frame_array, $e_phase;
                 push @end_frame_array, $e_end_phase;
                 push @rank_array, $rank;
                 push @exon_id_array, $e_stable_id;
                 $deleted_region = $deleted_region + ($e_end - $e_start) +1 ;
                 $hit = 1;
              }
              ### end of flanking deletion
		
              # Internal deletion
              elsif( $crispr_min > $e_start && $crispr_max < $e_end){ 
                 $crispr_stat = "Deletion within the exon";
                 $outofframe = "Deletion within the exon";
                 push @rank_array, $rank;
                 $deleted_region = $deleted_region + ($crispr_max - $crispr_min) +1;
                 $hit = 1;
              }
              ### end of internal deletion
              
	      # 5' part of exon deleted
                 elsif( $crispr_min < $e_start && $crispr_max > $e_start && $crispr_max < $e_end){ 
                 $crispr_stat = "Partial exon deletion";
                 $outofframe = "Partial exon deletion" ;
                 $deleted_region = $deleted_region + ($crispr_max - $e_start) +1 ;
                 print "tr_cds_region_start= $tr_cds_region_start , tr_cds_region_end=  $tr_cds_region_end\n";
                 push @rank_array, $rank;
                 $hit = 1;
              }
              ### end of 5' end of exon deleted
			  
	      # 3' part of exon deleted
              elsif( $crispr_min > $e_start && $crispr_max > $e_end && $crispr_min < $e_end){ 
                 print "3 part of exon\n";
                 $crispr_stat = "Partial exon deletion"; 
                 $outofframe = "Partial exon deletion";
                 $deleted_region = $deleted_region + ($e_end  - $crispr_min) +1 ;
                 push @rank_array, $rank; 
                 $hit = 1;
              }
              ### end of 3' end of exon deleted

              $temp_start = $e_start; #save start of last exon for next round. Redundant I think - should now happen in the array
              $temp_stop = $e_end; #save end of last exon for next round. ditto
           }
           ### end of plus strain
           else {
             ### negative strand. Not sure this makes a difference in the way Ensembl handles exon coordinates, but I'm too scared to take it out now.
             #print "looking on neg strand\n";
             if (($crispr_max >  $e_end)  && $crispr_min < $e_start){
                $crispr_stat = "Flanking";
                print "loop -1 ";
                push @start_frame_array, $e_phase;
                push @end_frame_array, $e_end_phase;     
                push @exon_id_array, $e_stable_id;
                push @rank_array, $rank;
                $deleted_region = $deleted_region + ($e_end - $e_start) +1 ;
                 $hit = 1;
             }
             elsif( $crispr_max < $e_end && $crispr_min > $e_start){
                $crispr_stat = "Deletion within the exon";
                $outofframe = "Deletion within the exon";
                push @rank_array, $rank;
                $deleted_region = $deleted_region + ($crispr_max - $crispr_min) +1;
                 $hit = 1;
             }
             elsif( $crispr_min < $e_start && $crispr_max > $e_start && $crispr_max < $e_end){
                $crispr_stat = "Partial exon deletion"; 
                $outofframe = "Partial exon deletion";
                push @rank_array, $rank; 
                $deleted_region = $deleted_region + ($crispr_max - $e_start) +1;
                 $hit = 1;
             }
             elsif( $crispr_min > $e_start && $crispr_max > $e_end && $crispr_min < $e_end){  
                $crispr_stat = "Partial exon deletion";
                $outofframe = "Partial exon deletion";
                push @rank_array, $rank;
                $deleted_region = $deleted_region + ($e_end  - $crispr_min) +1;
                 $hit = 1;
             }
             $temp_start_neg = $e_stop; #save end of last exon for next round. 
             $temp_stop = $e_start; #save end of last exon for next round.  
          } 
   
          $estring = sprintf( "%s\t%s\t%d\t%d\t(%+d)\t%s\t%s", $e_stable_id, $e_seq_region, $e_start, $e_end, $e_strand, $e_phase, $e_end_phase);
          print (OUT  "$in_gene\t$tr_id\t$estring\t$rank\t$crispr_min\t$crispr_max\t$crispr_stat\n");
          #print "$in_gene\t$tr_id\t$tr_ens\t$estring\t$rank\t$crispr_stat\n";
          if ($crispr_stat){
             push @crispr_stat_array, $crispr_stat; # make a summary of all affected exons.
          }
          $crispr_stat = "";

       } # end of each exon

       ### sort out frame issues
       $frame_array_size = @start_frame_array;
       #p@start_frame_array;
       #p@end_frame_array;
       print "crispr_start_frame = $crispr_start_frame, crispr_end_frame = $crispr_end_frame\n";
       if ($frame_array_size > 0){
          $crispr_start_frame = $start_frame_array[0];
          $crispr_end_frame = $end_frame_array[-1];
          $crispr_end_id = $exon_id_array[-1];
          print "crispr_start_frame = $crispr_start_frame, crispr_end_frame = $crispr_end_frame\n";                   
          if (($crispr_start_frame > -1 && $crispr_start_frame != $crispr_end_frame) && $crispr_end_id ne $last_exon){
             $outofframe = "Out of frame whole exon deletion";
          }
          elsif (($crispr_start_frame == $crispr_end_frame)&& $crispr_end_id ne $last_exon ){
             $outofframe = "In frame whole exon deletion";
          }
          elsif ($crispr_start_frame == -1){
             $outofframe = "first coding exon deleted";
          }
          elsif ($crispr_end_id eq $last_exon){
             $outofframe = "last coding exon deleted";
          }
          else {
             $outofframe = "check";
          }
       }
       ###

       ### Pull together all the above into a summary file
       print "gene $in_gene $status =  $outofframe\n";
       print "deleted size =  $deleted_region\n";
       if ($tr_cds_length > 0){
            #print "tr_loop"; 
            $percent_cds_del = int(($cds_deleted/$tr_cds_length)*100);
       }
       else {
            $percent_cds_del = "";
            #print "tr_loop2";
       }
       my $del_size = $crispr_max - $crispr_min;
       my $percent_del = int(($deleted_region/$tr_length)*100);
       print "del_size=$del_size, tr_length=$tr_length, percent_del=$percent_del, percent_cds_del= $percent_cds_del, deleted CDS length = $deleted_length, tr_cds_length=$tr_cds_length\n";
       $slice = $slice_adaptor->fetch_by_region( 'chromosome', $tr_chr,  $crispr_min, $crispr_max); # get deletion DNA sequence from ensembl.org
       my $del_sequence = $slice->seq();
       $del_sequence = uc($del_sequence);
	   
       ### Add to .summary file
       print (OUT2 "$project_name	$in_gene	$tr_id	$tr_ens	$tr_length	$tr_type	$tr_strand	$tr_source	$tr_chr	$crispr_min	$crispr_max	$del_size	".$rank_array[0]." to ".$rank_array[-1]."	$crispr_start_frame	$crispr_end_frame	$deleted_region	$percent_del	$tr_cds_length	$cds_deleted	$percent_cds_del	");
       foreach my $ii (@crispr_stat_array) {
          print (OUT2 "$ii\;");
       }
       print (OUT2 "	$outofframe	$cds_ATG	$crispr_source	");
       foreach my $h (@domain_summary) { 
          print (OUT2 "$h\;");
       }
       print (OUT2 "	$del_sequence	$pep_coord_start to $pep_coord_end	$pep_sequence	$pep_deleted");
       ###
       print "outofframe = $outofframe\n";
       if ($hit == 1){
          $hit_yes = "Yes";
          if ($outofframe eq "Out of frame whole exon deletion"){
             $summary = "Exon deleted and out of frame";
             print "OOF loop 1\n";
          }
          elsif ($percent_cds_del > 50 && $has_domain == 1){
             $summary = "More than 50% CDS deleted and domain disruption";
          }
          elsif ($has_domain == 1 && $percent_cds_del <= 50 && $cds_ATG eq "yes"){
             $summary = "Domain disruption and ATG removed";
          }
          elsif ($has_domain == 1 && $outofframe eq "Partial exon deletion" && $percent_cds_del <= 50){
             $summary = "Domain disruption and exon partially deleted"; 
          }
          elsif ($has_domain == 1 && $outofframe eq "In frame whole exon deletion" && $percent_cds_del <= 50){ 
             $summary = "Domain disruption and in-frame exon deletion";
          }
          elsif ($has_domain == 1 && $outofframe eq "Deletion within the exon" && $percent_cds_del <= 50){
             $summary = "Domain disruption and deletion within the exon";         
          }
          elsif ($has_domain == 1 && $percent_cds_del <= 50 && $cds_ATG eq "no"){
             $summary = "Domain disruption";
          }
          else { 
            $summary = "Less that 50% CDS deleted and no domain disruption";
          }
       }
       else {
         $hit_yes = "No";
         $summary = "";
       }
       print "summary = $summary\n";
       print (OUT2 "	$hit_yes	$summary	$url\n");

       ### A few conclusions to print to the user.
       $exon_count += scalar(@t_exons);
       print $tr->display_id." has ".scalar(@t_exons)." Exons\n\n\n\n";
	   ###
	   	   
  #} # end of CDDS loop
  } # end of transcript loop
    

} # end of gene loop



