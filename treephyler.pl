#!/usr/bin/perl -w
########################################################################
# Script name  :    treephyler.pl
#
# Date created :    February 2011
#
# Author       :    Fabian Schreiber <fab.schreiber@gmail.com>
# 
# . 
#
#  Starter file for analysis using Treephyler
#
########################################################################
## INCLUDES
####################
use Getopt::Std;
use Getopt::Long;
use File::Basename;
#use Bio::Perl;
use Pod::Usage;
#use Bio::SearchIO;
#use Bio::Tools::BPlite;
#use Bio::Tools::Run::StandAloneBlast;
use Data::Dumper;
#use Bio::Index::Fasta;
use Pod::Usage;
use File::Temp qw/ tempfile tempdir /;
use strict;
#use warnings;
no warnings;
#require 'Iostuff2.pm';
require 'Configuration.pm';
require 'Treephy.pm';
#require 'Stockholm.pm';

####################
# VARIABLES
####################
#my $options_file = "options.txt";
my $tmp_dir = "/tmp/";
my $dir = tempdir( CLEANUP => 1 );
my $start =  0;
my $end = 10000000000; 
my $execution_option = q();
my $infile = q();
my $all_sequences = q();
my $max_size = q();
my $all_predictions = q();
my $pfam_alignment_file = q();
my $pfam_taxonomy_file = q();
my $pfam_hmm_file = q();
my $pfam_fasta_file = q();
my $comet_format = q();
my %options = GetOptions ("mode|m=s" => \$execution_option,
											"assig|a=s" => \$all_predictions,
											"in|i=s" => \$infile,
											"msize=s" => \$max_size,
											"comet" => \$comet_format,
											"pfamA|p=s" => \$pfam_alignment_file,
											"pfamT|q=s" => \$pfam_taxonomy_file,
											"pfamH|r=s" => \$pfam_hmm_file,
											"pfamF|f=s" => \$pfam_fasta_file);
my $root_directory = $Configuration::root_directory || die "Please specify root directory\n";

Treephy::printhelp() if !$execution_option;

############## PARAMETERS START ##############
### INPUT FILES
my $all_sequences_converted = $infile."_converted";
my $all_predictions_converted = $all_predictions."_converted";
### OUTPUT_DIRECTORY
my $output_directory    = $root_directory."/".$Configuration::project_name;
### PFAM DIRECTORY
#my $pfam_directory = "/c1/scratch/fabian/db_meta/db_pfam/MetaPfammer/pfam";
my $pfam_directory = $root_directory."/pfam";
## SEQUENCE TYPE: Nucleotide or Amino Acid
my $sequence_type = $Configuration::sequence_type || die qq(Please specify sequence_type\n);
### EVALUE HMMSEARCH
my $evalue_hmmsearch = $Configuration::evalue_hmmsearch || die qq(Please specify evalue_hmmsearch\n);
### PERFORM COMPUTATION ON COMPUTER CLUSTER
my $computation_mode = $Configuration::computation_mode;
my $number_of_cores = $Configuration::number_of_cores;
### NUMBER OF PARALLEL JOBS ON COMPUTER CLUSTER
my $number_of_cluster_jobs = $Configuration::number_of_cluster_jobs;
### HMMER VERSION
my $hmmer_version = $Configuration::hmmer_version;
### MIINIMUM LENGTH FOR EGT SEQUENCE
my $minimum_length = $Configuration::minimum_length;
my %cluster_groups;
my $tree_computation_log = "tree_computation_log.txt";
############## PARAMETERS END ##############
# CHECK DEFINEDNESS
 Treephy::printhelp() and die  if(!defined $all_sequences || ! defined $all_predictions);

# CHECK EXISTENCE
BEGIN:{
require 'Treephy.pm';

Treephy::print_start();
}

####################
# MAIN
####################
    
MAIN: {

      #	Treephy::print_start();
      Treephy::printhelp() if $execution_option eq "help";
      if($execution_option ne "analyse" 
      && $execution_option ne "statistics"  
      && $execution_option ne "prepare" 
      && $execution_option ne "split" 
      && $execution_option ne "translate" 
      && $execution_option ne "eval"){
            print "\tUnrecognized option: \"$execution_option\".\n";
            exit();
      }
      ## make 6-frame translations
      if($execution_option eq "translate"){
            ## CHECK INPUT 
            if(! $infile ){
                  print "\tPlease specify an input file for this analysis\n";
                  Treephy::printhelp();
            }
            if(!-e $infile && ! -s $infile){
                  print "\tInput file $infile does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            print "\tStart translating sequences\n";
            my $translated_output_file = $infile."_translated";
            Treephy::translate_bioperl({all_egts_file=> $infile, egt_file =>$translated_output_file });
            print "\tPlease find the translated sequences in $translated_output_file\n";

      }
      if($execution_option eq "split"){
            ## CHECK INPUT 
            if(! $infile || ! $max_size){
                  print "\tPlease specify an input file or number of splits for this analysis\n";
                  Treephy::printhelp();
            }
            if(!-e $infile && ! -s $infile){
                  print "\tInput file $infile does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            print "\t[Splitting File] in progress.....\n";
            my $split_result = Treephy::split_file({split_file=> $infile,max_size => $max_size});
            if($split_result){
                  print "\tPlease find the splitted files here...\n";
            }

      }

      ## make statistics
      if($execution_option eq "statistics"){
            ## CHECK INPUT 
            if(! $infile){
                  print "\tPlease specify an input file for this analysis\n";
                  Treephy::printhelp();
            }
            if(!-e $infile && ! -s $infile){
                  print "\tInput file $infile does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            print "\t[MAKING STATISTICS] in progress.....\n";
            my $evaluation_result = Treephy::evaluate_statistics({new_assignments=> $infile});
            if($evaluation_result){
                  print "\tPlease find the statistics in $evaluation_result\n";
            }
      }
      if($execution_option eq "eval"){
            my $evaluation_result = Treephy::evaluate_all({input_dir=> $infile});
            if($evaluation_result){
                  print "\tPlease find the statistics in $evaluation_result\n";
            }
      }
      if($execution_option eq "prepare"){
            ## CHECK INPUT 
            if(!-e $pfam_alignment_file && ! -s $pfam_alignment_file){
                  print "\tPfam alignemnt file $pfam_alignment_file does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            if(!-e $pfam_fasta_file && ! -s $pfam_fasta_file){
                  print "\tFasta input file $pfam_fasta_file does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            if(!-e $pfam_taxonomy_file && ! -s $pfam_taxonomy_file){
                  print "\tTaxonomy input file $pfam_taxonomy_file does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            my $pfam_alignment_repository = "$pfam_directory/alignments";
            my $pfam_hmm_repository = "$pfam_directory/hmms";
            my $pfam_taxonomy_repository = "$pfam_directory/taxonomy";
            my %pfam_accessions_hash = ();
            mkdir($pfam_directory) if !-e $pfam_directory;
            mkdir($pfam_alignment_repository) if !-e $pfam_alignment_repository;
            mkdir($pfam_hmm_repository) if !-e $pfam_hmm_repository;
            mkdir($pfam_taxonomy_repository) if !-e $pfam_taxonomy_repository;
            if(! -e $pfam_directory || ! -e $pfam_alignment_repository || ! -e $pfam_hmm_repository || ! -e $pfam_taxonomy_repository){
                  print "\tThere was a problem creating the directories for PFAM alignments, etc. Exiting...\n";
                  exit;
            }
            
            ### FULL ALIGNMENT-FILE -> SINGLE FILES
            print "\tConverting PFAM alignments....\n";
            if (scalar <$pfam_alignment_repository/*>){
                  print "\t\talignments already exist. Skip this step\n";
            }else{
                  Treephy::convert_pfam_alignments({output_directory=> $pfam_alignment_repository, pfam_alignment_file=>$pfam_alignment_file});
            }
            print "\tFetching HMM records....\n";
            ### FETCH HMMS RECORDS 
            #Treephy::collect_pfam_accessionnumbers({ pfam_hmm_file=>$pfam_hmm_file, pfam_accessions_hashref =>\%pfam_accessions_hash });
            #Treephy::small_hmm_fetcher({output_directory=> $pfam_hmm_repository, pfam_hmm_file=>$pfam_hmm_file, pfam_accessions_hashref =>\%pfam_accessions_hash});
            
            Treephy::small_hmm_builder({output_directory=> $pfam_hmm_repository, alignment_repository=>$pfam_alignment_repository});
            
            ### FETCH TAXONOMY RECORDS 
            print "\tFetching Taxonomy records....\n";
            Treephy::collect_taxonomy_records({output_directory=> $pfam_taxonomy_repository, pfam_fasta_file => $pfam_fasta_file, taxonomy_file => $pfam_taxonomy_file});
      }
      if($execution_option eq "analyse"){
            ## CHECK INPUT 
            if(! $infile){
                  print "\tPlease specify an input file for this analysis\n";
                  Treephy::printhelp();
            }
            if(!-e $infile && ! -s $infile){
                  print "\tInput file $infile does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }
            if(!-e $all_predictions && ! -s $all_predictions){
                  print "\tPrediction file $infile does not exist or is empty.Aborting...\n";
                  Treephy::printhelp();
            }


            print "\t[PREPARATION] Cleaning old files before starting...";
            unlink($tree_computation_log);
            my @curr_pfams = glob("$output_directory/*");	
            my $sequences_assigned = 0;
            foreach(@curr_pfams){
                  unlink("$_/in_progress.txt");
                  #	   	`rm $_/*.tree`;
                  #	   	`rm $_/*.final`;
                  #	   	`rm $_/*.tmp`;
                  #		   	`rm $_/*.egt*`;
                  if(glob("$_/*.egt-all")){
                        $sequences_assigned = 1;
                  }
                  if(glob("$_/*.nucl")){
                        $sequences_assigned = 1;
                  }

            }
            print "done\n";
            if(!$sequences_assigned){
                  if(! -e $output_directory){ Treephy::create_single_directory($output_directory);}
                  my %curr_pfam_families = ();
                  my %assignment_by_acc_hash = ();
                  my $assignment_by_acc_hashref = \%assignment_by_acc_hash;
                  my %assignment_by_pfam = ();
                  my $assignment_by_pfam_hashref = \%assignment_by_pfam;

                  if(!Treephy::check_sequences_and_predictions({assignment_file=> $all_predictions, sequence_file =>$infile})){
                        #  		die "\tThere were headers in the sequence file that were not found in the prediction file\n";
                  }
                  if($comet_format){
                        ### Convert Comet output to UFO output
                        $assignment_by_acc_hashref =  Treephy::change_output_format({all_sequences=>$infile, 
                              old_assignments => $all_predictions,
                              new_assignments => $all_predictions_converted,
                              });
                              $all_predictions = $all_predictions_converted;
                  }
                        #    print "\tnew assignemnts in $all_predictions_converted, value of assignments: $all_predictions\n";
                        #   exit;
                        $assignment_by_acc_hashref =  Treephy::convert_and_collect({all_sequences_file=>$infile, 
                              assignment_by_acc_hashref => \%assignment_by_acc_hash,
                              all_sequences_file_converted => $all_sequences_converted,
                              all_predictions => $all_predictions,
                              });

                              ### ASSIGN EGTS TO PFAM-FAMILY FOLDERS
                              if(!keys(%{$assignment_by_acc_hashref})){
                                    die "\tThere was a problem reading the assignment and sequence file. Make sure the format is correct.Aborting..\n";
                              }
                              my $all_sequence_hashref = Treephy::preparation_step_slurp({ 
                                    egt_file=>$all_sequences_converted, 
                                    minimum_length => $minimum_length, 
                                    output_directory=> $output_directory, 
                                    assignment_hash_by_acc=>$assignment_by_acc_hashref, 
                                    sequence_type=>$sequence_type
                              });
                                    ### CHECK IF PRINTING WAS SUCCESSFUL
                                    print "\t[WRITE FILES] Successfully wrote files\n";
            }

            if($computation_mode eq "single"){
                  my $arg = "time perl process.pl  -a $all_predictions -i $infile  1> $output_directory/log.txt 2>$output_directory/warnings.txt";
                  #	print $arg."\n";
                  print "\t[START] Starting analysis on a single computer...($arg)\n";
                  `$arg`;
            }
            elsif($computation_mode eq "multi"){
                  for(my $jobs_counter = 1; $jobs_counter <= $number_of_cores; $jobs_counter ++){
                        my $arg = "time perl process.pl  -b $jobs_counter -a $all_predictions -i $infile 1> $output_directory/log$jobs_counter.txt 2>$output_directory/warnings$jobs_counter.txt&";
                        #	      print "$arg\n";
                        print "\t[START] Starting analysis on a single multi-core computer...process $jobs_counter\n";
                        ## START ON Multi-Core
                        system($arg);
                  }
            }
            elsif($computation_mode eq "cluster"){
                  ## Prepare bash script and start in cluster
                  ## Process i will analyse all sequences assigned to process i as given in the load_balance_file
                  for(my $jobs_counter = 1; $jobs_counter <= $number_of_cluster_jobs; $jobs_counter ++){
                        my $arg = "time perl process.pl -b $jobs_counter -a $all_predictions -i $infile ;";
                        my $bash_script = "bash_script.sh";
                        Treephy::write_bash_script($bash_script,$arg,"TTMPfam_$jobs_counter");
                        ## START ON CLUSTER
                        print "\t[START] Starting analysis on a computer cluster...process $jobs_counter\n";
                        `qsub -o $output_directory -e $output_directory $bash_script`;
                  }
            }
      }
      if($execution_option eq "t"){
            print "\t[Calculating time]\n";
            my $time_file = "time_file_$output_directory";
            open my $EGT_FILE, '<', $time_file || die "Could not open $time_file\n";

            open my $TIME_FILE, '>', "time_file_$output_directory" || die "Could not open time_file\n";
            ## grep real time
            print "grep real $output_directory/TTMPfam_*.e*\n";
            my @real_array = `grep real $output_directory/TTMPfam_*.e*`;
            my $real_value = 0;
            if(!scalar(@real_array)){
                  warn "\tCould not find real time output\n";
            }
            else{
                  foreach(@real_array){
                        /:real\s+(\d+)m/;
                        $real_value += $1;
                        print  "$1\n";
                  }
                  $real_value = int($real_value / scalar(@real_array));
                  ## grep user time
            }
            my @user_array = `grep user $output_directory/TTMPfam_*.e*`;
            my $user_value = 0;
            foreach(@user_array){
                  /:user\s+(\d+)m/;
                  $user_value += $1;
                  #			print {$TIME_FILE} "user $1 mins\n";
            }
            $user_value = int($user_value / scalar(@user_array));
            ## grep sys time
            my @sys_array = `grep sys $output_directory/TTMPfam_*.e*`;
            my $sys_value = 0;
            foreach(@sys_array){
                  /:sys\s+(\d+)m/;
                  $sys_value += $1;
                  #				print {$TIME_FILE} "sys $1 mins\n";

            }
            $sys_value = int($sys_value / scalar(@sys_array));
            #	print "real: $real_value\n\tuser: $user_value\n\tsys: $sys_value\n";
            close $TIME_FILE || die "Could not close TIME_FILE\n";

      }
}
####################
# END
####################
END{
#	Treephy::print_end();
}



