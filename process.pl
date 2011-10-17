#!/usr/bin/perl -w
########################################################################
# Script name  :    process.pl
#
# Date created :    February 2011
#
# Author       :    Fabian Schreiber <fab.schreiber@gmail.com>
# 
# . 
#
#  Manages the workflow of the treephyler analysis.
########################################################################
## INCLUDES
####################
use Getopt::Std;
use Getopt::Long;
use File::Basename;
#require 'sys/syscall.ph';
use Bio::Perl;
use Pod::Usage;
use Switch;
use Bio::SearchIO;
use Data::Dumper;
use File::Copy;


use Pod::Usage;
use File::Temp qw/ tempfile tempdir /;
use strict;
#use warnings;
no warnings;
	require 'Configuration.pm';
    require 'Treephy.pm';
####################
# VARIABLES
####################
my $options_file =  "options.txt";
my $start =  0;
my $end = 10000000000; 
my $incrementor = q{};
my $all_sequences = q();
my $all_predictions = q();
our  $own_process_id = q{};

### COMMAND LINE OPTIONS
my %options = GetOptions (
			"inc|b=s" => \$own_process_id,
			"assig|a=s" => \$all_predictions,
			"in|i=s" => \$all_sequences);
our $start_process_run = time();
############## PARAMETERS START ##############
my $root_directory = $Configuration::root_directory || die "Please specify root directory\n";
### INPUT FILES
### OUTPUT_DIRECTORY
my $output_directory    = $root_directory."/".$Configuration::project_name;
### PFAM DIRECTORY
my $pfam_directory = $root_directory."/pfam";
#my $pfam_directory = "/c1/scratch/fabian/db_meta/db_pfam/MetaPfammer/pfam";
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
my $minimum_length = $Configuration::minimun_length;
my %cluster_groups;
my $tree_computation_log = "tree_computation_log.txt";
############## PARAMETERS END ##############

### PFAM ALIGNMENT REPOSITORY
my $pfam_alignment_repository = "$pfam_directory/alignments";
my $pfam_hmm_repository = "$pfam_directory/hmms";
my $pfam_taxonomy_repository = "$pfam_directory/taxonomy";

my $all_sequences_converted = $all_sequences."_converted";
my $all_predictions_converted = $all_predictions."_converted";
my $gblocking = $Configuration::gblocking || "0";

### COUNTS NUMBER OF FASTA ENTRIES
my $no_current_fasta_entry      = 1;


BEGIN:{
	require 'Treephy.pm';
	Treephy::print_start();
}
####################
# MAIN
####################
    
MAIN: {
			if(! -e $output_directory){ Treephy::create_single_directory($output_directory);}
			my %curr_pfam_families = ();
			my ($pfams_to_take, $taxonomy_to_take, $trees_to_take,$hmms_to_take) = q{};
			print "[Main Analysis] Iteratively analyse PFAM families\n";
		### DIRECTORIES WITH ALIGNMENT, TAXONOMY, AND HMM FILES
					$pfams_to_take = $pfam_alignment_repository;
					$taxonomy_to_take = $pfam_taxonomy_repository;
					$hmms_to_take = $pfam_hmm_repository;
				my $number_of_curr_family = 1;
				my $file_type  = ($sequence_type eq "p")? "egt-all" : "nucl";
#####
## LIST ALL PFAM FAMILIES	
#####
my @curr_pfams = glob("$output_directory/*");
	if(!@curr_pfams){
		die "\tDirectory $output_directory does not contain any folder. Make sure the right project is selected in the configuration file\n";
	}
	my $number_of_families = scalar(@curr_pfams);
	PFAM_FAMILY:    
	foreach my $curr_domain(sort @curr_pfams){
		next if -f $curr_domain;
		#Skip the dirname part
		$curr_domain = basename($curr_domain);
		my $lockfile = "$output_directory/$curr_domain/in_progress.txt";
		if(!Treephy::create_lockfile({lockfile => $lockfile})){
			print "skipping $curr_domain. Is already being computed (lockfile)\n";
			next PFAM_FAMILY;
		}
		if(!-e $lockfile){
			print "\tcould not get lock for pfam $curr_domain\n";
			next PFAM_FAMILY;
		}
		print $curr_domain." (".$number_of_curr_family++." / $number_of_families)\n";
		if(!-e "$output_directory/$curr_domain/$curr_domain.$file_type"){
			warn "\tThere are no sequences for $curr_domain\n";
			next PFAM_FAMILY;
		}
		my $number_of_query_sequences_before_cleaning = `grep ">" $output_directory/$curr_domain/$curr_domain.$file_type | wc -l`;
		chomp($number_of_query_sequences_before_cleaning);
		print "\t[STATS] Group contains $number_of_query_sequences_before_cleaning query sequences\n";
### TRANSLATE SEQUENCES
			my $translate_status = 1;
			if(!($sequence_type eq "p")){
				$translate_status = Treephy::translate_egts({curr_domain=> $curr_domain, 
															output_directory => $output_directory,
								#							only_estscan => $only_estscan,
															pfam_alignment_repository => $pfams_to_take});
			}
			next if !$translate_status;
### REMOVE INSIGNIFICANT HITS FROM EGTS
			my $remove_status = 0;
			if($translate_status){
					print "\t[CLEAN SEQUENCES] Removing insignificant hits for $curr_domain...\n";
				$remove_status = Treephy::remove_insignificant_hits({ curr_pfam_family => $curr_domain, 
																	output_directory=> $output_directory,
																	pfam_alignment_repository => $pfams_to_take, 
																	evalue_hmmsearch=>$evalue_hmmsearch,
																	pfam_hmm_repository => $hmms_to_take,
																	hmmer_version => $hmmer_version});
			}
			
			next if !$remove_status;
			my $number_of_query_sequences = `grep ">" $output_directory/$curr_domain/$curr_domain.egt | wc -l`;
			chomp($number_of_query_sequences);
			print "\t[STATS] Group now contains $number_of_query_sequences query sequences\n";

### COMPUTE HMMALIGN ALIGNMENTS
			my $align_status = 0;
			if($remove_status){
				$align_status = Treephy::align_sequences({ curr_pfam_family => $curr_domain, 
															output_directory=> $output_directory, 
															pfam_alignment_repository => $pfams_to_take,
															pfam_hmm_repository => $hmms_to_take,
															hmmer_version => $hmmer_version});
			}
			next if !$align_status;
### COMPUTE TREES
			my $compute_status = 0;
			if($align_status){
				my $number_of_known_sequences = `grep ">" $output_directory/$curr_domain/$curr_domain.final | wc -l`;
				chomp($number_of_known_sequences);
		### GET ALIGNMENT LENGTH
		my $alignment_length = `grep LENG $hmms_to_take/$curr_domain.hmm`;
		$alignment_length =~ /LENG\s*(\d+)/;
		$alignment_length = $1;

### REJECT ALL WITH < 10 NEW SEQUENCES
		#		if($number_of_query_sequences < 10){
		#				print "\t[TREE BUILDING] Group contains too few ($number_of_query_sequences <10) query sequences.Skipping....\n";
		#				next PFAM_FAMILY;
		#		}
					if($number_of_known_sequences > 3000){
						print "\t[TREE BUILDING] Group contains too many ($number_of_known_sequences>3000) sequences.Skipping....\n";
						next PFAM_FAMILY;
					}
### MEASURE TIME FOR TREE BUILDING
				our $start_run = time();
			my $tree_method_used = q{};
### SELECTION OF DIFFERENT TREE BUILDING MODI	
							print "\t[TREE BUILDING] Computing tree using FastTree: Alignment contains ".($number_of_query_sequences + $number_of_known_sequences)." taxa\n";
								$compute_status = Treephy::compute_tree({ curr_pfam_family => $curr_domain, 
																		    output_directory=> $output_directory, 
																		    pseudocounts => 0,
																		    standard => 1,
																		    gblocking =>$gblocking });
								$tree_method_used = "FastTree";
				my $end_run = time();
				my $run_time = $end_run -  $start_run;
			my $current_hmm_alignment_file = $hmms_to_take."/".$curr_domain.".hmm";
		my $grep_length_of_alignment = `grep "LENG" $current_hmm_alignment_file`;
		$grep_length_of_alignment =~ /LENG\s+(\d+)/;
		$grep_length_of_alignment = $1;
		if(!defined $grep_length_of_alignment){
			"\t\t[STATISTICS] Could not grep alignment length\n";
		}
		
## LOG COMPUTATION TIME OF TREE
			my $tree_computation = "$run_time domain:$curr_domain met:$tree_method_used len:$grep_length_of_alignment seq:$number_of_known_sequences ($number_of_query_sequences)";
			if(!$compute_status){
				$tree_computation .= " (failed)";
			}
			my $tree_computation_log = "$output_directory/tree_computation_log.txt";
			open  my $TREE_COMPUTATION_FH, '>>', $tree_computation_log or warn "Could not open tree_computation_log: $!";
			flock $TREE_COMPUTATION_FH, 2;                 # try to lock the file
			print {$TREE_COMPUTATION_FH} $tree_computation."\n" or warn "Could not write tree_computation_log: $!";
			flock $TREE_COMPUTATION_FH, 8;
			close $TREE_COMPUTATION_FH or warn "Could not close tree_computation_log: $!";


			}
			### REMOVE ERROR FILES;
			if(glob("core.*")){
			`rm core.*`;
			}
			next if !$compute_status;
### PARSE TREES
			if($compute_status){
				print "\t[TREE PARSING] Parsing tree\n";
				Treephy::parse_tree({ curr_pfam_family => $curr_domain, 
												output_directory=> $output_directory, 
												taxonomy_repository => $taxonomy_to_take,
												id => $own_process_id});
			}
	    }


}
####################
# END
####################

END{
			Treephy::print_end();
			my $end_run = time();
			my $run_time = $end_run -  $start_process_run;
			my $time_computation = "$own_process_id $run_time";
			my $time_computation_log = "$output_directory/time_logger.txt";
			open  my $time_COMPUTATION_FH, '>>', $time_computation_log or warn "Could not open time_computation_log: $!";
			flock $time_COMPUTATION_FH, 2;                 # try to lock the file
			print {$time_COMPUTATION_FH} $time_computation."\n" or warn "Could not write time_computation_log: $!";
			flock $time_COMPUTATION_FH, 8;
			close $time_COMPUTATION_FH or warn "Could not close time_computation_log: $!";
}

####################
### LIST OF MAIN METHODS
####################



	


