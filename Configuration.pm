########################################################################
# Script name  :    Configuration.pm
#
# Date created :    February 2011
#
# Author       :    Fabian Schreiber <fab.schreiber@gmail.com>
#  
# This is the configuration file for Treephyler
# See the User manual for detailed descriptions
# NOTE: Paths must end with "/"
########################################################################
package Configuration;
{
####
# SYSTEM PARAMETERS
#
# ENTER HERE THE PATH OF THE PROGRAMS
# Make sure you do not accidently delete the semicolon at the end of the lines
#####
## HMMER BIN DIRECTORY
        $hmmer = "/Users/bin/hmmer/";
	# hmmer version 3: "3"
	# hmmer version 2: "2"
	$hmmer_version = "3";
## FILENAME OF FASTTREE BINARY
	$fasttree = "/Users/bin/FastTree";
#####
# PROJECT PARAMETERS
#####

## PROJECT NAME
	$project_name = "project_1";
## ROOT DIRECTORY FOR ANALYSIS (absolute pathname required)
	$root_directory = "/user/work/treephyler/";
## SEQUENCE TYPE
	# PROTEINS = "p"
	# NUCLEOTIDES = "n"
	$sequence_type = "p";
## COMPUTATION ON CLUSTER	
	# "single" - single core 
	# "multi" - multi-core
	# "cluster" - (SGE) cluster
	$computation_mode = "single";
	## multi-core
	$number_of_cores = "1";
## NUMBER OF CLUSTER JOBS
	$number_of_cluster_jobs = "20";
## MINIMAL LENGTH OF EGT	
	$minimum_length = "100";
### EVALUE FOR HMMSEARCH [SIGNIFICANT HITS]
	$evalue_hmmsearch = "0.01";
}

1;
