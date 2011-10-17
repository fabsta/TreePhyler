package Treephy;
########################################################################
# Script name  :    Treephy
#
# Date created :    February 2011
#
# Author       :    Fabian Schreiber <fab.schreiber@gmail.com>
# 
#
# Library of functions for treephyler
########################################################################
use warnings;
use Bio::Perl;
use File::Basename;
#use File::Dirname;
use Data::Dumper;
use Bio::AlignIO;
use Switch;
use File::Copy;
use Bio::SimpleAlign;
use Cwd;
use Getopt::Std; 
#use Bio::Tools::Run::StandAloneBlast;
use Bio::TreeIO;
use Bio::SeqIO;
use IO::Handle;
#use List::MoreUtils qw(any all);
use File::Temp qw/ tempfile tempdir /;
use strict;
#require 'Stockholm.pm';
require 'Configuration.pm';
#use Time::HiRes qw/usleep/;

our $LINE_LENGTH = 64;

####################################################################
# PARAMETERS
# all_sequences - File containing all input sequences
# new_assignments - of egt to be further considered
# old_assignments - for results of the analysis
####################################################################
sub change_output_format(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $all_sequences = $arg_ref->{all_sequences};
my $new_assignments = $arg_ref->{new_assignments};
my $old_assignments = $arg_ref->{old_assignments};

print "\tConvert assignment file\n";
## READ ALL FASTA_HEADERS INTO ARRAY
my @all_fasta_headers = `grep ">" $all_sequences`;
print "\tFound ".scalar(@all_fasta_headers)." fasta header\n";

## READ THE ASSIGNMENT FILE INT
### Slurping File
my $fh = new IO::Handle;
open $fh, '<', $old_assignments;
my $text = q{};
sysread $fh, $text, -s $fh;
$fh->flush;
close $fh;
#print "\tAll assignments from $old_assignments read into memory\n";
#exit;
my @whole_lines = split(/\n/,$text);
if(!scalar(@whole_lines)){
  print "\t[Convert] no elements in whole_lines array\n";
  print "\tmistake\n";
  return();
}
	open my $NEW_ASSIGNMENT_FH, '>', $new_assignments or warn "Can't open ".$new_assignments.": $!";

#print "\twriting new assignment file $new_assignments\n";
my %multiple_assignments_hash = ();
my $counter = 0;
my $wrote_first_line = 0;
foreach my $line(@whole_lines){
  $line =~ /\s*(\d+)\s*,\s*(PF\d+)/;
  my ($seq_number,$pfam_family) = ($1,$2);
#	print "$seq_number \t $pfam_family\n";
#last if $counter > 10;
  if(!defined $seq_number || !defined $pfam_family){
  print "\t[Convert] no elements in $line\n";
    return();
  }
#  while(length($pfam_family) < 5){
 #   $pfam_family = "0".$pfam_family;
 # }
  my $corr_header = $all_fasta_headers[$seq_number-1];
  chomp($corr_header);
  $corr_header =~ s/>//;
 # print "\tcorr header is $corr_header -> $pfam_family\n";
  #exit;
  if(exists($multiple_assignments_hash{$corr_header})){
	  print $NEW_ASSIGNMENT_FH "$pfam_family\n";
  }
  else{
	print 	$NEW_ASSIGNMENT_FH "\n" if $wrote_first_line;	
  print $NEW_ASSIGNMENT_FH $corr_header."\n$pfam_family\n";
  $wrote_first_line = 1;
  }
  $multiple_assignments_hash{$corr_header} = 1;
#print $all_fasta_headers[$seq_number]."\nPF$pfam_family";
}

	close $NEW_ASSIGNMENT_FH or warn "Can't close ".$new_assignments.": $!";
print "\t[Done]\n";
}

####################################################################
# PARAMETERS
# egt_file - File containing all input sequences
# minimum_length - of egt to be further considered
# output_directory - for results of the analysis
# assignment_hash_by_acc - contains the sequences to write as a reference to a hash
# sequence_type - protein or nucleotide type 
sub preparation_step_slurp(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
#	print Dumper $arg_ref;
#	exit;
## File containing all input sequences
	my $egt_file = $arg_ref->{egt_file};
## Minimum length for sequence
	my $minimum_length = $arg_ref->{minimum_length};
## output directory for analysis
	my $output_directory = $arg_ref->{output_directory};
	#my $pfam2ufo_mapping_hashref = $arg_ref->{pfam2ufo_mapping_hashref};
	my $sequence_folder = $arg_ref->{sequence_folder};
## sequences to save
	my $assignment_hashref_by_acc = $arg_ref->{assignment_hash_by_acc};
## protein or nucleotide sequences
	my $sequence_type = $arg_ref->{sequence_type};
#### OTHER VARIABLES
	my $temp_alignment_file_in = q{};
	my $template_file = q{};
	my $hmm_fetch_file = q{};
	my $hmm_emit_file = q{};
	my $egt_counter = 1;
	my $text = q{};
	my $found_egt_counter = 0;
	my %all_sequences_hash = ();

## CHECK DEFINEDNESS
if(!defined $egt_file || !defined $minimum_length || !defined $output_directory || !defined $assignment_hashref_by_acc || !defined $sequence_type){
		warn "\t\t[WRITE FILES] Problem with parameters for function preparation_step_slurp\n";
		return();
}																
																

	print "\t[WRITE FILES] This process will analyse ".keys(%$assignment_hashref_by_acc)." sequences\n";

## CHECK EXISTENCE
	# INPUT FILE EXISTS?
	if(!-e $egt_file || ! -s $egt_file){
		warn "\t\t[WRITE FILES] File containing all sequences does not exist or is empty ($egt_file)\n";
		return 0;
	}
	## Output directory exists?
	die "$output_directory does not exist" if (!-e $output_directory);

	### Slurping File into memory (very fast)
my $fh = new IO::Handle;
open $fh, '<', $egt_file;
sysread $fh, $text, -s $fh;
$fh->flush;
close $fh;

if(!keys(%$assignment_hashref_by_acc)){
	warn "\t[] No accs to keep\n" and return undef;
}

print "\t[WRITE FILES] All Sequences read into memory!\n";
my @whole_lines = split(/\n/,$text);
my $curr_sequence =q{};
my $debug_counter = 0;
my $debug_counter2 = 0;
my $curr_header = q{};
my $curr_domain = q{};
my %acc2taxname_hash = ();

print "\t[WRITE FILES] Process all sequences in $egt_file\n";
#foreach my $taxon(<$sequence_folder/*>){
#foreach my $curr_line(@whole_lines){
#print "\tHas to process ".scalar(@whole_lines)." lines\n";
foreach my $line (@whole_lines){
#	print "\tline: $line\n";	
## Iterate over all lines from input file
		## all sequences to find already found? 
		if($line =~ />(QSQ\d+)/){
		#	$curr_header = $1;
			if($curr_sequence){
		#			print "\t\tseq exists (header: \"$curr_header\")\n";
						if(exists $$assignment_hashref_by_acc{$curr_header}){
								($curr_domain = $$assignment_hashref_by_acc{$curr_header});
		### NOW CHECK IF PFAM FAMILY IS A GOOD ONE
		### Sequence has the defined minimal length
								if(length($curr_sequence) >= $minimum_length){
		### In case of multiple assignments for a sequence
			## CONVERT SEQ_HEADER: 
			#Acropor_mi_DMPC1875133AT --> *.fa
				#								print "\t\t\t$curr_domain\n";
												if($curr_domain =~ /#/){
												### Save each assignment
													my @tmp_split = split(/#/,$curr_domain);
													foreach my $splitted_domain(@tmp_split){
														next if $splitted_domain eq q{};
														$all_sequences_hash{$splitted_domain} .= ">$curr_header\n$curr_sequence\n";
						#								print "\t\t\t\t$curr_domain: $curr_header\n".substr($curr_sequence,0,10)."\n";
													}
												}
												else{
					#								print "\t\t\t\t$curr_domain: $curr_header\n".substr($curr_sequence,0,10)."\n";
													$all_sequences_hash{$curr_domain} .= ">$curr_header\n$curr_sequence\n";
												}
								}
								else{
							### Log too short sequences
									open TOO_SHORT, '>>', "too_short.txt";
									print TOO_SHORT "$curr_header, length: ".length($curr_sequence)."\n";
									close TOO_SHORT;	
								}
						}
					else{
			#			print "\tWarn: Could not find $curr_header in data structure\n";
					}
			}		
			$line =~ />(QSQ\d+)/;
			$curr_header = $1;
#			print "\t\theader to $curr_header\n";
				#last if $found_egt_counter > 100;
			$found_egt_counter++;
	#		print "$found_egt_counter / $number_of_egts_to_find = ".(($found_egt_counter * 100) / $number_of_egts_to_find)." % found \n";
			$curr_sequence = q{};
		}
		else{
			next if $line =~ /^$/;
			$curr_sequence .= $line;
		}
	}	

	if(!keys(%all_sequences_hash)){
		die "\t[WRITE FILES] Could not read input sequences in $egt_file\n";
	}
#}
### WRITING CONVERSION FILE "header_conversions.txt"
				my $header_conversion_file = "header_conversions.txt";
				unlink($header_conversion_file) if -e $header_conversion_file;
				open my $HEADER_CONVERSION_FILE, '>', $header_conversion_file || die "Could not open $header_conversion_file\n";
				foreach(keys %acc2taxname_hash){
					print {$HEADER_CONVERSION_FILE} $_."-->".$acc2taxname_hash{$_}."\n";
				}
				close $HEADER_CONVERSION_FILE || die "Could not close $header_conversion_file\n";

	print "\t[WRITE FILES] Writing ".(keys (%all_sequences_hash))." Files\n\n";
		foreach my $curr_domain (sort keys %all_sequences_hash){
#				print "\t[WRITE FILES] writing sequences for Pfam family ".$$pfam2ufo_mapping_hashref{$curr_domain}."\n";
#			my $assignments_dir = $output_directory."/".$$pfam2ufo_mapping_hashref{$curr_domain}."/assignments/";
	#		foreach my $curr_taxon (sort keys(%{$all_sequences_hash{$curr_domain}})){
			my $template_file = $output_directory."/".$curr_domain."/".$curr_domain;
	### Nucleotid or protein sequences
		# nucleotide sequences will be temporarily saved in .nucl. The translated sequences then will be found in .egt-all
		# protein sequences will be saved in .egt-all
				$template_file .= $sequence_type eq "p"? ".egt-all":".nucl";
				&create_single_directory($output_directory."/".$curr_domain);
			#	&create_single_directory($output_directory."/".$$pfam2ufo_mapping_hashref{$curr_domain}."/assignments");
				open my $FASTA_FILE, '>>', $template_file || die "Could not open $template_file\n";
				print {$FASTA_FILE} $all_sequences_hash{$curr_domain};
				close $FASTA_FILE || die "Could not close $template_file\n";
			#}
		#	exit;
	}
return 1;
}
####################################################################
# PARAMETERS
# output_directory
# curr_pfam_family
# pfam_alignment_repository
# pfam_hmm_repository
sub align_sequences(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $output_directory = $arg_ref->{output_directory};
	my $curr_pfam_family = $arg_ref->{curr_pfam_family};
	my $pfam_alignment_repository = $arg_ref->{pfam_alignment_repository};
	my $pfam_hmm_repository = $arg_ref->{pfam_hmm_repository};
	my $hmmer_version = $arg_ref->{hmmer_version};
#### OTHER VARIABLES
	my $curr_path = $output_directory."/".$curr_pfam_family;
	my $current_hmm_alignment_file = $pfam_hmm_repository."/".$curr_pfam_family.".hmm";
#	my $current_seed_alignment_file = "/c1/scratch/fabian/db_meta/db_pfam/Treephyler/hmms_full/$curr_pfam_family.full";
	my $current_seed_alignment_file = $pfam_alignment_repository."/$curr_pfam_family.full";
	my $current_tmp_alignment_file = $curr_path."/".$curr_pfam_family.".tmp";
	my $current_tmp2_alignment_file = $curr_path."/".$curr_pfam_family.".tmp2";
	my $current_final_alignment_file = $curr_path."/".$curr_pfam_family.".final";
	my $current_egt_file = $curr_path."/".$curr_pfam_family.".egt";
		
		
	# CHECK DEFINEDNESS
if(!defined $output_directory || !defined $curr_pfam_family || !defined $pfam_alignment_repository){
		warn "\t\t[HMMALIGN] Problem with parameters for function align_sequences\n";
		return();
}	
																										
	# CHECK EXISTENCS
		if(!-e $current_hmm_alignment_file){
			warn "\t[HMMALIGN] No HMM for $curr_pfam_family found\n";
			return 0;
		}
		if(!-e $current_seed_alignment_file){
			warn "\t[HMMALIGN] No Alignment file for $curr_pfam_family found $current_seed_alignment_file\n";
			return 0;
		}
		if(!-e $current_egt_file){
			warn "\t[HMMALIGN] There were no significant hits for $curr_pfam_family\n";
			return 0;
		}
		print "\t[HMMALIGN] Align new sequences to existing alignment\n";

		### HMMER 2
		
		use IPC::Open3;
		use File::Spec;
		use Symbol qw(gensym);
		my $cmd;
		$cmd = $Configuration::hmmer."/hmmalign --trim --mapali $current_seed_alignment_file -o $current_tmp_alignment_file  $current_hmm_alignment_file $current_egt_file";
	#	print $cmd."\n";
		open(NULL, ">", File::Spec->devnull);
#		print "\tperforming "
		my $pid = open3(gensym, \*PH, ">&NULL", "$cmd");
		while( <PH> ) { }
		waitpid($pid, 0);
		
#		`$Configuration::hmmer/hmmalign -m --outformat clustal --withali $current_seed_alignment_file -o $current_tmp_alignment_file  $current_hmm_alignment_file $current_egt_file 2> /dev/null`;
	## ALIGNING WAS SUCCESSFULL?
		if(!-e $current_tmp_alignment_file){
			warn "\t[HMMALIGN] There was a problem usign hmmalign for $curr_pfam_family\n";
			return();
		}
		print "\t[FORMAT CONVERSION] Converting Stockholm --> FASTA\n";
		my $stockholm2fasta_call = "perl stockholm2fasta.pl -g $current_tmp_alignment_file > $current_final_alignment_file";
		`$stockholm2fasta_call`;
	#	print $stockholm2fasta_call."\n";
		## CONVERSION WAS SUCCESSFULL?		
		if(!-e $current_final_alignment_file && ! -s $current_final_alignment_file){
			warn "\t\tNo tmp_fasta_alignment for $curr_pfam_family found. Skipping this family\n";
			return();
		}
	## ALIGNMENT CONTAINS QUERY SEQUENCES?
		my @query_sequences = `grep ">QSQ" $current_final_alignment_file`;
		
		if(!@query_sequences){
			warn "\t\tNo query sequences in current_final_alignment_file for $curr_pfam_family found.Skipping...\n";
			return();
		}
		return 1;
	
}

####################################################################
# PARAMETERS
# curr_pfam_family
# output_directory
# pseudocounts
# gblocking
# standard
sub compute_tree(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $output_directory = $arg_ref->{output_directory};
	my $curr_pfam_family = $arg_ref->{curr_pfam_family};
	my $pseudocounts = $arg_ref->{pseudocounts};
	my $gblocking = $arg_ref->{gblocking};
	my $standard = $arg_ref->{standard};
	
	
#### OTHER VARIABLES
	my $final_alignment_file = "$output_directory/$curr_pfam_family/$curr_pfam_family.final";
	my $tree_file = "$output_directory/$curr_pfam_family/$curr_pfam_family.tree";
	
# CHECK DEFINEDNESS
	if(!defined $output_directory || !defined $curr_pfam_family || !defined $pseudocounts){
		warn "\t\t[TREE BUILDING] Problem with parameters for function compute_tree\n";
		return();
}	

## CHECK EXISTENCE
	if (! -e $final_alignment_file){
		warn "\t\t[TREE BUILDING] Alignment file for $curr_pfam_family does not exist\n";
		return 0;
	}
	## REMOVE EXISTING TREE FILE
	`rm $tree_file` if -e $tree_file;
	### CALL FASTTREE
	if($standard){
#		print $Configuration::fasttree." -quiet  $final_alignment_file> $tree_file\n";
		`$Configuration::fasttree -quiet  $final_alignment_file> $tree_file`;
		if (! -e $tree_file){
			warn "\t\t[TREE BUILDING] Tree file for $curr_pfam_family does not exist\n";
			return();
		}
		return 1;
	}
	if($pseudocounts){
#		print "\tpseudocount\n";
		`$Configuration::fasttree -quiet -pseudo -fastest $final_alignment_file> $tree_file`;
	}
	else{
#		`$Configuration::fasttree -quiet -fastest $final_alignment_file> $tree_file`;
#	print $Configuration::fasttree." -quiet -fastest -noml -boot 0 $final_alignment_file> $tree_file\n";
		`$Configuration::fasttree -quiet -fastest -noml -boot 0 $final_alignment_file> $tree_file`;
	}
	if (! -e $tree_file){
		warn "\t\t[TREE BUILDING] Tree file for $curr_pfam_family does not exist\n";
		return();
	}
	return 1;
}

####################################################################
# PARAMETERS
# curr_pfam_family
# output_directory
sub parse_tree(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $output_directory = $arg_ref->{output_directory};
	my $curr_pfam_family = $arg_ref->{curr_pfam_family};
	my $taxonomy_repository = $arg_ref->{taxonomy_repository};
	my $id = $arg_ref->{id};
	### Supress warnings
	no warnings 'recursion';
	
#### OTHER VARIABLES
	my $curr_tree_file = "$output_directory/$curr_pfam_family/$curr_pfam_family.tree";
	my $curr_tax_file = "$taxonomy_repository/$curr_pfam_family.tax";
	my $curr_egt_file = "$output_directory/$curr_pfam_family/$curr_pfam_family.egt";
	my $annotation_file = $output_directory."/taxonomic_assignments$id.txt";
	my %conversion_hash = ();
	my $conversion_file = "conversion.txt";
## CHECK DEFINEDNESS
	if(!defined $output_directory || !defined $curr_pfam_family || !defined $taxonomy_repository){
		warn "\t\t[PARSE TREE] Problem with parameters for function parse_tree\n";
		return();
}	
																										
																										
## CHECK EXISTENCE
	if(! -e $curr_tree_file || ! -s $curr_tree_file){
		warn "\t\t[PARSE TREE] No Tree file (curr_tree_file) for $curr_pfam_family\n";
		return();
	}

	if(! -e $curr_tax_file || ! -s $curr_tax_file){
		warn "\t\t[PARSE TREE] No Taxonomy file ($curr_tax_file) for $curr_pfam_family\n";
		return();
	}

	if(! -e $curr_egt_file || ! -s $curr_egt_file){
		warn "\t\t[PARSE TREE] No EGT file ($curr_egt_file) for $curr_pfam_family\n";
		return();
	}
	### Check existence of annotation and conversion file
	if(! -e $conversion_file || ! -s $conversion_file){
		warn "\tConversion file $conversion_file does not exist\n";
		return 0;
	}
### READ CONVERSIONS e.g. tr|A2U4H6|A2U4H6_BACCO Catalase OS=Bacillus coagulans 36D1 GN=BcoaDRAFT_2314-->QS23
		read_conversion_file({hash_ref => \%conversion_hash});
	#print "\tCollecting header from fasta file\n";
	### Collect the new headers
	my @new_headers_array = `grep ">" $curr_egt_file`;
	my %new_header_hash = ();
	foreach(@new_headers_array){
		my $tmp = $_;
		$tmp =~ s/\|/_/g;
		$tmp =~ s/>|\n//g;
		$new_header_hash{$tmp} = 1;
	}
	if(!%new_header_hash){
		print "\tCould not grep headers from\n";
		return();
	}
#	print "\tReading tree file\n";
## READ TREE
	my $treeio = Bio::TreeIO->new(-format => 'newick',
			      -file => $curr_tree_file,
			      -internal_node_id => 'bootstrap');
	my $tree = q{};
	$tree = $treeio->next_tree;
	my @cleaned_nodes = (); 
	if(!$tree->get_nodes){
		print "\tCould not extract nodes from tree. Skipping..\n";
	}
#	print "\tChecking single nodes\n";
	foreach($tree->get_nodes){
		next if !defined $_->id;
#		print $_->id."\n" if $_->id =~ /QS/;
		next if $_->id =~ /^\d/; # QUERY nodes start with QSQ
		next if $_->id eq q{};
		#print $_->id."\n" if $_->id =~ /tr/;
		push(@cleaned_nodes,$_); # ONLY LOOK AT THE QUERY NODES
		#if exists $new_header_hash{$_->id};
	}
	if(!@cleaned_nodes){
		print "\tCould not find query sequences in final tree.Skipping...\n";
		return();
	}
#	print "\titerating over nodes\n";

	CLEANED_NODE:
	foreach(@cleaned_nodes){
		next if $_ eq '';
		if($_->id =~ /^QSQ/){
			my $curr_node = $_;
			my $temp_node = $curr_node;
			my $taxonomy_string = q();
	#		print "\t\tDetermining phylogenetic position for for ".$curr_node->id."\t";
			
			my $enough_information = 0;
			while(my $curr_ancestor = $temp_node->ancestor){
			#	print "\t\t\tCurr ancestor is: ".$curr_ancestor->id."\t";
				my @clade_nodes = $curr_ancestor->get_all_Descendents;
				#print "with ".scalar(@clade_nodes)." descendents\n";
			#	print "\t\t\t\t";
				my %fast_hash = ();
				foreach(@clade_nodes){
					next if !defined $_->id;
					next if $_->id =~ /QSQ/;
					next if $_->id =~ /^\d/;
					#print $_->id."\t";
					my $tmp_header = $_->id;
				#	print "\tmp_header is $tmp_header\n";
					$tmp_header =~ s/\/.*//;
					#print $1."\n";
					$fast_hash{$tmp_header} = 1;
				}
				#print Dumper %fast_hash;
				#print "with ".keys(%fast_hash)." descendents\n";
				if(keys %fast_hash < 3){
					$temp_node = $curr_ancestor;
					next;
				}
				#exit;
				$taxonomy_string = &determine_taxonomic_information_for_clade({clade_nodes => \%fast_hash, curr_tax_file => $curr_tax_file});	
				#print $taxonomy_string." with return\n";
				#exit;
				last if $taxonomy_string ne q{};
				$temp_node = $curr_ancestor;
				#print "\n";
			}
		#	print "\t\t".$curr_node->id.": \t".($taxonomy_string eq ""? "spezies unknown": $taxonomy_string)."\n";
			
			my $query_id = $curr_node->id;
		#	print $query_id."\t";
			$query_id =~ /(QSQ\d+)/;
	#		print $1."\n";
			if($1){
				if(exists $conversion_hash{$1}){
				#print "\tconverted ..\n"; 
				$query_id = $conversion_hash{$1};}
				#else{ print "\t$1 does not exist in conversion.txt\n";}
			}
			
			my $annotation_arg = $query_id."\t".($taxonomy_string eq q{}? "spezies unknown": $taxonomy_string);
#			print $annotation_arg."\n";
			
			open  my $annotation_fh, '>>', $annotation_file or warn "Could not open: $!";
			flock $annotation_fh, 2 or warn "locking taxonomy file failed.Skipping...\n" and next CLEANED_NODE;                 # try to lock the file
			print {$annotation_fh} $annotation_arg."\n" or warn "Could not write: $!";
			flock $annotation_fh, 8;
			close $annotation_fh or warn "Could not close: $!";
		}
	}
#	print Dumper @taxa;
#exit;
	
return 1;
}

####################################################################
# PARAMETERS
# acc2taxonomy_mapping_hashref
# output_directory
# curr_pfam_family
sub determine_taxonomic_information_for_clade(){
####################################################################
	my ($arg_ref) = @_;
	my $clade_nodes_hashref = $arg_ref->{clade_nodes};
	my %current_headers_hash = ();
	my $taxonomy_output_file = $arg_ref->{curr_tax_file};
	if(! -e $taxonomy_output_file){
		warn "Taxonomy file $taxonomy_output_file does not exist\n";
		return 0;
	}
	my @taxo_array = ();
	my @debug_array = ();
	
	open my $INFILE, '<',$taxonomy_output_file or die "\t\t\tCouldn't open $taxonomy_output_file\n";
	while(<$INFILE>){
		next if /^$/;
		if(/(\w*_\w*)\s*(.*)/){
			if(exists $$clade_nodes_hashref{$1}){
				push @taxo_array, [ split(/;\s?/,$2) ];
				push @debug_array, $2;
				my @arr_split = split /;/, $2;
				$current_headers_hash{$1} = @arr_split;
			}
		}
	}
	close $INFILE or die "Could not close $taxonomy_output_file\n";
	if(scalar(@taxo_array) == 0){
		warn "Current Pfam family has no sequences (Determining taxonomic information)\n";
		return "spezies unknown";
	}
	my %taxonomic_hash = ();
	my $taxonomy_level_counter = 0;
	my $taxonomy_string = "";
	for($taxonomy_level_counter = 0; $taxonomy_level_counter < 10; $taxonomy_level_counter++){
		my %temp_hash = ();
		foreach(@taxo_array){
			my @tmp_arr = @{$_};
			next if not defined $tmp_arr[$taxonomy_level_counter];
			my $curr_taxo = $tmp_arr[$taxonomy_level_counter];
			$curr_taxo =~ s/\.//g;
			$temp_hash{$curr_taxo}++;
		}
		if(keys %temp_hash > 1){
			$taxonomy_string = ($taxonomy_string eq "")? "spezies unknown": $taxonomy_string;
			return $taxonomy_string;
		}
		my $good_hit = 0;
		foreach my $value (sort {$temp_hash{$a} cmp $temp_hash{$b} } keys %temp_hash){
			$taxonomy_string .= $value."; " if $temp_hash{$value} > 2;
		### STOP HERE IF MORE THAN TWO WITH > 2 HITS
			$good_hit = 1 if $temp_hash{$value} > 2;
		}
		last if !$good_hit;
	}
	return $taxonomy_string;
}

####################################################################
# PARAMETERS
# acc2taxonomy_mapping_hashref
# output_directory
# curr_pfam_family
sub determine_taxonomic_information(){
####################################################################
	my ($arg_ref) = @_;
	my $acc2taxonomy_mapping_hashref = $arg_ref->{acc2taxonomy_mapping_hashref};
	my $output_directory = $arg_ref->{output_directory};
	my $curr_pfam_family = $arg_ref->{curr_pfam_family};
	my %current_headers_hash = ();
	#my $current_headers_hashref = \%current_headers_hash;
	
	my $taxonomy_output_file = "$output_directory/$curr_pfam_family/$curr_pfam_family.tax";
#	print "\t\tOpening $taxonomy_output_file\n";
	
	if(! -e $taxonomy_output_file){
		warn "Taxonomy output for $curr_pfam_family doese not exist\n";
		return 0;
	}
	my @taxo_array = ();
	open my $INFILE, '<',$taxonomy_output_file or die "\t\t\tCouldn't open $taxonomy_output_file\n";
	while(<$INFILE>){
		next if /^$/;
			if(/(\w*_\w*)\s*(.*)/){
				    push @taxo_array, [ split(/;/,$2) ];
				my @arr_split = split /;/, $2;
				$current_headers_hash{$1} = @arr_split;
		}
	}
	close $INFILE or die "Could not close $taxonomy_output_file\n";
	if(scalar(@taxo_array) == 0){
		warn "Current Pfam family has no sequences (Determining taxonomic information)\n";
		return 0;
	}
	my %taxonomic_hash = ();
	my $taxonomy_level_counter = 0;
	for($taxonomy_level_counter = 0; $taxonomy_level_counter < 10; $taxonomy_level_counter++){
		my %temp_hash = ();
		foreach(@taxo_array){
			my @tmp_arr = @{$_};
			next if not defined $tmp_arr[$taxonomy_level_counter];
			my $curr_taxo = $tmp_arr[$taxonomy_level_counter];
			$curr_taxo =~ s/\.//g;
			$temp_hash{$curr_taxo}++;
		}
		my $good_hit = 0;
		foreach my $value (sort {$temp_hash{$a} cmp $temp_hash{$b} } keys %temp_hash){
			$good_hit = 1 if $temp_hash{$value} > 2;
		}
		last if !$good_hit;
	}
	print "\t\tTaxonomic resolution until level: $taxonomy_level_counter\n";
return $taxonomy_level_counter;
}


####################################################################
# PARAMETERS:
# pfam_alignment_repository
# pfam_hmm_repository
# curr_pfam_family
# output_directory
# evalue_hmmsearch
sub remove_insignificant_hits(){
####################################################################
#	die "False number of arguments for function hmmsearch\n" if @_ != 2;
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $pfam_alignment_repository = $arg_ref->{pfam_alignment_repository};
	my $pfam_hmm_repository = $arg_ref->{pfam_hmm_repository};
	my $curr_pfam_family = $arg_ref->{curr_pfam_family};
	my $output_directory = $arg_ref->{output_directory};
	my $evalue_hmmsearch = $arg_ref->{evalue_hmmsearch};
	my $hmmer_version = $arg_ref->{hmmer_version};
#### OTHER VARIABLES
	my $searchio;
	my %accession_number_to_keep_hash = ();
	my $count_significant_hits = 0;
	my %tmp_sequence_hash = ();
	my $pfam_hmm_file = "$pfam_hmm_repository/$curr_pfam_family.hmm";
#	print "\thmm is $pfam_hmm_file\n";
	my $egt_file = "$output_directory/$curr_pfam_family/$curr_pfam_family.egt";
	my $egt_file_all = "$output_directory/$curr_pfam_family/$curr_pfam_family.egt-all";
	my $output_template = "$output_directory/$curr_pfam_family/HMMSEARCH";

## CHECK DEFINEDNESS
	if(!defined $pfam_alignment_repository || !defined $curr_pfam_family || !defined $output_directory || !defined $evalue_hmmsearch){
		warn "\t\t[CLEAN SEQUENCES] Problem with parameters for function parse_tree\n";
		return();
}	


	
## CHECK EXISTENCE	
	if(! -e $egt_file_all || ! -s $egt_file_all){
		warn "\t\t[CLEAN SEQUENCES] File containing egts does not exist [remove_insignificant_hits].Skipping this Pfam family....\n";
		return 0;
	}
	if(! -e $pfam_hmm_file || ! -s $pfam_hmm_file){
		warn "\t\t[CLEAN SEQUENCES] File containing the HMM does not exist or is empty [remove_insignificant_hits].Skipping this Pfam family....\n";
		return 0;
	}
	### CLEAN OLD FILES
	unlink($egt_file) if -e $egt_file;

	my $hmm_nc_score = `grep NC $pfam_hmm_file`;
	my $number_of_sequences_to_check = `grep ">" $egt_file_all | wc -l`;
	chomp($number_of_sequences_to_check);
	my ($fh, $tmp_filename) = tempfile();
	$hmm_nc_score =~ /NC\s+(-?\d*\.?\d*)/;
	$hmm_nc_score = $1;
	## COPY EGT --> EGT-all
	if($evalue_hmmsearch){
		#	print $Configuration::hmmer."/hmmsearch -E $evalue_hmmsearch $pfam_hmm_file $egt_file_all\n";
	#		print "$Configuration::hmmer/hmmsearch -E $evalue_hmmsearch $pfam_hmm_file $egt_file_all\n";
			$searchio = `$Configuration::hmmer/hmmsearch -E $evalue_hmmsearch $pfam_hmm_file $egt_file_all 2> /dev/null`;
	}
		else{
	#		print "$Configuration::hmmer/hmmsearch -E $evalue_hmmsearch $pfam_hmm_file $egt_file_all\n";
			$searchio = `$Configuration::hmmer/hmmsearch --cut_ga $pfam_hmm_file $egt_file_all 2> /dev/null`;
		}
		my $hmm_out = "$egt_file_all.hmm";
		open my $INFILE, '>',$hmm_out or die "\t\t\tCouldn't open $hmm_out\n";
		print {$INFILE} $searchio;
		close $INFILE;
		my @tmp_arr = split(/\n/, $searchio);
		my @foo = grep(/Total hits:\s+(\d)/, @tmp_arr);    # weed out comments
		my $hit_found = 0;
		### HMMER 2
		my @found_hits;
		my %keep_hash = ();
		if($hmmer_version eq "2"){
			@found_hits = grep(/^QSQ\d+[-|+]?\d?/, @tmp_arr);    # weed out comments
			if(!scalar(@found_hits)){
				print "\t\t[CLEAN SEQUENCES] No significant hits for $curr_pfam_family. Skipping this Pfam family.... \n";
				return();
			}
			foreach(@found_hits){
				/^(QSQ\d+[-|+]?\d?)/;
				if(!defined $1){
					warn "\t\t[CLEAN SEQUENCES] Could not parse line from hmmsearch output. Skipping this Pfam family....\n";
				}
				$accession_number_to_keep_hash{$1} = 1;
				$count_significant_hits++;
			}
		}
		### HMMER 3
		else{
#			print "using hmmer3\n";
			@found_hits = grep(/^>>\s+QSQ\d+/, @tmp_arr);    # weed out comments
			if(!scalar(@found_hits)){
				print "\t\t[CLEAN SEQUENCES] No significant hits for $curr_pfam_family. Skipping this Pfam family.... \n";
				return();
			}
			foreach(@found_hits){
				/^>>\s+(QSQ\d+)/;
				if(!defined $1){
					warn "\t\t[CLEAN SEQUENCES] Could not parse line from hmmsearch output. Skipping this Pfam family....\n";
				}
				$accession_number_to_keep_hash{$1} = 1;
				$count_significant_hits++;
			}
		}

	### Copy only those sequences with close hits
	my $in  = Bio::SeqIO->new(-file => $egt_file_all , '-format' => 'Fasta');
	my $out  = Bio::SeqIO->new(-file => ">$egt_file" , '-format' => 'Fasta');
			
	$count_significant_hits = 0;
	my $count_all_hits = 0;
	while ( my $seq = $in->next_seq() ) {
		if(exists($accession_number_to_keep_hash{$seq->id})){
	## remove 6-frame stuff
			my $tmp_header = $seq->id;
			$tmp_header =~ /(QSQ\d+[-|+]?\d?)/;
			if(!$1){
				warn "\tProblem removing 6-frame translation header\n";
			}
			$tmp_header = $1;
			$seq->id($tmp_header);
			$out->write_seq($seq);
			$count_significant_hits++;
		}
	}
	
	## REMOVE 6-frame translation endings
		if($evalue_hmmsearch){	
			print "\t\t[CLEAN SEQUENCES] $count_significant_hits / $number_of_sequences_to_check significant Hits (Evalue: $evalue_hmmsearch)\n";
		}
		else{
			print "\t\t[CLEAN SEQUENCES] $count_significant_hits / $number_of_sequences_to_check significant Hits (Cutoff: gathering treshold)\n";
		}
		return 1;
}


####################################################################
sub create_single_directory($){
####################################################################
#### PARAMETER VARIABLES
my $name = shift;
## CHECK DEFINEDNESS
		if(!defined $name){
		warn "\t\t[CLEAN SEQUENCES] Problem to create directory\n";
		return();
}	

	#if ( $name  eq q{}) {
	#	warn "\t"
	#}
if ( !-e $name ) {
	mkdir($name, 0777 ) or die "could not create '$name': $!";
    }
    return 1;
}

####################################################################
sub write_bash_script($$$){
####################################################################
#### PARAMETER VARIABLES
my $file = shift;
my $arg = shift;
my $identifier = shift || q{};

open my $BASH_SCRIPT,">",$file || die "Couldn't open '$file': $!";
print {$BASH_SCRIPT} "#!/bin/bash\n";
print {$BASH_SCRIPT} "#\$ -S /bin/bash\n";
print {$BASH_SCRIPT} "#\$ -cwd\n";
print {$BASH_SCRIPT} "#\$ -N ".$identifier."\n";
print {$BASH_SCRIPT} "#\$ -V\n";
print {$BASH_SCRIPT} $arg."\n";
close($BASH_SCRIPT)|| warn "Couldn't close '$file': $!";
}


####################################################################
## CONVERTS EGT-FILE, WRITES CONVERSION.TXT and SAVES ASSIGNMENTS IN NICE HASH
sub convert_and_collect(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $all_sequences_file = $arg_ref->{all_sequences_file};
my $all_sequences_file_converted = $arg_ref->{all_sequences_file_converted};
#my $assignment_by_acc_hashref => $arg_ref->{assignment_by_acc_hashref};
my %assignment_by_acc_hash = ();
my $all_predictions = $arg_ref->{all_predictions};
#### OTHER VARIABLES

# CHECK DEFINEDNESS
		if(!defined $all_sequences_file || !defined $all_predictions){
		warn "\t\t[CONVERT SEQUENCES] problem with parameters for function convert_egt_file\n";
		return();
}	

## FORMAT OF ALL_PREDICTIONS
# 9,534,0.9887
# 18,3793,0.6311
# 28,4909,0.9515
# CAN HAVE MULTIPLE ASSIGNMENTS
# 776,2803,0.9878
# 776,2801,0.7971

if(!-e $all_sequences_file && !-s $all_sequences_file && !-e $all_predictions && !-s $all_predictions){
#if(-e $assignment_outfile &&  -e $outfile ){
	print "$all_predictions\t$all_sequences_file_converted\n";
	print "\t[CONVERT SEQUENCES] Important files are missing. Skipping\n";
	return 1;
}

my %assignment_pfam_hash = ();

### CONVERT ALL SEQUENCES
	my $conversion_file = "conversion.txt";
	my %converter_hash = ();
	my %assignment_conversion_hash = ();
	my $curr_sequence = q{};
	my $curr_seq_counter = 0;
	my $curr_header = q{};
	my $curr_domain = q{};
### Slurping EGT File
	my $fh = new IO::Handle;
	open $fh, '<', $all_sequences_file;
	my $text = q{};
	sysread $fh, $text, -s $fh;
	$fh->flush;
	close $fh;
	my @whole_lines = split(/\n/,$text);
	open my $OUTFILE_FH, '>', $all_sequences_file_converted or warn "Can't open ".$all_sequences_file_converted.": $!";
	open my $CONVERSION_FH, '>', $conversion_file or warn "Can't open ".$conversion_file.": $!";
	my $sequence_counter = 1;
	print "\tConverting sequence file $all_sequences_file\n";
	foreach(@whole_lines){
		if(/>(.*)/){
			$curr_header = $1;
			$curr_seq_counter++;
			my $converted_header = ">QSQ$sequence_counter";
			print {$OUTFILE_FH} ">QSQ$sequence_counter\n";
			print {$CONVERSION_FH} "$curr_header-->QSQ$sequence_counter\n";
			$converter_hash{$curr_header} = "QSQ$sequence_counter";
			$sequence_counter++;
	#		print "\t$curr_header --> $converted_header\n" if $sequence_counter < 10;
			### CONVERTING ASSIGNMENTS
			## FIND ALL ASSIGNMENTS
		}
		else{
			print {$OUTFILE_FH} $_."\n";
		}
	}
	print "\t[CONVERT SEQUENCES] Done converting headers of fasta file\n";
	#exit;
	close $OUTFILE_FH or warn "Can't close ".$all_sequences_file_converted.": $!";
	close $CONVERSION_FH or warn "Can't close ".$conversion_file.": $!";

###READ ASSIGNMENT FILE
### Slurping ASSIGNMENT File
	my $fh_assignment = new IO::Handle;
	open $fh_assignment, '<', $all_predictions;
	my $text_assignment = q{};
	sysread $fh_assignment, $text_assignment, -s $fh_assignment;
	$fh_assignment->flush;
	close $fh_assignment;
	my @whole_lines_assignment = split(/\n/,$text_assignment);
	#print "\tread ".scalar(@whole_lines_assignment)." assignments into memory\n";
	my $curr_pfam = q();
	my $lame_counter = 0;
	print "\tReading and converting prediction file $all_predictions\n";

	foreach my $line(@whole_lines_assignment){
		next if $line=~ /^$/;
		next if $line =~ /No assignments/;
		if($line =~ /^(PF\d+)/){
			my $grepped_pfam = $1;
			if($curr_pfam eq ""){
				$curr_pfam = $grepped_pfam;
			}
			else {
				$curr_pfam .= "#$grepped_pfam";
			}
		}
		else{
			if($curr_pfam ne ""){
				## save last entry
				if(exists $converter_hash{$curr_header}){
					$assignment_by_acc_hash{$converter_hash{$curr_header}} = $curr_pfam;
				}
				else{
				}
			}
			$curr_header = $line;
			$curr_pfam = "";
		}
}
return \%assignment_by_acc_hash;
}

####################################################################
sub check_sequences_and_predictions(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $assignment_file = $arg_ref->{assignment_file};
my $sequence_file = $arg_ref->{sequence_file};
my %headers_in_assignment_file = ();
my $no_missing_headers = 0;

open my $ASSIGN_FH, '<', $assignment_file or warn "Can't open ".$assignment_file.": $!";
while(<$ASSIGN_FH>){
	next if /^$/;
	next if /^PF/;
	$headers_in_assignment_file{$_} = 1;
}
my $co= 0;

close $ASSIGN_FH || warn "\tCould not close $assignment_file\n";
my @all_fasta_headers = `grep ">" $sequence_file`;
foreach(@all_fasta_headers){
	/>(.*)/;
	my $curr_header = $1;
			if(! exists $headers_in_assignment_file{$curr_header}){
				$no_missing_headers++;
			}
}
if($no_missing_headers){
#	print "\tThere are $no_missing_headers in the sequence but not in the predictions file\n";
;}
return (!$no_missing_headers)?1:undef;
}


####################################################################
# Converts an alignment in clustal format in an alignment in fasta format
# PARAMETERS
# infile - sequences in clustal format
# outfile - sequences in fasta format
sub clustalw2fasta(){    
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $infile = $arg_ref->{infile};
	my $outfile = $arg_ref->{outfile};
#### OTHER VARIABLES
	my %sequences_hash =();
## CHECK DEFINEDNESS
		if(!defined $infile || !defined $outfile){
		warn "\t\t[CLUSTAL2FASTA] problem with parameters for function read_conversion_file\n";
		return();
}	

## CHECK EXISTENCE
	if(! -e $infile || ! -s $infile){
		warn "\t\t[CLUSTAL2FASTA] Input file for format conversion does not exist\n";
		return();
	}

## Open file in clustal format
	open my $INFILE_FH, '<', $infile or warn "Can't open ".$infile.": $!";	
	while(my $line = <$INFILE_FH>){
	#	print $line."\n";
		next if lc($line) =~ /clustal/;
		next if $line =~ /^$/;
		#mathes tr|A3I0L9|A3I0L9_9SPHI *ISFTLV-.FM.LGLIAIAQA.QD--RVVKGVV.TADGL-------P.MP
		$line =~ /^\s*(.*) (.*)/;
		my ($header, $sequence) = ($1,$2);
		if(! defined $header || ! defined $sequence){
			warn "\t\t[CLUSTAL2FASTA] Warning in conversion ($line)\n"; 
		}
		else{
#		print "\theader is $header\n";
			chomp($sequence);
	### CONVERT TO unique GAP SYMBOLS
			$sequence =~ s/\*/-/g;
			$sequence =~ s/\./-/g;
			$sequence =~ s/X/-/g;
			$sequences_hash{$header} .= $sequence;
		}
	}
	### CHECK IF SEQUENCES HAVE SAME LENGTH
	## DISCARD SEQUENCES WITH DIFFERENT LENGTH
	my $seq_length = q{};
	my $counter = 0;
	foreach my $curr_sequence(keys(%sequences_hash)){
		my $current_length = length($sequences_hash{$curr_sequence});
		if(!$seq_length){
			$seq_length = $current_length;
		}
		if($current_length != $seq_length){
			print "\tDelete $curr_sequence ($current_length != $seq_length)\n";
			delete($sequences_hash{$curr_sequence});
		}
	}
	close $INFILE_FH or warn "Can't close ".$infile.": $!";
	## Write in fasta format
	open my $OUTFILE_FH, '>', $outfile or warn "Can't open ".$outfile.": $!";	
	foreach my $seq(sort keys(%sequences_hash)){
		### CHANGE HEADERS
		print {$OUTFILE_FH} ">$seq\n".$sequences_hash{$seq}."\n";
	}
	close $OUTFILE_FH or warn "Can't close ".$outfile.": $!";
	
	## CHECK EXISTENCE
	if(! -e $outfile || ! -s $outfile){
		warn "\t\t[CLUSTAL2FASTA] Output file for format conversion does not exist\n";
		return();
	}
	return 1;
}


####################################################################
# Convert alignment and tree file so that RAxML can deal with it
# PARAMETERS
# sequence_file
# raxml_output_dir
# output_directory
# converted_Short2Orig_hash
# converted_Orig2Short_hash
# reference_tree
sub convert_fasta_and_tree_file(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $fasta_file = $arg_ref->{sequence_file};
my $raxml_output_dir  = $arg_ref->{raxml_output_dir};
my $converted_Short2Orig_hashref = $arg_ref->{converted_Short2Orig_hash};
my $converted_Orig2Short_hashref = $arg_ref->{converted_Orig2Short_hash};
my $reference_tree = $arg_ref->{reference_tree};

#### OTHER VARIABLES
my $reference_tree_converted = "$raxml_output_dir/".basename($reference_tree)."_converted";
my $outfile = $fasta_file."_converted";
my %converter_hash = ();

# CHECK DEFINEDNESS
		if(!defined $fasta_file || !defined $raxml_output_dir || !defined $converted_Short2Orig_hashref || !defined $converted_Orig2Short_hashref || !defined $reference_tree){
		warn "\t\t[CONVERT TREE FILE] problem with parameters for function convert_fasta_and_tree_file\n";
		return();
}	



### CONVERSION OF ALIGNMENT
	open my $FASTA_FILE, '<', $fasta_file or warn "Can't open ".$fasta_file.": $!";
	open my $OUTFILE_FH, '>', $outfile or warn "Can't open ".$outfile.": $!";
	my $sequence_counter = 1;
	foreach(<$FASTA_FILE>){
		if(/>QS\d+/){
			print {$OUTFILE_FH} $_;
			next;
		}
		if(/>(.*)/){
			my $curr_header = $1;
			print {$OUTFILE_FH} ">KS$sequence_counter\n";
			$$converted_Short2Orig_hashref{"KS$sequence_counter"} = $curr_header;
			$$converted_Orig2Short_hashref{$curr_header} = "KS$sequence_counter";
			$sequence_counter++;
		}
		else{
			print {$OUTFILE_FH} $_;
		}
	}
	close $OUTFILE_FH or warn "Can't close ".$outfile.": $!";
	close $FASTA_FILE or warn "Can't close ".$fasta_file.": $!";
	
### CONVERSION OF TREE
	local $/=undef;
	open my $TREE_FILE_IN, '<', $reference_tree or warn "Can't open ".$reference_tree.": $!";
	open my $TREE_FILE_OUT, '>', $reference_tree_converted or warn "Can't open ".$reference_tree_converted.": $!";
	my $tree_string = <$TREE_FILE_IN>;
	close $TREE_FILE_IN or warn "Can't close ".$reference_tree.": $!";

	foreach my $orig_taxon(keys(%$converted_Orig2Short_hashref)){
		$tree_string =~ s/$orig_taxon/$$converted_Orig2Short_hashref{$orig_taxon}/;
	}
	$tree_string =~ s/\n//g;
	$tree_string .= "\n";
	print {$TREE_FILE_OUT} $tree_string;
	close $TREE_FILE_OUT or warn "Can't close ".$reference_tree_converted.": $!";
	
}

####################################################################
# PARAMETERS
# curr_pfam_family - Current PFam family
# raxml_output_dir - Directory of RAxML output files
# converted_Short2Orig_hash - 
# raxml_labelled_tree - RAxML output tree
sub reconvert_fasta_and_tree_file(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $curr_pfam_family = $arg_ref->{curr_pfam_family};
	my $raxml_output_dir  = $arg_ref->{raxml_output_dir};
	my $converted_Short2Orig_hash = $arg_ref->{converted_Short2Orig_hash};
	my $raxml_labelled_tree = $arg_ref->{raxml_labelled_tree};
#### OTHER VARIABLES
	my $reference_tree_converted = "$raxml_output_dir/".$curr_pfam_family.".tree";
	my %converter_hash = ();

	# CHECK DEFINEDNESS
		if(!defined $curr_pfam_family || !defined $raxml_output_dir || !defined $converted_Short2Orig_hash || !defined $raxml_labelled_tree){
		warn "\t\t[RECONVERT TREE FILE] problem with parameters for function reconvert_fasta_and_tree_file\n";
		return();
}	


## Check if tree exists
if(! -e $raxml_labelled_tree || ! -s $raxml_labelled_tree){
	warn "\tRAxML tree does not exist\n";
	return 0;
}
## Check hash is non-empty
if(! keys(%$converted_Short2Orig_hash)){
	warn "\tHash with recoding information is empty\n";
	return 0;
}
## CONVERSION OF ALIGNMENT
	local $/=undef;
	open my $TREE_FILE_IN, '<', $raxml_labelled_tree or warn "Can't open ".$raxml_labelled_tree.": $!";
	open my $TREE_FILE_OUT, '>', $reference_tree_converted or warn "Can't open ".$reference_tree_converted.": $!";
	my $tree_string = <$TREE_FILE_IN>;
	close $TREE_FILE_IN or warn "Can't close ".$raxml_labelled_tree.": $!";
	foreach my $orig_taxon(keys(%$converted_Short2Orig_hash)){
		$tree_string =~ s/$orig_taxon/$$converted_Short2Orig_hash{$orig_taxon}/;
	}
	$tree_string =~ s/\n//g;
	## RECODE RAxML specific stuff
	$tree_string =~ s/\[I\d+\]//g;
	$tree_string =~ s/QUERY___//;       
	$tree_string =~ s/___\d+//;       
	$tree_string .= "\n";
	print {$TREE_FILE_OUT} $tree_string;
	print "\t\tLabeled tree written to $reference_tree_converted\n";
	close $TREE_FILE_OUT or warn "Can't close ".$reference_tree_converted.": $!";

	return 1;	
}

####################################################################
# Converts the sequence identifier from the temporary format (e.g. "QSQ8475") to the original format
# PARAMETERS
# output_directory - containing the file "taxonomic_assignments.txt"
sub reconvert_taxonomic_file(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $output_directory = $arg_ref->{output_directory};
#### OTHER VARIABLES
	my $conversion_file = "conversion.txt";
	my $annotation_file = $output_directory."/taxonomic_assignments.txt";
	my $annotation_file_out = $output_directory."/taxonomic_assignments_converted.txt";
	my %conversion_hash = ();

	# CHECK DEFINEDNESS
		if(!defined $output_directory){
		warn "\t\t[RECONVERT TAXONOMIC FILE] problem with parameters for function reconvert_taxonomic_file\n";
		return();
}	

	###
	# Put all predictions in one file
	###
#	print "cat $output_directory/taxonomic_assignments* > $annotation_file\n";
#	exit;
if(! -e $annotation_file || ! -s $annotation_file){
	`cat $output_directory/taxonomic_assignments* > $annotation_file`;
}
if(! -e $annotation_file || ! -s $annotation_file){
	warn "\tConversion file $annotation_file does not exist\n";
	return 0;
}
### Check existence of annotation and conversion file
if(! -e $conversion_file || ! -s $conversion_file){
	warn "\tConversion file $conversion_file does not exist\n";
	return 0;
}
### READ CONVERSIONS e.g. tr|A2U4H6|A2U4H6_BACCO Catalase OS=Bacillus coagulans 36D1 GN=BcoaDRAFT_2314-->QS23
	read_conversion_file({hash_ref => \%conversion_hash});
## READ Annotations e.g. QS1876  Bacteria; Proteobacteria; Gammaproteobacteria;
#	after conversion e.g. tr|A2U4H6|A2U4H6_BACCO Catalase OS=Bacillus coagulans 36D1 GN=BcoaDRAFT_2314  Bacteria; Proteobacteria; Gammaproteobacteria;
##
my $number_of_converted_entries = 0;
#print "\topening $annotation_file\n";
#exit;
open my $ANNOTATION_FILE_IN, '<', $annotation_file or warn "Can't open ".$annotation_file.": $!";
open my $ANNOTATION_FILE_OUT, '>', $annotation_file_out or warn "Can't open ".$annotation_file_out.": $!";
while(<$ANNOTATION_FILE_IN>){
#print $_."\t";
	/(QSQ\d+)\s(.*)/;
	my ($qs,$taxonomic_classification) = ($1,$2);
#	print $qs."\t$taxonomic_classification\n";
	if(! defined $qs || ! defined $taxonomic_classification){
		warn "\t\tCould not read entry from $conversion_file.Skipping...\n";
		next;
	}
	if(! exists $conversion_hash{$qs}){
		warn "\t\tCould not convert entry: $qs.Skipping.....\n";
		next;
	}
#	print "\t $qs -->  ".$conversion_hash{$qs}."\n";
	$number_of_converted_entries++;
	print {$ANNOTATION_FILE_OUT} $conversion_hash{$qs}."-->".$taxonomic_classification."\n";
}
print "\tcould convert $number_of_converted_entries entries\n";

close $ANNOTATION_FILE_IN or warn "Can't close ".$annotation_file.": $!";
close $ANNOTATION_FILE_OUT or warn "Can't close ".$annotation_file_out.": $!";

return 1;

}


sub stockholm2fasta(){
####################################################################
my ($arg_ref) = @_;
#print Dumper $arg_ref;
my $curr_fasta_file = $arg_ref->{fasta_file};
my $curr_stockholm_file = $arg_ref->{stockholm_file};
## CHECK DEFINEDNESS
	die "One/several of the parameters for 'fasta2stockholm' was/were not defined.\n"
                  if(!defined $curr_fasta_file || !defined $curr_stockholm_file);
# loop through FASTA files
# read FASTA file
    my %seq;
    my @name;
    my $name;
    my $columns = 50;
    my $gapped = 0;

    open STOCKHOLM, "<$curr_stockholm_file" or die "Couldn't open '$curr_stockholm_file': $!";
    while (<STOCKHOLM>) {
     next unless /\S/;
     next if /^\s*\#/;
          if (/^\s*\/\//) { &printseq() }
         else {
	      chomp;
	      my ($name, $seq) = split;
	      #$seq =~ s/[\.\-]//g unless $gapped;
	      $seq{$name} .= $seq;
     }

	&printseq();
    close STOCKHOLM;

    sub printseq(){
        open FASTA, ">>", $curr_fasta_file or die "Couldn't open '$curr_fasta_file': $!";

            while (my ($name, $seq) = each %seq) {
    		  	print FASTA ">$name\n";
	       	  	for (my $i = 0; $i < length $seq; $i += $columns) {
	       			print FASTA substr ($seq, $i, $columns), "\n";
			     }
   	 	      }
            close FASTA;
        }
    if(-e $curr_fasta_file && -s $curr_fasta_file){
        return 1;
    }
    else{
        return 0;
    }
}
}


####################################################################
# PARAMETERS
# hash_ref - reference to hash that will contain "ORIGINAL_FASTA_HEADER-->QSQXYZ"
sub read_conversion_file(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $hash_ref = $arg_ref->{hash_ref};
#### OTHER VARIABLES
	my $conversion_file = "conversion.txt";
	my %conversion_hash = ();

	# CHECK DEFINEDNESS
		if(!defined $hash_ref){
		warn "\t\t[READ CONVERSION FILE] problem with parameters for function read_conversion_file\n";
		return();
}	


## Check if file exists
	if(! -e $conversion_file || ! -s $conversion_file){
		warn "\tConversion file $conversion_file does not exist\n";
		return 0;
	}
	### READ CONVERSIONS e.g. tr|A2U4H6|A2U4H6_BACCO Catalase OS=Bacillus coagulans 36D1 GN=BcoaDRAFT_2314-->QSQ23
	open my $CONVERSION_FILE_IN, '<', $conversion_file or warn "Can't open ".$conversion_file.": $!";
	while(<$CONVERSION_FILE_IN>){
		/(.*)-->(QSQ\d+)/;
		my ($fasta_header,$key) = ($1,$2);
		if(! defined $fasta_header || ! defined $key){
			warn "\t\tCould not read entry from $conversion_file.Skipping...\n";
			next;
		}
		$$hash_ref{$key} = $fasta_header;
	}
	close $CONVERSION_FILE_IN or warn "Can't close ".$conversion_file.": $!";
	if(!keys(%{$hash_ref})){
		warn "\t\tNo entries read from conversion file\n";
		return 0;
	}
return 1;	
}

####################################################################
# PARAMETERS
# curr_domain
# output_directory
# seq
# pfam_accessions_hashref
# only_estscan
sub translate_egts(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $curr_domain = $arg_ref->{curr_domain}; 
	my $output_directory = $arg_ref->{output_directory};
	my $pfam_alignment_repository = $arg_ref->{pfam_alignment_repository};
	
	
	# CHECK DEFINEDNESS
		if(!defined $curr_domain || !defined $output_directory || !defined $pfam_alignment_repository){
		warn "\t\t[PROTEIN TRANSLATION] Problem with parameters for function translate_egts\n";
		return();
}	


	
	$curr_domain =~ s/\s*//g;
					# Fetched HMM FROM HMM-DATABASE
					my $hmm_fetch_file = $output_directory."/".$curr_domain."/".$curr_domain.".hmm";
					# Emitted SEQ FROM HMM
					my $hmm_emit_file = $output_directory."/".$curr_domain."/".$curr_domain.".seq";
					my $egt_file = $output_directory."/".$curr_domain."/".$curr_domain.".egt-all";
					my $temp_egt_file = $output_directory."/".$curr_domain."/".$curr_domain.".tmp";
					my $all_egts_file = $output_directory."/".$curr_domain."/".$curr_domain.".nucl";
					my $all_translated_egts_file = $output_directory."/".$curr_domain."/".$curr_domain.".trans";
			if(!-e $all_egts_file){
				warn "There are no sequences for $curr_domain [Translate EGT]\n";
				return 0;
			}

			print "\t[PROTEIN TRANSLATION] using Bioperl 6-frame...\n";
			## CALL FUNCTION USING SEQOBJECT AS PARAMETER
			if(!translate_bioperl({all_egts_file => $all_egts_file, egt_file => $egt_file})){
						warn "\t[PROTEIN TRANSLATION] There was an error. File with translated sequences empty\n";
					return();
			}
	## CHECK IF RESULTING FILE exists and is non-empty
	if(! -e $egt_file|| ! -s $egt_file){
		warn "\t[PROTEIN TRANSLATION] There was an error. File with translated sequences empty\n";
		return();
	}

	return 1;
}


####################################################################
sub translate_bioperl{
####################################################################

my ($arg_ref) = @_;
my $all_egts_file = $arg_ref->{all_egts_file};
my $egt_file = $arg_ref->{egt_file};


	my $seq_in = Bio::SeqIO->new('-file' => "<$all_egts_file",
                                      '-format' => "fasta");
         my $seq_out = Bio::SeqIO->new('-file' => ">$egt_file",
                                       '-format' => "fasta");

         # write each entry in the input file to the output file
         while (my $seq = $seq_in->next_seq) {
         my %translated_hash = ();
		     my $len  = $seq->length();
		     for(my $frame=1; $frame < 4; $frame++){
		    	 my $seqobj = Bio::Seq->new(
                        -seq        => $seq->subseq($frame,$seq->length),
                        -display_id => "",
                        -alphabet   => 'dna'
                    );
		     	$translated_hash{$seq->id}{seq} .= $seqobj->translate('X', 'X', '0', '1', '0', '0', '0', '0')->seq; 
     		}
    
		    my $rev = $seq->revcom;
		     for(my $frame=1; $frame < 4; $frame++){
		    	 my $seqobj = Bio::Seq->new(
                        -seq        => $rev->subseq($frame,$rev->length),
                        -display_id => "najaa",
                        -alphabet   => 'dna'
                    );
		     		$translated_hash{$seq->id}{seq} .= $seqobj->translate('X', 'X', '0', '1', '0', '0', '0', '0')->seq; 
     			}
     			foreach(sort keys %translated_hash){
     				 my $seqobj = Bio::Seq->new(
                        -seq        => $translated_hash{$_}{seq},
                        -display_id => $_,
                        -alphabet   => 'protein'
                    );
     			$seq_out->write_seq($seqobj);
     			}
     	}
return 1;
}


####################################################################
# PARAMETERS
# egt_file
# minimum_length
# output_directory
# assignment_hash_by_acc
sub evaluate_statistics(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $new_assignments_file = $arg_ref->{new_assignments};
my $statistics_file = "statistics_file";

#### OTHER VARIABLES
my (%new_assignments_hash,%reference_assignments_hash) =();
if(! -e $new_assignments_file|| ! -s $new_assignments_file){
	warn "\tFile with new assignments missing\n";
	return();
}
#if(! -e $reference_assignments_file|| ! -s $reference_assignments_file){
#	warn "\tFile with reference assignments missing\n";
#	return();
#}

my %superkingdoms = ("Bacteria" => 1,
									"Eukaryota" => 1,
									"Archaea" => 1
);


my %phyla = ("Acidobacteria" => 1,
"Actinobacteria" => 1,
"Alphaproteobacteria" => 1,
"Aquificae" => 1,
"Bacteroidetes" => 1,
"Betaproteobacteria" => 1,
"Chlamydiae" => 1,
"Chlorobi" => 1,
"Chloroflexi" => 1,
"Cyanobacteria" => 1,
"Deinococcus-Thermus" => 1,
"Deltaproteobacteria" => 1,
"Epsilonproteobacteria" => 1,
"Firmicutes" => 1,
"Gammaproteobacteria" => 1,
"Lentisphaerae" => 1,
"Planctomycetes" => 1,
"Spirochaetes" => 1,
"Thermotogae" => 1,
"Zetaproteobacteria" => 1);


##############################
##############################
### READ REFERENCE ASSIGNMENTS
##############################
##############################
my %layer_taxa_statistics_hash = ();

my $total_number_of_entries = `cat $new_assignments_file | wc -l`;
my $number_of_unknown_entries = `grep unknown $new_assignments_file -c`;
my $total_number_of_assignments = $total_number_of_entries - $number_of_unknown_entries;

foreach my $phylum(keys %phyla){
#print "grep $phylum $new_assignments_file\n";
	my @hits = `grep $phylum $new_assignments_file`;
		$layer_taxa_statistics_hash{$phylum} = scalar(@hits);
		next;
	foreach(@hits){
#		/(\w*-?\d+_-?\d+)\s*(.*)/;
		/(QSQ\d+[+|-]?\d)\s+(.*)/;
		my ($acc_name,$taxonomy_string) = ($1,$2);
		if(!defined $acc_name || !defined $taxonomy_string){
			warn "\tCould not parse file\n";
		}
		$layer_taxa_statistics_hash{$phylum}++;
	}
}
### WRITE STATISTICS
open my $CARMA_STAT_FH, '>', $statistics_file or warn "Can't open ".$statistics_file.": $!";
		foreach my $taxon(sort keys(%layer_taxa_statistics_hash)){
			my $percentage = sprintf("%02.2f", (100* $layer_taxa_statistics_hash{$taxon}/$total_number_of_assignments));
			print {$CARMA_STAT_FH} $taxon."\t".$layer_taxa_statistics_hash{$taxon}."\t$percentage\n";
}
close $CARMA_STAT_FH || warn "Could not close $statistics_file\n";
return (-e $statistics_file && -s $statistics_file)? $statistics_file : undef;
}

####################################################################
# PARAMETERS
# egt_file
# minimum_length
# output_directory
# assignment_hash_by_acc
sub evaluate_all(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $input_dir = $arg_ref->{input_dir};
my $statistics_file = "statistics_file";

#### OTHER VARIABLES
my (%new_assignments_hash,%reference_assignments_hash) =();
#if(! -e $input_dir|| ! -s $input_dir){
#	warn "\tFile with new assignments missing\n";
#	return();
#}
#if(! -e $reference_assignments_file|| ! -s $reference_assignments_file){
#	warn "\tFile with reference assignments missing\n";
#	return();
#}

my %superkingdoms = ("Bacteria" => 1,
									"Eukaryota" => 1,
									"Archaea" => 1
);
my %phyla = (
"Acidobacteria" => 1,
"Actinobacteria" => 1,
"Alphaproteobacteria" => 1,
"Bacteroidetes" => 1,
"Betaproteobacteria" => 1,
"Chlamydiae" => 1,
"Chlorobi" => 1,
"Chloroflexi" => 1,
"Cyanobacteria" => 1,
"Deinococcus-Thermus" => 1,
"Deltaproteobacteria" => 1,
"Epsilonproteobacteria" => 1,
"Firmicutes" => 1,
"Gammaproteobacteria" => 1,
"Planctomycetes" => 1,
"Spirochaetes" => 1,
);
my %sixtenS = (
"Acidobacteria" => "1",
"Actinobacteria" => "10,5",
"Alphaproteobacteria" => "8,5",
"Bacteroidetes" => "30",
"Betaproteobacteria" => "41",
"Chlamydiae" => 0,
"Chlorobi" => 0,
"Chloroflexi" => "1",
"Cyanobacteria" => "1",
"Deinococcus-Thermus" => 0,
"Deltaproteobacteria" => "2",
"Epsilonproteobacteria" => 0,
"Firmicutes" => 0,
"Gammaproteobacteria" => "5",
"Planctomycetes" => 0,
"Spirochaetes" => 0);

##############################
##############################
### READ REFERENCE ASSIGNMENTS
##############################
##############################
my %layer_taxa_statistics_hash = ();

print "\tgrep $input_dir/*\n";
foreach my $phylum(keys %phyla){
print "grep $phylum $input_dir/*\n";
	my @hits = `grep $phylum $input_dir/*`;
		#$layer_taxa_statistics_hash{$phylum} = scalar(@hits);
	#	next;
	foreach(@hits){
		/(.*):\w+-?\w*\s+(\d+)\s+(\d+\.?\d*)/;
		my ($file,$total_count,$perc_count) = ($1,$2,$3);
		if(!defined $file || !defined $total_count || !defined $perc_count){
			warn "\tCould not parse file\n";
		}
		$file = basename($file);
		$file =~ s/\.txt//g;
	#	$layer =~ s/0//;
		$layer_taxa_statistics_hash{$file}{$phylum}{"total"} = $total_count;
		$layer_taxa_statistics_hash{$file}{$phylum}{"perc"} = $perc_count;
	}
}
#print Dumper %layer_taxa_statistics_hash;
#exit;
### WRITE STATISTICS
open my $CARMA_STAT_FH, '>', $statistics_file or warn "Can't open ".$statistics_file.": $!";
				print {$CARMA_STAT_FH} (join("\t",sort keys(%phyla)))."\n";
		foreach my $mode(sort keys(%layer_taxa_statistics_hash)){
		
			#	print {$CARMA_STAT_FH} $mode."\t";
			foreach my $phylum(sort keys(%{$layer_taxa_statistics_hash{$mode}})){
				print {$CARMA_STAT_FH} $layer_taxa_statistics_hash{$mode}{$phylum}{"perc"}."\t";
			}
				print {$CARMA_STAT_FH} "\n";
		}
			foreach my $phylum(sort keys(%phyla)){
					print {$CARMA_STAT_FH} $sixtenS{$phylum}."\t";
			}
				print {$CARMA_STAT_FH} "\n";
close $CARMA_STAT_FH || warn "Could not close $statistics_file\n";
return (-e $statistics_file && -s $statistics_file)? $statistics_file : undef;

}


####################################################################
# Converts tree file to replace sequence names with taxonomy string
# PARAMETERS
# output_directory - output directory of analysis
# curr_pfam_family - current Pfam family
# taxonomy_repository - directory with taxonomy information
sub create_lockfile(){
####################################################################
#### PARAMETER VARIABLES
	my ($arg_ref) = @_;
	my $lockfile = $arg_ref->{lockfile};
	
	# CHECK DEFINEDNESS
if(!defined $lockfile){
		warn "\t\t[PREPARATION] Problem with parameters for function create_lockfile\n";
		return();
}	
# CHECK EXISTENCE
		return() if(-e $lockfile);
## TRY TO OPEN FILE
		open  my $LOCK_FH, '>', $lockfile or warn "Trying to open locked file.Skipping...\n" and return();
		flock $LOCK_FH, 2  or warn "locking file failed.Skipping...\n" and return();                 # try to lock the file
		print {$LOCK_FH} "0\n" or warn "Could not write: $!";
		flock $LOCK_FH, 8;
		close $LOCK_FH or warn "Could not close: $!";
# CHECK EXISTENCE
	(-e $lockfile)? return 1: return ;
}

####################################################################
# PARAMETERS
# pfam_hmm_file
sub collect_pfam_accessionnumbers{
####################################################################
	my ($arg_ref) = @_;
	my $pfam_hmm_file = $arg_ref->{pfam_hmm_file}; 
	my $pfam_accessions_hashref = $arg_ref->{pfam_accessions_hashref}; 
	if(!-e $pfam_hmm_file){
		warn "File containing all PFAM accessionnumbers does not exist\n";
		return 0;
	}
	my @grepped_accs = `grep ACC $pfam_hmm_file`;
#	print scalar(@grepped_accs)."\n";
	if(!scalar(@grepped_accs)){
		warn "There was a problem getting the accessionnumbers from file $pfam_hmm_file\n";
		return 0;
	}
	foreach(@grepped_accs){
		if(/ACC\s+((PF\d*)\.\d+)/){
			my($acc_long, $acc_short) = ($1,$2);
			if(!defined $acc_long || !defined $acc_short){
				warn "There was a problem with the accessionnumber format in file $pfam_hmm_file\n";
			}
	#		print "$acc_long --> $acc_short \n";
			$$pfam_accessions_hashref{$acc_short} = $acc_long;
		}
	#exit();
	}
	return 1;
}

####################################################################
# PARAMETERS
# output_directory
# pfam_alignment_file
# assignment_hash_by_pfam
sub convert_pfam_alignments(){
####################################################################
my ($arg_ref) = @_;
my $output_directory = $arg_ref->{output_directory};
my $pfam_alignment_file = $arg_ref->{pfam_alignment_file};

print "Read all alignments into memory\n";
### Slurping File
my $fh = new IO::Handle;
open $fh, '<', $pfam_alignment_file;
print "Process all alignments....\n";
my $curr_sequence ="";

(my $number_of_alignments = `grep "GF AC" -c $pfam_alignment_file`) =~ s/\s//g;

my $ali_counter = 0;
my $curr_pfam = "";
	while(<$fh>){
		if(/\#\s+STOCKHOLM/){
				### SAVE
				my $seed_file = "$output_directory/$curr_pfam.full";
				if(!$curr_sequence eq ""){
					open my $seed_fh, '>', $seed_file || die "Could not open $seed_file\n";
					print $seed_fh $curr_sequence;
					close $seed_fh || die "Could not close $seed_file\n";
				### HMMMCALIBRATE
				#	print "\tCalibrate HMM\n";
				#	`hmmcalibrate --seed 0 $hmm_file`;
				}
				#}
	#		">$output_directory/$output_pfam_family/$output_pfam_family.seed"
	#	}
			$curr_sequence = "";
		}
		if(/#=GF AC   (\w*)\.\d+/){
		#	print $1."\n";
			$curr_pfam = $1;
			$ali_counter++;
			#last if $ali_counter > 5;
	#		print "$ali_counter / $number_of_alignments\n";
		}
#			$curr_sequence .= $_."\n";
			$curr_sequence .= $_;
	}
close $fh;

}

####################################################################
# PARAMETERS
# output_directory
# pfam_sequence_file_short
sub collect_taxonomy_records(){
####################################################################
my ($arg_ref) = @_;
my $output_directory = $arg_ref->{output_directory};
my $pfam_sequence_file_short = $arg_ref->{pfam_sequence_file_short};
my $taxonomy_file = $arg_ref->{taxonomy_file};
my $pfam_fasta_file = $arg_ref->{pfam_fasta_file};
#my %assignment_hash_by_pfam = %{$arg_ref->{assignment_hash_by_pfam}};
my %acc2taxonomy_mapping_hash = ();
my %interesting_headers = ();
my $interesting_headers_hashref = \%interesting_headers;
print "\tRead all taxonomy information into memory\n";


my %acc2pfamfamily_mapping_hash = ();
my @all_fasta_headers = `grep ">" $pfam_fasta_file`;
foreach(@all_fasta_headers){
	/>(.*)\/.*(PF\d+)/;
	$acc2pfamfamily_mapping_hash{$1} = $2;
}
print "\tRead acc2pfam mappings for ".keys(%acc2pfamfamily_mapping_hash)." accs \n";
my $co = 0;
my %ass2pfam_taxonomy_hash = ();
#$output_directory = "/c1/scratch/fabian/db_meta/db_pfam/Treephyler/taxonomy_files/full";
open my $TAXIN_FH, '<',$taxonomy_file or die "\t\t\tCouldn't open $taxonomy_file\n";
while(<$TAXIN_FH>){
	my @splitted_line = split(/'/,$_);
	my ($curr_id, $curr_taxonomy) = ($splitted_line[3],$splitted_line[21]);
	if(!exists $acc2pfamfamily_mapping_hash{$curr_id}){
		#print "\tCould not find $curr_id in Pfam-A.fasta file\n";
		next;
	}
	my $curr_pfam = $acc2pfamfamily_mapping_hash{$curr_id};
	$ass2pfam_taxonomy_hash{$curr_pfam} .= "$curr_id\t$curr_taxonomy\n";
}
close $TAXIN_FH || warn "\tCould not close taxonomy file\n";
print "\tread information for ".keys(%ass2pfam_taxonomy_hash)." Pfam families\n";
#print Dumper %ass2pfam_taxonomy_hash;
#exit;

foreach my $curr_pfam_family(keys(%ass2pfam_taxonomy_hash)){
	my $taxonomy_output_file = "$output_directory/$curr_pfam_family.tax";
	open my $INFILE, '>',$taxonomy_output_file or die "\t\t\tCouldn't open $taxonomy_output_file\n";
		print {$INFILE} $ass2pfam_taxonomy_hash{$curr_pfam_family};
	close $INFILE or die "\t\t\tCould not close taxonomy file for $curr_pfam_family\n";
}
	return \%acc2taxonomy_mapping_hash;

}


####################################################################
sub small_hmm_builder(){
####################################################################
	my ($arg_ref) = @_;
	#print Dumper $arg_ref;
	my $output_directory = $arg_ref->{output_directory};
	my $alignment_repository = $arg_ref->{alignment_repository};
	my @files = <$alignment_repository/*>;
#	print "Iterating over $alignment_repository\n";
#      exit;
	### HMMINDEX
	foreach my $curr_alignment(@files){
		my $curr_pfam = basename($curr_alignment);
		$curr_pfam =~ s/\.full//;
		print "Computing $curr_pfam\n";
		my $hmm_output_file = $output_directory."/".$curr_pfam.".hmm";
#		next if -e $hmm_output_file;
		my $cmd = $Configuration::hmmer."/hmmbuild $hmm_output_file $curr_alignment";
#		print "$cmd\n";
#		exit;
		`$cmd`;
		if(!-e $hmm_output_file){
			die "Could not compute hmm for $curr_pfam\n";
		}
		
	}
	
}


####################################################################
sub small_hmm_fetcher(){
####################################################################
	my ($arg_ref) = @_;
	my $output_directory = $arg_ref->{output_directory};
	my $pfam_hmm_file = $arg_ref->{pfam_hmm_file};
	my $pfam_accessions_hashref = $arg_ref->{pfam_accessions_hashref};
	my @files = <$output_directory/*>;

	### HMMINDEX
	my $hmmindex_file = $pfam_hmm_file.".ssi";	
	if(!-e $hmmindex_file && ! -s $hmmindex_file){
		print "\tIndexing Pfam database...(this can take a moment)\n";
		`$Configuration::hmmer/hmmfetch --index $pfam_hmm_file`;
	}
	foreach(sort keys(%{$pfam_accessions_hashref})){
		my $curr_domain=$_;
		
#		$curr_domain =~ s/\..*//;
		my $hmm_fetch_file = $output_directory."/".$curr_domain.".hmm";
#		$curr_domain = basename($curr_domain);
		#print $curr_domain."\t $hmm_fetch_file\n";
		#print $curr_domain."\n";
		my $curr_pfamm_acc = $$pfam_accessions_hashref{$curr_domain};
		if(! defined $curr_pfamm_acc ){
			warn "Couldn't grep accession number from Pfam's Hmm file\n";
			exit;
		}
	#	print "\t$curr_pfamm_acc\n";
		#next;
		`hmmfetch $pfam_hmm_file $curr_pfamm_acc > $hmm_fetch_file`;
		if(!-e $hmm_fetch_file){
			die "Could not fetch hmm for $curr_pfamm_acc from $pfam_hmm_file\n";
		}
	}
}
####################################################################
sub split_file(){
####################################################################
my ($arg_ref) = @_;
my $infile	= $arg_ref->{split_file};
my $max_size	=  $arg_ref->{max_size};
my $i		= 0;
my $j		= 1;
$max_size = $max_size -1;
$max_size = $max_size * 100;
my $number_of_entries = 0;
my $curr_file_name = "$infile"."_$j.fasta";
my $current_fasta_file_string = q();
### OBJECT INSTANTIATION
my $in	= Bio::SeqIO->new(
			-file	=> $infile,
			-format	=> 'Fasta',
			);

### PROCESS FILES
while ( my $seq = $in->next_seq() ) {

	use bytes;
	my $byte_size = length($current_fasta_file_string);
	## Convert to MB
	$byte_size = $byte_size / 10000;
	if ($byte_size >= $max_size) {
		print "\t\twriting file $j\n";
		$curr_file_name = "$infile"."_$j.fasta";
				open my $INFILE, '>',$curr_file_name or die "\t\t\tCouldn't open $curr_file_name\n";
				print {$INFILE} $current_fasta_file_string;
				close $INFILE or die "\t\t\tCould not close fasta file $curr_file_name\n";
		$j++;
		$current_fasta_file_string = q();
	}
	$current_fasta_file_string .= ">".$seq->id."\n".$seq->seq."\n";
	$number_of_entries++;
	}
		$curr_file_name = "$infile"."_$j.fasta";
		print "\t\twriting file $j\n";
	open my $INFILE, '>',$curr_file_name or die "\t\t\tCouldn't open $curr_file_name\n";
	print {$INFILE} $current_fasta_file_string;
	close $INFILE or die "\t\t\tCould not close fasta file $curr_file_name\n";
return();
}


####################################################################
sub print_start(){
####################################################################
print STDOUT <<START;
------------------------------------------------------------------------------------------------------------
					Welcome to TreePhyler!
------------------------------------------------------------------------------------------------------------
START
}


####################################################################
sub print_end(){
####################################################################
print STDOUT <<END;
------------------------------------------------------------------------------------------------------------
					Thanks for using TreePhyler!
------------------------------------------------------------------------------------------------------------      
END
}


####################################################################
sub printhelp {
####################################################################
print STDERR <<ENDHELP;
Parameters of $0:
perl treephyler.pl -m ("statistics"|"analyse"|"prepare"|"split"|translate) -i Fastafile  -a Assignment_file [-nsplits|-pfamA|-pfamT|-pfamH|-pfamF]

Required parameters
	-i : Contains the query sequences in Fasta format
	-a: Contains PFAM predictions in UFO format
	-comet: if predictions are from the comet server. Internally converts them into UFO format
	-m: modus, available options are:
			"statistics" - Compute statistics 
				 e.g. perl treephyler.pl -m statistics -i assignments.fa
			"prepare" - Extracts information about alignments, taxonomy, and hmms from PFAM files and saves them in the subfolder "pfam"
				e.g. perl treephyler.pl -m prepare -p Pfam-A.full -pfamT pfamseq.txt -pfamH Pfam-A.hmm -pfamF Pfam-A.fasta
			"analyse" - Start analysis
				e.g. perl treephyler.pl -m analyse -i gletscher.fas -a Ufo.out (-comet)
			"split" - Splits a file into several smaller chunks
				e.g. perl treephyler.pl -m split -i input_sequences.fa -msize 50
			"translate" - Translates sequences in all six reading frames
				e.g. perl treephyler.pl -m translate -i input_sequences.fa 
			"help" - print this help message             
Example call:
	perl treephyler.pl -m analyse -i gletscher.fas -a Ufo.out
ENDHELP
exit;
}

1
