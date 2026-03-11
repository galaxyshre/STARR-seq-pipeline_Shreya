#!/usr/bin/perl

#
# A script to assign TSS positions according to STARR-seq template-switch data
# @author Teemu Kivioja 
#

use strict;
use warnings;

use AppConfig ':argcount';
use AppConfig::Getopt;
use AppConfig::File;

use File::Basename;
use File::Temp qw/ tempfile tempdir /;

use Bio::SeqIO;

use FindBin qw($Bin);
use lib "$Bin";

use TemplateSwitchPosition;

#
# Here we assign the TSS positions in the uncertain cases in such a way that it is compatible with the certain cases
# i.e. leads to consistent overall probabilities 

# Changes to previous versions
# June 4, 2019: changed cutadapt allowed error rate from the default 0.10 to 0.05 to make sure sequence without splice site does 
# not match by chance
# June 5, rewriting the matching to input to use join instead of seqkit grep because did not scale
# June 7, added the UMI counting and file that contains the RNA-INPUT id pairs (not utilized yet)

# August 21, 2019,
# - use index to match only to the input R1 reads that had the matching R2 start
# - new primer design, longer UMI, putting data from the new run also as the default set

# October 16, 2019
# - skip input information when making UTR file to get a larger set

# November 28, 2019
# - new version that also has the no input UTR and annotation parts from the 1907 code
# - matching simplified: use the 20 base suffix of R1 to find the matching candidates instead of R2 start

#
# December 16, 2019: Add to the annotation file a line for those that do not have the full tss flank
#


my $DEBUG = 0;

my $conf;
configure();

# keep the results the same between runs
srand(42);

my $min_len = 20;
my $UMI_LEN = $conf->UMI_length();;


my $r1_rna_fastq = $conf->R1_rna_fastq_file();
my $r2_rna_fastq = $conf->R2_rna_fastq_file();

my ($base_name, $dirs, $suffix) = fileparse($r1_rna_fastq, (".fastq.gz"));
print STDERR "Script start file: $base_name, $dirs, $suffix\n";

my $r1_trimmed_fastq = $conf->output_dir.$base_name.".intron_trimmed.fastq.gz";
my $r1_trim_info = $conf->output_dir.$base_name.".intron_info.txt";

# RNA R1 random region suffix before constant collected for input matching
my $rna_rand_suffix_file = $conf->output_dir.$base_name.".rand_suffix.tsv";

# output files, before TSS, after TSS, and both, and finally a count file
my $promoter_output_file = $conf->output_dir.$base_name.".promoters.fasta";
my $utr_output_file = $conf->output_dir.$base_name.".utrs.fasta";
my $both_flanks_output_file = $conf->output_dir.$base_name.".tss_flanks.fasta";
my $umi_count_file = $conf->output_dir.$base_name.".umi_counts.tsv";

my $anno_file = $conf->output_dir.$base_name.".tss_flank_anno.tsv";

#
# Step 1: trim R1 at the spliced intron and keep only those that have it
#         - importantly produce the info file that contains the position of the match
if (1) {
    print "Step 1: Recognizing and trimming the spliced intron from R1...\n";
    my $min_overlap = $conf->min_intron_overlap();
    my $spliced_intron_seq = $conf->spliced_intron_seq();
    
    my $command = "cutadapt -j 8 --discard-untrimmed -e 0.05 --overlap=".$min_overlap." --minimum-length=".$min_len;
    $command .= " -a $spliced_intron_seq";
    $command .= " -o $r1_trimmed_fastq $r1_rna_fastq";
    $command .= " --info-file $r1_trim_info";
    
    print STDERR "\tcommand: $command\n"; 
    system($command) == 0 or die "Error: trimming command $command failed: $?";    
}


#
# Step 2, NEW: match the suffix of the RNA random region to the suffix of the input R1 read
#
# make RNAs already 
my $rnas = parse_template_switched_from_trim_info_file($r1_trim_info, $min_len, $rna_rand_suffix_file);

my $rna_input_matches_file = $conf->output_dir.$base_name.".rna_input_matches.tsv.gz";
my $rna_input_matches_id_file = $conf->output_dir.$base_name.".rna_input_matches_ids.txt";
my $rna_in_input_r1_file = $conf->output_dir.$base_name.".rna_in_input_R1.fastq"; # keep uncompressed for the next step (for now)

if (1) {
    # match the last 20 bases to the input and keep the matching sequence and both ids, need to sort the RNA file 
    my $command2 = "join -1 2 -2 2 <(zcat ".$conf->R1_input_tab_sorted_file.") <(sort -k2 -b ".$rna_rand_suffix_file.") | ";
    $command2 .= 'cut -f 1,2,3 -d " " | perl -ple \'s/\s/\t/g\' | gzip -c >'.$rna_input_matches_file;

    print STDERR "\tcommand2: $command2\n";
    my @bash_args2 = ("bash", "-c", $command2);
    system(@bash_args2) == 0 or die "Error: command2 $command2 failed: $?";     

    # make a version that has only the input id for fishing out the r1 sequences, take a way the duplicates
    my $command3 = "zcat ".$rna_input_matches_file." | cut -f 2 | sort -u > ".$rna_input_matches_id_file;
    system($command3) == 0 or die "Error:command3 $command3 failed: $?";
}
    
if (1) {
    # collect the R1s of the found
    # - this is now slow because goes through the huge input file 
    #
    
    my $command3 =  "seqkit grep -f ".$rna_input_matches_id_file." -o ".$rna_in_input_r1_file." ";
    $command3 .= $conf->R1_input_fastq_file;
    
    my @bash_args3 = ("bash", "-c", $command3);
    print STDERR "\tcommand3: $command3\n"; 
    system(@bash_args3) == 0 or die "Error: command3 $command3 failed: $?";
}


#
# Step 3: get the sequence upstream of TSS
#
if (1) {
    print "Step 3: Extracting the sequences around TSS ...\n";

    my $rnas = estimate_tss_probs($rnas, $rna_input_matches_file, $rna_in_input_r1_file, $min_len);
    get_tss_seqs($rnas, $promoter_output_file, $utr_output_file, $both_flanks_output_file, $umi_count_file, $anno_file,
    		 $min_len, $conf->promoter_length, $conf->utr_length);
    
}

#
# Step 4: additional stats etc.
# 
if ($conf->run_revcoms_and_stats) {

    print "Step 4: Reverse complement and collect sequence stats\n";
    my $revcom_file = $conf->output_dir.$base_name.".tss_flanks.revcom.fasta";
    my $command = "seqkit seq -r -p ".$both_flanks_output_file." -o ".$revcom_file;

    print STDERR "\tcommand: $command\n"; 
    system($command) == 0 or die "Error: command $command failed: $?";
    
    my $out_nuc_dist_file = $conf->output_dir.$base_name.".tss_flanks_nuc_dist.txt";
    my $command2 = "cat ".$both_flanks_output_file." | fasta-get-markov -norc -m 5 > ".$out_nuc_dist_file; 

    print STDERR "\tcommand2: $command2\n"; 
    system($command2) == 0 or die "Error: command2 $command2 failed: $?";

    my $prom_revcom_file = $conf->output_dir.$base_name.".promoters.revcom.fasta";
    my $command3 = "seqkit seq -r -p ".$promoter_output_file." -o ".$prom_revcom_file;

    print STDERR "\tcommand: $command3\n"; 
    system($command3) == 0 or die "Error: command $command3 failed: $?";
    
}


print STDERR "All done!\n";


#
# Generate an index of the input R1 sequences matching transcripts 
# - index according to the matching RNA id
sub generate_r1_match_index {
    my ($rna_input_matches_file, $rna_in_input_r1_file) = @_;

    my %r1_match_index = ();
    
    # read the id pairs
    my %input_to_rna_ids = ();

    # compressed file
    open(IDFH, "zcat $rna_input_matches_file |") or die "Could not zcat read match id file $rna_input_matches_file: $!";
    while (my $line = <IDFH>) {
	chomp $line;
	# note: now matching is based on R1
	my ($r1_match_seq, $input_id, $rna_id) = split /\t/, $line;

	# note: several RNA:s can map to the same input sequence
	push @{ $input_to_rna_ids{$input_id} }, $rna_id; 
    }

    # read the corresponding R1 input sequences and put to the index
    my $r1_in = Bio::SeqIO->new(-file => $rna_in_input_r1_file, 				 
				-format => 'fastq');
    
    # store, can take a lot of memory
    my $num_indexed = 0;
    while (my $seq_obj = $r1_in->next_seq) {
	my $input_id = $seq_obj->display_id;
	if (defined($input_to_rna_ids{$input_id})) {
	    for my $rna_id (@{ $input_to_rna_ids{$input_id} }) {
		# can be many
		push @{ $r1_match_index{$rna_id} }, $seq_obj;
		$num_indexed++;
	    }
	} else {
	    die "Error: could not find matching RNA ids for input id $input_id";  
	}
    }

    print STDERR "Indexed $num_indexed RNA input sequence pairs\n";
    return \%r1_match_index;
}


#
# Match indexed 
# 
sub match_to_index {
    my ($pattern, $id, $match_index) =  @_;
    
    if (defined($match_index->{$id})) {
	my $found = 0;
	my $match_pos = -1;
	my $match_seq_obj;
	for my $seq_obj (@{ $match_index->{$id}}) {
	    my $str = $seq_obj->seq;
	    my $pos = index($str, $pattern);
	    if ($pos != -1) {
	        $found++;
		$match_pos = $pos;
		$match_seq_obj = $seq_obj;
	    } 
	}
	if ($found == 1) {
	    return ($match_pos, $match_seq_obj);
	} else {
	    return (-1);
	}   
	
    } else {
	##print STDERR "Fail 2\n";
	return (-1);
    }


}


#
# match pattern exactly to all in a list 
# using the basic index function as could not figure out aindex
# return pos and index if unique match to one only
#
sub match_to_all {
    my ($pattern, $str_list) = @_;

    my $found = 0;
    my $match_pos = -1;
    my $match_i = -1;
    for my $i (0..scalar(@$str_list)-1) {
	my $str = $str_list->[$i];
	my $pos = index($str, $pattern);
	if ($pos != -1) {
	    $found++;
	    $match_pos = $pos;
	    $match_i = $i;
	}
    }
    if ($found == 1) {
	return ($match_pos, $match_i);
    } else {
	return (-1, -1);
    }
}


#
# Estimate TSS probabilitities from the cases without G:s in the input
# - as a side effect add the input sequence information to the RNAs so that we do not need find them again
# - modified so that rnas come as parameter
sub estimate_tss_probs {
    my ($rnas, $rna_input_matches_file, $input_file, $trim_min_len) = @_;

    my $tssPos = TemplateSwitchPosition->new();
    
    # prefixes of the transcripts starting from the TSS, together with their UMI labels
    # modified: get as 
    #my $rnas = parse_template_switched_from_trim_info_file($trim_info_file, $trim_min_len);

    my $rna2inputseq = generate_r1_match_index($rna_input_matches_file, $input_file);

    my $r1_in = Bio::SeqIO->new(-file => $input_file, 				 
				-format => 'fastq');
    
    my @input_seq_list = ();
    my @input_id_list = ();
    my @input_seq_obj_list = ();
    # store to list, can take a lot of memory
    while (my $seq_obj = $r1_in->next_seq) {
	push @input_seq_list, $seq_obj->seq;
	push @input_id_list, $seq_obj->display_id;
	push @input_seq_obj_list, $seq_obj;
    }
    
    my $num_matching = 0;
    my ($num_not_matching_G, $num_not_matching) = (0, 0);
    
    for my $rna_id (keys %{ $rnas }) {
	my $rna_seq = $rnas->{$rna_id}->{seq};

	# TODO
	# Can modify this part to use the last 20 bases?
	
	# match few bases downstream in case some extra G:s
	my $match_seq = substr($rna_seq, 2, $min_len);
	
	# match to index
	my ($match_pos, $input_seq_obj) = match_to_index($match_seq, $rna_id, $rna2inputseq);
     
		
	if ($match_pos >= 0) {
	    $match_pos -= 2;
	    $num_matching++;
	    #print STDERR "Found unique match for $rna_id, rna length ", length($rna_seq), "\n";

	    my $input_seq = $input_seq_obj->seq;
	    my $input_id = $input_seq_obj->display_id;

	    $rnas->{$rna_id}->{input_seq} = $input_seq;
	    $rnas->{$rna_id}->{input_id} = $input_id;
	    $rnas->{$rna_id}->{match_pos} = $match_pos;

	    my $input_match_seq = substr($input_seq, $match_pos, $trim_min_len);
	    
	    my $pair = $tssPos->add_pair($rna_seq, $input_match_seq);

	    if ($DEBUG) {
		print STDERR "$pair\n";
		print STDERR "rna $rna_seq\n";
		print STDERR "in  $input_match_seq\n";
	    }
	} else {
	    my $rna_seq = $rnas->{$rna_id}->{seq};
	    my $rna_b1 = substr($rna_seq, 0, 1);
	    $num_not_matching++;
	    $num_not_matching_G++ if ($rna_b1 eq "G")
	    #print STDERR "Did not find unique match for $rna_id rna length ", length($rna_seq), "\n";
	}	
    }
    my $num_tot = scalar(keys %{ $rnas });
    print STDERR "Number of matching $num_matching, total $num_tot\n";
    my $pG_utr = $num_not_matching_G / $num_not_matching; 
    print STDERR "Number of non_matching $num_not_matching, of which starts with G $num_not_matching_G, ", sprintf("%.3f", $pG_utr), "\n";
    $tssPos->print_counts();
    $tssPos->estimate_1base_probabilities();

    
    # go through the RNAs again and assign probabilities according to the estimated probabilities
    my ($at, $ae, $a0, $a1, $a2) = (0, 0, 0, 0, 0);
    my ($at_utr,  $a0_utr, $a1_utr, $g_utr) = (0, 0, 0, 0);
    for my $rna_id (keys %{ $rnas }) {
	my $rna_seq = $rnas->{$rna_id}->{seq};

	# continue only with those for which input match was found
	if (defined($rnas->{$rna_id}->{input_seq})){
	    my $input_seq = $rnas->{$rna_id}->{input_seq};
	    my $match_pos = $rnas->{$rna_id}->{match_pos};
	    my $input_match_seq = substr($input_seq, $match_pos, $trim_min_len);
	    
	    my $tss_shift = $tssPos->assign_tss_shift_max_one($rna_seq, $input_match_seq);
	    #print STDERR "shift: $tss_shift\n";
	    if ($tss_shift == 0) {
		$rnas->{$rna_id}->{tss_shift} = 0;
		$a0++;
		$at++;
	    } elsif ($tss_shift == 1) {
		$rnas->{$rna_id}->{tss_shift} = 1;
		$a1++;
		$at++;
	    } elsif ($tss_shift == -1) {
		$rnas->{$rna_id}->{tss_shift} = -1;
		$ae++;
	    } else {
		die "Error: could not handle TSS shift $tss_shift\n";
	    }
	} else {
	    # mark as not matched, note that both inconsistency with input and no input are marked with negative value
	    $rnas->{$rna_id}->{tss_shift} = -2;

	    # as we want to print the UTR anyway, make a separate call for these without input
	    # - still need to decide the position
	    my $tss_shift = $tssPos->assign_tss_without_input_shift_max_one($rna_seq, $pG_utr);
	    my $rna_b1 = substr($rna_seq, 0, 1);
	    $g_utr++ if ($rna_b1 eq "G");
	    if ($tss_shift == 0) {
		$a0_utr++;
		$at_utr++;
		$rnas->{$rna_id}->{no_input_tss_shift} = 0;
	    } elsif ($tss_shift == 1) {
		$a1_utr++;
		$at_utr++;
		$rnas->{$rna_id}->{no_input_tss_shift} = 1;
	    } else {
		die "Error: could not handle no input TSS shift $tss_shift\n"
	    }
	    
	}
    }
    printf(STDERR "Assigned probabilities:  p(0) = %.3f, p(1) = %.3f, assigned %d, rejected %d\n", 
	   $a0 / $at, $a1 / $at, $at, $ae);
     printf(STDERR "Assigned probabilities for those lacking input:  p(0) = %.3f, p(1) = %.3f, assigned %d, p(G) = %.3f\n", 
	   $a0_utr / $at_utr, $a1_utr / $at_utr, $at_utr, $g_utr / $at_utr);
    return $rnas;
}



#
# match each RNA sequence to the INPUT sequences, if match found return the sequences around
# i.e. the promoter, UTR and both
sub get_tss_seqs {
    my ($rnas, $prom_out_file, $tss_out_file, $prom_tss_out_file, $umi_counts_out_file, $anno_file,
	$trim_min_len, $prom_len, $tss_len) = @_;
    
    my $prom_out = Bio::SeqIO->new(-file   => ">$prom_out_file",
				   -format => 'Fasta');
    my $tss_out = Bio::SeqIO->new(-file   => ">$tss_out_file",
 				  -format => 'Fasta');
    my $both_out = Bio::SeqIO->new(-file   => ">$prom_tss_out_file",
 				   -format => 'Fasta');

    open(AFH, ">$anno_file") or die("Could not open TSS annotation file for writing");
    printf(AFH join("\t", "input.id", "rna.id", "tss.shift", "tss.pos", "numg", "tss.seq", "tss.flank.seq", "rna.seq", "inputstr", "umi", "\n"));

    # match each rna to all input, can be slow
    my ($num_matching, $num_not_matching) = (0, 0);
    my %umi_counts = ();
    my $num_dup = 0;
    my $num_written = 0;

    for my $rna_id (keys %{ $rnas }) {
	my $rna_seq = $rnas->{$rna_id}->{seq};
	my $umi_label = $rnas->{$rna_id}->{umilabel};
	my $tss_shift = $rnas->{$rna_id}->{tss_shift};

	if ($tss_shift >= 0) {
	    $num_matching++;
	    
	    # extract the sequence before TSS and output to a file
	    my $input_str = $rnas->{$rna_id}->{input_seq};
	    my $input_id =  $rnas->{$rna_id}->{input_id};
	    my $tss_pos = $rnas->{$rna_id}->{match_pos} + $tss_shift;
	    
	    my $start_pos = $tss_pos - $prom_len; # 0-based coordinates
	    if ($start_pos < 0) {
		#print STDERR "trying to take too long promoter, start pos $start_pos in $input_id\n";
	    } else {
		# only one sequence per input promoter written out, now the first encountered, 
                # TODO choose the most frequent pos
		if (!defined($umi_counts{$input_id})) {   
		    # write out the promoter sequence
		    my $prom_str = substr($input_str, $start_pos, $prom_len);
		    my $prom_id = $input_id.":pos_".$tss_pos;   
		    my $prom_seq_obj = Bio::Seq->new(-seq        => $prom_str,
						     -display_id => $prom_id);		
		    $prom_out->write_seq($prom_seq_obj);
		    $num_written++;
		    
		    # write out the sequence starting at TSS if long enough (and requested in the first place)
		    # and also the prom + utr
		    my $tss_str = substr($input_str, $tss_pos, $tss_len);
		    if (($tss_len > 0) && (length($tss_str) >= $tss_len)) {
			my $tss_id = $input_id.":pos_".$tss_pos;
			my $tss_seq_obj = Bio::Seq->new(-seq        => $tss_str,
							-display_id => $tss_id);		
			$tss_out->write_seq($tss_seq_obj);

			my $prom_tss_str = $prom_str.$tss_str;
			my $prom_tss_seq_obj = Bio::Seq->new(-seq        => $prom_tss_str,
					 		     -display_id => $tss_id);
			$both_out->write_seq($prom_tss_seq_obj);

			# annotate the choice TSS
			my $numg = $rnas->{$rna_id}->{numg};
			printf(AFH join("\t", $input_id, $rna_id, $tss_shift, $tss_pos, $numg, 
					$tss_str, $prom_tss_str, $rna_seq, $input_str, $umi_label, "\n")); 
			
		    } else {
			# write the annotation anyway, put . to missing values
			my $numg = $rnas->{$rna_id}->{numg};
			printf(AFH join("\t", $input_id, $rna_id, $tss_shift, $tss_pos, $numg, 
					'.', '.', $rna_seq, $input_str, $umi_label, "\n")); 
		    }
		    
		} else {
		    $num_dup++;
		}
		
		# update the UMI counts
		#print STDERR "Observed TSS $input_id : $match_pos : $umi_label\n";
		if (!defined($umi_counts{$input_id}->{$tss_pos}->{$umi_label})) {
		    $umi_counts{$input_id}->{$tss_pos}->{$umi_label} = 1;
		} else {
		    $umi_counts{$input_id}->{$tss_pos}->{$umi_label}++;
		}		
	    } 
	    
	} else {
	    $num_not_matching++;	    
	}
	    
    }
    
    # write out the UMI counts
    open(UFH, ">$umi_counts_out_file") or die("Could not open output file for UMI counts: $!");
    for my $id (sort keys(%umi_counts)) {
	for my $pos (sort {$a <=> $b} keys %{ $umi_counts{$id} }) {
	    for my $umi_label (sort keys %{ $umi_counts{$id}->{$pos} }) {
		my $count =  $umi_counts{$id}->{$pos}->{$umi_label};
		print UFH "$id\t$pos\t$umi_label\t$count\n";
	    }
	}
    }
    close(UFH);
    
    print STDERR "\tnumber of matched sequences $num_matching, number of sequences not matched (or inconsistent) $num_not_matching\n";
    print STDERR "\twrote in total ", $num_written, " sequences\n";
    print STDERR "\tNumber of RNA:s matching to the same input: $num_dup\n";
}







#
# parse the cutadapt file that contains the spliced intron matches and importantly the preceding sequence
# - also check and remove the UMI and template switch Gs
#    Read name
#    Number of errors, -1 means that no match was found
#    0-based start coordinate of the adapter match
#    0-based end coordinate of the adapter match
#    Sequence of the read to the left of the adapter match (can be empty)
#    Sequence of the read that was matched to the adapter
#    Sequence of the read to the right of the adapter match (can be empty)
#    Name of the found adapter.
#    Quality values corresponding to sequence left of the adapter match (can be empty)
#    Quality values corresponding to sequence matched to the adapter (can be empty)
#    Quality values corresponding to sequence to the right of the adapter match (can be empty)
#
# addition November 21, 2019: writ also the suffix of the random region for input matching
sub parse_template_switched_from_trim_info_file {
    my ($filename, $min_rna_len, $suffix_file) = @_;

    open(IFH, "<$filename") or die("Could not open trim info file $filename: $!");
    
    my $num_match = 0;
    my $num_nogs = 0;
    my $num3g = 0;
    my $num4g = 0;
    my $num5g = 0;

    open(OF, ">$suffix_file") or die("Could not open suffix file $suffix_file for writing: $!");
    
    my %rnas = ();
    while (my $line = <IFH>) {
	chomp $line;

	my ($id, $ne, $pos, $end_pos, $left_flank_seq, $match_seq, $right_flank_seq, $adapter_name, $qleft, $qmatch, $qright) =
	    split /\t/, $line;
	
	if ($ne >= 0 && length($left_flank_seq) > $min_rna_len) {
	    #print STDERR "\tflank: $left_flank_seq\n";
	    
	    # remove UMI
	    my $umi_label = substr($left_flank_seq, 0, $UMI_LEN);
	    my $no_umi_seq = substr($left_flank_seq, $UMI_LEN);
	    #print STDERR "\tno UMI: $no_umi_seq\n";

	    # remove the part after space from id
	    $id =~ s/\s[\w:]+//;
	    
	    # count template switch Gs that should now be in the beginning
	    # require at least 3, can be also 4 or 5
	    my $rna_seq = "";
	    my $numg = 0;
	    if (substr($no_umi_seq, 0, 5) eq "GGGGG") {
		$num5g++;
		$numg = 5;
	    } elsif (substr($no_umi_seq, 0, 4) eq "GGGG") {
		$num4g++;
		$numg = 4;
	    } elsif (substr($no_umi_seq, 0, 3) eq "GGG") {
		$num3g++;
		$numg = 3;
	    } else {
		$num_nogs++;
	    }
	    if ($numg > 0) {
		$rna_seq = substr($no_umi_seq, 3);
	    }
	    #print STDERR substr($no_umi_seq, 0, 10), ", ", substr($rna_seq, 0, 10), "\n";
	    
	    if (length($rna_seq) > $min_rna_len) {
				
		$rnas{$id}->{seq} = $rna_seq;
		$rnas{$id}->{umilabel} = $umi_label;
		$rnas{$id}->{numg} = $numg;
  
		my $rand_region_suffix = substr($left_flank_seq, -1 * $min_rna_len);
		$rnas{$id}->{rand_suffix} = $rand_region_suffix;
		print(OF "$id\t$rand_region_suffix\n"); 
		$num_match++;

		#print STDERR "\tRNA:$rna_seq\n";
	    }
	}
    }
    print STDERR "\tparsing trimmed: starting with 5 Gs $num5g, 4 Gs $num4g; 3 Gs $num3g; no 3 G $num_nogs\n";
    print STDERR "\tparsing trimmed: parsed $num_match sequences flanking the trimmed sequence and at least $min_rna_len long\n";
    return \%rnas;
}



#-----------------------------------------------------
# Read option values
# defaults here are the NextSeq run with new primers
#-----------------------------------------------------
sub configure {

    $conf  = AppConfig::State->new();
    my $getopt = AppConfig::Getopt->new($conf);

    # new longer one
    $conf->define("UMI_length", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => 10
		  });

    
    # RNA sequence files
    $conf->define("R1_rna_fastq_file", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => ""
		  });


    # fastq files
    $conf->define("R1_input_fastq_file", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => ""
		  });



    # input R1 first 20 bases sorted text file, note that has only two columns, qualities omitted
    $conf->define("R1_input_tab_sorted_file", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => ""

		  });
    
    # output 
    $conf->define("output_dir", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => ""
		  });


    # how many bases to output before TSS
    $conf->define("promoter_length", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => 40
		  });
    # how many bases to output after TSS
    $conf->define("utr_length", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => 40
		  });   



   
    # the sequence that should be in the 3' end RNA derived sequence 
    $conf->define("spliced_intron_seq", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => "AAGCTTCTGCCTTCTCCCTCCTGTGAGTTTGGTTGGTGTACAGTAGCTTCCACC"
		  });
   
    # the minimum length we should see to be be sure that we have the splicing
    # the splice site is after 31 bases
    # AAGCTTCTGCCTTCTCCCTCCTGTGAGTTTG
    $conf->define("min_intron_overlap", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => "40"
		  });

    # the minimum length we should see to be be sure that we have the splicing
    # the splice site is after 31 bases
    # AAGCTTCTGCCTTCTCCCTCCTGTGAGTTTG
    $conf->define("run_revcoms_and_stats", 
		  {
		      ARGCOUNT => ARGCOUNT_ONE,
		      DEFAULT => "0"
		  });

     
    $conf->define("config_file", {
	              ARGCOUNT => ARGCOUNT_ONE
                  });

    
    # Read command line arguments
    $getopt->parse();
    
    if (defined($conf->config_file)) {
	printf(STDERR "Random Starr-seq: reading configuration file %s\n", 
	       $conf->config_file);
  
       # Read configuration file
       my $cfgfile = AppConfig::File->new($conf, $conf->config_file);
    }
}


