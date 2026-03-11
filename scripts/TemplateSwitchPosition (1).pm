#
# package to estimate and store template switch position information 
# @author Teemu Kivioja
#

package TemplateSwitchPosition;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);


# the classification of the template switch patterns
# C: certain TSS position 
# U: uncertain
# E: error, observed RNA and input DNA are conflicting under our assumptions
our %ts_class = (
    "HH-HH" => "C",
    "HH-HG" => "E",
    "HH-GH" => "E",
    "HH-GG" => "E",  
    "HG-HH" => "E",
    "HG-HG" => "U",
    "HG-GH" => "E",
    "HG-GG" => "E",    
    "GH-HH" => "C",
    "GH-HG" => "E",
    "GH-GH" => "U",
    "GH-GG" => "E",    
    "GG-HH" => "C",
    "GG-HG" => "U",
    "GG-GH" => "U",
    "GG-GG" => "U"    
);


# collapse to four classes according to the first base
our %ts_first_class = (
    "HH-HH" => "H-H",
    "HH-HG" => "H-H",
    "HH-GH" => "H-G",
    "HH-GG" => "H-G",  
    "HG-HH" => "H-H",
    "HG-HG" => "H-H",
    "HG-GH" => "H-G",
    "HG-GG" => "H-G",    
    "GH-HH" => "G-H",
    "GH-HG" => "G-H",
    "GH-GH" => "G-G",
    "GH-GG" => "G-G",    
    "GG-HH" => "G-H",
    "GG-HG" => "G-H",
    "GG-GH" => "G-G",
    "GG-GG" => "G-G"       
    );

our %ts_first_class_type = (
    "H-H" => "C",
    "H-G" => "E",
    "G-H" => "C",
    "G-G" => "U"
    );



#--------------
# STATIC
#--------------


#
# Get the class of the first base comparison from the two we have stored
#
sub get_ts_first_pattern {
    my $pattern_pair = shift;

    return $ts_first_class{$pattern_pair};

}


# 
# dinucleotide pattern, is G or not G i.e H
#
sub ts_pattern {
    my ($b1, $b2) = @_;

    if ($b1 ne "G" && $b2 ne "G") {
	return "HH";
    } elsif ($b1 ne "G" && $b2 eq "G") {
	return "HG";
    } elsif ($b1 eq "G" && $b2 ne "G") {
	return "GH";	
    } elsif ($b1 eq "G" && $b2 eq "G") {
	return "GG";
    } else {
	die "Could not handle dinucleotide $b1 $b2";
    }
}


#
# nucleotide pattern in RNA and input at the template switch position  
#
sub ts_pattern_pair {
    my ($rna_seq, $input_seq) = @_;

    # the first two bases to be compared
    my $rna_b1 = substr($rna_seq, 0, 1);
    my $rna_b2 = substr($rna_seq, 1, 1);
    my $input_b1 = substr($input_seq, 0, 1);
    my $input_b2 = substr($input_seq, 1, 1);

    my $rp = ts_pattern($rna_b1, $rna_b2);
    my $ip = ts_pattern($input_b1, $input_b2);
   
    my $ri_pair = join("-", ($rp, $ip));
    #print STDERR "pair: $rp $ip $ri_pair\n";

    return $ri_pair;
}
    




#
# NEW object
#
sub new {
    my ($class) = @_;
    
    my $self = {};
    bless($self, $class);

    # initialize observed counts
    $self->{ts_counts} = {};

    for my $pat (keys %ts_class) {
	$self->{ts_counts}->{$pat} = 0;
    }

    $self->{ts_first_counts} = {};

    for my $pat (keys %ts_first_class) {
	$self->{ts_first_counts}->{$pat} = 0;
    }


    
    # keep the results the same between runs - moved up to the script level
    #srand(42);
    
    return $self;
}


sub print_counts {
    my $self = shift;

    for my $pat (sort keys %ts_class) {
	print $pat, "\t", $ts_class{$pat}, "\t", $self->{ts_counts}->{$pat}, "\n";
    }
    
}


#
# add TSS template switch pattern seen one RNA/INPUT sequence pair to the counts
#
sub add_pair {
    my ($self, $rna_seq, $input_seq) = @_;
    
    my $ts_pair = ts_pattern_pair($rna_seq, $input_seq);
    #print STDERR "ts pair $ts_pair\n";

    if (!defined($self->{ts_counts}->{$ts_pair})) {
	die "Error: unknown TS pair $ts_pair";
    } else {
	$self->{ts_counts}->{$ts_pair}++;
    }
    return $ts_pair;
}


#
# sum the 2 base counts to get the first base counts
#
sub sum_ts_pair_counts {
    my $self = shift;

    for my $pair_pat (keys %{ $self->{ts_counts} }) {
	my $pat = $ts_first_class{$pair_pat};
	$self->{ts_first_counts}->{$pat} += $self->{ts_counts}->{$pair_pat};
    }

}



#
# use the counts from the cases without uncertainty in the template switch position 
# to estimate the probabilities for the rest, uncertain ones  
# - this one does for two bases
sub estimate_2base_probabilities {
    my ($self) = @_;
 
    my $ts_counts = $self->{ts_counts};
      
    # class total counts
    my ($nc, $nu, $ne, $nt)  = 0;
    for my $pat (keys %ts_class) {
	if ($ts_class{$pat} eq "C") {
	    $nc += $ts_counts->{$pat};
	} elsif ($ts_class{$pat} eq "U") {
	    $nu += $ts_counts->{$pat};
	} else {
	    $ne += $ts_counts->{$pat};
	}
	$nt += $ts_counts->{$pat};
    }
    # check
    die "Error: class totals do not add up: $nc + $nu + $ne not equal to $nt" if (($nc + $nu + $ne) != $nt);
    
    # probabilities for the TSS positions 0, 1 and 2 estimated from those without any Gs in the input 
    my $p0 = $ts_counts->{"HH-HH"} / $nc;
    my $p1 = $ts_counts->{"GH-HH"} / $nc;
    my $p2 = $ts_counts->{"GG-HH"} / $nc;
    print STDERR "Observed probabilities: p(0) = $p0, p(1) = $p1, p(2) = $p2\n";
    
    # check
    die("Error: estimated TSS probabilities do not sum to one") if (abs(($p0 + $p1 + $p2) - 1) > 1e-3);
    
    # probabilities of the patterns inside the two main classes
    my %p_uts = ();
    my %p_cts = ();
    my $pt_uts = 0;
    my $pt_cts = 0;
    for my $pat (sort keys %ts_class) {
	if ($ts_class{$pat} eq "U") {
	    $p_uts{$pat} = $ts_counts->{$pat} / $nu;
	    $pt_uts += $p_uts{$pat};
	    print STDERR "U: $pat prob ", $p_uts{$pat}, "\n";
	} elsif ($ts_class{$pat} eq "C") {
	    $p_cts{$pat} = $ts_counts->{$pat} / $nc;
	    $pt_cts += $p_cts{$pat};
	    print STDERR "C: $pat prob ", $p_cts{$pat}, "\n";
	}
    }
    print STDERR "sum of uncertain classes: $pt_uts\n";
    print STDERR "sum of uncertain classes: $pt_cts\n";
    
    # estimate the unknowns needed
    my $a1 = $p1 / $p0;
    my $a2 = $p2 / $p1;
    
    my $x = 1 / ($a1 + 1);
    my $y = 1 / ($a2 + 1);
    
    my $z = ($p1 - $p_uts{"GH-GH"} * (1 - $x) - $p_uts{"GG-HG"} * $y) / $p_uts{"GG-GG"};
    

    print STDERR "x = $x, y = $y, z = $z\n";

    # the conditional probabilities of interest
    my %p_cond = ();
 
    $p_cond{"GH-GH"}->{0} = $x;
    $p_cond{"GG-GG"}->{0} = ($p0 - $p_uts{"HG-HG"} - $p_uts{"GH-GH"} * $x) / $p_uts{"GG-GG"};

    $p_cond{"GH-GH"}->{1} = 1 - $x;
    $p_cond{"GG-HG"}->{1} = $y;
    $p_cond{"GG-GG"}->{1} = $z;

    $p_cond{"GG-HG"}->{2} = 1 - $y;
    $p_cond{"GG-GG"}->{2} = ($p2 - $p_uts{"GG-GH"} - $p_uts{"GG-HG"} * (1 - $y)) / $p_uts{"GG-GG"};
    
    my $z2 = 1 - $p_cond{"GG-GG"}->{0} - $p_cond{"GG-GG"}->{1};
    
    print STDERR "z2: $z2\n";
    
    
    # check that all probabilities are between 0 and 1 and sum to one
    my $pct = 0;
    for my $pat (sort keys %p_cond) {
	for my $pos (sort keys %{ $p_cond{$pat} }) {
	    my $pc = $p_cond{$pat}->{$pos};
	    print STDERR "cond prob for $pat at $pos is $pc\n";
	    #die "Error: conditional probability for $pat at $i not between 0 and 1: $pc"
		#if ($pc > 1 || $pc < 0);
	    $pct += $pc;
	}
    }
    print STDERR "Conditional probabilities add to: $pct\n"; 
    
    $self->{pcond} = \%p_cond;

    # the cumulative conditional probabilities i.e. the TSS is here or before it
    # - easier to use for assignment
    my %p_cum = ();
    
    # the certain ones
    # 0 
    $p_cum{"HH-HH"}->{0} = 1;
    $p_cum{"HG-HG"}->{0} = 1;

    # 1
    $p_cum{"GH-HH"}->{0} = 0;
    $p_cum{"GH-HH"}->{1} = 1;

    # 2
    $p_cum{"GG-HH"}->{0} = 0;
    $p_cum{"GG-HH"}->{1} = 0;
    $p_cum{"GG-HH"}->{2} = 1;

    $p_cum{"GG-GH"}->{0} = 0;
    $p_cum{"GG-GH"}->{1} = 0;
    $p_cum{"GG-GH"}->{2} = 1;
    
    # the uncertain ones
    $p_cum{"GH-GH"}->{0} = $p_cond{"GH-GH"}->{0};
    $p_cum{"GH-GH"}->{1} = 1;

    $p_cum{"GG-HG"}->{0} = 0;
    $p_cum{"GG-HG"}->{1} = $p_cond{"GG-HG"}->{1};
    $p_cum{"GG-HG"}->{2} = 1;

    $p_cum{"GG-GG"}->{0} = $p_cond{"GG-GG"}->{0};
    $p_cum{"GG-GG"}->{1} = $p_cond{"GG-GG"}->{0} + $p_cond{"GG-GG"}->{1};
    $p_cum{"GG-GG"}->{2} = 1;

    $self->{pcum} = \%p_cum;
    
    return \%p_cum;
}


#
# use the counts from the cases without uncertainty in the template switch position 
# to estimate the probabilities for the rest, uncertain ones  
# - the first position after GGG only
sub estimate_1base_probabilities {
    my ($self) = @_;

    # move from 2 base counts to only the first one
    $self->sum_ts_pair_counts();  
	
    my $ts_counts = $self->{ts_first_counts};
      
    # class total counts
    my ($nc, $nu, $ne, $nt)  = 0;
    for my $pat (keys %ts_first_class_type) {
	if ($ts_first_class_type{$pat} eq "C") {
	    $nc += $ts_counts->{$pat};
	} elsif ($ts_first_class_type{$pat} eq "U") {
	    $nu += $ts_counts->{$pat};
	} else {
	    $ne += $ts_counts->{$pat};
	}
	$nt += $ts_counts->{$pat};
    }
    # check
    die "Error: class totals do not add up: $nc + $nu + $ne not equal to $nt" if (($nc + $nu + $ne) != $nt);
    
    # probabilities for the TSS positions 0 and 1 estimated from those without any Gs in the input 
    my $p0 = $ts_counts->{"H-H"} / $nc;
    my $p1 = $ts_counts->{"G-H"} / $nc;
    print STDERR "Observed probabilities: p(0) = $p0, p(1) = $p1\n";

    $self->{est_p0} = $p0;
    $self->{est_p1} = $p1;
    
    # check
    die("Error: estimated TSS probabilities do not sum to one") if (abs(($p0 + $p1) - 1) > 1e-3);
    
    # probabilities of the patterns inside the two main classes
    my %p_uts = ();
    my %p_cts = ();
    my $pt_uts = 0;
    my $pt_cts = 0;
    for my $pat (sort keys %ts_first_class_type) {
	
	if ($ts_first_class_type{$pat} eq "U") {
	    $p_uts{$pat} = $ts_counts->{$pat} / $nu;
	    $pt_uts += $p_uts{$pat};
	    print STDERR "U: $pat prob ", $p_uts{$pat}, "\n";
	    
	} elsif ($ts_first_class_type{$pat} eq "C") {
	    $p_cts{$pat} = $ts_counts->{$pat} / $nc;
	    $pt_cts += $p_cts{$pat};
	    print STDERR "C: $pat prob ", $p_cts{$pat}, "\n";
	}
    }
    print STDERR "sum of uncertain classes: $pt_uts\n";
    print STDERR "sum of certain classes: $pt_cts\n";
    
    # estimate the unknowns needed
    my $a = $p0 / $p1;
    my $x = 1 / ($a + 1);
        
    #print STDERR "x = $x, y = $y, z = $z\n";

    # the conditional probabilities of interest
    my %p_cond = ();
 
    $p_cond{"G-G"}->{0} = 1 - $x;
    $p_cond{"G-G"}->{1} = $x;
    
    # check that all probabilities are between 0 and 1 and sum to one
    my $pct = 0;
    for my $pat (sort keys %p_cond) {
	for my $pos (sort keys %{ $p_cond{$pat} }) {
	    my $pc = $p_cond{$pat}->{$pos};
	    print STDERR "cond prob for $pat at $pos is $pc\n";
	    #die "Error: conditional probability for $pat at $i not between 0 and 1: $pc"
		#if ($pc > 1 || $pc < 0);
	    $pct += $pc;
	}
    }
    print STDERR "Conditional probabilities add to: $pct\n"; 
    
    $self->{pcond_first} = \%p_cond;

    # the cumulative conditional probabilities i.e. the TSS is here or before it
    # - easier to use for assignment
    my %p_cum = ();

    # 0 
    $p_cum{"H-H"}->{0} = 1;

    # 1
    $p_cum{"G-H"}->{0} = 0;
    $p_cum{"G-H"}->{1} = 1;

    # uncertain
    $p_cum{"G-G"}->{0} = $p_cond{"G-G"}->{0};
    $p_cum{"G-G"}->{1} = 1;

    $self->{pcum_first} = \%p_cum;

    
    return \%p_cum;
}
 



#
# assign TSS position according to the rules and probabilities
#
sub assign_tss_shift_max_two {
    my ($self, $rna_seq, $input_seq) = @_;

    my $pat = ts_pattern_pair($rna_seq, $input_seq);
    #print STDERR "pat $pat\n";

    if (($ts_class{$pat} eq "C") || ($ts_class{$pat} eq "U")) {
	my $r = rand();
	if ($r < $self->{pcum}->{$pat}->{0}) {
	    return 0;
	} elsif ($r < $self->{pcum}->{$pat}->{1}) {
	    return 1;
	} elsif ($self->{pcum}->{$pat}->{2} = 1) {
	    return 2;
	} else {
	    die "Error: could not assign TSS to pattern $pat";
	}
    } elsif ($ts_class{$pat} eq "E") {
	return -1; # one of the cases where input and rna are not consistent
    } else {
	die "Error: Could not assign TSS to pattern $pat";
    }
}


#
# assign TSS position according to the rules and probabilities
# - shift by max one
sub assign_tss_shift_max_one {
    my ($self, $rna_seq, $input_seq) = @_;

    my $pair_pat = ts_pattern_pair($rna_seq, $input_seq);
    #print STDERR "pat $pat\n";

    if (($ts_class{$pair_pat} eq "C") || ($ts_class{$pair_pat} eq "U")) {

	my $pat = get_ts_first_pattern($pair_pat);
	
	my $r = rand();
	if ($r < $self->{pcum_first}->{$pat}->{0}) {
	    return 0;
	} elsif ($r < $self->{pcum_first}->{$pat}->{1}) {
	    return 1;
	} else {
	    die "Error: could not assign TSS to pattern $pat";
	}
    } elsif ($ts_class{$pair_pat} eq "E") {
	return -1; # one of the cases where input and rna are not consistent
    } else {
	die "Error: Could not assign TSS to pattern $pair_pat";
    }
}


#
# assign TSS position according to the rules and probabilities
# - shift by max one
# case without input sequence
sub assign_tss_without_input_shift_max_one {
    my ($self, $rna_seq, $pG) = @_;

    # first check if there is G in RNA or not
    my $rna_b1 = substr($rna_seq, 0, 1);
    if ($rna_b1 eq "G") {

	# G, use observed probability of G and the zero prob estimated
	# p(0|G) = (P(0) - P(H)) / P(G)
	my $p0 = $self->{est_p0};
	my $pH = 1 - $pG;
	my $p0G = ($p0 - $pH) / $pG; 
	
	my $r = rand();	
	if ($r < $p0G) {
	    return 0;
	} else {
	    return 1;
	}
    } else {
	# no G, no shift
	return 0;
    }
}



#
# MAIN - for testing only
# - This code tests the estimation of 2 bases after GGG
if (0) {
    # simulated probabilities
    my ($p0, $p1, $p2) = (0.5, 0.4, 0.1);
   

    print STDERR "testing ...\n";
    my $testSwitch = TemplateSwitchPosition->new();
    
    my $test_file = "../tss_seqs/TS-rGPE_combined_1906/TS-rGPE_1906_R1.promoter_controls_noN.txt";
  
    # read input sequences and add some G:s to simulate template switch
    open(FH, "<$test_file") or die "Could not open test file $test_file: $!";
    my @input_seqs = ();
    my @rna_seqs = ();
    my ($g0, $g1, $g2) = (0, 0, 0);
    while (my $input = <FH>) {
	chomp($input);

	# simulate some bias
	if (rand() < 0.5 && (substr($input, 0, 1) ne "G")) {
	    #print STDERR "Skipping $input\n";
	    next;
        }
	    
	push @input_seqs, $input;
	my $rna = $input;
	my $r = rand();
	if ($r < $p2) {
	    $g2++;
	    substr($rna, 0, 2, "GG"); 
	} elsif ($r < $p2 + $p1) {
	    $g1++;
	    substr($rna, 0, 1, "G");
	} else {
	    $g0++;
	}
	push @rna_seqs, $rna;
    }
    my $n = $g0 + $g1 + $g2;
    printf(STDERR "Added Gs: 0 = %.2f, 1 = %.2f, 2 = %.2f\n", $g0 / $n, $g1 / $n, $g2 / $n);
	    
    # process and see if we get the same probabilities
    
    for my $i (0..scalar(@rna_seqs)-1) {
	#print STDERR "adding ", $rna_seqs[$i], ", ", $input_seqs[$i], " ..\n";
	$testSwitch->add_pair($rna_seqs[$i], $input_seqs[$i]);
    }
    $testSwitch->print_counts();
    $testSwitch->estimate_2base_probabilities();

    # go through the sequences again and assign each a TSS position
    my ($at, $ae, $a0, $a1, $a2) = (0, 0, 0, 0, 0);
    for my $i (0..scalar(@input_seqs)-1) {
	my $tss = $testSwitch->assign_tss_shift_max_two($rna_seqs[$i], $input_seqs[$i]);
	#print STDERR "assigned TSS pos: $tss\n";
	$at++;
	if ($tss == 0) {
	    $a0++;
	} elsif ($tss == 1) {
	    $a1++;	    
	} elsif ($tss == 2) {
	    $a2++;
	} else {
	    $ae++;
	}
    }
    printf(STDERR "Simulated probabilities: p(0) = %.3f, p(1) = %.3f, p(2) = %.3f\n", $p0, $p1, $p2);
    printf(STDERR "Assigned probabilities:  p(0) = %.3f, p(1) = %.3f, p(2) = %.3f\n", $a0 / $at, $a1 / $at, $a2 / $at);


    
}



#
# MAIN - for testing only
# - This code tests the estimation of only 1 after GGG
if (0) {
    # simulated probabilities
    my ($p0, $p1, $p2) = (0.6, 0.4, 0.0);   

    print STDERR "testing ...\n";
    my $testSwitch = TemplateSwitchPosition->new();
    
    my $test_file = "../tss_seqs/TS-rGPE_combined_1906/TS-rGPE_1906_R1.promoter_controls_noN.txt";
  
    # read input sequences and add some G:s to simulate template switch
    open(FH, "<$test_file") or die "Could not open test file $test_file: $!";
    my @input_seqs = ();
    my @rna_seqs = ();
    my ($g0, $g1, $g2) = (0, 0, 0);
    while (my $input = <FH>) {
	chomp($input);

	# simulate some bias
	if (rand() < 0.5 && (substr($input, 0, 1) ne "G")) {
	    #print STDERR "Skipping $input\n";
	    next;
        }
	    
	push @input_seqs, $input;
	my $rna = $input;
	my $r = rand();
	if ($r < $p2) {
	    $g2++;
	    substr($rna, 0, 2, "GG"); 
	} elsif ($r < $p2 + $p1) {
	    $g1++;
	    substr($rna, 0, 1, "G");
	} else {
	    $g0++;
	}
	push @rna_seqs, $rna;
    }
    my $n = $g0 + $g1 + $g2;
    printf(STDERR "Added Gs: 0 = %.2f, 1 = %.2f, 2 = %.2f\n", $g0 / $n, $g1 / $n, $g2 / $n);
	    
    # process and see if we get the same probabilities
    
    for my $i (0..scalar(@rna_seqs)-1) {
	#print STDERR "adding ", $rna_seqs[$i], ", ", $input_seqs[$i], " ..\n";
	$testSwitch->add_pair($rna_seqs[$i], $input_seqs[$i]);
    }
    $testSwitch->print_counts();
    $testSwitch->estimate_1base_probabilities();

    # go through the sequences again and assign each a TSS position
    my ($at, $ae, $a0, $a1) = (0, 0, 0, 0);
    for my $i (0..scalar(@input_seqs)-1) {
	my $tss = $testSwitch->assign_tss_shift_max_one($rna_seqs[$i], $input_seqs[$i]);
	#print STDERR "assigned TSS pos: $tss\n";
	$at++;
	if ($tss == 0) {
	    $a0++;
	} elsif ($tss == 1) {
	    $a1++;	    
	} else {
	    $ae++;
	}
    }
    printf(STDERR "Simulated probabilities: p(0) = %.3f, p(1) = %.3f, p(2) = %.3f\n", $p0, $p1, $p2);
    printf(STDERR "Assigned probabilities:  p(0) = %.3f, p(1) = %.3f\n", $a0 / $at, $a1 / $at);

    
}


1;
