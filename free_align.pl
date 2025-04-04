#!/usr/bin/perl

##========================================================================
##
## Usage: free_align.pl [options] seq1 seq2
##    
##    free_align takes two RNA sequences and attempts to find the
##    minimum free energy binding between the two (the most stable binding)
##
##    seq1 can be a string or a file, in FASTA format, and is assumed to be 
##      written in 5'->3' orientation.
##    seq2 can be a string or a file, in FASTA format, and is assumed to be 
##      written in 3'->5' orientation.
##
##    NOTES ON OUTPUT:
##    
##    | = pair is a member of a helix    
##    0 = pair is a member of a loop
##    P = bulge in top strand
##    b = bulge in bottom strand
##
##    Any and all lowercase bases in the output were not used to calculate
##    the best alignment, but are there to show that the end of the helix is
##    not terminal.
##
## Author: Joshua Starmer <jdstarme@unity.ncsu.edu>, 2004
## 
## Copyright (C) 2004, Joshua Starmer
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
## 
##========================================================================

use warnings;
use strict;

use Getopt::Std; 
use FreeAlign;


##========================================================================
##
## DEFINE GLOBAL CONSTANTS
##
##========================================================================

use constant TRUE => 1;
use constant FALSE => 0;



##========================================================================
##
## DEFINE GLOBAL VARIABLES
##
##========================================================================

my $gEnergyOnly = FALSE;
my $gDebugMode = FALSE;
my $gForceBind = FALSE;
my $gHelixOnly = FALSE;
my $gDoubletID = "XIA_MATHEWS";
my $gLoopID = "JAEGER";
my $gTerminalBulgeRule = "ALLOWED";
my $gInternalBulgeRule = "ALLOWED";
my $gLoopRule = "ALLOWED";
my $gTemperature = 37;
my $gTest = FALSE;

##========================================================================
##
## DEFAULT HELP MESSAGE
##
##========================================================================
my $sHelp = <<END_OF_HELP;

Usage: free_align.pl [options] seq1_5'->3' seq2_3'->5'
    
    free_align takes two RNA sequences and attempts to find the
    minimum free energy binding between the two (the most stable binding)

    seq1 can be a string or a file, in FASTA format, and is assumed to be 
      written in 5'->3' orientation.
    seq2 can be a string or a file, in FASTA format, and is assumed to be 
      written in 3'->5' orientation.

    example:
    ./free_align.pl aug uac
    Doublet Binding Parameters: XIA_MATHEWS
    Internal Loop Parameters: JAEGER
    Terminal Bulges: ALLOWED
    Internal Bulges: ALLOWED
    Internal Loops: ALLOWED
    Binding Temperature: 37 C

    Delta-G for best pairing: -3.19596 + 4.075225 = 0.879264999999997

    Length of bound region between sequences: 3

    seq1 binding begins on base index: 0
    seq2 binding begins on base index: 0

    | = pair is a member of a helix
    0 = pair is a member of a loop
    P = bulge in top strand
    b = bulge in bottom strand

    seq1 0: AUG
          : |||
    seq2 0: UAC


    NOTES ON OUTPUT:

    Any and all lowercase bases in the output were not used to calculate
    the best alignment, but are there to show that the end of the helix is
    not terminal.

Options

-h
   print this help message.

-f
   Force an alignment between the two strands of equal length, assuming that 
   they both pair together from their first bases on.  This method ignores 
   both loops and bulges as possibilities.  Only bases in this forced 
   alignment that can form helices are scored.  
   "-f" takes precedence over "-o".
   This method scores doublets as both internal and terminal, depending on
   where they are in the helix.

-o
   Helix only.  This is similar to the '-f' option, execpt that the helix
   does not have to begin formation from the first bases in both sequence.
   In other words, this is the same as '-f' except there is no penalty for
   terminal bulges (gaps).  "-f" takes precedence over "-o".

-B [0|1]  -- NOT FULLY TESTED
   Choose between allowing and disallowing both terminal and interanl bulges 
   in the pairing between strands.  1, the default, allows terminal bulges.  
   For the time being, this option forces all doublets to be scored as
   "internal" doublets.

-b [0|1]
   Choose between allowing and disallowing internal bulges in the pairing 
   between strands.  1, the default, allows internal bulges.  0 does not allow 
   them.

-L [0|1]  -- NOT FULLY IMPLEMENTED!!!!
   Choose between allowing and disallowing terminal loops in the pairing 
   between strands.  1, the default, allows terminal loops.  0 sets the penalty
   for terminal loops to the value of BIG_NUM (see FreeAlign.pm)

-l [0|1]
   Choose between allowing and disallowing internal loops in the pairing 
   between strands.  1, the default, allows internal loops.  0 does not allow 
   them.

-e
   Only print out the energy given off (the Delta-G value ) by the most
   stable pairing of the two sequences (and do not print out anything else).

-p FREIER|SANTALUCIA|XIA_MATHEWS
   Allows you to determine which set of parameters are used to simulate binding
   between the to strands of nucleotides.  The default value is XIA_MATHEWS.
   FREIER - RNA binding parameters from 1986
   SANTALUCIA - DNA binding parameters from 1998
   XIA_MATHEWS - RNA binding parameters from 1998

-t temperature
   The temperature, in Celsius, at which the binding is supposed to take place.

-d
   Debug Mode.  Prints out all sorts of debugging information...

END_OF_HELP


##========================================================================
##
## here is the "main"
##
##========================================================================
{
    
    # check to see if anything is coming in to the program (via a pipe or
    # on the command line).  If there is no input at all, print out the help
    if (-t STDIN && !@ARGV) {
        print STDERR $sHelp;
        exit 1;
    }
    
    # process the command line arguments
    my %opts; # a hash table to store file names passed in as agruments
    getopts('hdeB:b:L:l:fop:t:T', \%opts);
     
    if ($opts{'h'}) { # print help and exit
        print STDERR $sHelp;
        exit 1;
    }

    if ($opts{'e'}) { # restrict output to delta-g only
	$gEnergyOnly = TRUE;
    }

    if ($opts{'d'}) { # turn on debugging mode!
	$gDebugMode = TRUE;
    }

    # now that we know we are not just going to exit the program, create
    # an Align object...
    my $align_obj = new FreeAlign();

    if (defined($opts{'B'})) {
	if ($opts{'B'} == 0) {
	    $align_obj->no_bulge_possible(TRUE);
	    $gTerminalBulgeRule = "NOT ALLOWED";
	    $gInternalBulgeRule = "NOT ALLOWED";
	}
    }

    if (defined($opts{'b'})) {
	if ($opts{'b'} == 0) {
	    $align_obj->internal_bulge_possible(FALSE);
	    $gInternalBulgeRule = "NOT ALLOWED";
	}
    }

    if (defined($opts{'L'})) {
	if ($opts{'L'} == 0) {
#	    $gTerminalLoopPenalty = $align_obj->BIG_NUM;
#	    $align_obj->terminal_loop_penalty($gTerminalLoopPenalty);
	}
    }

    if (defined($opts{'l'})) {
	if ($opts{'l'} == 0) {
	    $align_obj->loop_possible(FALSE);
	    $gLoopRule = "NOT ALLOWED";
	}
    }

    if (defined($opts{'t'})) {
	$gTemperature = $opts{'t'};
	$align_obj->set_temp($gTemperature);
    }

    if (defined($opts{'p'})) {
	$gDoubletID = $opts{'p'};
    }

    if ($gDoubletID eq "FREIER") {
	$align_obj->select_parameters($align_obj->FREIER);
    } elsif ($gDoubletID eq "SANTALUCIA") {
	$align_obj->select_parameters($align_obj->SANTALUCIA);
    } elsif ($gDoubletID eq "XIA_MATHEWS") {
	$align_obj->select_parameters($align_obj->XIA_MATHEWS);
    } else {
	print STDERR "ERROR: gDoubletID : ".$gDoubletID." undefined\n";
	exit;
    }

    if (defined($opts{'f'})) {
	# the user is "forcing" an alignment between the two strands...
	$gForceBind = TRUE;
    }

    if (defined($opts{'o'})) {
	# the user wants a helix only alignment between the two strands...
	# force bind has precedence over helix only...
	if (!$gForceBind) {
	    $gHelixOnly = TRUE;
	}
    }

    if (defined($opts{'T'})) {
	$gTest = TRUE;
    }

    # get the first sequence
    my $seq1_file = shift(@ARGV);
    if (!defined($seq1_file)) {
	print STDERR "\nERROR: you must specify two sequences\n\n";
        print STDERR $sHelp;
	exit;
    }
    my $seq1;
    my $seq1_accession;
    my $seq1_comments;
    # check to see if seq1 is a file, and if so, extract the sequence...
    if (-f $seq1_file) {
	($seq1, $seq1_accession, $seq1_comments) = 
	    $align_obj->load_fasta($seq1_file);
    } else {
	$seq1 = $seq1_file;
    }
    $align_obj->dna2rna(\$seq1);
    $seq1 = uc($seq1);
    
    # get the second sequence
    my $seq2_file = shift(@ARGV);
    if (!defined($seq2_file)) {
	print STDERR "\nERROR: you must specify two sequences\n\n";
        print STDERR $sHelp;
	exit;
    }
    my $seq2;
    my $seq2_accession;
    my $seq2_comments;
    # check to see if seq2 is a file, and if so, extract the sequence...
    if (-f $seq2_file) {
	($seq2, $seq2_accession, $seq2_comments) = 
	    $align_obj->load_fasta($seq2_file);
    } else {
	$seq2 = $seq2_file;
    }
    $align_obj->dna2rna(\$seq2);
    $seq2 = uc($seq2);

    my $init_penalty = $align_obj->init_penalty();

    if (!$gEnergyOnly) {
	# start printing out header stuff...
	print "Doublet Binding Parameters: ".$gDoubletID."\n";
	print "Internal Loop Parameters: ".$gLoopID."\n";
	if ($gForceBind) {
	    print "Force Bind: ON\n";
	} elsif ($gHelixOnly) {
	    print "Helix Only: ON\n";
	} else {
	    print "Terminal Bulges: ".$gTerminalBulgeRule."\n";
	    print "Internal Bulges: ".$gInternalBulgeRule."\n";
#	    print "Terminal Loop Penalty: ".$gTerminalLoopPenalty."\n";
	    print "Internal Loops: ".$gLoopRule."\n";
	}
	print "Binding Temperature: ".$gTemperature." C\n";

	if ($gTest) {
	    print "Test: ON\n";
	}
    }

    my @helix_matrix;
    my @loop_matrix;
    my @loop_length_matrix;
    my @bulge_x_matrix;
    my @bulge_x_length_matrix;
    my @bulge_y_matrix;
    my @bulge_y_length_matrix;

    my $lowest_energy;
    my $seq1_pos;
    my $seq2_pos;
    my $matrix_id;

    my $helix_length;
    my $helix_start_pos;	

    if ($gHelixOnly) {
	($lowest_energy, $seq1_pos, $seq2_pos) =
	    $align_obj->helix_only($seq1, $seq2, \@helix_matrix);
	
    } elsif ($gTest) {
	($lowest_energy) =
	    $align_obj->test($seq1, $seq2);
	

	my $free_energy = $align_obj->total_free_energy();
	if ($gEnergyOnly) {
	    if ($free_energy < 0) {
		print $free_energy."\n";
	    } else {
		print "0\n";
	    }
	    exit;
	}
	
	print "\nDelta-G for best pairing: ".$lowest_energy." + ".$init_penalty." = ".$free_energy."\n";
	
	exit;
		
    } elsif ($gForceBind) {
	($lowest_energy, $helix_start_pos, $helix_length) =
	    $align_obj->force_bind($seq1, $seq2);
	
	if (!defined($lowest_energy)) {
	    print STDERR "free_align->force_bind() failed...\n";
	    exit;
	}

	my $seq1_bind = substr($seq1, $helix_start_pos, $helix_length);
       	my $seq2_bind = substr($seq2, $helix_start_pos, $helix_length);

	#print STDERR "lowest_energy: $lowest_energy\n";
	#print STDERR "helix_start_pos: $helix_start_pos\n";
	#print STDERR "helix_length: $helix_length\n";
	#print STDERR "seq1_bind: $seq1_bind\n";
	#print STDERR "seq2_bind: $seq2_bind\n";

	$seq1_pos = $helix_start_pos + $helix_length;
	$seq2_pos = $seq1_pos;

    } else {
    
	# fill the dynamic programming matrices to find the pairing of
	# the two strands that gives off the most free energy.
	($lowest_energy, $seq1_pos, $seq2_pos, $matrix_id) = 
	    $align_obj->fill_matrices($seq1, $seq2,
				      \@helix_matrix, 
				      \@loop_matrix, 
				      \@loop_length_matrix,
				      \@bulge_x_matrix, 
				      \@bulge_x_length_matrix,
				      \@bulge_y_matrix, 
				      \@bulge_y_length_matrix);
	
    }

    if ($gDebugMode) {
	
	my $matrix_seq1 = " ".$seq1;
	my $seq1_mlength = length($matrix_seq1);
	my @seq1_array = split(//, $matrix_seq1);
	
	my $matrix_seq2 = " ".$seq2;
	my $seq2_mlength = length($matrix_seq2);
	my @seq2_array = split(//, $matrix_seq2);
	
	print STDERR "helix_matrix:\n";
	$align_obj->print_matrix(\@helix_matrix, 
				 $seq1_mlength, $seq2_mlength, 
				 \@seq1_array, \@seq2_array, TRUE);
	
	print STDERR "loop_matrix:\n";
	$align_obj->print_matrix(\@loop_matrix, $seq1_mlength, 
				 $seq2_mlength, 
				 \@seq1_array, \@seq2_array, TRUE);

	print STDERR "loop_length_matrix:\n";
	$align_obj->print_matrix(\@loop_length_matrix, $seq1_mlength, 
				 $seq2_mlength, 
				 \@seq1_array, \@seq2_array, TRUE);
	
	print STDERR "bulge_x_matrix:\n";
	$align_obj->print_matrix(\@bulge_x_matrix, $seq1_mlength, 
				 $seq2_mlength,
				 \@seq1_array, \@seq2_array, TRUE);
	
#	print STDERR "bulge_x_length_matrix:\n";
#	$align_obj->print_matrix(\@bulge_x_length_matrix, $seq1_mlength, 
#				 $seq2_mlength,
#				 \@seq1_array, \@seq2_array, TRUE);

	print STDERR "bulge_y_matrix:\n";
	$align_obj->print_matrix(\@bulge_y_matrix, $seq1_mlength, 
				 $seq2_mlength,
				 \@seq1_array, \@seq2_array, TRUE);
	    
#	print STDERR "bulge_y_length_matrix:\n";
#	$align_obj->print_matrix(\@bulge_y_length_matrix, $seq1_mlength, 
#				 $seq2_mlength,
#				 \@seq1_array, \@seq2_array, TRUE);
	
    }	

    # get the score for the best binding...
    my $free_energy = $align_obj->total_free_energy();
    if ($gEnergyOnly) {
	if ($free_energy < 0) {
	    print $free_energy."\n";
	} else {
	    print "0\n";
	}
	exit;
    }

    print "\nDelta-G for best pairing: ".$lowest_energy." + ".$init_penalty." = ".$free_energy."\n";

    my $seq1_bind;
    my $seq2_bind;
    my $pairs;

    if ($gHelixOnly) {
	($seq1_bind, $seq2_bind, $pairs) =
	    $align_obj->helix_only_trace_route(substr($seq1, 0, $seq1_pos), 
					       substr($seq2, 0, $seq2_pos),
					       \@helix_matrix);

    } elsif($gForceBind) {
	$seq1_bind = substr($seq1, $helix_start_pos, $helix_length);
       	$seq2_bind = substr($seq2, $helix_start_pos, $helix_length);
	$pairs = "|"x$helix_length;
	
    } else {
	
	($seq1_bind, $seq2_bind, $pairs) = 
	    $align_obj->trace_route(substr($seq1, 0, $seq1_pos), 
				    substr($seq2, 0, $seq2_pos),
				    \@helix_matrix,
				    \@loop_matrix, \@loop_length_matrix,
				    \@bulge_x_matrix, \@bulge_x_length_matrix,
				    \@bulge_y_matrix, \@bulge_y_length_matrix,
				    $matrix_id);
    }

    my $length1 = $seq1_bind;
    $length1 =~ s/[\s\-]//g;
    $length1 = length($length1);

    my $length2 = $seq2_bind;
    $length2 =~ s/[\s\-]//g;
    $length2 = length($length2);

    my $seq1_start_pos = $seq1_pos - $length1;
    my $seq2_start_pos = $seq2_pos - $length2;

    my $binding_length = $pairs;

    # trim off leading and trailing white space
    $binding_length =~ s/\A\s*//;
    $binding_length =~ s/\s*\Z//;
    $binding_length = length($binding_length);

    print "\nLength of bound region between sequences: ".$binding_length."\n";

    print "\nseq1 binding begins on base index: ".$seq1_start_pos."\n";
    print "seq2 binding begins on base index: ".$seq2_start_pos."\n";

    print "\n| = pair is a member of a helix\n";
    print "0 = pair is a member of a loop\n";
    print "P = bulge in top strand\n";
    print "b = bulge in bottom strand\n\n";


    ################################################
    #  from here on, all we are doing is tinkering with the output a little
    #  bit by adding additional characters to the beginning or the end of
    #  the alignments so that they are a little more clear
    ################################################
    my $seq1_length = length($seq1);
    my $seq2_length = length($seq2);    

    # any characters in the second half of seq1 not involved in binding
    # should be in lower case...
    my $seq1_added = 0;
    if ($seq1_pos < $seq1_length) {

	$seq1_bind = $seq1_bind.lc(substr($seq1, $seq1_pos));
	$seq1_added = $seq1_length - $seq1_pos;
    }

    # any characters in the second half of seq2 not involved in binding
    # should be in lower case...
    my $seq2_added = 0;
    if ($seq2_pos < $seq2_length) {
	$seq2_bind = $seq2_bind.lc(substr($seq2, $seq2_pos));
	$seq2_added = $seq2_length - $seq2_pos;
    }

    # pad the ends of alignment strings with blank spaces so that they will
    # all have the same length and can be printed out correctly.
    if ($seq1_added > $seq2_added) {
	my $difference = $seq1_added - $seq2_added;
	$seq2_bind = $seq2_bind." "x$difference;
	$pairs = $pairs." "x($difference+$seq2_added);
    } elsif ($seq2_added > $seq1_added) {
	my $difference = $seq2_added - $seq1_added;
	$seq1_bind = $seq1_bind." "x$difference;
	$pairs = $pairs." "x($difference+$seq1_added);
    } else {
	$pairs = $pairs." "x$seq1_added;
    }

    my $just_bases_1 = $seq1_bind;
    $just_bases_1 =~ s/[\s\-]//g;

    my $just_bases_2 = $seq2_bind;
    $just_bases_2 =~ s/[\s\-]//g;

    my $bind1_length = length($just_bases_1);
    my $bind2_length = length($just_bases_2);

    # any characters in the first half of seq1 not involved in binding
    # should be in lower case...
    $seq1_added = 0;
    if ($bind1_length < $seq1_length) {
	$seq1_added = $seq1_length - $bind1_length;
	$seq1_bind = lc(substr($seq1, 0, $seq1_added)).$seq1_bind;
    }

    # any characters in the first half of seq2 not involved in binding
    # should be in lower case...
    $seq2_added = 0;
    if ($bind2_length < $seq2_length) {
	$seq2_added = $seq2_length - $bind2_length;
	$seq2_bind = lc(substr($seq2, 0, $seq2_added)).$seq2_bind;
    }

    # pad the beginings of the alignment sequences with white spaces 
    # so that everything will be aligned correctly when printed out.
    if ($seq1_added > $seq2_added) {
	my $difference = $seq1_added - $seq2_added;
	$seq2_bind = " "x$difference.$seq2_bind;
	$pairs = " "x($difference+$seq2_added).$pairs;
    } elsif ($seq2_added > $seq1_added) {
	my $difference = $seq2_added - $seq1_added;
	$seq1_bind = " "x$difference.$seq1_bind;
	$pairs = " "x($difference+$seq1_added).$pairs;
    } else {
	$pairs = " "x$seq1_added.$pairs;
    }

    $align_obj->print_binding($seq1_bind, "seq1", 
			      $seq2_bind, "seq2", $pairs);

    exit 1;
}
