package FreeAlign;

##========================================================================
##
## This is a perl module or class that implements some common 
## subroutines that free energy calculations by dynamic programming may want
## to use
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
## Values for the various parameters are from:
##
## Freier et al., PNAS (1986) Vol. 83, pp. 9373-9377
##
## Jaeger et al., PNAS (1989) Vol. 86, pp. 7706-7710
##
## SantaLucia, John, PNAS (1998) Vol. 95, pp. 1460-1465
##
## Xia et al., Biochemistry (1998), Vol. 37, pp. 14719-14735
##
## Mathews et al., J. Mol. Biol. (1999), Vol. 288, pp. 911-940
##
## Hopefully each parameter used in the program will have a comment that
## will help determine exactly which publication the value came from.  If there
## is no comment, you can infer that the value is arbitrary and made up.
## 
##========================================================================

use warnings;
use strict;


##========================================================================
#
# DEFINE GLOBAL CONSTANTS
#
##========================================================================

use constant TRUE => 1;
use constant FALSE => 0;
use constant FAIL => -1;


# id numbers for various matrices...
use constant HELIX_ID => 0;
use constant LOOP_ID => 1;
use constant BULGE_X_ID => 2;
use constant BULGE_Y_ID => 3;

use constant LEFT_DANGLE => '-2';
use constant RIGHT_DANGLE => '-1';


#use constant DEBUG => TRUE;
use constant DEBUG => FALSE;

use constant BIG_NUM => 30000;
use constant SMALL_NUM => -30000;

use constant HELIX_CHAR => "|";
use constant BULGE_X_CHAR => "P";
use constant BULGE_Y_CHAR => "b";
use constant LOOP_CHAR => "O";


#
# CONSTANTS FOR PRETTY PRINTING VARIOUS THINGS
#
use constant MATRIX_CELL_WIDTH => 8; # good for floating point numbers...
use constant LINE_LENGTH => 77; # good for xterms
#use constant LINE_LENGTH => 60; # good for email


#
# CONSTANTS FOR SCORING THE BINDING/ALIGNMENT:
#

use constant FREIER => 1; # RNA parameters (1986)
use constant SANTALUCIA => 2; # DNA (NOT RNA) parameters (1998)
use constant XIA_MATHEWS => 3; # RNA parameters(1998)

# NOTE: FREIER_INIT_PENALTY is a one time penalty for forming any bonds between
# two different RNA strands (Lewin, Genes VI, page 102)
use constant FREIER_INIT_PENALTY => 3.4; # Freier (1986)

# NOTE" HELIX_START_GC and HELIX_START_AT are both for DNA hybridization
# and should be added on for both ends.  i.e. If the helix starts with 
# a G/C pair and ends with a G/C pair, then you should use HELIX_START_GC
# twice.
use constant HELIX_START_GC => 0.98; # SantaLucia (1998)
use constant HELIX_START_AT => 1.03; # SantaLucia (1998)

# single_stranded_end corresponds to a terminal gap penalty in alignment
# terminology.
use constant SINGLE_STRANDED_END => 0;


use constant MIN_LOOP_STEP => 2;


##========================================================================
#
# DEFINE GLOBAL VARIABLES
#
##========================================================================

# these global variables are needed to keep track of the sequence lengths
# for the "trace_route" routine
my $gSeqXLength = 0;
my $gSeqYLength = 0;

# this keeps track of the best score for binding two sequences together.
my $gBestScore = 0;

# keep track of whether or not bulges are allowed
my $gInternalBulgePossible = TRUE; # by default we will allow them
my $gNoBulgePossible = FALSE; # by default we will allow (terminal/internal)

# keep track of whether or not loops are allowed
my $gLoopPossible = TRUE; # by default we will allow them

# keep track of which set of binding parameter values are being used...
#my $gParameterID = FREIER;
my $gParameterID;

#my $gInitPenalty = FREIER_INIT_PENALTY; # the initial helix formation penalty
my $gInitPenalty; # the initial helix formation penalty

my $gTemperature = 37; # (Celsius) The temperature at which the binding is 
# supposed to be taking place.  Only parameters based on delta_H and
# delta_S take this into account.                


##========================================================================
##
## SUBROUTINE new() - a constructor for a "FreeAlign" object.
##
##========================================================================
sub new {
    my $class = shift;
    my $self = {};

    bless($self, $class);

    return $self;
}


##========================================================================
##
## SUBROUTINE select_parameters() - enables you specify which set of
##   binding parameters that you wish to use.
##
##   ARGUMENTS:
##   $parameter_id is the idenifier for the set of parameters you wish to
##     use.  These identifiers are defined as constants.
##
##========================================================================
sub select_parameters {
    my ($self, $parameter_id) = @_;
    if ($parameter_id == FREIER) {
	$gParameterID = FREIER;
	$gInitPenalty = FREIER_INIT_PENALTY;
    } elsif ($parameter_id == SANTALUCIA) {
	$gParameterID = SANTALUCIA;
    } elsif ($parameter_id == XIA_MATHEWS) {
	$gParameterID = XIA_MATHEWS;

	# set the helix initiation penalty...
	my $deltaH = 3.61; # Xia (1998)
	my $deltaS = -1.5; # Xia (1998)
	
	$gInitPenalty = $deltaH - ($gTemperature + 273.15)*($deltaS/1000);
    } else {
	print STDERR "ERROR: FreeAlign::select_parameters()\n";
	print STDERR "       There is no set of parameters specified with\n";
	print STDERR "       parameter_id: $parameter_id\n";
    }
}


##========================================================================
##
## SUBROUTINE intneral_bulge_possible() - enables you to allow/disallow 
##   internal bulges in optimal bindings/alignments.
##
##   ARGUMENTS:
##   bool - a boolean value that is eather True (!= 0) or False (= 0)
##
##========================================================================
sub internal_bulge_possible {
    my ($self, $bool) = @_;
    
    $gInternalBulgePossible = $bool;
}


##========================================================================
##
## SUBROUTINE no_bulge_possible() - enables you to allow/disallow 
##   terminal and internal bulges inoptimal bindings/alignments.
##
##   ARGUMENTS:
##   bool - a boolean value that is eather True (!= 0) or False (= 0)
##
##========================================================================
sub no_bulge_possible {
    my ($self, $bool) = @_;
    
    $gNoBulgePossible = $bool;
    $gInternalBulgePossible = !($bool);
}


##========================================================================
##
## SUBROUTINE loop_possible() - enables you to allow/disallow loops in 
##   optimal bindings/alignments.
##
##   ARGUMENTS:
##   bool - a boolean value that is eather True (!= 0) or False (= 0)
##
##========================================================================
sub loop_possible {
    my ($self, $bool) = @_;

    $gLoopPossible = $bool;
}


##========================================================================
##
## SUBROUTINE terminal_loop_penalty() - enables you to set the value for
##   terminal loops (unmatched 5' or 3' ends).  The default value is 0.
##
##   ARGUMENTS:
##   penalty - the penalty for terminal loops.
##
##========================================================================
sub terminal_loop_penalty {
    my ($self, $penalty) = @_;
    
#    $gTerminalLoopPenalty = $penalty;
}


##========================================================================
##
## SUBROUTINE init_penalty() returns the initial penalty for forming a helix
##
##   RETURN VALUES:
##   - the intial penalty for forming a helix between to strands of nucleotides
##
##========================================================================
sub init_penalty {
    return $gInitPenalty;
}


##========================================================================
##
## SUBROUTINE total_free_energy() returns the total amount of free energy
##   given off by allowing two sequences of nucleotides to bind to each other.
##
##   RETURN VALUES:
##   - the total free energy given off by allowing two sequences to bind
##
##========================================================================
sub total_free_energy {
    if (($gParameterID == FREIER) || ($gParameterID == XIA_MATHEWS)) {
	return ($gBestScore + $gInitPenalty);
    } else {
	print STDERR "Warning: FreeAlign::total_free_energy()\n";
	print STDERR "         This function has not yet been modified to\n";
	print STDERR "         function with gParameterID: $gParameterID\n";
    }
}


##========================================================================
##
## SUBROUTINE set_seq_lengths() sets $gSeqXLength and $gSeqYLength.
##
##   ARGUMENTS:
##     x_length - the length of a sequence that defines the x-axis of the
##       dynamic programming matrices.
##
##     y_length - the length of a sequence that defines the y-axis of the
#        dynamic programming matrices.
##
##========================================================================
sub set_seq_lengths {
    my ($self, $x_length, $y_length) = @_;
    
    $gSeqXLength = $x_length;
    $gSeqYLength = $y_length;
}

##========================================================================
##
## SUBROUTINE set_temp() - enables you to set the temperature at which
##   binding parameters are calculated (only used for binding parameters
##   that are calculated on the fly...).
##
##   ARGUMENTS:
##   temp - the temperature, in Celsius.
##
##========================================================================
sub set_temp {
    my ($self, $temp) = @_;
    
    $gTemperature = $temp;
}


##========================================================================
##
## SUBROUTINE max() computes and returns the maximum of any number of 
##   parameters
##
##   ARGUMENTS:
##   - a comma separated list of values
##
##   RETURN VALUES:
##   most - the maximum value of the arguments.
##
##========================================================================
sub max {
    my ($self, $most, @rest) = @_;
    
    foreach my $x (@rest) {
	if ($x > $most) {
	    $most = $x;
	}
    }
    return $most;
}


##========================================================================
##
## SUBROUTINE min() computes and returns the minimum of any number of
##   parameters
##
##   ARGUMENTS:
##   - a comma separated list of values
##
##   RETURN VALUES:
##   least - the minimum value of the arguments.
##
##========================================================================
sub min {
    my ($self, $least, @rest) = @_;
    
    foreach my $x (@rest) {
	if ($x < $least) {
	    $least = $x;
	}
    }
    return $least;
}


##========================================================================
##
## SUBROUTINE complement(): takes a reference to a sequence (a string) and
##   returns the complement of it.  If string is DNA, then you get the 
##   following substitutions: A->T, T->A, G->C and C->G (if the string is
##   RNA, U is used instead of T)
##
##   ARGUMENTS:
##   seq_ref - a reference to a sequence (a reference to a string)
##
##   is_rna - a boolean that decides whether T->A or U->A etc.
##
##========================================================================
sub complement {
    my ($self, $seq_ref, $is_rna) = @_;
    
    if (!$is_rna) {
	$$seq_ref =~ tr/atcgATCG/tagcTAGC/;
    } else {
	$$seq_ref =~ tr/aucgAUCG/uagcUAGC/;
    }
}


##========================================================================
##
## SUBROUTINE dna2rna(): takes a reference to a DNA sequence (a string) and
##   returns the RNA version of it.
##
##   ARGUMENTS:
##   seq_ref - a reference to a DNA sequence (a reference to a string)
##
##========================================================================
sub dna2rna {
    my ($self, $seq_ref) = @_;

    $$seq_ref =~ tr/tT/uU/;
}


##========================================================================
##
## SUBROUTINE internal_loop(): calculates and returns the free energy for an
##   internal loop.  This is similar to a mismatch penalty in a normal
##   alignment.
##
##
##   ARGUMENTS:
##   length - the number of nucleotides in the loop
##
##   RETURN VALUES:
##   - a free energy value (a score) for the loop.
##
##========================================================================
sub internal_loop {
    my ($self, $length) = @_;

    # are we allowed to consider loops in an optimal binding?    
    if (!$gLoopPossible) {
	return BIG_NUM;
    }

    if ($length == 0) {
	# if the loop is length 0, then there is no loop to score...
	return 0;
    }

    if ($length == 1) {
	# if there is only one base in the loop, then this is really a single
	# base bulge...
	print STDERR "WARNING: FreeAlign::internal_loop()\n";
	print STDERR "         loop length is 1 (so it's not a real loop)\n";
	return BIG_NUM;
    }

    if ($gParameterID == FREIER) {
	return $self->loop_jaeger($length);
    } elsif ($gParameterID == SANTALUCIA) {
	return $self->loop_santalucia($length);
    } elsif ($gParameterID == XIA_MATHEWS) {
	return $self->loop_jaeger($length);
	#return $self->loop_xia_mathews($length);
    } else {
	print STDERR "ERROR: FreeAlign::internal_loop()\n";
	print STDERR "       gParameterID: $gParameterID, undefined\n";
	return BIG_NUM;
    }
}


##========================================================================
##
## SUBROUTINE loop_jaeger(), parameters values are from Jaeger (1989),
##
##   See "internal_loop" for more detials.
##
##   ARGUMENTS:
##   length - the number of nucleotides in the loop
##
##   RETURN VALUES:
##   - a free energy value (a score) for the loop.
##
##========================================================================
sub loop_jaeger {
    my ($self, $length) = @_;


    # NOTE: THESE SHORT LOOPS (with fewer than 10 bases) ARE SCORED
    # INDEPENDENTLY OF THE TEMPERATURE....  PERHAPS WE SHOULD PUT A
    # SWITCH IN THAT MAKES EXCLUSIVE USE OF THE FORMULA AT THE BOTTOM
    # WHEN THE TEMPERATURE IS NOT SET TO 37.

    if ($length == 2) {
	return 4.1; # Jaeger (1989) - experimentally derived
    } elsif ($length == 3) {
	return 4.5; # Jaeger (1989) - experimentally derived
    } elsif ($length == 4) {
	return 4.9; # Jaeger (1989) - experimentally derived
    } elsif ($length == 5) {
	return 5.1; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 6) {
	return 5.7; # Jaeger (1989) - experimentally derived
    } elsif ($length == 7) {
	return 5.9; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 8) {
	return 6.0; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 9) {
	return 6.1; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 10) {
	return 6.3; # Jaeger (1989) - not experimentally derived
    } else {
	# Formula from Jaeger (1989): NOTE in the formula in the publication
	# needs to be converted to return Kcal per mol (thus, the division by
	# 1000 is included here).
	return 5.7 + (1.75 * 1.987/1000 * (273.15 + $gTemperature) 
		      * log($length/6)); 
    }
}


##========================================================================
##
## SUBROUTINE loop_mathews(), parameters values are from Mathews (1999),
##
##   See "internal_loop" for more detials.
##
##   ARGUMENTS:
##   length - the number of nucleotides in the loop
##
##   RETURN VALUES:
##   - a free energy value (a score) for the loop.
##
##========================================================================
sub loop_mathews {
    my ($self, $n1, $n2, 
	$open_t5, $open_b3, $close_t3, $close_b5,
	$first_mm_t5, $first_mm_b3, $last_mm_t3, $last_mm_b5) = @_;

    my $length = $n1 + $n2;

    my $Delta_G_init;

    my $asymm_diff = abs($n1 - $n2);
    my $asymm_penalty = $asymm_diff * 0.48;
    
    my $au_gu_penalty = 0;
    if ((($open_t5 eq "G") && ($open_b3 eq "U")) ||	
	(($open_t5 eq "U") && ($open_b3 eq "G")) ||
	(($open_t5 eq "A") && ($open_b3 eq "U")) ||
	(($open_t5 eq "U") && ($open_b3 eq "A"))) {
	$au_gu_penalty = $au_gu_penalty + 0.2;
    }

    my $ga_ag_bonus = 0;
    if ((($first_mm_t5 eq "G") && ($first_mm_b3 eq "A")) ||
	(($first_mm_t5 eq "A") && ($first_mm_b3 eq "G"))) {
	$ga_ag_bonus = $ga_ag_bonus + (-1.7);
    }    

    my $uu_bonus = 0;
    if (($first_mm_t5 eq "U") && ($first_mm_b3 eq "U")) {
	$uu_bonus = $uu_bonus + (-0.7);
    }

    if (($n1 > 1) && ($n2 > 1)) {
	if ((($close_t3 eq "G") && ($close_b5 eq "U")) ||	
	    (($close_t3 eq "U") && ($close_b5 eq "G")) ||
	    (($close_t3 eq "A") && ($close_b5 eq "U")) ||
	    (($close_t3 eq "U") && ($close_b5 eq "A"))) {
	    $au_gu_penalty = $au_gu_penalty + 0.2;
	}

	if ((($last_mm_t3 eq "G") && ($last_mm_b5 eq "A")) ||
	    (($last_mm_t3 eq "A") && ($last_mm_b5 eq "G"))) {
	    $ga_ag_bonus = $ga_ag_bonus + (-1.7);
	}
	
	if (($last_mm_t3 eq "U") && ($last_mm_b5 eq "U")) {
	    $uu_bonus = $uu_bonus + (-0.7);
	}
    }

    
    if ($length == 4) {
	$Delta_G_init = 1.7;
    } elsif ($length == 5) {
	$Delta_G_init = 1.8;
    } elsif ($length == 6) {
	$Delta_G_init = 2.0;
    } else {
	# NOTE: this formula is essentially the same one used in Jaeger (1989)
	$Delta_G_init = 2.0 + (1.75 * 1.987/1000 * (273.15 + $gTemperature) 
			       * log($length/6));
    }

    my $loop_penalty = $Delta_G_init + $asymm_penalty + $au_gu_penalty  
	+ $ga_ag_bonus + $uu_bonus;

}


##========================================================================
##
## SUBROUTINE bulge_energy(): calculates and returns the free energy for a 
##   bulge.  This is similar to a gap penalty in a normal alignment.
##
##   NOTE:  When scoring by hand, remember that the bases in the bulge
##   are ignored and you calculate an "doublet" for the surrounding
##   bases.
##   example:           +3.9
##                        |    
##            G   G   G   A   G   G   G
##            |   |   |   P   |   |   |
##            C   C   C   -   C   C   C
##             |    |     |     |    |
##           -2.9 -2.9  -2.9  -2.9 -2.9
##
##                       ||
##                       
##                        A
##                   G G G G G G
##                   | | | | | |
##                   C C C C C C
##
##  all of the GG doublets are treated as if in a contiguous helix (which we
##             CC
##
##  giving a score of -2.9 for each doublet (ignoring terminal ends...))
##
##   Total is (-14.5) + 3.9 = -10.6 (not counting the helix init penalty)
##
##   ARGUMENTS:
##   length - the number of nucleotides in the loop
##
##   RETURN VALUES:
##   - a free energy value (a score) for the bulge.
##
##========================================================================
sub bulge_energy {
    my ($self, $length) = @_;

    # are we allowed to consider bulges in an optimal binding?
    if (!$gInternalBulgePossible) {
	return BIG_NUM;
    }

    if ($length == 0) {
	return 0;
    }

    if ($length == 1) {
	return 3.9; # Jaeger (1989) - experimentally derived
    } elsif ($length == 2) {
	return 3.1; # Jaeger (1989) - experimentally derived
    } elsif ($length == 3) {
	return 3.5; # Jaeger (1989) - experimentally derived
    } elsif ($length == 4) {
	return 4.2; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 5) {
	return 4.8; # Jaeger (1989) - experimentally derived
    } elsif ($length == 6) {
	return 5.0; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 7) {
	return 5.2; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 8) {
	return 5.3; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 9) {
	return 5.4; # Jaeger (1989) - not experimentally derived
    } elsif ($length == 10) {
	return 5.5; # Jaeger (1989) - not experimentally derived
    } else {
	# Formula from Jaeger (1989): NOTE in the formula in the publication
	# needs to be converted to return Kcal per mol (thus, the division by
	# 1000 is included here).
	return 4.8 + (1.75 * 1.987/1000 * 310.15 * log($length/5)); 
    }
}


##========================================================================
##
## SUBROUTINE internal_doublet() scores a "doublet" of RNA residues based on 
##   how well they would bind, assuming that the doublet is found inside a
##   helix structure (thus, this only considers Watson/Crick paris plus G/U
##   pairs, see NOTE) and returns the score.   A "doublet" consists of two 
##   pairs of RNA residues.
##
##   NOTE: This only scores matches or internal (to a helix structure, i.e.
##     not in a loop or bulge) G/U mismatches.  To score internal mismatches 
##     then you are actually scoring a loop or bulge...  To score terminal
##     mismatches, use "terminal_doublet()".
##
##   EXAMPLES: GC would get one score and AU would get another...
##             CG                         UG
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub internal_doublet {
    my ($self, $t5, $t3, $b3, $b5, $t55, $t33, $b33, $b55) = @_;
    
    if ($gParameterID == FREIER) {
	return $self->doublet_freier($t5, $t3, $b3, $b5);
    } elsif ($gParameterID == SANTALUCIA) {
	return $self->doublet_santalucia($t5, $t3, $b3, $b5);
    } elsif ($gParameterID == XIA_MATHEWS) {
	return $self->doublet_xia_mathews($t5, $t3, 
					  $b3, $b5,
					  $t55, $t33,
					  $b33, $b55);
    } else {
	print STDERR "ERROR: FreeAlign::internal_doublet()\n";
	print STDERR "       gParameterID: $gParameterID, undefined\n";
	return BIG_NUM;
    }
}


##========================================================================
##
## SUBROUTINE doublet_freier(), parameter values are from Freier
##   (1986).
##
##   See "doublet" for more details.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub doublet_freier {
    my ($self, $t5, $t3, $b3, $b5) = @_;

    if (!($self->valid_pair($t5, $b3) && $self->valid_pair($t3, $b5))) {
	return BIG_NUM;
    }

    if (($t5 eq "A") && ($b3 eq "U")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # AA
	    return -0.9; # Freier (1986)      UU	
	} 
	if (($t3 eq "U") && ($b5 eq "A")) { # AU
	    return -0.9; # Freier (1986)      UA	
	} 
	if (($t3 eq "G") && ($b5 eq "C")) { # AG
	    return -1.7; # Freier (1986)      UC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # AC
	    return -2.1; # Freier (1986)      UG
	}

	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # AG
	    return -0.5; # Freier (1986)      UU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # AU
	    return -0.7; # Freier (1986)      UG
	}
    }

    if (($t5 eq "U") && ($b3 eq "A")) {
	# watson/crick matches...
	if (($t3 eq "U") && ($b5 eq "A")) { # UU
	    return -0.9; # Freier (1986)      AA
	}
	if (($t3 eq "A") && ($b5 eq "U")) { # UA
	    return -1.1; # Freier (1986)      AU
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # UG
	    return -1.8; # Freier (1986)      AC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # UC
	    return -2.3; # Freier (1986)      AG
	}

	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # UG
	    return -0.7; # Freier (1986)      AU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # UU
	    return -0.5; # Freier (1986)      AG - not experimentally derived
	}
    }

    if (($t5 eq "C") && ($b3 eq "G")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # CA
	    return -1.8; # Freier (1986)      GU
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # CU
	    return -1.7; # Freier (1986)      GA
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # CC
	    return -2.9; # Freier (1986)      GG
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # CG
	    return -2.0; # Freier (1986)      GC
	}

	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # CG
	    return -1.5; # Freier (1986)      GU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # CU
	    return -1.5; # Freier (1986)      GG
	}
    }
    
    if (($t5 eq "G") && ($b3 eq "C")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # GA
	    return -2.3; # Freier (1986)      CU
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # GU
	    return -2.1; # Freier (1986)      CA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    return -2.9; # Freier (1986)      CC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    return -3.4; # Freier (1986)      CG
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # GG
	    return -1.3; # Freier (1986)      CU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # GU
	    return -1.9; # Freier (1986)      CG
	}
    }

    # double G/U mismatches...
    if (($t5 eq "G") && ($b3 eq "U")) {	
	if (($t3 eq "G") && ($b5 eq "U")) { # GG
	    return -0.5; # Freier (1986)      UU - not experimentally derived
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # GU
	    return -0.5; # Freier (1986)      UG - not experimentally derived
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # GU
	    return -0.7; # Freier (1986)      UA
	}
	if (($t3 eq "A") && ($b5 eq "U")) { # GA
	    return -0.5; # Freier (1986)      UU - not experimentally derived
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    return -1.5; # Freier (1986)      UC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    return -1.9; # Freier (1986)      UG
	}
    }

    # more double G/U mismatches...
    if (($t5 eq "U") && ($b3 eq "G")) {	
	if (($t3 eq "G") && ($b5 eq "U")) { # UG
	    return -0.6; # Freier (1986)      GU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # UU
	    return -0.5; # Freier (1986)      GG - not experimentally derived
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # UU
	    return -0.5; # Freier (1986)      GA - not experimentally derived
	}
	if (($t3 eq "A") && ($b5 eq "U")) { # UA
	    return -0.7; # Freier (1986)      GU
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # UG
	    return -1.5; # Freier (1986)      GC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # UC
	    return -1.3; # Freier (1986)      GG
	}
    }
    print STDERR "ERROR: FreeAlign::doublet_freier()\n";
    print STDERR "ERROR: doublet $t5$t3\n";
    print STDERR "               $b3$b5 not scored yet...\n";
    return BIG_NUM;
}


##========================================================================
##
## SUBROUTINE doublet_santalucia(), parameter values are from 
##   SantaLucia (1998)
##
##   NOTE: THESE PARAMETERS ARE FOR DNA ONLY!
##
##   See "doublet" for more details.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub doublet_santalucia {
    my ($self, $t5, $t3, $b3, $b5) = @_;
    
    ## since we are scoring DNA, only score WC pairs...
    if (!($self->wc_pair($t5, $b3) && $self->wc_pair($t3, $b5))) {
	return BIG_NUM;
    }
    
    if (($t5 eq "A") && ($b3 eq "T")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "T")) { # AA
	    return -1;    # SantaLucia (1998) TT
	}
	if (($t3 eq "T") && ($b5 eq "A")) { # AT
	    return -0.88; # SantaLucia (1998) TA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # AG
	    return -1.28; # SantaLucia (1998) TC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # AC
	    return -1.44; # SantaLucia (1998) TG
	}
    }
    
    if (($t5 eq "T") && ($b3 eq "A")) {
	# watson/crick matches...
	if (($t3 eq "T") && ($b5 eq "A")) { # TT
	    return -1;    # SantaLucia (1998) AA
	}
	if (($t3 eq "A") && ($b5 eq "T")) { # TA
	    return -0.58; # SantaLucia (1998) AT
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # TG
	    return -1.45; # SantaLucia (1998) AC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # TC
	    return -1.30; # SantaLucia (1998) AG
	}
    }
    
    if (($t5 eq "C") && ($b3 eq "G")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "T")) { # CA
	    return -1.45; # SantaLucia (1998) GT
	}
	if (($t3 eq "T") && ($b5 eq "A")) { # CT
	    return -1.28; # SantaLucia (1998) GA
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # CC
	    return -1.42; # SantaLucia (1998) GG
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # CG
	    return -2.17; # SantaLucia (1998) GC
	}
    }
    
    if (($t5 eq "G") && ($b3 eq "C")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "T")) { # GA
	    return -1.30; # SantaLucia (1998) CT
	}
	if (($t3 eq "T") && ($b5 eq "A")) { # GT
	    return -1.44; # SantaLucia (1998) CA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    return -1.42; # SantaLucia (1998) CC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    return -1.84; # SantaLucia (1998) CG
	}       
    }
    
    print STDERR "ERROR: FreeAlign::doublet_santalucia()\n";
    print STDERR "ERROR: doublet $t5$t3\n";
    print STDERR "               $b3$b5 not scored yet...\n";
    return BIG_NUM;
    
}


##========================================================================
##
## SUBROUTINE doublet_xia_matthews(), parameter values are from
##   Xia (1998) (for Watson-Crick pairs) and Mathews (1999) (for G/U 
##   mismatches).
##
##   NOTE: These parameter values are generated on the fly...
##
##   See "doublet" for more details.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub doublet_xia_mathews {
    my ($self, $t5, $t3, $b3, $b5, $t55, $t33, $b33, $b55) = @_;
    
    if (!($self->valid_pair($t5, $b3) && $self->valid_pair($t3, $b5))) {
	return BIG_NUM;
    }
    
    my $deltaH;
    my $deltaS;
    
    if (($t5 eq "A") && ($b3 eq "U")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # AA
	    $deltaH = -6.82;  # Xia (1998)    UU
	    $deltaS = -19.0;  # Xia (1998)
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # AU
	    $deltaH = -9.38;  # Xia (1998)    UA
	    $deltaS = -26.7;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # AG
	    $deltaH = -10.48; # Xia (1998)    UC
	    $deltaS = -27.1;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # AC
	    $deltaH = -11.40; # Xia (1998)    UG
	    $deltaS = -29.5;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  AG
	    $deltaH = -3.21; # Mathews (1999)  UU
	    $deltaS = -8.6;  # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  AU
	    $deltaH = -8.81; # Mathews (1999)  UG
	    $deltaS = -24.0; # Mathews (1999)
	}
    }
    
    elsif (($t5 eq "U") && ($b3 eq "A")) {
	# watson/crick matches...
	if (($t3 eq "U") && ($b5 eq "A")) { # UU
	    $deltaH = -6.82;  # Xia (1998)    AA
	    $deltaS = -19.0;  # Xia (1998)
	}
	if (($t3 eq "A") && ($b5 eq "U")) { # UA
	    $deltaH = -7.69;  # Xia (1998)    AU
	    $deltaS = -20.5;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # UG
	    $deltaH = -10.44; # Xia (1998)    AC
	    $deltaS = -26.9;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # UC
	    $deltaH = -12.44; # Xia (1998)    AG
	    $deltaS = -32.5;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  UG
	    $deltaH = -6.99; # Mathews (1999)  AU
	    $deltaS = -19.3; # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  UU
	    $deltaH = -12.83; # Mathews (1999) AG 
	    $deltaS = -37.3;  # Mathews (1999)
	}
    }
    
    elsif (($t5 eq "C") && ($b3 eq "G")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # CA
	    $deltaH = -10.44; # Xia (1998)    GU
	    $deltaS = -26.9;  # Xia (1998)
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # CU
	    $deltaH = -10.48; # Xia (1998)    GA
	    $deltaS = -27.1;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # CC
	    $deltaH = -13.39; # Xia (1998)    GG
	    $deltaS = -32.7;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # CG
	    $deltaH = -10.64; # Xia (1998)    GC
	    $deltaS = -26.7;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  CG
	    $deltaH = -5.61; # Mathews (1999)  GU
	    $deltaS = -13.5; # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  CU
	    $deltaH = -12.11; # Mathews (1999) GG
	    $deltaS = -32.2;  # Mathews (1999)
	}
    }
    
    elsif (($t5 eq "G") && ($b3 eq "C")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # GA
	    $deltaH = -12.44; # Xia (1998)    CU
	    $deltaS = -32.5;  # Xia (1998)
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # GU
	    $deltaH = -11.40; # Xia (1998)    CA
	    $deltaS = -29.5;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    $deltaH = -13.39; # Xia (1998)    CC
	    $deltaS = -32.7;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    $deltaH = -14.88; # Xia (1998)    CG
	    $deltaS = -36.9;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  GG
	    $deltaH = -8.33; # Mathews (1999)  CU
	    $deltaS = -21.9; # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  GU
	    $deltaH = -12.59; # Mathews (1999) CG
	    $deltaS = -32.5;  # Mathews (1999)
	}
    } 

    elsif (($t5 eq "G") && ($b3 eq "U")) {	
	if (($t3 eq "G") && ($b5 eq "U")) { #  GG
	    $deltaH = -13.47; # Mathews (1999) UU
	    $deltaS = -41.82; # Based on Mathews (1999) - see note in Table 4
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  GU
	    $deltaH = -14.59; # Mathews (1999) UG
	    $deltaS = -51.2;  # Mathews (1999)
	    
	    if (defined($t55) && defined($t33) &&     # only scoring middle two
		defined($b33) && defined($b55)) {     #  /\
		if (($t55 eq "G") && ($b33 eq "C") && # GGUC
		    ($t33 eq "C") && ($b55 eq "G")) { # CUGG
		    
		    # Mathews (1999) scores the entire GGUC/CUGG unit with
		    # deltaH = -30.80 and deltaS = -86.0.  The values below
		    # are calculated by... 
		    # deltaH = -30.80-(2*-8.33), -8.33 is the deltaH for GG/CU
		    # deltaS = -86-(2*-21.9), -21.9 is the deltaS for GG/CU
		    
		    $deltaH = -14.14; # Based on Mathews (1999)
		    $deltaS = -42.2;  # Based on Mathews (1999)
		}
	    }	
	}
	if (($t3 eq "U") && ($b5 eq "A")) { #  GU
	    $deltaH = -8.81;  # Mathews (1999) UA
	    $deltaS = -24.0;  # Mathews (1999)	
	}
	if (($t3 eq "A") && ($b5 eq "U")) { #  GA
	    $deltaH = -12.83; # Mathews (1999) UU
	    $deltaS = -37.3;  # Mathews (1999)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { #  GG
	    $deltaH = -12.11; # Mathews (1999) UC
	    $deltaS = -32.2;  # Mathews (1999)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { #  GC
	    $deltaH = -12.59; # Mathews (1999) UG
	    $deltaS = -32.5;  # Mathews (1999)
	}
    }    

    elsif (($t5 eq "U") && ($b3 eq "G")) {	
	if (($t3 eq "G") && ($b5 eq "U")) { #  UG
	    $deltaH = -9.26;  # Mathews (1999) GU
	    $deltaS = -30.8;  # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  UU
	    $deltaH = -13.47; # Mathews (1999) UU
	    $deltaS = -41.82; # Based on Mathews (1999) - see note in Table 4
	}
	if (($t3 eq "U") && ($b5 eq "A")) { #  UU
	    $deltaH = -3.21;  # Mathews (1999) GA
	    $deltaS = -8.6;   # Mathews (1999)
	}
	if (($t3 eq "A") && ($b5 eq "U")) { #  UA
	    $deltaH = -6.99;  # Mathews (1999) GU
	    $deltaS = -19.3;  # Mathews (1999)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { #  UG
	    $deltaH = -5.61;  # Mathews (1999) GC
	    $deltaS = -13.5;  # Mathews (1999)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { #  UC
	    $deltaH = -8.33;  # Mathews (1999) GG
	    $deltaS = -21.9;  # Mathews (1999)
	}
    }

    if (!defined($deltaH) || !defined($deltaS)) {
	print STDERR "ERROR: FreeAlign::doublet_xia_mathews()\n";
	print STDERR "ERROR: doublet $t5$t3\n";
	print STDERR "               $b3$b5 not scored yet...\n";
	return BIG_NUM;
    }

#    print STDERR "deltaH: $deltaH\n";
#    print STDERR "deltaS: $deltaS\n";

    return $deltaH - ($gTemperature + 273.15)*($deltaS/1000);
}



##========================================================================
##
## SUBROUTINE terminal_doublet() scores a "doublet" of RNA residues based on 
##   how well they would bind, assuming that the doublet is found at a
##   terminal end of a helix structure and returns the score.
##   A "doublet" consists of two pairs of RNA residues.
##
##   NOTE: TERMINAL PAIRS MUST BE AT THE VERY ENDS OF BOTH STRANDS
##
##   EXAMPLES: GC would get one score and AU would get another...
##             CG                         UG
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   left_side - is the terminal end of the helix on the left side,
##     or on the right side?
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_doublet {
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;

    if ($gParameterID == FREIER) {
	return $self->terminal_doublet_freier($t5, $t3, $b3, $b5, $left_side);
    } elsif ($gParameterID == SANTALUCIA) {
	return $self->terminal_doublet_santalucia($t5, $t3, 
						  $b3, $b5, $left_side);
    } elsif ($gParameterID == XIA_MATHEWS) {
	return $self->terminal_doublet_xia_mathews($t5, $t3, 
						   $b3, $b5, $left_side);
    } else {
	print STDERR "ERROR: FreeAlign::terminal_doublet()\n";
	print STDERR "       gParameterID: $gParameterID, undefined\n";
	return BIG_NUM;
    }
}

##========================================================================
##
## SUBROUTINE terminal_doublet_xia_mathews(), parameter values are from
##   Xia (1998) and Mathews (1999)
##
##   See "terminal_dobulet" for more detials.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   left_side - is the terminal end of the helix on the left side,
##     or on the right side?
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_doublet_xia_mathews {
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;
    
    my $doublet_score = $self->doublet_xia_mathews($t5, $t3, $b3, $b5);

    my $terminal_penalty;
    if ($left_side) {
	$terminal_penalty = $self->terminal_pair_xia_mathews($t5, $b3);
    } else {
	$terminal_penalty = $self->terminal_pair_xia_mathews($t3, $b5);
    }

    return $doublet_score + $terminal_penalty;
}


##========================================================================
##
## SUBROUTINE terminal_pair_xia_mathews() scores a pair of RNA residues based 
##   on how well they would bind, assuming that the pair is found at a
##   terminal end of a helix structure and returns the score.
##
##   NOTE: TERMINAL PAIRS MUST BE AT THE VERY ENDS OF BOTH STRANDS
##
##   EXAMPLES: G would get one score and G would get another...
##             C                         U
##
##   ARGUMENTS:
##   t - top strand
##     
##   b - bottom strand
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_pair_xia_mathews {
    my ($self, $t, $b) = @_;

    my $deltaH = 0;
    my $deltaS = 0;

    if ((($t eq "A") && ($b eq "U")) ||
	(($t eq "U") && ($b eq "A")) ||
	(($t eq "G") && ($b eq "U")) ||
	(($t eq "U") && ($b eq "G"))) {
	$deltaH = 3.72; # Xia (1998) and Mathews (1999)
	$deltaS = 10.5; # Xia (1998) and Mathews (1999) 
    }
    
    return $deltaH - ($gTemperature + 273.15)*($deltaS/1000);

}


##========================================================================
##
## SUBROUTINE terminal_doublet_freier(), parameters values from Freier
##   (1986).
##
##   See "terminal_doublet" for more details.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   left_side - is the terminal end of the helix on the left side,
##     or on the right side?
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_doublet_freier {
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;

    #   NOTE: TERMINAL PAIRS MUST BE AT THE VERY ENDS OF BOTH STRANDS

    if ((defined($left_side)) && $left_side == TRUE) {
	# if we are looking at the left side of the helix, then we will need
	# to flip the bases around so that we will score them correctly.       
	my $temp = $t5;
	$t5 = $b5;
	$b5 = $temp;

	$temp = $t3;
	$t3 = $b3;
	$b3 = $temp;
    }

    # we only need the first pair to match... (and it can't be G/U mismatch)
    if (!($self->wc_pair($t5, $b3))) {
	return BIG_NUM;	
    }

    # now score all possible terminating 
    if (($t5 eq "A") && ($b3 eq "U")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # AA
	    return -0.9; # Freier (1986)      UU
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # AU
	    return -0.9; # Freier (1986)      UA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # AG
	    return -1.7; # Freier (1986)      UC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # AC
	    return -2.1; # Freier (1986)      UG
	}

	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # AG
	    return -0.9; # Freier (1986)      UU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # AU
	    return -0.9; # Freier (1986)      UG
	}

	# misatches..
	if (($t3 eq "A") && ($b5 eq "A")) { # AA
	    return -0.8; # Freier (1986)      UA
	}
	if (($t3 eq "C") && ($b5 eq "C")) { # AC
	    return -0.7; # Freier (1986)      UC
	}
	if (($t3 eq "G") && ($b5 eq "G")) { # AG
	    return -1.0; # Freier (1986)      UG
	}
	if (($t3 eq "U") && ($b5 eq "U")) { # AU
	    return -0.8; # Freier (1986)      UU
	}
	
	if (($t3 eq "A") && ($b5 eq "C")) { # AA
	    return -1.0; # Freier (1986)      UC
	}
	if (($t3 eq "A") && ($b5 eq "G")) { # AA
	    return -1.0; # Freier (1986)      UG
	}
	if (($t3 eq "C") && ($b5 eq "A")) { # AC
	    return -0.7; # Freier (1986)      UA
	}
	if (($t3 eq "G") && ($b5 eq "A")) { # AG
	    return -0.8; # Freier (1986)      UA
	}

	if (($t3 eq "U") && ($b5 eq "C")) { # AU
	    return -0.8; # Freier (1986)      UC
	}
	if (($t3 eq "C") && ($b5 eq "U")) { # AC
	    return -0.7; # Freier (1986)      UU
	}
    }

    if (($t5 eq "U") && ($b3 eq "A")) {
	# watson/crick matches...
	if (($t3 eq "U") && ($b5 eq "A")) { # UU
	    return -0.9; # Freier (1986)      AA
	}
	if (($t3 eq "A") && ($b5 eq "U")) { # UA
	    return -1.1; # Freier (1986)      AU
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # UG
	    return -1.8; # Freier (1986)      AC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # UC
	    return -2.3; # Freier (1986)      AG
	}

	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # UG
	    return -0.9; # Freier (1986)      AU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # UU
	    return -1.0; # Freier (1986)      AG
	}

	# mismatches...
	if (($t3 eq "A") && ($b5 eq "A")) { # UA
	    return -1.0; # Freier (1986)      AA
	}
	if (($t3 eq "C") && ($b5 eq "C")) { # UC
	    return -0.6; # Freier (1986)      AC
	}
	if (($t3 eq "G") && ($b5 eq "G")) { # UG
	    return -1.2; # Freier (1986)      AG
	}
	if (($t3 eq "U") && ($b5 eq "U")) { # UU
	    return -0.5; # Freier (1986)      AU
	}

	if (($t3 eq "A") && ($b5 eq "C")) { # UA
	    return -0.8; # Freier (1986)      AC
	}
	if (($t3 eq "A") && ($b5 eq "G")) { # UA
	    return -1.1; # Freier (1986)      AG
	}
	if (($t3 eq "C") && ($b5 eq "A")) { # UC
	    return -0.7; # Freier (1986)      AA
	}
	if (($t3 eq "G") && ($b5 eq "A")) { # UG
	    return -1.1; # Freier (1986)      AA
	}

	if (($t3 eq "U") && ($b5 eq "C")) { # UU
	    return -0.6; # Freier (1986)      AC
	}
	if (($t3 eq "C") && ($b5 eq "U")) { # UC
	    return -0.5; # Freier (1986)      AU
	}
    }

    if (($t5 eq "C") && ($b3 eq "G")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # CA
	    return -1.8; # Freier (1986)      GU
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # CU
	    return -1.7; # Freier (1986)      GA
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # CC
	    return -2.9; # Freier (1986)      GG
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # CG
	    return -2.0; # Freier (1986)      GC
	}

	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { # CG
	    return -1.6; # Freier (1986)      GU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # CU
	    return -1.9; # Freier (1986)      GG
	}

	# mismatches...
	if (($t3 eq "A") && ($b5 eq "A")) { # CA
	    return -1.9; # Freier (1986)      GA
	}
	if (($t3 eq "C") && ($b5 eq "C")) { # CC
	    return -1.1; # Freier (1986)      GC
	}
	if (($t3 eq "G") && ($b5 eq "G")) { # CG
	    return -1.9; # Freier (1986)      GG
	}
	if (($t3 eq "U") && ($b5 eq "U")) { # CU
	    return -1.2; # Freier (1986)      GU
	}

	if (($t3 eq "A") && ($b5 eq "C")) { # CA
	    return -2.0; # Freier (1986)      GC
	}
	if (($t3 eq "A") && ($b5 eq "G")) { # CA
	    return -1.9; # Freier (1986)      GG
	}
	if (($t3 eq "C") && ($b5 eq "A")) { # CC
	    return -1.0; # Freier (1986)      GA
	}
	if (($t3 eq "G") && ($b5 eq "A")) { # CG
	    return -1.9; # Freier (1986)      GA
	}

	if (($t3 eq "U") && ($b5 eq "C")) { # CU
	    return -1.5; # Freier (1986)      GC
	}
	if (($t3 eq "C") && ($b5 eq "U")) { # CC
	    return -0.8; # Freier (1986)      GU
	}
    }
    
    if (($t5 eq "G") && ($b3 eq "C")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # GA
	    return -2.3; # Freier (1986)      CU
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # GU
	    return -2.1; # Freier (1986)      CA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    return -2.9; # Freier (1986)      CC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    return -3.4; # Freier (1986)      CG
	}
	
	# G/U mismatches...
 	if (($t3 eq "G") && ($b5 eq "U")) { # GG
	    return -1.4; # Freier (1986)      CU
	}
	if (($t3 eq "U") && ($b5 eq "G")) { # GU
	    return -2.3; # Freier (1986)      CG
	}

	# mismatches...              
 	if (($t3 eq "A") && ($b5 eq "A")) { # GA
	    return -1.1; # Freier (1986)      CA
	}
 	if (($t3 eq "C") && ($b5 eq "C")) { # GC
	    return -0.6; # Freier (1986)      CC
	}
 	if (($t3 eq "G") && ($b5 eq "G")) { # GG
	    return -1.4; # Freier (1986)      CG
	}
 	if (($t3 eq "U") && ($b5 eq "U")) { # GU
	    return -0.7; # Freier (1986)      CU
	}

 	if (($t3 eq "A") && ($b5 eq "C")) { # GA
	    return -1.3; # Freier (1986)      CC
	}
 	if (($t3 eq "A") && ($b5 eq "G")) { # GA
	    return -1.3; # Freier (1986)      CG
	}
 	if (($t3 eq "C") && ($b5 eq "A")) { # GC
	    return -1.1; # Freier (1986)      CA
	}
 	if (($t3 eq "G") && ($b5 eq "A")) { # GG
	    return -1.6; # Freier (1986)      CA
	}

 	if (($t3 eq "U") && ($b5 eq "C")) { # GU
	    return -0.8; # Freier (1986)      CC
	}
 	if (($t3 eq "C") && ($b5 eq "U")) { # GC
	    return -0.5; # Freier (1986)      CU
	}	
    }

    print STDERR "ERROR: FreeAlign::terminal_doublet_freier()\n";
    print STDERR "ERROR: doublet $t5$t3\n";
    print STDERR "               $b3$b5 not scored yet...\n";	
    return BIG_NUM;
}


##========================================================================
##
## SUBROUTINE dangling_3prime() scores 3' dangling RNA residues at the ends of
##   helices based on how well they would bind and returns the score.
##
##   EXAMPLES: GC would get one score and AU would get another...
##             C                          U
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, paired with t5.
##
##   RETURN VALUES:
##   - a free energy value (a score) for the unpaired terminal nucleotide
##
##========================================================================
sub dangling_3prime {
    my ($self, $t5, $t3, $b3) = @_;

    if ($gParameterID == FREIER) {
	return $self->dangling_3prime_freier($t5, $t3, $b3);
    } elsif ($gParameterID == SANTALUCIA) {
	return $self->dangling_3prime_santalucia($t5, $t3, $b3);
    } elsif ($gParameterID == XIA_MATHEWS) {
	# CURRENTLY USING THE FREIER PARAMETERS... SHOULD USE PARAMETERS
	# FOUND IN...
	# Serra and Turner, 1995 - Predicting thermodyamic properties of RNA.
	#  Methods in Enzymology 259, 242-261
	return $self->dangling_3prime_freier($t5, $t3, $b3);
    } else {
	print STDERR "ERROR: FreeAlign::dangling_3prime()\n";
	print STDERR "       gParameterID: $gParameterID, undefined\n";
	return BIG_NUM;
    }
}

##========================================================================
##
## SUBROUTINE dangling_3prime_freier() scores 3' dangling RNA residues at the 
##   ends of helices based on how well they would bind and returns the score.
##
##   EXAMPLES: GC would get one score and AU would get another...
##             C                          U
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, paired with t5.
##
##   RETURN VALUES:
##   - a free energy value (a score) for the unpaired terminal nucleotide
##
##========================================================================
sub dangling_3prime_freier {
    my ($self, $t5, $t3, $b3) = @_;

    # we need the pair to match (and no G/U mismatches)...
    if (!($self->wc_pair($t5, $b3))) {
	return BIG_NUM;
    }
    
    # are we allowing terminal bulges?
    if ($gNoBulgePossible) {
	return BIG_NUM;
    }

    if (($t5 eq 'A') && ($b3 eq 'U')) {
	if ($t3 eq 'A') {                # AA
	    return -0.8; # Freier (1986)   U
	}
	if ($t3 eq 'C') {                # AC
	    return -0.5; # Freier (1986)   U
	}
	if ($t3 eq 'G') {                # AG
	    return -0.8; # Freier (1986)   U
	}
	if ($t3 eq 'U') {                # AU
	    return -0.6; # Freier (1986)   U
	}
    }

    if (($t5 eq 'C') && ($b3 eq 'G')) {
	if ($t3 eq 'A') {                # CA
	    return -1.7; # Freier (1986)   G
	}
	if ($t3 eq 'C') {                # CC
	    return -0.8; # Freier (1986)   G
	}
	if ($t3 eq 'G') {                # CG
	    return -1.7; # Freier (1986)   G
	}
	if ($t3 eq 'U') {                # CU
	    return -1.2; # Freier (1986)   G
	}
    }

    if (($t5 eq 'G') && ($b3 eq 'C')) {
	if ($t3 eq 'A') {                # GA
	    return -1.1; # Freier (1986)   C
	}
	if ($t3 eq 'C') {                # GC
	    return -0.4; # Freier (1986)   C
	}
	if ($t3 eq 'G') {                # GG
	    return -1.3; # Freier (1986)   C
	}
	if ($t3 eq 'U') {                # GU
	    return -0.6; # Freier (1986)   C
	}
    }

    if (($t5 eq 'U') && ($b3 eq 'A')) {
	if ($t3 eq 'A') {                # UA
	    return -0.7; # Freier (1986)   A
	}
	if ($t3 eq 'C') {                # UC
	    return -0.1; # Freier (1986)   A
	}
	if ($t3 eq 'G') {                # UG
	    return -0.7; # Freier (1986)   A
	}
	if ($t3 eq 'U') {                # UU
	    return -0.1; # Freier (1986)   A
	}
    }

    print STDERR "ERROR: FreeAlign::dangling_3prime()\n";
    print STDERR "ERROR: 3' dangle $t5$t3\n";
    print STDERR "                 $b3    not scored yet...\n";	    
    return BIG_NUM;
}


##========================================================================
##
## SUBROUTINE dangling_5prime() scores 5' dangling RNA residues at the ends of
##   helices based on how well they would bind and returns the score.
##
##   EXAMPLES: GC would get one score and AU would get another...
##              C                          U
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b5 - bottom strand, paired with t3.
##
##   RETURN VALUES:
##   - a free energy value (a score) for the unpaired terminal nucleotide.
##
##========================================================================
sub dangling_5prime {
    my ($self, $t5, $t3, $b5) = @_;

    # we need the pair to match (and no G/U mismatches)...
    if (!($self->wc_pair($t3, $b5))) {
	return BIG_NUM;
    }    
    
    # are we allowing terminal bulges?
    if ($gNoBulgePossible) {
	return BIG_NUM;
    }
    
    if (($t3 eq 'A') && ($b5 eq 'U')) {
	if ($t5 eq 'A') {                # AA
	    return -0.3; # Freier (1986)    U
	}
	if ($t5 eq 'C') {                # CA
	    return -0.3; # Freier (1986)    U
	}
	if ($t5 eq 'G') {                # GA
	    return -0.4; # Freier (1986)    U
	}
	if ($t5 eq 'U') {                # UA
	    return -0.2; # Freier (1986)    U
	}
    }

    if (($t3 eq 'C') && ($b5 eq 'G')) {
	if ($t5 eq 'A') {                # AC
	    return -0.5; # Freier (1986)    G
	}
	if ($t5 eq 'C') {                # CC
	    return -0.2; # Freier (1986)    G
	}
	if ($t5 eq 'G') {                # GC
	    return -0.2; # Freier (1986)    G
	}
	if ($t5 eq 'U') {                # UC
	    return -0.1; # Freier (1986)    G
	}
    }

    if (($t3 eq 'G') && ($b5 eq 'C')) {
	if ($t5 eq 'A') {                # AG
	    return -0.2; # Freier (1986)    C
	}
	if ($t5 eq 'C') {                # CG
	    return -0.3; # Freier (1986)    C
	}
	if ($t5 eq 'G') {                # GG
	    return -0.0; # Freier (1986)    C
	}
	if ($t5 eq 'U') {                # UG
	    return -0.0; # Freier (1986)    C
	}
    }

    if (($t3 eq 'U') && ($b5 eq 'A')) {
	if ($t5 eq 'A') {                # AU
	    return -0.3; # Freier (1986)    A
	}
	if ($t5 eq 'C') {                # CU
	    return -0.2; # Freier (1986)    A
	}
	if ($t5 eq 'G') {                # GU
	    return -0.2; # Freier (1986)    A
	}
	if ($t5 eq 'U') {                # UU
	    return -0.2; # Freier (1986)    A
	}
    }

    print STDERR "ERROR: FreeAlign::dangling_5prime()\n";
    print STDERR "ERROR: 5' dangle $t5$t3\n";
    print STDERR "                    $b5 not scored yet...\n";	    
    return BIG_NUM;
}


##========================================================================
##
## SUBROUTINE valid_pair(): takes two charactres and returns true if they 
##   are a form a valid helix pair (either standard Watson-Crick pair or
##   a G-U mismatch)
##
##   ARGUMENTS:
##   char_a - a DNA or RNA nucleotide (ACGTU)
##
##   char_b - a DNA or RNA nucleotide (ACGTU)
##
##   RETURN VALUES:
##   - TRUE or FALSE
##
##========================================================================
sub valid_pair {
    my ($self, $char_a, $char_b) = @_;

    if (($char_a eq "G") && (($char_b eq "C") || ($char_b eq "U"))) {
	return TRUE;
    }
    if ((($char_a eq "C") || ($char_a eq "U")) && ($char_b eq "G")) {
	return TRUE;
    }


    if (($char_a eq "A") && (($char_b eq "U") || ($char_b eq "T"))) {
	return TRUE;
    }
    if ((($char_a eq "U") || ($char_a eq "T")) && ($char_b eq "A")) {
	return TRUE;
    }
    
    return FALSE;
}


##========================================================================
##
## SUBROUTINE wc_pair(): takes two characters and returns true if they are
##   a standard Watson-Crick pair.
##
##   ARGUMENTS:
##   char_a - a DNA or RNA nucleotide
##
##   char_b - a DNA or RNA nucleotide
##
##   RETURN VALUES:
##   - TRUE or FALSE
##
##========================================================================
sub wc_pair {
    my ($self, $char_a, $char_b) = @_;

    if (!defined($char_a) || !defined($char_b)) {
	print STDERR "FreeAlign:wc_pair() - missing argument(s)\n";
	return;
    }

    if (($char_a eq "G") && ($char_b eq "C")) {
	return TRUE;
    }
    if (($char_a eq "C") && ($char_b eq "G")) {
	return TRUE;
    }


    if (($char_a eq "A") && (($char_b eq "U") || ($char_b eq "T"))) {
	return TRUE;
    }
    if ((($char_a eq "U") || ($char_a eq "T")) && ($char_b eq "A")) {
	return TRUE;
    }
    
    return FALSE;
}


##========================================================================
##
## SUBROUTINE force_bind() forces an alignment between the two strands, 
##   assuming that they both pair together from their first bases on.  
##   This method ignores both loops and bulges as possibilities.  Only bases 
##   in this forced alignment that can form helices are scored, and only the
##   score of the best helix is returned (that is to say, if there are two
##   or more helices formed, separated by gaps, then only the score of the
##   lowest scoring helix is returned.)
##
##   This is roughly
##   equivelent to assuming that both strands are pulled tight (imagine
##   holding two pieces of rope of equal length together by their ends and
##   pulling the ends apart from each other as far as you can).
##
##   NOTE: The sequences that are paired together are assumed to have the
##   same length.
##
##   ARGUMENTS:
##   seq_x - one of two sequences to be paired together.
##
##   seq_y - one of two sequences to be paired together.
##
##   RETURN VALUES:
##   best_score - The score for the lowest scoring helix in the binding.
##     For example, if the binding is...
##       GGGGCCCCCggAAUU
##       CCCCGGGGGggUUAA
##     Then only the score for the first helix (the one containing Gs and Cs)
##     is returned.
##     NOTE: This score DOES NOT INCLUDE the one time penalty for initiating
##     a helix.  This penalty can be subtracted later using the 
##     gInitPenalty value.
##
##   best_start - the position in the sequences where the best sub-helix
##     begins.
##
##   helix_length - the length of the best sub-helix.
##
##========================================================================
sub force_bind {
    my ($self, $seq_x, $seq_y) = @_;    

    $gBestScore = 0;

    my $x_length = length($seq_x);
    my $y_length = length($seq_y);

    if ($x_length != $y_length) {
	print STDERR "FreeAlign::force_bind() - unequal seqquence lengths.\n";
	return;
    }

    if (($x_length < 2) || ($y_length < 2)) {
	# you need at least two bases in each sequence in order to form
	# a structure between them.
	print STDERR "FreeAlign::force_bind() - sequences are too short.\n";
	return;	
    }


    my $current_score = 0;
    my $best_score = 0;

    my $helix_start_pos = 0;
    my $best_start;

    my $end_of_helix = 0;   

    for (my $i=0; $i<($x_length-1); $i++) {
	my $x_5 = substr($seq_x, $i, 1);
	my $x_3 = substr($seq_x, ($i+1), 1);
	
	my $y_3 = substr($seq_y, $i, 1);
	my $y_5 = substr($seq_y, ($i+1), 1);
	
	if ($current_score == 0) {
	    # start each helix structure with a "terminal doublet"
	    $current_score = $self->terminal_doublet($x_5, $x_3, 
						     $y_3, $y_5, TRUE);
	} else {
	    my $x_55 = "";
	    my $x_33 = "";
	    my $y_33 = "";
	    my $y_55 = "";

	    if ($i < $x_length-2) {
		# What is going on here is that some internal "doublets" need 
		# more context in order to be scored correctly.  See the 
		# scoring for GU in doublet_xia_mathews() for more details.
		#             UG
	        $x_55 = substr($seq_x, ($i-1), 1);
		$x_33 = substr($seq_x, ($i+2), 1);

		$y_33 = substr($seq_y, ($i-1), 1);
		$y_55 = substr($seq_y, ($i+2), 1);
	    }

	    $current_score = $current_score + 
		$self->internal_doublet($x_5, $x_3, 
					$y_3, $y_5, 
					$x_55, $x_33,
					$y_33, $y_55);
	}
	
	if ($current_score > 0) {
	    $helix_start_pos = $i + 1;
	    $current_score = 0;
	} elsif ($current_score < $best_score) {
	    $end_of_helix = $i;
	    $best_score = $current_score;
	    $best_start = $helix_start_pos;
	}
    }

    if ($best_score < 0) {
	# now swap out the last internal doublet score for a terminal doublet
	# score so that the helix starts and ends with a terminal doublet.
	# We have to do this now, since we can not know in advance
	# where this terminal doublet is going to be...

	my $x_5 = substr($seq_x, $end_of_helix, 1);
	my $x_3 = substr($seq_x, ($end_of_helix+1), 1);
	
	my $y_3 = substr($seq_y, $end_of_helix, 1);
	my $y_5 = substr($seq_y, ($end_of_helix+1), 1);
	
	my $internal_score = $self->internal_doublet($x_5, $x_3, $y_3, $y_5);
	
	my $terminal_score = $self->terminal_doublet($x_5, $x_3, 
						     $y_3, $y_5, FALSE);
    
	$best_score = $best_score - $internal_score + $terminal_score;

	# Here is where we should check for "self symmetry", however, given
	# the assumptions made with "force_bind", I'm not entirely sure it
	# is appropriate...  With the force_bind, we are assuming that the
	# strands are stretched out, and thus, not stuck to themselves forming
	# their own hairpin loop.  My understanding is that the "self symmetry"
	# penalty results from the individual strands being stuck together 
	# in their own hairpins.
    }

    $gBestScore = $best_score;

    my $helix_length = 0;
    if (defined($best_start)) {
	$helix_length = $end_of_helix - $best_start + 2;
    }

    return ($best_score, $best_start, $helix_length);
}


##========================================================================
##
## SUBROUTINE test() is very similar to force_bind() with the exception that
##   loops (mismatches) are allowed.
##
##   NOTE: The sequences that are paired together are assumed to have the
##   same length.
##
##   ARGUMENTS:
##   seq_x - one of two sequences to be paired together.
##
##   seq_y - one of two sequences to be paired together.
##
##   RETURN VALUES:
##   score - The score for the optimal pairing between two strands of
##     RNA with the assumption that both strands begin pairing from their 
##     first bases on.  Only bases that form a helix are scored.  Loops are
##     allowed.
##     NOTE: This score DOES NOT INCLUDE the one time penalty for initiating
##     a helix.  This penalty can be subtracted later using the 
##     gInitPenalty value.
##
##========================================================================
sub test {
    my ($self, $seq_x, $seq_y) = @_;    

    $gBestScore = 0;

    my $x_length = length($seq_x);
    my $y_length = length($seq_y);

    if ($x_length != $y_length) {
	print STDERR "FreeAlign::test() - unequal seqquence lengths.\n";
	return;
    }

    if (($x_length < 2) || ($y_length < 2)) {
	# you need at least two bases in each sequence in order to form
	# a structure between them.
	print STDERR "FreeAlign::test() - sequences are too short.\n";
	return;	
    }

    my $scored_left_terminal = FALSE;
    my @helix_scores;
    my $helix_total = 0;
    my $first_helix_start_pos = 0;

    my $current_score = 0;
    my $best_score = 0;

    my $helix_length = 1; # remember, the minimum length for a helix is 2
    my $helix_start_pos = 0;
    my $best_start;

    my $end_of_helix = 0;   
    for (my $i=0; $i<($x_length-1); $i++) {
	my $x_5 = substr($seq_x, $i, 1);
	my $x_3 = substr($seq_x, ($i+1), 1);
	
	my $y_3 = substr($seq_y, $i, 1);
	my $y_5 = substr($seq_y, ($i+1), 1);
	
	if (!$scored_left_terminal) {
	    # start the first helix structure with a "terminal doublet"
	    $current_score = $self->terminal_doublet($x_5, $x_3, 
						     $y_3, $y_5, TRUE);

	    if ($current_score < 0) {
		$scored_left_terminal = TRUE;
		$first_helix_start_pos = $i;
	    }

	} else {
	    my $x_55 = "";
	    my $x_33 = "";
	    my $y_33 = "";
	    my $y_55 = "";

	    if ($i < $x_length-2) {
		# What is going on here is that some internal "doublets" need 
		# more context in order to be scored correctly.  See the 
		# scoring for GU in doublet_xia_mathews() for more details.
		#             UG
	        $x_55 = substr($seq_x, ($i-1), 1);
		$x_33 = substr($seq_x, ($i+2), 1);

		$y_33 = substr($seq_y, ($i-1), 1);
		$y_55 = substr($seq_y, ($i+2), 1);
	    }

	    $current_score = $current_score + 
		$self->internal_doublet($x_5, $x_3, 
					$y_3, $y_5, 
					$x_55, $x_33,
					$y_33, $y_55);
	}

	if ($current_score > 0) {
	    $helix_start_pos = $i + 1;
	    push(@helix_scores, $helix_total);
	    $current_score = 0;
	    $helix_total = 0;
#	} elsif ($current_score < $best_score) {
#	    $helix_length = $helix_length + 1;
#	    $end_of_helix = $i;
#	    $best_score = $current_score;
#	    $best_start = $helix_start_pos;
#	}
	} elsif ($current_score < 0) {
	    $helix_total = $current_score;
	    $end_of_helix = $i;	    
	}

    }

    if ($best_score < 0) {
	# now swap out the last internal doublet score for a terminal doublet
	# score so that the helix starts and ends with a terminal doublet.
	# We have to do this now, since we can not know in advance
	# where this terminal doublet is going to be...

	my $x_5 = substr($seq_x, $end_of_helix, 1);
	my $x_3 = substr($seq_x, ($end_of_helix+1), 1);
	
	my $y_3 = substr($seq_y, $end_of_helix, 1);
	my $y_5 = substr($seq_y, ($end_of_helix+1), 1);
	
	my $internal_score = $self->internal_doublet($x_5, $x_3, $y_3, $y_5);
	
	my $terminal_score = $self->terminal_doublet($x_5, $x_3, 
						     $y_3, $y_5, FALSE);
    
	$best_score = $best_score - $internal_score + $terminal_score;

	# Here is where we should check for "self symmetry", however, given
	# the assumptions made with "force_bind", I'm not entirely sure it
	# is appropriate...  With the force_bind, we are assuming that the
	# strands are stretched out, and thus, not stuck to themselves forming
	# their own hairpin loop.  My understanding is that the "self symmetry"
	# penalty results from the individual strands being stuck together 
	# in their own hairpins.
    }

    $gBestScore = $best_score;

    return ($best_score, $best_start, $helix_length);
}


##========================================================================
##
## SUBROUTINE helix_only() this is similar to force_bind() execpt that
##   the alignment does not have to begin on the first base of both sequences
##   Another way to phrase it is to say that this is the same as force_bind()
##   with the exception that terminal bulges (gaps) are allowed at no penalty.
##   Another significant difference between helix_only() and force_bind()
##   is that all bases are scored using internal_doublet().
##
##   NOTE: helix_only() allows the sequences to have different lengths and
##   it could be possible that the opening gap on one sequence has a different
##   length than the other.  For example, you could get an alignment that looks
##   like this:
##
##   seq1  0: ccccAAAA
##          :     ||||
##   seq2  0:   ccUUUUcc
##
##
##   ARGUMENTS:
##   seq_x - one of two sequences to be paired together (x axis in matrix)
##
##   seq_y - one of two sequences to be paired together (y axix in matrix)
##
##   helix_matrix - a pointer to an array.
##     helix_matrix assumes that the current bases, i, j, (where i is an index
##     into seq_x and j is an index into seq_y) should be paired together in a
##     helix structure.
##
##   RETURN VALUES:
##   best_score - The score for the optimal pairing between two strands of
##     RNA with the assumption that there are no loops or bulges.  This means
##     that it looks for the best single helix structure and returns the score
##     for it alone, without considering other possible helices that might
##     not have such a good score. 
##     NOTE: This score DOES NOT INCLUDE the one time penalty for initiating
##     a helix.  This penalty can be subtracted later using the 
##     gInitPenalty value.
##
##   best_x_pos - the x coordinate for the location of the best score.
##
##   best_y_pos - the y coordinate for the location of the best score.
##
##========================================================================
sub helix_only {
    my ($self, $seq_x, $seq_y, $helix_matrix) = @_;    

    $gBestScore = 0;

    my $x_length = length($seq_x);
    my $y_length = length($seq_y);
    
    if (($x_length < 2) || ($y_length < 2)) {
	# you need at least two bases in each sequence in order to form
	# a structure between them.
	print STDERR "FreeAlign:helix_only() - sequences are too short.\n";
	return;	
    }

    $self->set_seq_lengths($x_length, $y_length);

    my $terminal_bulge_pen = 0;
    my $internal_bulge_pen = BIG_NUM;
    
    my $loop_pen = BIG_NUM;

    # initialize the upper left hand corner of the matrix
    $$helix_matrix[0][0] = 0;

    # initialize the first two rows of the matrix;
    for (my $x = 1; $x <= $x_length; $x++) {
	$$helix_matrix[$x][0] = $terminal_bulge_pen;
	$$helix_matrix[$x][1] = $terminal_bulge_pen;
    }

    # initialize the first two columns of the matrix;
    for (my $y = 1; $y <= $y_length; $y++) {
	$$helix_matrix[0][$y] = $terminal_bulge_pen;
	$$helix_matrix[1][$y] = $terminal_bulge_pen;
    }
    
    my $best_score = 0;
    my $best_x_pos = 0;
    my $best_y_pos = 0;

    for (my $x = 2; $x <= $x_length; $x++) {
	for (my $y = 2; $y <= $y_length; $y++) {
	    # x_seq reads 5'->3'
	    my $x_char5 = substr($seq_x, $x-2, 1);
	    my $x_char3 = substr($seq_x, $x-1, 1);
	    
	    # y_seq reads 3'<-5'
	    my $y_char5 = substr($seq_y, $y-1, 1);
	    my $y_char3 = substr($seq_y, $y-2, 1);

	    my $score = $self->internal_doublet($x_char5, $x_char3, 
						$y_char3, $y_char5);
	    
	    $$helix_matrix[$x][$y] = 
		$self->min($$helix_matrix[$x-1][$y-1] + $score,
			   $$helix_matrix[$x-1][$y] + $internal_bulge_pen,
			   $$helix_matrix[$x][$y-1] + $internal_bulge_pen,
			   0);

	    if ($$helix_matrix[$x][$y] < $best_score) {
		$best_score = $$helix_matrix[$x][$y];
		$best_x_pos = $x;
		$best_y_pos = $y;
	    }
	}
    }

    $gBestScore = $best_score;

    return ($best_score, $best_x_pos, $best_y_pos);
}


##========================================================================
##
## SUBROUTINE helix_only_trace_route() this is similar to trace_route() but 
##   that it is pared down to work the the results generated by helix_only().
##   Basically, it takes a "helix_matrix" and start coordinates within that
##   matrix for the best alignment and figures out what that alighment looks
##   like.
##
##   ARGUMENTS:
##   seq_x - one of two sequences to be paired together (x axis in matrix)
##
##   best_x_pos - the x coordinate for the location of the best score.
##
##   seq_y - one of two sequences to be paired together (y axix in matrix)
##
##   best_y_pos - the y coordinate for the location of the best score.
##
##   helix_matrix - a pointer to an array.
##     helix_matrix assumes that the current bases, i, j, (where i is an index
##     into seq_x and j is an index into seq_y) should be paired together in a
##     helix structure.
##
##   RETURN VALUES:
##   x_route - a string of RNA characters to reflect how
##     $seq_x binds with $seq_y
##   
##   y_route - a string of RNA characters to reflect how
##     $seq_x binds with $seq_y
##
##   pairs - a string that shows how x_route and y_route are paired
##
##========================================================================
sub helix_only_trace_route {
    my ($self, $seq_x, $seq_y, $helix_matrix) = @_;

    my $x = length($seq_x);
    my $y = length($seq_y);

    if (($x == 0) || ($y == 0)) {
	# since we are only performing a local alignment, when we reach the
	# end of either sequence, we can stop trying to extend things...
	return ("", "", "");
    }

    my $x_char = substr($seq_x, $x-1, 1);
    my $y_char = substr($seq_y, $y-1, 1);
    

    if ($$helix_matrix[$x][$y] == 0) {
	return ($x_char, $y_char, HELIX_CHAR);
    } else {
	my ($x_bind, $y_bind, $pairs) =
	    $self->helix_only_trace_route(substr($seq_x, 0, -1),
					  substr($seq_y, 0, -1),
					  $helix_matrix);

	return ($x_bind.$x_char, $y_bind.$y_char, $pairs.HELIX_CHAR);	
    }	
}


##========================================================================
##
## SUBROUTINE fill_matrices() fills matrices for determining the minimum free
##   energy binding configuration between two RNA sequences using a
##   dynamic programming algorithm...
##
##   Currently, the idea is to do this in a way similar to how you do normal 
##   local sequence alignments (Smither-Watterman algorithm) using affine 
##   gap penalties (separate penalties for starting a gap and extending a gap).
##   The big difference is that we'd be treating multiple things (bulges,
##   internal loops and helices) as things that can start and be extended in
##   different ways.
##
##   TERMINOLOGY MAP...
##   bulge = gap
##   internal loop = mismatch
##   helix = match
##
##   NOTE: We are only looking at pairing between the two strands and are 
##   completely ignoring the possibility that there could be pairing within
##   one of the strands.  Thus, we can ignore hairpin loops for now.  This is
##   justified by the fact that we are interested in the lowest free energy 
##   binding between the two strands, and not just the most optimized
##   conformation (which may be two separate hairpin loops).
##
##   ARGUMENTS:
##   seq_x - one of two sequences to be paired together in an optimal way
##     (x-axis)
##
##   seq_y - one of two sequences to be paired together in an optimal way.
##     (y-axis)
##
##   helix_matrix - a pointer to an array.
##     helix_matrix assumes that the current bases, i, j, (where i is an index
##     into seq_x and j is an index into seq_y) should be paired together in a
##     helix structure.
##
##   loop_matrix - a pointer to an array.
##     loop_matrix assumes that the current bases, i, j, should be part of
##     a loop.
##
##   loop_length_matrix - a pointer to an array.
##     this matrix keeps track of the lengths of loops being scored.  This
##     enables us to determine exaclty which charaters in the strings should
##     be compared.
##
##   bulge_x_matrix - a pointer to an array.
##     bulge_x_matrix assumes that the base given by i (an index into seq_x) is
##     part of a buldge.
##
##   bulge_x_length_matrix - a pointer to an array.
##     this matrix keeps track of the lengths of bulges being scored.  This
##     enables us to determine exaclty which charaters in the strings should
##     be compared.
##
##   bulge_y_matrix - a pointer to an array.
##     bulge_x_matrix assumes that the base given by j (an index into seq_y) is
##     part of a buldge.
##
##   bulge_y_length_matrix - a pointer to an array.
##     this matrix keeps track of the lengths of bulges being scored.  This
##     enables us to determine exaclty which charaters in the strings should
##     be compared.
##
##   RETURN VALUES:
##   best_score - The score for the optimal pairing between two strands of
##     RNA.  This score DOES NOT INCLUDE the one time penalty for initiating
##     a helix.  This penalty can be subtracted later using the 
##     gInitPenalty value.
##
##   best_x_pos - the x coordinate for the location of the best score.
##
##   best_y_pos - the y coordinate for the location of the best score.
##
##   best_matrix - the ID for the matrix that contains the best score.
##
##========================================================================
sub fill_matrices {
    my ($self, 
	$seq_x, $seq_y, 
	$helix_matrix, 
	$loop_matrix, $loop_length_matrix,
	$bulge_x_matrix, $bulge_x_length_matrix,
	$bulge_y_matrix, $bulge_y_length_matrix) = @_;    

    $gBestScore = 0;

    my $x_length = length($seq_x);
    my $y_length = length($seq_y);
    
    if (($x_length < 2) || ($y_length < 2)) {
	# you need at least two bases in each sequence in order to form
	# a structure between them.
	print STDERR "FreeAlign:fill_matrices() - sequences are too short.\n";
	return;	
    }

    $self->set_seq_lengths($x_length, $y_length);

    # the upper left hand corner which represents the empty string matching
    # the empty string.
    $$helix_matrix[0][0] = 0;
    $$loop_matrix[0][0] = 0; # terminal "mismatch"
    $$loop_length_matrix[0][0] = 0;
    if (!$gNoBulgePossible) {
	$$bulge_x_matrix[0][0] = 0; # terminal gap
	$$bulge_y_matrix[0][0] = 0; # terminal gap
    } else {
	$$bulge_x_matrix[0][0] = BIG_NUM; # terminal gap
	$$bulge_y_matrix[0][0] = BIG_NUM; # terminal gap
    }
    $$bulge_x_length_matrix[0][0] = 0;
    $$bulge_y_length_matrix[0][0] = 0;

    #
    # initialize the top row of cells in the matrices
    #
    # x_seq reads 5'->3'
    my $x_char5;
    my $x_char3;

    # y_seq reads 3'<-5'
    my $y_char5 = substr($seq_y, 0, 1);
    my $y_char3;

    for (my $x=1; $x<=$x_length; $x++) {
	$$helix_matrix[$x][0] = 0;
	if (!$gNoBulgePossible) {
	    $$loop_matrix[$x][0] = 0; # terminal "mismatch"
	    $$bulge_x_matrix[$x][0] = 0; # terminal gap
	    $$bulge_y_matrix[$x][0] = 0; # terminal gap
	} else {
	    $$loop_matrix[$x][0] = BIG_NUM; # terminal "mismatch"
	    $$bulge_x_matrix[$x][0] = BIG_NUM; # terminal gap
	    $$bulge_y_matrix[$x][0] = BIG_NUM; # terminal gap
	}
	$$loop_length_matrix[$x][0] = 0;
	$$bulge_x_length_matrix[$x][0] = 0;
	$$bulge_y_length_matrix[$x][0] = 0;

	my $dangle_score = 0;
	if ($x > 1) {	    
	    $x_char5 = substr($seq_x, $x-2, 1);
	    $x_char3 = substr($seq_x, $x-1, 1);
	    
	    $dangle_score = $self->dangling_5prime($x_char5, $x_char3, 
						   $y_char5);
	}

	if (!$gNoBulgePossible) {
	    $$helix_matrix[$x][1] = 0;	
	    $$loop_matrix[$x][1] = 0; # terminal "mismatch"
	    $$bulge_x_matrix[$x][1] = $dangle_score; # terminal dangle
	    $$bulge_x_length_matrix[$x][1] = LEFT_DANGLE;
	    $$bulge_y_matrix[$x][1] = 0; # terminal gap
	} else {
	    if ($x == 1) {
		$$helix_matrix[$x][1] = 0;	
		$$loop_matrix[$x][1] = 0; # terminal "mismatch"
	    } else {
		$$helix_matrix[$x][1] = BIG_NUM;	
		$$loop_matrix[$x][1] = BIG_NUM; # terminal "mismatch"
	    }
	    $$bulge_x_matrix[$x][1] = BIG_NUM; # terminal dangle
	    $$bulge_x_length_matrix[$x][1] = 0;
	    $$bulge_y_matrix[$x][1] = BIG_NUM; # terminal gap
	}
	$$loop_length_matrix[$x][1] = 0;
	$$bulge_y_length_matrix[$x][1] = 0;	
    }

    #
    # initialize the left column of cells in the matrices
    #

    $x_char5 = substr($seq_x, 0, 1);

    for (my $y=1; $y<=$y_length; $y++) {
	$$helix_matrix[0][$y] = 0;
	if (!$gNoBulgePossible) {
	    $$loop_matrix[0][$y] = 0; # terminal "mismatch"
	    $$bulge_x_matrix[0][$y] = 0; # terminal gap
	    $$bulge_y_matrix[0][$y] = 0; # terminal gap
	} else {
	    $$loop_matrix[0][$y] = BIG_NUM; # terminal "mismatch"
	    $$bulge_x_matrix[0][$y] = BIG_NUM; # terminal gap
	    $$bulge_y_matrix[0][$y] = BIG_NUM; # terminal gap
	}
	$$loop_length_matrix[0][$y] = 0;
	$$bulge_x_length_matrix[0][$y] = 0;
	$$bulge_y_length_matrix[0][$y] = 0;

	my $dangle_score = 0;
	if ($y > 1) {
	    $y_char5 = substr($seq_y, $y-1, 1);
	    $y_char3 = substr($seq_y, $y-2, 1);
	    
	    $dangle_score = $self->dangling_3prime($y_char5, $y_char3, 
						   $x_char5);
	}

	if (!$gNoBulgePossible) {
	    $$helix_matrix[1][$y] = 0;
	    $$loop_matrix[1][$y] = 0; # terminal "mismatch"
	    $$bulge_x_matrix[1][$y] = 0; # terminal gap
	    $$bulge_y_matrix[1][$y] = $dangle_score; # terminal dangle
	    $$bulge_y_length_matrix[1][$y] = LEFT_DANGLE;
	} else {
	    if ($y == 1) {
		$$helix_matrix[1][$y] = 0;
		$$loop_matrix[1][$y] = 0; # terminal "mismatch"
	    } else {
		$$helix_matrix[1][$y] = BIG_NUM;
		$$loop_matrix[1][$y] = BIG_NUM; # terminal "mismatch"
	    }
	    $$bulge_x_matrix[1][$y] = BIG_NUM; # terminal gap
	    $$bulge_y_matrix[1][$y] = BIG_NUM; # terminal dangle
	    $$bulge_y_length_matrix[1][$y] = 0;
	}
	$$loop_length_matrix[1][$y] = 0;
	$$bulge_x_length_matrix[1][$y] = 0;
    }

    # fill the remaining cells in the matrices
    my $best_score = 0;
    my $best_x_pos = 0;
    my $best_y_pos = 0;
    my $best_matrix = HELIX_ID;

    my $bulge_start = $self->bulge_energy(1); # bulge of length 1
    my $loop_start = $self->internal_loop(MIN_LOOP_STEP);
    
    for (my $x=2; $x<=$x_length; $x++) {
	for (my $y=2; $y<=$y_length; $y++) {

	    if($gNoBulgePossible && ($x != $y)) {
		# What would speed things up would be if we didn't have to
		# fill up the remaining portions of the matrices and could
		# just skip to the next index.  However, to do this we would
		# have to check to see if a position is defined before
		# extending or begining a bulge from it when $x == $y.

		$$helix_matrix[$x][$y] = BIG_NUM;
		$$loop_matrix[$x][$y] = BIG_NUM; # terminal "mismatch"
		$$loop_length_matrix[$x][$y] = 0;
		$$bulge_x_matrix[$x][$y] = BIG_NUM; # terminal gap
		$$bulge_y_matrix[$x][$y] = BIG_NUM; # terminal gap
		$$bulge_x_length_matrix[$x][$y] = 0;
		$$bulge_y_length_matrix[$x][$y] = 0;

		next;
	    }

	    # x_seq reads 5'->3'
	    my $x_char5 = substr($seq_x, $x-2, 1);
	    my $x_char3 = substr($seq_x, $x-1, 1);

	    # y_seq reads 3'<-5'
	    my $y_char5 = substr($seq_y, $y-1, 1);
	    my $y_char3 = substr($seq_y, $y-2, 1);


	    ################################################
	    #
	    # score the three (3) ways the current bases can be in a helix
	    #
	    ################################################
	    
	    # 1) score for extending the existing helix, or returning to
	    # a helix after a loop	    
	    
	    my $score;

	    if (($x == 2) && ($y == 2)) {
		# are we at the start of both of the strands?
		#
		# example:
		#    GAGGGGetc.
		#    CTCCCCetc.
		#
		# We want A & T (the second pair) to match, but it's OK 
		# if the first pair (G & C) do not match.

		$score = $self->terminal_doublet($x_char5, $x_char3, 
						 $y_char3, $y_char5, TRUE);

	    } elsif (($x == $x_length) && ($y == $y_length)) {
		# are we at the end of both of the strands?

		$score = $self->terminal_doublet($x_char5, $x_char3, 
						 $y_char3, $y_char5);
		
	    } else { # we are not at then ends of both strands...

		$score = $self->internal_doublet($x_char5, $x_char3, 
						 $y_char3, $y_char5);
	    }
	    
	    # 2) score for ending a bulge in seq x and returning to a helix
	    my $bulge_x_score = $score;
	    my $bulge_x_length = $$bulge_x_length_matrix[$x-1][$y-1];
	    
	    if ($bulge_x_length == LEFT_DANGLE) {
		$bulge_x_score = $self->internal_doublet($x_char5, $x_char3, 
							 $y_char3, $y_char5);
	    } elsif ($bulge_x_length > 0) {
		my $bulge_x_char5 = substr($seq_x, $x-$bulge_x_length-2, 1);
		
		$bulge_x_score = $self->internal_doublet($bulge_x_char5, 
							 $x_char3, 
							 $y_char3, $y_char5);
	    }

	    # 3) score for ending a bulge in seq y and returning to a helix
	    my $bulge_y_score = $score;
	    my $bulge_y_length = $$bulge_y_length_matrix[$x-1][$y-1];
	    if ($bulge_y_length == LEFT_DANGLE) {
		$bulge_x_score = $self->internal_doublet($x_char5, $x_char3, 
							 $y_char3, $y_char5);
	    } elsif ($bulge_y_length > 0) {
		my $bulge_y_char3 = substr($seq_y, $y-$bulge_y_length-2, 1);
		
		$bulge_y_score = $self->internal_doublet($x_char5, $x_char3, 
							 $bulge_y_char3, 
							 $y_char5);
	    }

	    $$helix_matrix[$x][$y] = 
		$self->min($$helix_matrix[$x-1][$y-1] + $score,
			   $$loop_matrix[$x-2][$y-2] + $score,
			   $$bulge_x_matrix[$x-1][$y-1] + $bulge_x_score,
			   $$bulge_y_matrix[$x-1][$y-1] + $bulge_y_score,
			   0);

	    ################################################
	    #
	    # score the ways to be in a loop (starting or extending)
	    #
	    ################################################

	    my $loop_extend =
		$$loop_matrix[$x-1][$y-1] + 
		($self->internal_loop($$loop_length_matrix[$x-1][$y-1]
				      + MIN_LOOP_STEP)
		 -
		 $self->internal_loop($$loop_length_matrix[$x-1][$y-1]));

	    $$loop_matrix[$x][$y] = 
		$self->min($$helix_matrix[$x-1][$y-1] + $loop_start,
			   $loop_extend,
			   $$bulge_x_matrix[$x-1][$y-1] + $loop_start,
			   $$bulge_y_matrix[$x-1][$y-1] + $loop_start,
			   0);


	    if ($$loop_matrix[$x][$y] < 0) {
		if ($$loop_matrix[$x][$y] == $loop_extend) {
		    $$loop_length_matrix[$x][$y] 
			= $$loop_length_matrix[$x-1][$y-1] + MIN_LOOP_STEP;
		} else {
		    $$loop_length_matrix[$x][$y] = MIN_LOOP_STEP;
		}
	    } else {
		$$loop_length_matrix[$x][$y] = 0;
	    }	    

	    ################################################
	    #
	    # score the ways to be in a bulge (starting or extending) on X
	    #
	    ################################################


	    if (($x <= $x_length) && ($y == $y_length)) {
		# do we have a dangling 3 prime end? 
		# (in other words, do we have a terminal bulge?)
		
		# example
		#   ....AAAUU
		#   ....UUU
		
		my $dangle_score = $self->dangling_3prime($x_char5, $x_char3, 
							  $y_char5);

		$$bulge_x_matrix[$x][$y] = $dangle_score 
		    + $$helix_matrix[$x-1][$y];
		$$bulge_x_length_matrix[$x][$y] = RIGHT_DANGLE;
	    } else {

		# calculate the penalty for extending a bulge on the x strand.
		my $x_bulge_extend = 
		    $$bulge_x_matrix[$x-1][$y] + 
		    ($self->bulge_energy($$bulge_x_length_matrix[$x-1][$y]+1)
		     -
		     $self->bulge_energy($$bulge_x_length_matrix[$x-1][$y]));
		
		$$bulge_x_matrix[$x][$y] =
		    $self->min($$helix_matrix[$x-1][$y] + $bulge_start,
			       $$loop_matrix[$x-1][$y] + $bulge_start,
			       $x_bulge_extend,
			       $$bulge_y_matrix[$x-1][$y] + $bulge_start,
			       0);

		
		# keep track of the length of the bulge as it grows...
		if ($$bulge_x_matrix[$x][$y] < 0) {
		    if ($$bulge_x_matrix[$x][$y] == $x_bulge_extend) {
			$$bulge_x_length_matrix[$x][$y] =
			    $$bulge_x_length_matrix[$x-1][$y] + 1;
		    } else {
			$$bulge_x_length_matrix[$x][$y] = 1;
		    }
		} else {
		    $$bulge_x_length_matrix[$x][$y] = 0;
		}
	    }


	    ################################################
	    #
	    # score the ways to be in a bulge (starting or extending) on Y
	    #
	    ################################################

	    if (($x == $x_length) && ($y <= $y_length)) {
		# do we have a dangling 5' end? 
		# (in other words, do we have a terminal bulge?)
		
		# example
		#   ....AAA
		#   ....UUUUU
		
		my $dangle_score = $self->dangling_5prime($y_char5, $y_char3, 
							            $x_char3);

		$$bulge_y_matrix[$x][$y] = $dangle_score 
		    + $$helix_matrix[$x][$y-1];
		$$bulge_y_length_matrix[$x][$y] = RIGHT_DANGLE;

	    } else {
		my $y_bulge_extend = 
		    $$bulge_y_matrix[$x][$y-1] +
		    ($self->bulge_energy($$bulge_y_length_matrix[$x][$y-1]+1)
		     -
		     $self->bulge_energy($$bulge_y_length_matrix[$x][$y-1]));
				
		$$bulge_y_matrix[$x][$y] =
		    $self->min($$helix_matrix[$x][$y-1] + $bulge_start,
			       $$loop_matrix[$x][$y-1] + $bulge_start,
			       $$bulge_x_matrix[$x][$y-1] + $bulge_start,
			       $y_bulge_extend,
			       0);

		# keep track of the length of the bulge as it grows...
		if ($$bulge_y_matrix[$x][$y] < 0) {
		    if ($$bulge_y_matrix[$x][$y] == $y_bulge_extend) {
			$$bulge_y_length_matrix[$x][$y] =
			    $$bulge_y_length_matrix[$x][$y-1] + 1;
		    } else {
			$$bulge_y_length_matrix[$x][$y] = 1;
		    }
		} else {
		    $$bulge_y_length_matrix[$x][$y] = 0;
		}
	    }
	    
	    ################################################
	    #
	    # keep track of the best way for the bases to bond...
	    #
	    ################################################

	    if ($$helix_matrix[$x][$y] < $best_score) {
		$best_score = $$helix_matrix[$x][$y];
		$best_x_pos = $x;
		$best_y_pos = $y;
		$best_matrix = HELIX_ID;
	    }

	    if ($$bulge_x_matrix[$x][$y] < $best_score) {
		$best_score = $$bulge_x_matrix[$x][$y];
		$best_x_pos = $x;
		$best_y_pos = $y;
		$best_matrix = BULGE_X_ID;
	    }

	    if ($$bulge_y_matrix[$x][$y] < $best_score) {
		$best_score = $$bulge_y_matrix[$x][$y];
		$best_x_pos = $x;
		$best_y_pos = $y;
		$best_matrix = BULGE_Y_ID;
	    }


	}      
    }
    
    $gBestScore = $best_score;

    return ($best_score, $best_x_pos, $best_y_pos, $best_matrix);    
}


##========================================================================
##
## SUBROUTINE trace_route() looks at the last characters in $seq_x and $seq_y,
##   and, given the matrix specified by $current_matrix, looks to see if these
##   characters are starting/continuing a helix, loop or bulge.
##   trace_route then recurses on the ramaining characters in $seq_x and 
##   $seq_y, generating a string that relfects the minimum free energy binding 
##   between $seq_x and $seq_y.
##
##   NOTE: since we are looking at the last characters in $seq_x and $seq_x and
##   recursing on the remaining characters, we are essentially looking at
##   the pairing/alignment from right to left.
##
##   ARGUMENTS:
##   seq_x - a substring of the sequence used to fill the matrices (x axis)
##
##   seq_y - a substring of the sequence used to fill the matrices (y axis)
##
##   helix_matrix - a pointer to an array.  See fill_matrices() for details.
##
##   loop_matrix - a pointer to an array.  See fill_matrices() for details.
##
##   loop_length_matrix - a pointer to an array.  See fill_matrices()...
##
##   bulge_x_matrix - a pointer to an array.  See fill_matrices() for details.
##
##   bulge_x_length_matrix - a pointer to an array.  See fill_matrices()...
##
##   bulge_y_matrix - a pointer to an array.  See fill_matrices() for details.
##
##   bulge_y_length_matrix - a pointer to an array.  See fill_matrices()...
##
##   current_matrix - Used to determine how to look at the matrices...
##     current_matrix should assume values defined by the contants:
##     HELIX_ID, LOOP_ID, BULGE_X_ID or BULGE_Y_ID
##
##   RETURN VALUES:
##   x_route - a string of RNA characters plus others to reflect how
##     $seq_x binds with $seq_y
##   
##   y_route - a string of RNA characters plus others to reflect how
##     $seq_x binds with $seq_y
##
##   pairs - a string that shows how x_route and y_route are paired
##
##========================================================================
sub trace_route {
    my ($self, 
	$seq_x, $seq_y, 
	$helix_matrix, 
	$loop_matrix, $loop_length_matrix,
	$bulge_x_matrix, $bulge_x_length_matrix,
	$bulge_y_matrix, $bulge_y_length_matrix,
	$current_matrix) = @_;

    no warnings qw(recursion);
    
    my $x = length($seq_x);
    my $y = length($seq_y);
        
    if (($x == 0) || ($y == 0)) {
	# since we are only performing a local alignment, when we reach the
	# end of either sequence, we can stop trying to extend things...
	return ("", "", "");
    }

    if ((($x > 1) && ($y == 1)) && ($current_matrix == HELIX_ID)) {    
	# this means we have a dangling 5' thing...
	my $x_char5 = substr($seq_x, $x-2, 1);
	my $x_char3 = substr($seq_x, $x-1, 1);
	
	my $y_char5 = substr($seq_y, $y-1, 1);

	return ($x_char5.$x_char3, " ".$y_char5, " ".HELIX_CHAR);	
    }


    if ((($y > 1) && ($x == 1)) && ($current_matrix == HELIX_ID)) {    
	# this means we have a dangling 3' thing...
	my $y_char5 = substr($seq_y, $y-1, 1);
	my $y_char3 = substr($seq_y, $y-2, 1);
	
	my $x_char3 = substr($seq_x, $x-1, 1);

	return (" ".$x_char3, $y_char3.$y_char5, " ".HELIX_CHAR);	
    }

    if (($current_matrix == HELIX_ID) &&
	(($x == 2) && ($y == 2))) {
	# this means that we have a terminal end beginning to our helix, thus 
	# we allow a mismatch in the first pair if it exists...
	my $x_char5 = substr($seq_x, $x-2, 1);
	my $y_char3 = substr($seq_y, $y-2, 1);
	my $first_valid = $self->valid_pair($x_char5, $y_char3);

	my $x_char3 = substr($seq_x, $x-1, 1);
	my $y_char5 = substr($seq_y, $y-1, 1);

	my $return_string = " ";
	if ($first_valid) {
	    $return_string = HELIX_CHAR.HELIX_CHAR;
	} else {
	    $return_string = " ".HELIX_CHAR;
	}

	return ($x_char5.$x_char3, $y_char3.$y_char5, $return_string);	
    }
    

    if (($current_matrix == HELIX_ID) && ($$helix_matrix[$x][$y] == 0)) {
	# we have traced back our local alignment as far as we need to go
	# now we just need to see if the current characters should be
	# bound together or not...

	my $x_char = substr($seq_x, $x-1, 1);
	my $y_char = substr($seq_y, $y-1, 1);
	
	if ($self->valid_pair($x_char, $y_char)) {
	    return ($x_char, $y_char, HELIX_CHAR);
	} else {
	    # NOTE: this case can happen as a result of some sort of 
	    # dangling end...
	    return ($x_char, $y_char, " ");
	}	
    }

    if ((($current_matrix == LOOP_ID) && ($$loop_matrix[$x][$y] == 0)) ||
	(($current_matrix == BULGE_X_ID) && ($$bulge_x_matrix[$x][$y] == 0)) ||
	(($current_matrix == BULGE_Y_ID) && ($$bulge_y_matrix[$x][$y] == 0))) {
	# we have traced back our local alignment as far as we need to go
	return ("", "", "");
    }
    
    
    ####################################################################
    #
    # extract the last characters from the seq_x and seq_y strings...
    #
    ####################################################################

    # x_seq reads 5'->3'
    my $x_char5;
    my $x_char5_bulge;
    my $x_char3;
    
    $x_char3 = substr($seq_x, $x-1, 1);
    if ($x > 1) {
	$x_char5 = substr($seq_x, $x-2, 1);
    } else {
	$x_char5 = "";
    }
    $x_char5_bulge = $x_char5;

    # y_seq reads 3'<-5'
    my $y_char3;
    my $y_char3_bulge;
    my $y_char5;

    $y_char5 = substr($seq_y, $y-1, 1);
    if ($y > 1) {
	$y_char3 = substr($seq_y, $y-2, 1);
    } else {
	$y_char3 = "";
    }
    $y_char3_bulge = $y_char3;


    ####################################################################
    #
    # score the three (3) ways the current bases can be in a helix
    #
    ####################################################################

    # 1) score for extending the helix...
    my $score = $self->internal_doublet($x_char5, $x_char3, 
			       $y_char3, $y_char5);
    
    # 2) score for ending a bulge in seq x...
    my $bulge_x_score = $score;
    my $end_bulge_x_length = $$bulge_x_length_matrix[$x-1][$y-1];
    if ($end_bulge_x_length > 0) {
	my $seq_x_length = length($seq_x);
	$x_char5_bulge = substr($seq_x, 
				$seq_x_length-$end_bulge_x_length-2, 1);
	
	$bulge_x_score = $self->internal_doublet($x_char5_bulge, $x_char3, 
					$y_char3, $y_char5);
    }       

    # 3) score for ending a bulge in seq y...
    my $bulge_y_score = $score;
    my $end_bulge_y_length = $$bulge_y_length_matrix[$x-1][$y-1];
    if ($end_bulge_y_length > 0) {
	my $seq_y_length = length($seq_y);
	$y_char3_bulge = substr($seq_y, 
				$seq_y_length-$end_bulge_y_length-2, 1);

	$bulge_y_score = $self->internal_doublet($x_char5, $x_char3, 
					$y_char3_bulge, $y_char5);
    }

    ####################################################################
    #
    # calculate the various penalties for loops and bulges...
    #
    ####################################################################

    # bulge_start/loop_start are the penalties for having a bulge or loop 
    # of length 1
    my $bulge_start = $self->bulge_energy(1);
    my $loop_start = $self->internal_loop(MIN_LOOP_STEP);

    # 1) calculate the loop penalty...
    my $loop_extend = 0;
    my $loop_length = 0;
    if (($x != $gSeqXLength) && ($y != $gSeqYLength)) {
	$loop_length = $$loop_length_matrix[$x][$y];
	if ($loop_length > 0) {
	    $loop_extend = $self->internal_loop($loop_length)
		- $self->internal_loop($loop_length - MIN_LOOP_STEP);
	}
    }

    # 2) calculate the bulge penalty (unpaired bases in X)
    my $bulge_x_penalty = 0;
    my $bulge_x_bind_char;
    my $bulge_x_pairs_char;
    my $bulge_x_length = $$bulge_x_length_matrix[$x][$y];
    if ($y != $gSeqYLength) {
	if ($bulge_x_length > 0) {
	    $bulge_x_penalty = $self->bulge_energy($bulge_x_length)
		- $self->bulge_energy($bulge_x_length - 1);
	}
	$bulge_x_bind_char = "-";
	$bulge_x_pairs_char = BULGE_X_CHAR;
    } else {
	$bulge_x_penalty = 0;
	$bulge_x_bind_char = " ";
	$bulge_x_pairs_char = " ";
    }	

    # 2) calculate the bulge penalty (unpaired bases in Y)
    my $bulge_y_penalty = 0;
    my $bulge_y_bind_char;
    my $bulge_y_pairs_char;
    my $bulge_y_length = $$bulge_y_length_matrix[$x][$y];    
    if  ($x != $gSeqXLength) {
	if ($bulge_y_length > 0) {
	    $bulge_y_penalty = $self->bulge_energy($bulge_y_length)
		- $self->bulge_energy($bulge_y_length - 1);
	}
	$bulge_y_bind_char = "-";
	$bulge_y_pairs_char = BULGE_Y_CHAR;
    } else {
	$bulge_y_penalty = 0;
	$bulge_y_bind_char = " ";
	$bulge_y_pairs_char = " ";
    }	    
    
    # x_bind is the string characters that describes how the bases in 
    # seqeucne X are bound to sequence Y.  For example, if there are
    # bulges in sequence Y (extra bases in Y that correspond to gaps in
    # X) then x_bind might look something like this: AC--GGGA where '-'
    # represents the bulge or gap.
    my $x_bind;

    # y_bind is the string characters that describes how the bases in 
    # seqeucne Y are bound to sequence X.
    my $y_bind;

    # pairs is a string that relates x_bind to y_bind, showing complementary
    # pairing with '|', or bulges with "b" or "P" etc.
    my $pairs;

    # next_matrix stores the ID for the matrix that we should look at next
    # (after we finish looking at "current_matrix")
    my $next_matrix;

    # first check to see which matrix we are in and then figure out which 
    # matrix we will go to
    # REMEMBER: since we are looking at the last characters in $seq_x and 
    # $seq_x and recursing on the remaining characters, we are essentially 
    # looking at the pairing/alignment from right to left.
    if (!defined($current_matrix)) {
	print STDERR "WARNING: FreeAlign::trace_route()\n";
	print STDERR "WARNING: current_matrix is not defined!\n";
	print STDERR "x: $x, y: $y\n";
    }

    if ($current_matrix == HELIX_ID) {
	if (($x == $gSeqXLength) && ($y == $gSeqYLength)) {
	    # this means that we have a terminal ending to our helix, thus 
	    # we allow a mismatch in the first pair if it exists...

	    my $second_valid = $self->valid_pair($x_char3, $y_char5);
	    my $return_string = " ";

	    if ($second_valid) {
		$return_string = HELIX_CHAR;
	    } else {
		$return_string = " ";
	    }

	    $next_matrix = HELIX_ID;

	    ($x_bind, $y_bind, $pairs) =
		$self->trace_route(substr($seq_x, 0, -1), 
				   substr($seq_y, 0, -1),
				   $helix_matrix, 
				   $loop_matrix, $loop_length_matrix,
				   $bulge_x_matrix, $bulge_x_length_matrix,
				   $bulge_y_matrix, $bulge_y_length_matrix,
				   $next_matrix);

	    return ($x_bind.$x_char3, 
		    $y_bind.$y_char5, $pairs.$return_string);

	} elsif ($$helix_matrix[$x][$y] == 
		 $$loop_matrix[$x-2][$y-2] + $score) {
	    # transitions from a helix to a loop require us to gobble up 
	    # additional characters from the ends of seq_x and seq_y due to 
	    # the fact that we are scoring doublets or triplets instead of 
	    # each pair in the helix.

	    # here we are making a transition from a helix to a loop
	    # REMEMBER: we are looking from right to left...

	    $next_matrix = LOOP_ID;

	    ($x_bind, $y_bind, $pairs) =
		$self->trace_route(substr($seq_x, 0, -2), 
				   substr($seq_y, 0, -2),
				   $helix_matrix, 
				   $loop_matrix, $loop_length_matrix,
				   $bulge_x_matrix, $bulge_x_length_matrix,
				   $bulge_y_matrix, $bulge_y_length_matrix,
				   $next_matrix);
	    return ($x_bind.$x_char5.$x_char3, 
		    $y_bind.$y_char3.$y_char5, $pairs.HELIX_CHAR.HELIX_CHAR);
	} else {
	    if ($$helix_matrix[$x][$y] == 
		$$helix_matrix[$x-1][$y-1] + $score) {

		# we are just continuing an existing helix...
		$next_matrix = HELIX_ID;
	    } elsif ($$helix_matrix[$x][$y] == 
		     $$bulge_x_matrix[$x-1][$y-1] + $bulge_x_score) {
		# a transition from a helix to a bulge in seq_x
		$next_matrix = BULGE_X_ID;
		$x_char5 = $x_char5_bulge;
	    } elsif ($$helix_matrix[$x][$y] == 
		     $$bulge_y_matrix[$x-1][$y-1] + $bulge_y_score) {
		# a transition from a helix to a bulge in seq_y
		$next_matrix = BULGE_Y_ID;
		$y_char3 = $y_char3_bulge;
	    } else {
		print STDERR "HELIX_ID: oops!\n";
		print STDERR "  x: $x, y: $y\n";
		print STDERR "attempted to score: $x_char5 $x_char3\n";
		if ($y_char3 eq "") {
		    $y_char3 = " ";
		}
		print STDERR "                    $y_char3 $y_char5\n";
		print STDERR "  score: $score\n";
		print STDERR "  bulge_x_score: $bulge_x_score\n";
		print STDERR "  bulge_y_score: $bulge_y_score\n";
		exit;
	    }
	    ($x_bind, $y_bind, $pairs) =
		$self->trace_route(substr($seq_x, 0, -1), 
				   substr($seq_y, 0, -1),
				   $helix_matrix, 
				   $loop_matrix, $loop_length_matrix,
				   $bulge_x_matrix, $bulge_x_length_matrix,
				   $bulge_y_matrix, $bulge_y_length_matrix,
				   $next_matrix);
	    return ($x_bind.$x_char3, $y_bind.$y_char5, $pairs.HELIX_CHAR);
	}
    } elsif ($current_matrix == LOOP_ID) {
	if ($$loop_matrix[$x][$y] ==
	    $$helix_matrix[$x-1][$y-1] + $loop_start) {
	    # we are going from a loop to a helix 
	    # REMEMBER: we are looking from right to left...
	    $next_matrix = HELIX_ID;
	} elsif ($$loop_matrix[$x][$y] == 
		 $$loop_matrix[$x-1][$y-1] + $loop_extend) {
	    # here we are staying in a loop
	    $next_matrix = LOOP_ID;
	} elsif ($$loop_matrix[$x][$y] == 
		 $$bulge_x_matrix[$x-1][$y-1] + $loop_start) {
	    # here we are making a transition from a loop to a bulge in seq_x
	    $next_matrix = BULGE_X_ID;
	} elsif ($$loop_matrix[$x][$y] == 
		 $$bulge_y_matrix[$x-1][$y-1] + $loop_start) {
	    # here we are making a transition from a loop to a bulge in seq_y
	    $next_matrix = BULGE_Y_ID;
	} else {
	    print STDERR "LOOP_ID: there is something funny going on here\n";
	    print STDERR "  x: $x, y: $y\n";
	    print STDERR "  loop_length: $loop_length\n";
	    print STDERR "  loop_start: $loop_start\n";
	    print STDERR "  loop_extend: $loop_extend\n";
	    exit;
	}
	($x_bind, $y_bind, $pairs) =
	    $self->trace_route(substr($seq_x, 0, -1), substr($seq_y, 0, -1),
			       $helix_matrix, 
			       $loop_matrix, $loop_length_matrix,
			       $bulge_x_matrix, $bulge_x_length_matrix,
			       $bulge_y_matrix, $bulge_y_length_matrix,
			       $next_matrix);
	return ($x_bind.$x_char3, $y_bind.$y_char5, $pairs.LOOP_CHAR);
    } elsif ($current_matrix == BULGE_X_ID) {
	if ($$bulge_x_length_matrix[$x][$y] == LEFT_DANGLE) {
	    # this means that the first base in the x strand used in the 
	    # binding is a dangling 5' base (it also means that we are done!)
	    return ($x_char5.$x_char3, " ".$y_char5, " |");
	       
	} elsif ($$bulge_x_length_matrix[$x][$y] == RIGHT_DANGLE) {
	    # this means that the last base in the x strand used in the binding
	    # is a dangling 3' base...
	    $next_matrix = HELIX_ID;
	} elsif ($$bulge_x_matrix[$x][$y] ==
	    $$helix_matrix[$x-1][$y] + $bulge_start) {
	    # we are going from a bulge to a helix 
	    # REMEMBER: we are looking from right to left...
	    $next_matrix = HELIX_ID;
	} elsif ($$bulge_x_matrix[$x][$y] == 
		 $$loop_matrix[$x-1][$y] + $bulge_start) {
	    # here we are going from a bulge to a loop
	    $next_matrix = LOOP_ID;
	} elsif ($$bulge_x_matrix[$x][$y] == 
		 $$bulge_x_matrix[$x-1][$y] + $bulge_x_penalty) {
	    # here we are staying in a bulge in seq_x
	    $next_matrix = BULGE_X_ID;
	} elsif ($$bulge_x_matrix[$x][$y] == 
		 $$bulge_y_matrix[$x-1][$y] + $bulge_start) {
	    # here we are making a transition from a bulge in seq_x 
	    # to a bulge in seq_y
	    $next_matrix = BULGE_Y_ID;
	} else {
	    print STDERR "BULGE_X_ID: oops!\n";
	    print STDERR "  x: $x, y: $y\n";
	    print STDERR "  bulge_x_matrix[$x][$y]: ".$$bulge_x_matrix[$x][$y];
	    print STDERR "\n";
	    print STDERR "  bulge_x_length: $bulge_x_length\n";
	    print STDERR "  bulge_start: $bulge_start\n";
	    print STDERR "  bulge_x_penalty: $bulge_x_penalty\n";
	    exit;
	}		
	($x_bind, $y_bind, $pairs) =
	    $self->trace_route(substr($seq_x, 0, -1), $seq_y,
			       $helix_matrix, 
			       $loop_matrix, $loop_length_matrix,
			       $bulge_x_matrix, $bulge_x_length_matrix,
			       $bulge_y_matrix, $bulge_y_length_matrix,
			       $next_matrix);
	return ($x_bind.$x_char3, $y_bind.$bulge_x_bind_char, 
		$pairs.$bulge_x_pairs_char);	
    } elsif ($current_matrix == BULGE_Y_ID) {
	if ($$bulge_y_length_matrix[$x][$y] == LEFT_DANGLE) {
	    # this means that the first base in the x strand used in the 
	    # binding is a dangling 5' base (it also means that we are done!)
	    return (" ".$x_char3, $y_char3.$y_char5, " |");

	} elsif ($$bulge_y_length_matrix[$x][$y] == RIGHT_DANGLE) {
	    # this means that the last base in the y strand used in the binding
	    # is a dangling 5' base...
	    $next_matrix = HELIX_ID;
	} elsif ($$bulge_y_matrix[$x][$y] ==
	    $$helix_matrix[$x][$y-1] + $bulge_start) {
	    # we are going from a bulge to a helix
	    # REMEMBER: we are looking from right to left...
	    $next_matrix = HELIX_ID;
	} elsif ($$bulge_y_matrix[$x][$y] == 
		 $$loop_matrix[$x][$y-1] + $bulge_start) {
	    # here we are going from a bulge to a loop
	    $next_matrix = LOOP_ID;
	} elsif ($$bulge_y_matrix[$x][$y] == 
		 $$bulge_x_matrix[$x][$y-1] + $bulge_start) {
	    # here we are making a transition from a bulge in seq_y
	    # to a bulge in seq_x
	    $next_matrix = BULGE_X_ID;
	} elsif ($$bulge_y_matrix[$x][$y] == 
		 $$bulge_y_matrix[$x][$y-1] + $bulge_y_penalty) {
	    # here we are staying in a bulge in seq_y
	    $next_matrix = BULGE_Y_ID;
	} else {
	    print STDERR "BULGE_Y_ID: oops!\n";
	    print STDERR "  x: $x, y: $y\n";
	    print STDERR "  bulge_y_length: $bulge_y_length\n";
	    print STDERR "  bulge_start: $bulge_start\n";
	    print STDERR "  bulge_y_penalty: $bulge_y_penalty\n";
	    exit;
	}		
	($x_bind, $y_bind, $pairs) =
	    $self->trace_route($seq_x, substr($seq_y, 0, -1),
			       $helix_matrix, 
			       $loop_matrix, $loop_length_matrix,
			       $bulge_x_matrix, $bulge_x_length_matrix,
			       $bulge_y_matrix, $bulge_y_length_matrix,
			       $next_matrix);
	return ($x_bind.$bulge_y_bind_char, $y_bind.$y_char5, 
		$pairs.$bulge_y_pairs_char);
    }

    return;
}



##========================================================================
##
## SUBROUTINE print_matrix() prints out a matrix to STDERR
##
##   ARGUMENTS:
##   matrix - a reference to a matrix to print
##
##   width - the number of columns in the matrix
##
##   height - the number of rows in the matrix
##
##   row - a reference to an array of strings to label the rows
##
##   column - a reference to an array of strings to label the columns
##
##   floats - a boolean that indicates whether or not the matrix contains
##     floating point numbers or just integers.
##
##========================================================================
sub print_matrix {
    my ($self, $matrix, $width, $height, $row, $column, $floats) = @_;

    no warnings;

    my $print_string = "\%".MATRIX_CELL_WIDTH."d";

    if (defined($row)) {
	print STDERR " ";
	for (my $i=0; $i<$width; $i++) {
	    print STDERR " "x(MATRIX_CELL_WIDTH-1).$$row[$i];
	}
	print STDERR "\n";
    }

    for (my $y=0; $y<$height; $y++) {
	if (defined($column)) {
	    print STDERR $$column[$y];
	}

	for (my $x=0; $x<$width; $x++) {
	    if (defined($floats) && $floats) {
		my $temp = sprintf("\%.1f", $$matrix[$x][$y]);
		$print_string = " "x(MATRIX_CELL_WIDTH - length($temp)).$temp;
	    }
	    printf(STDERR $print_string, $$matrix[$x][$y]);
	}
	print STDERR "\n";
    }
}



##========================================================================
##
## SUBROUTINE print_binding():  pretty prints out the binding
##   between the two RNA strands.
##
##   ARGUMENTS:
##   seq_x - a sequence that has been bound to another sequence
##
##   x_label - the name of seq_x
##
##   seq_y - a sequence that has been bound to another sequence
##
##   y_label - the name of seq_y
##
##   match_string - an optional string that is intended to go between the
##     two sequences.  It is intended to show relationships between the two 
##     sequences (exact matches, close matches...)
##
##========================================================================
sub print_binding {
    my ($self, $seq_x, $x_label, $seq_y, $y_label, $match_string) = @_;

    my $label_length;
    my $x_label_length = length($x_label);
    my $y_label_length = length($y_label);
    if ($x_label_length > $y_label_length) {
	$label_length = $x_label_length;
	# pad $y_label with white space so that it is the same width as
	# $x_label...
	my $pad = $x_label_length - $y_label_length;
	$y_label = $y_label." "x$pad;
    } else {
	$label_length = $y_label_length;
	# pad $x_label with white space to make it the same width as 
	# $y_label...
	my $pad = $y_label_length - $x_label_length;
	$x_label = $x_label." "x$pad;	
    }

    my $alignment_length = length($seq_x);
    my $index_length = length($alignment_length);
    my $index_str = sprintf("%%%dd", $index_length);
    my $index_pad = " "x$index_length;

    my $line_length = LINE_LENGTH - $label_length - $index_length - 2;
    my $match_title = " "x$label_length;
    my $i=0;
    my $x_index = 0;
    my $y_index = 0;
    while ($i<$alignment_length) {
	my $x_sub = substr($seq_x, $i, $line_length);
	my $match_sub;
	if ($match_string) {
	    $match_sub = substr($match_string, $i, $line_length);
	}
	my $y_sub = substr($seq_y, $i, $line_length);	
	printf("%s $index_str: %s\n", $x_label, $x_index, $x_sub);
	if ($match_string) {	    
	    printf("%s %s: %s\n", $match_title, $index_pad, $match_sub);
	}
	printf("%s $index_str: %s\n\n", $y_label, $y_index, $y_sub);
	
	$x_index += ($x_sub =~ tr/a-zA-Z0-9//);
	$y_index += ($y_sub =~ tr/a-zA-Z0-9//);
	$i=$i+$line_length;
    }
}


##========================================================================
##
## SUBROUTINE extract_seq(): extract sequence from a string
##
##   ARGUMENTS:
##   filehandle - a reference to a filehandle to read the sequence from.  
##
##   RETURN VALUES:
##   seq_ref - a reference to a string of sequence data
##   
##   starts_and_stops - a string that specifies the location of the substring.
##     This string is very similar to the format used in Genbank files
##     for specifing CDS and gene locations.
##
##   leader_offset - the number of bases prior to the start, from 
##     starts_and_stops, that you want to extract.
##
##   trailer - the number of bases after the stop, from starts_and_stops, that
##     you want to extract.
##
##   num_bases - the number of bases beyond the start base, or the 
##     leader_offset that you want to extract.
##
##========================================================================
sub extract_seq {
    my ($self, $seq_ref, $starts_and_stops, 
	$leader, $trailer, $num_bases) = @_;

    if (!defined($leader)) {
	$leader = 0;
    }
    if (!defined($trailer)) {
	$trailer = 0;
    }

    my $seq_length = length($$seq_ref);

    my $isRNA = FALSE;
    # check to see if the sequence is RNA...
    if ($$seq_ref =~ /[uU]/) {
	$isRNA = TRUE;
    }

    # strip off any comments...
    $starts_and_stops =~ s/\#(.*)\Z//;
    my $comment = $1;

    # check for the 'complement' keyword...
    my $is_complement = FALSE;
    if ($starts_and_stops =~ /\Acomplement\(/) {
	$is_complement = TRUE;
	# strip off the "complement(" and the closing ")"
	$starts_and_stops =~ s/\Acomplement\(//;
	$starts_and_stops =~ s/\)\s*\Z//;
    }
    
    # strip off join and last ')' if there...
    $starts_and_stops =~ s/\Ajoin\(//;
    $starts_and_stops =~ s/\)\s*\Z//;
    
    my @starts_and_stops = split(',', $starts_and_stops);
    
    my $output_seq = "";
    my $num_starts_and_stops = scalar(@starts_and_stops);	
    for(my $i=0; $i < $num_starts_and_stops; $i++) {
	my $start_and_stop;
	
	if ((defined($num_bases)) && $is_complement) {
	    $start_and_stop = $starts_and_stops[-1];
	} else {
	    $start_and_stop = $starts_and_stops[$i];
	}

	my ($start, $stop) = split(/\.\./, $start_and_stop);
	
	if ($is_complement) {
	    my $temp = $start;
	    $start = $stop;
	    $stop = $temp;
	}
	

	my $seq_data;
	if ($start < $stop) {
	    if (($i == 0) && $leader) {
		$start = $start - $leader;
	    } 
	    if (($i == $num_starts_and_stops-1) && $trailer) {
		$stop = $stop + $trailer;
	    }

	    my $sub_seq_length;
	    if (defined($num_bases)) {
		$sub_seq_length = $num_bases;
		# $seq_data = substr($$seq_ref, $start-1, $num_bases);
	    } else {
		$sub_seq_length = ($stop-$start+1);
		# $seq_data = substr($$seq_ref, $start-1, ($stop-$start+1));
	    }

	    if ($sub_seq_length > $seq_length) {
		print STDERR "ERROR: You are trying to extract a\n";
		print STDERR " subsequence that is longer than the sequence\n";
		exit 0;
	    }

	    $seq_data = substr($$seq_ref, $start-1, $sub_seq_length);
	    my $seq_data_length = length($seq_data);

	    if ((($start <= 0) || ($stop > $seq_length)) 
		&& ($sub_seq_length > $seq_data_length)) {
		# in this case, we need to "wrap around"...
		$sub_seq_length = $sub_seq_length - $seq_data_length;
		my $extra_seq = substr($$seq_ref, 0, $sub_seq_length);
		$seq_data = $seq_data.$extra_seq;
	    }

	    $output_seq = $output_seq.$seq_data;
	    
	} else {
	    if (($i == 0) && $leader) {
		$start = $start + $leader;
	    } 
	    if (($i == $num_starts_and_stops-1) && $trailer) {
		$stop = $stop - $trailer;
	    }

	    my $sub_seq_length;
	    my $start_position;
	    if (defined($num_bases)) {
		$sub_seq_length = $num_bases;
		$start_position = ($start - $num_bases);
		#$seq_data = substr($$seq_ref, $start-$num_bases, $num_bases);
	    } else {
		$sub_seq_length = ($start - $stop + 1);
		$start_position = ($stop - 1);
		#$seq_data = substr($$seq_ref, $stop-1, ($start-$stop+1));
	    }
	    

	    if ($sub_seq_length > $seq_length) {
		print STDERR "ERROR: You are trying to extract a\n";
		print STDERR " subsequence that is longer than the sequence\n";
		exit 0;
	    }

	    $seq_data = substr($$seq_ref, $start_position, $sub_seq_length);
	    my $seq_data_length = length($seq_data);
	    
	    if ((($stop <= 0) || ($start > $seq_length))
		&& ($sub_seq_length > $seq_data_length)) {
		# in this case, we need to "wrap around"...
		$sub_seq_length = $sub_seq_length - $seq_data_length;
		my $extra_seq = substr($$seq_ref, 0, $sub_seq_length);
		$seq_data = $seq_data.$extra_seq;
	    }
	    
	    $self->complement(\$seq_data, $isRNA);
	    $seq_data = reverse($seq_data);
	    $output_seq = $seq_data.$output_seq;
	}
	if (defined($num_bases)) {
	    last;
	}
    }
    return ($output_seq);
}


##========================================================================
##
##
##========================================================================
sub extract_seq_swap {
    my ($self, $seq_ref, $starts_and_stops, 
	$leader, $num_bases, $swap_codon) = @_;

    my $seq_length = length($$seq_ref);

    my $isRNA = FALSE;
    # check to see if the sequence is RNA...
    if ($$seq_ref =~ /[uU]/) {
	$isRNA = TRUE;
    }

    # strip off any comments...
    $starts_and_stops =~ s/\#(.*)\Z//;
    my $comment = $1;

    # check for the 'complement' keyword...
    my $is_complement = FALSE;
    if ($starts_and_stops =~ /\Acomplement\(/) {
	$is_complement = TRUE;
	# strip off the "complement(" and the closing ")"
	$starts_and_stops =~ s/\Acomplement\(//;
	$starts_and_stops =~ s/\)\s*\Z//;
    }
    
    # strip off join and last ')' if there...
    $starts_and_stops =~ s/\Ajoin\(//;
    $starts_and_stops =~ s/\)\s*\Z//;
    
    my @starts_and_stops = split(',', $starts_and_stops);
    
    my $output_seq = "";
    my $num_starts_and_stops = scalar(@starts_and_stops);	
    my $start_and_stop;

    if ((defined($num_bases)) && $is_complement) {
	$start_and_stop = $starts_and_stops[-1];
    } else {
	$start_and_stop = $starts_and_stops[0];
    }
    
    my ($start, $stop) = split(/\.\./, $start_and_stop);
    
    if ($is_complement) {
	my $temp = $start;
	$start = $stop;
	$stop = $temp;
    }
    
    my $seq_data;
    if ($start < $stop) {
	$start = $start - $leader;

	my $sub_seq_length;
	$sub_seq_length = $num_bases;

	if ($sub_seq_length > $seq_length) {
	    print STDERR "ERROR: You are trying to extract a\n";
	    print STDERR " subsequence that is longer than the sequence\n";
	    exit 0;
	}
	
	$seq_data = substr($$seq_ref, $start-1, $sub_seq_length);
	my $seq_data_length = length($seq_data);

	if ((($start <= 0) || ($stop > $seq_length)) 
	    && ($sub_seq_length > $seq_data_length)) {
	    # in this case, we need to "wrap around"...
	    $sub_seq_length = $sub_seq_length - $seq_data_length;
	    my $extra_seq = substr($$seq_ref, 0, $sub_seq_length);
	    $seq_data = $seq_data.$extra_seq;
	}

    } else {
	$start = $start + $leader;

	my $sub_seq_length;
	my $start_position;
	$sub_seq_length = $num_bases;
	$start_position = ($start - $num_bases);
	
	
	if ($sub_seq_length > $seq_length) {
	    print STDERR "ERROR: You are trying to extract a\n";
	    print STDERR " subsequence that is longer than the sequence\n";
	    exit 0;
	}
	
	$seq_data = substr($$seq_ref, $start_position, $sub_seq_length);
	my $seq_data_length = length($seq_data);
	
	if ((($stop <= 0) || ($start > $seq_length))
	    && ($sub_seq_length > $seq_data_length)) {
	    # in this case, we need to "wrap around"...
	    $sub_seq_length = $sub_seq_length - $seq_data_length;
	    my $extra_seq = substr($$seq_ref, 0, $sub_seq_length);
	    $seq_data = $seq_data.$extra_seq;
	}
	
	$self->complement(\$seq_data, $isRNA);
	$seq_data = reverse($seq_data);
    }

    $output_seq 
	= substr($seq_data, 0, $leader).$swap_codon.substr($seq_data, ($leader+3));
    return ($output_seq);
}



##========================================================================
##
## SUBROUTINE read_fasta(): read a file in FASTA format.
##
##   NOTE: this function reads from wherever <$filehandle> last left off and
##     returns the next FASTA sequence.  To read all of the sequences in
##     a FASTA file, you may have to call this function multiple times.
##
##   ARGUMENTS:
##   filehandle - a reference to a filehandle to read the sequence from.  
##
##   RETURN VALUES:
##   seq - a sequence from the file.
##   
##   seq_accession - the accession number in the file.  NOTE: depending on
##     how the line of comments is formatted, this may or may not contain
##     the actual accession number.
##
##   seq_comments - the comments that come before the sequence.
##
##========================================================================
sub read_fasta {
    my ($self, $filehandle) = @_;

    # eliminate all of the filler until the first sequence
    my $input = "";
    my $done = FALSE;
    while ((!$done) && ($input = <$filehandle>)) {
	if ($input =~ /\A\x3e/) {  # this line starts with '>'
	    $done = TRUE;
	}
    }

    if (!$input) {
	# if there is nothing left in the file, we can just return 
	# everything undefined...
	return;
    }

    chomp($input);

    # anything on the same line as the first '>' character is a comment
    my $seq_comments = $input;
    # extract the accession # for the sequence from the first line of comments
    $input =~ /.*?\|.*?\|.*?\|(.*?)\|/;
    my $seq_accession = $1;

    # get the sequence
    my $seq;
    $done = FALSE;
    while((!$done) && ($input = <$filehandle>)) {
	if ($input =~ /\A\W/) { # this line starts with white space
	    $done = TRUE;
	} else {
	    chomp($input);
	    # eliminate any other whitespaces
	    $input =~ s/\s//;
	    $seq .= $input;
	}
    }
    return ($seq, $seq_accession, $seq_comments);
}


##========================================================================
##
## SUBROUTINE load_fasta(): opens and reads the first sequence in a FASTA file
##
##   NOTE: this function only returns the very first sequence in the FASTA
##     file.  If there are more than one sequence in the file and you want
##     to access them, you should pass the filehandle to "read_fasta()".
##
##   ARGUMENTS:
##   input_file - the name of the file to read the sequence from or it is a 
##     filehandle (like STDIN).  
##
##   RETURN VALUES:
##   seq - a sequence from the file.
##   
##   seq_accession - the accession number in the file.  NOTE: depending on
##     how the line of comments is formatted, this may or may not contain
##     the actual accession number.
##
##   seq_comments - the comments that come before the sequence.
##
##   filehandle - the filehandle for the newly opened file, or, if the 
##     input_file was a filehandle to begin with, input_file.
##
##========================================================================
sub load_fasta {
    my ($self, $input_file) = @_;
    
    my $input = "";
    my $done = FALSE;
    my $filehandle;
    if (-f $input_file) {
	open($filehandle, $input_file) || die("can't open $input_file: $!");
    } else {
	$filehandle = $input_file;
    }
    
    my ($seq, $seq_accession, $seq_comments) = $self->read_fasta($filehandle);
    
    return ($seq, $seq_accession, $seq_comments, $filehandle);
}


1;
