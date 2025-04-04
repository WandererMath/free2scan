#!/usr/bin/perl

##========================================================================
##
## Usage: extract_seq.pl [options] sequence_file.fasta [index_data]
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

use strict;
use warnings;

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

my $gOutputFasta = FALSE;


##========================================================================
##
## DEFAULT HELP MESSAGE
##
##========================================================================
my $sHelp = <<END_OF_HELP;

Usage: extract_seq.pl [options] sequence_file.fasta [index_data]
    
    extract_seq takes a DNA/RNA sequence file (sequence_file), in FASTA
    format, and extracts subsequences from it (prints it to STDOUT) specified
    by index_data.   
    
    If index_data is ommitted, than all of the sequence data is returned.

    index_data can be either a string that specifies start and stop 
    positions within sequence_file.fasta in the same format as used in Genbank
    files, or it can be a filename for a file that contains start and stop
    information for multiple regions.  In the case that index_data is a 
    filename, then there will be one region of sequence extracted per line 
    in the file.

    NOTE: If index_data is a string, make sure you put it in quotation marks
    in order to prevent the shell from interpreting the parentheses.  For
    example:

    shell> ./extract_seq.pl -f sequence_file.fasta  "complement(388..500)"

    The format for the index_data file is very similar to how sequence 
    locations are annotated in Genbank files.  Specifically:

    100..500 # all comments are optional...
    join(100..150,200..250,300..500)  # same as above but without introns
    500..100 # here start > stop, this will cause the extracted sequence
             # to be translated to the complement (a->t, g->c etc.)
             # and then reversed. 
    complement(100..500) # this does the same as the previous indices.
    complement(join(100..150,200..250)) # reverse strand without an intron
    100 # do this only with the '-w' option.
    100,100,100 # do this only with the '-w' option.
    etc.
    
Options:

-h
   print this help message.

-f
   Print the output in FASTA format.

-T trailer_size
   The amount of sequence that should go after the last stop position in
   index_data.  If index_data is a file, this is done for each line in
   that file.

-L leader_offset
   This parameter allows you to define an offset from the first start 
   position in index_data.  If index_data is a file, this is done for 
   each line in that file.

-n num_bases
   This option allows you to specify how many bases from the first start
   position in index_data (or, if index_data is a file, each line in the
   file), or the position defined by "-l",  you wish to extract.
   For example: "-l 10 -n 10", will extract the 10 elements prior to the
   first start positions defined by index_data.

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
    getopts('hw:L:T:fn:s:', \%opts);

    if ($opts{'h'}) { # print help and exit
        print STDERR $sHelp;
        exit 1;
    }

    my $swap_codon;
    if ($opts{'s'}) { # THIS IS AN UNDOCUMENTED FEATURE...
	$swap_codon = $opts{'s'};
    }

    if ($opts{'w'}) { # we are no longer supporting this...
	print STDERR "The window feature is not supported, use -L and -n\n";
	exit 1;
    }

    if ($opts{'f'}) {
	$gOutputFasta = TRUE;
    }

    my $leader_offset = 0;
    if ($opts{'L'}) {
	$leader_offset = $opts{'L'};
    }
    
    my $num_bases;
    if ($opts{'n'}) {
	$num_bases = $opts{'n'};
    }

    my $trailer = 0;
    if ($opts{'T'}) {
	$trailer = $opts{'T'};
    }
    
    # get the sequence file name
    my $seq_file = shift(@ARGV);
    if (!defined($seq_file)) {
	print STDERR "\nERROR: you must specify a sequence file\n\n";
        print STDERR $sHelp;
	exit;
    }
    
    my $align_obj = new FreeAlign();
    
    my ($seq) = $align_obj->load_fasta($seq_file);
    
    # get the index_data
    my $index_data = shift(@ARGV);
    my @index_data_array;
    if (!defined($index_data)) {
	if ($gOutputFasta) {
	    print "> SEQUENCE FILE: $seq_file\n";
	}
	print $seq."\n";
	exit; # all done!
    } elsif (!(-f $index_data)) {
	push(@index_data_array, $index_data);
    } else {
	open(CMD_FILE, $index_data)
	    || die("can't open $index_data: $!");

	@index_data_array = <CMD_FILE>
    }

    foreach my $starts_and_stops (@index_data_array) {
	chomp($starts_and_stops);	

	if ($starts_and_stops eq "") {
	    next;
	}

	# strip off any comments...
	$starts_and_stops =~ s/\#(.*)\Z//;
	my $comment = $1;

	if ($gOutputFasta) {
	    if (-f $index_data) {
		print "> INDEX DATA: $index_data, $starts_and_stops";
	    } else {
		print "> INDEX DATA: $starts_and_stops";
	    }
	    if (defined($comment) && ($comment ne "")) {
		$comment =~ /\A\s*(.*?)\s*\Z/;
		print " \# $comment";
	    }
	    print ", SEQUENCE FILE: $seq_file";

	    if ($leader_offset || $trailer) {
		print ", leader: $leader_offset, trailer, $trailer";
	    }

	    if (defined($num_bases)) {
		print ", num_bases: ".$num_bases;
	    }

	    print "\n";
	}

	my $output_seq = "";
	if (defined($swap_codon)) {
	    $output_seq = $align_obj->extract_seq_swap(\$seq, 
						       $starts_and_stops, 
						       $leader_offset,
						       $num_bases,
						       $swap_codon);
	} else {
	    $output_seq = $align_obj->extract_seq(\$seq, $starts_and_stops, 
						  $leader_offset, $trailer, 
						  $num_bases);
	}

	print $output_seq."\n";
    }
}
