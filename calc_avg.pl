#!/usr/bin/perl

##========================================================================
##
## Usage: calc_avg.pl [options] input_file
##
##    Calculates the average value for each "position" in the input file.  The
##    format for the input file is (and the comments are optional)
##    
##    position_1   # comment
##    position_2
##    position_3
##    .
##    .
##    position_m
##
##    position_1   # comment
##    position_2
##    .
##    .
##    position_n
##
##    position_1   # comment
##    etc.
##
##    options:
##
##    -h
##       print this help message.
##
##    -n
##       print out the number of elements used to calculate the average
##       in the next column.
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


##========================================================================
##
## DEFINE GLOBAL CONSTANTS
##
##========================================================================

use constant TRUE => 1;
use constant FALSE => 0;


##========================================================================
##
## DEFAULT HELP MESSAGE
##
##========================================================================
my $sHelp = <<END_OF_HELP;

Usage: calc_avg.pl [options] input_file

    Calculates the average value for each "position" in the input file.  The
    format for the input file is (and the comments are optional)
    
    position_1   # comment
    position_2
    position_3
    .
    .
    position_m

    position_1   # comment
    position_2
    .
    .
    position_n

    position_1   # comment
    etc.

    options:

    -h
       print this help message.

    -n
       print out the number of elements used to calculate the average
       in the next column.

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
    getopts('n', \%opts);

    my $print_count = FALSE;
    if ($opts{'n'}) {
	$print_count = TRUE; # print out the number of elements
    }


    my $index = 0;
    my @total_vec;
    my @count_vec;
    
    foreach my $line (<>) {
	
	# look to see if the line is nothing but a new line, if it is, reset
	# the $index etc...
	if ($line eq "\n") {
	    $index = 0;
	    next;
	}
	
	# get the number...
	$line =~ /\A(.*?)\s/;
	my $num = $1;
	
	if (!defined($num)) {
	    print STDERR "Error!, num was not defined!\n";
	    print STDERR $line."\n";
	    exit 0;
	}
	
	if (defined($total_vec[$index])) {
	    $total_vec[$index] = $total_vec[$index] + $num;
	    $count_vec[$index] = $count_vec[$index] + 1;
	} else {
	    $total_vec[$index] = $num;
	    $count_vec[$index] = 1;
	}
	
	$index++
	}
    
    for(my $i=0; $i < scalar(@total_vec); $i++) {
	my $avg = $total_vec[$i] / $count_vec[$i];
	
	print $avg;

	if ($print_count) {
	    print "\t".$count_vec[$i];
	}
	print "\n";
    }
}
