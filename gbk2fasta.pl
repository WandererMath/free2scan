#!/usr/bin/perl

##========================================================================
##
## Usage: gbk2fasta.pl genbank_file
##
##    gbk2fasta takes a genbank file and attempts to create a FASTA file
##    from it.
##
##    options: there are no options as of yet.
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

Usage: gbk2fasta.pl genbank_file

    gbk2fasta takes a genbank file and attempts to create a FASTA file
    from it.

    options:

    -h
       print this help message.

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
    getopts('h', \%opts);

    if ($opts{'h'}) { # print help and exit
        print STDERR $sHelp;
        exit 1;
    }

    my $filename = $ARGV[0];

    open(GBK, $filename) || die ("can't open $filename: $!");

    # look for the sequence definition data.
    my $done = FALSE;
    my $definition;
    my $accession;
    my $input = <GBK>;
    while (!$done && $input) {
	if ($input =~ /\ADEFINITION/) { # save the definition
	    chomp($input);
	    $input = substr($input, 12);
	    $definition = $input;
	} elsif ($input =~ /\AACCESSION/) { # save the accession number
	    chomp($input);
	    $input = substr($input, 12);
	    $accession = $input;
	} elsif ($input =~ /\AORIGIN/) { # ORIGIN occurs right before sequence
	    $done = TRUE;
	} else {
	    $input = <GBK>;
	}
    }
    
    # print out the sequence accession number and definition data as the title 
    # for the sequence
    if (!defined($accession)) {
	$accession = "sequence data from ".$filename;
    }
    if (!defined($definition)) {
	$definition = "";
    }
    print ">gi|".$accession."|".$definition."\n";
    
    $input = <GBK>;
    while (($input) && !($input =~ /\A\/\//)) {
	$input = substr($input, 10);
	$input =~ s/\W//g; # remove all white-spaces
	print $input."\n";;
	
	$input = <GBK>;
    }
}
