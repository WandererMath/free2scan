#!/usr/bin/perl

##========================================================================
##
## Usage: extract_indices.pl [options] genbank_file > index_data
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

use constant TRUE => 1;
use constant FALSE => 0;

use constant CDS => 1;
use constant GENE => 2;
use constant MRNA => 3;


##========================================================================
#
# DEFINE GLOBAL VARIABLES
#
##========================================================================


my $gRelaxed = FALSE;
my $gDebug = FALSE;
my $gExtractFeature = CDS;


##========================================================================
##
## DEFAULT HELP MESSAGE
##
##========================================================================
my $sHelp = <<END_OF_HELP;

 Usage: extract_indices.pl [options] genbank_file > index_data

    extract_indices takes a genbank file and prints index data 
    (start/stop locations, in genbank style format) for 
    the each gene that has a "/gene=" tag, ingnoring those that also contain
    one or more of the following tags:

    /transposon=
    /pseudo
    /note=\"Derived by\"
    
    You can choose to extract the indices for the entire gene, its first 
    "CDS" or first "mRNA" (see options below). 
    

    options:

    -h
      print this help message

    -c
      extract the CDS indices (THIS IS THE DEFAULT)

    -g
      extract the indices for the entire gene.

    -m
      extract the mRNA indices.

    -r
      This puts the program into "relaxed" mode, where it is not so picky about
      whether the gene has a "/gene=" tag or a "/locus_tag" tag.

    -d
      Turns on debug mode.  Lots of stuff gets printed out.  Brace yourself.

    -n
      Refrain from printing the name of the gene as a comment, i.e.
      " # gene_name " at the end of each line.  
 
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
    getopts('hcgmndr', \%opts);
    
    if ($opts{'h'}) { # print help and exit
        print STDERR $sHelp;
        exit 1;
    }

    if ($opts{'c'}) {
	$gExtractFeature = CDS;
    }

    if ($opts{'g'}) {
	$gExtractFeature = GENE;
    }

    if ($opts{'m'}) {
	$gExtractFeature = MRNA;
    }
    
    if ($opts{'r'}) {
	$gRelaxed = TRUE;
    }

    if ($opts{'d'}) {
	$gDebug = TRUE;
	
	print STDERR "\nDEBUG MODE IS ON\n";
	print STDERR "gExtractFeature is set to: $gExtractFeature\n";
    }

    my $gbk_file = $ARGV[0];
    open(GBK_FILE, $gbk_file) 
	|| die("can't open $gbk_file: $!");

    my $line;
    my $gene_name = "";
    my %gene_names;
    my $locus_tag = "";
    my $output_line = "";
    while ((defined($line = <GBK_FILE>)) && (!($line =~ /\AORIGIN/))) {
	chomp($line);

	# extract the name of the gene...
	if ($line =~ /\A     gene/) {	
	    $output_line = "";
	    $gene_name = "";

	    if ($gExtractFeature == GENE) {
		$output_line = substr($line, 21);
	    }

	    my $done = FALSE;
	    # loop until we find another gene, snRNA, mRNA or some other thing
	    while ((!$done) 
		   && (defined($line = <GBK_FILE>)) 
		   && ($line !~ /\A     \w/)) {

		chomp($line);
		if ($line =~ /\A                     \/gene=/) {
		    $gene_name = substr($line, 27);
		    # strip off qutation marks
		    $gene_name =~ s/\"//g;
		    if ($gDebug) {
			print STDERR $gene_name." -- examining...\n";
		    }
		} elsif ($gRelaxed && 
			 ($line =~ /\A                     \/locus\_tag=/)) {
		    $locus_tag = substr($line, 32);
		    # strip off qutation marks
		    $locus_tag =~ s/\"//g;
		    if ($gDebug) {
			print STDERR $locus_tag." -- examining...\n";
		    }

		} 
		
		if (($line =~ /\A                     \/transposon=/) ||
		    ($line =~ /\A                     \/pseudo/) || 
		    ($line =~ /\A                     \/note=\"Derived by/)) {
		    
		    if ($gDebug) {
			print STDERR $gene_name." --- not verified as a true gene\n";
		    }
		    
		    $gene_name = "";
		    $done = TRUE;
		}
	    }
	}
	
	if (!($gene_name eq "") || !($locus_tag eq "")) { 
	    
	    if ($gDebug) {
		print STDERR "looking for ".$gExtractFeature." in : ".$line;
	    }

	    if ((($gExtractFeature == CDS) && ($line =~ /\A     CDS/)) ||
		(($gExtractFeature == MRNA) && ($line =~ /\A     mRNA/))) {
		
		if ($gDebug) {
		    print STDERR "found what we were looking for...\n";
		}

		# strip off "     CDS             " or
		# strip off "     mRNA            " or		
		$line = substr($line, 21);

		# if the annotation is vague (contains '<' or '>' characters)
		# then we will just skip the gene...
		if ($line =~ /[<>]/) {
		    next;
		}

		my $done_with_indices = FALSE;
		do {
		    chomp($line);
		    $output_line = $output_line.$line;
		    
		    $line = <GBK_FILE>;
		    chomp($line);
		    $line = substr($line, 21);
		    if (!($line =~ /\A[0-9]/)) {
			# if $line does not begin with a number, 
			# then we are done.
			$done_with_indices = TRUE;	    
		    }
		    
		} while (!$done_with_indices);
	    }
	    

	    if (!($output_line eq "")) {
		if ($opts{'n'}) {
		    print $output_line."\n";
		} elsif ($gene_name eq "") {
		    print $output_line." \# ".$locus_tag."\n";		   
		} else {
		    if ($gene_names{$gene_name}) {
			my $unique_name = FALSE;
			my $counter = 2;
			while(!$unique_name) {
			    if (!$gene_names{$gene_name."-".$counter}) {
				$gene_name = $gene_name."-".$counter;
				$unique_name = TRUE;
			    } else {
				$counter++;
			    }
			}
		    }
		    $gene_names{$gene_name} = TRUE;
		    print $output_line." \# ".$gene_name."\n";
		}

		$output_line = "";
		$gene_name = ""; # just look at one version of the gene for now
		$locus_tag = "";
	    }
	}
    }
}
