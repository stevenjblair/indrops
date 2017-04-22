#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $help_flag;

my $min_reads = 750;

my $usage = <<__EOUSAGE__;

#######################################
#
#  --bam <string>        bam file
#
#  optional:
#
#  --min_reads <int>     default: $min_reads
#
#######################################

__EOUSAGE__

    ;


my $bam_file;


&GetOptions ( 'h' => \$help_flag,
              'bam=s' => \$bam_file,
              'min_reads=i' => \$min_reads,
    );

unless ($bam_file) {
    die $usage;
}

main: {
     
    
    my $fh;
    if ($bam_file =~ /\.sam$/) {
        open($fh, $bam_file) or die "Error, cannot open file $bam_file";
    }
    elsif ($bam_file =~ /\.bam$/) {
        open($fh, "samtools view $bam_file |") or die "Error, cannot samtools view $bam_file";
    }

    my %read_to_target;

    {
        print STDERR "-parsing bam file\n";
        my $line_counter = 0;
        while (<$fh>) {
            chomp;
            if (/^\@/) {
                # header line
                next;
            }
            $line_counter++;
            if ($line_counter % 100000 == 0) {
                print STDERR "\r[$line_counter]      ";
            }
            
            my @x = split(/\t/);
            my $read = $x[0];
            my $hit_acc = $x[2];
            if ($hit_acc eq '*') { next; }
            
            $read_to_target{$read}->{$hit_acc} = 1;
        }
        close $fh;
    }

    
    print STDERR "\n\nBuilding count matrix\n";

    my %gene_sample_counter;
    my %cells;
    
    {
        foreach my $read (keys %read_to_target) {
            my @hits = keys %{$read_to_target{$read}};

            my ($sample, $barcode, @rest) = split(/-/, $read);
            my $sample_barcode = "${sample}_${barcode}";
            foreach my $hit (@hits) {
                $gene_sample_counter{$hit}->{$sample_barcode}++;
            }
            $cells{$sample_barcode} += 1;
        }
    }

    # remove cells w/ fewer than min_reads
    {
        my @remove_cells;

        my $total_cells = 0;
        foreach my $cell (keys %cells) {
            $total_cells++;
            my $num_reads = $cells{$cell};
            if ($num_reads < $min_reads) {
                push (@remove_cells, $cell);
            }
        }
        my $num_removed = 0;
        foreach my $cell (@remove_cells) {
            delete $cells{$cell};
            $num_removed++;
        }

        print STDERR "-removed $num_removed / $total_cells as below $min_reads read threshold\n";
    }
    

    print STDERR "\n\nPrinting count matrix\n";
    # output matrix
    {
        my @cells = sort keys %cells;

        print "\t" . join("\t", @cells) . "\n";
        foreach my $gene (sort keys %gene_sample_counter) {
            my @vals = ($gene);
            
            foreach my $cell (@cells) {
                my $count = $gene_sample_counter{$gene}->{$cell} || 0;
                push (@vals, $count);
            }

            print join("\t", @vals) . "\n";
        }
    }


    exit(0);
}

