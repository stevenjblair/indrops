#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $help_flag;

my $min_reads = 750;

my $max_top_cells = 5000;
my $min_cells_per_gene = 3;
my $max_read_mappings = 10;

my $usage = <<__EOUSAGE__;

#######################################
#
#  --bam <string>        bam file
#
#  optional:
#
#  --min_reads <int>     default: $min_reads
#
#  --min_cells_per_gene <int>    default: $min_cells_per_gene
#
#  --max_top_cells <int>   restrict to top read count cells (0=no restriction) (default: $max_top_cells)
#
#  --max_read_mappings <int>      default: $max_read_mappings
#
#######################################




__EOUSAGE__

    ;


my $bam_file;


&GetOptions ( 'h' => \$help_flag,
              'bam=s' => \$bam_file,
              'min_reads=i' => \$min_reads,
              'min_cells_per_gene=i' => \$min_cells_per_gene,
              'max_top_cells=i' => \$max_top_cells,
              'max_read_mappings=i' => \$max_read_mappings,
    );

if ($help_flag) {
    die $usage;
}

unless ($bam_file) {
    die $usage;
}

main: {
     
    
    my $fh;
    if ($bam_file =~ /\.sam$/) {
        open($fh, $bam_file) or die "Error, cannot open file $bam_file";
    }
    elsif ($bam_file =~ /\.sam\.gz$/) {
        open($fh, "gunzip -c $bam_file | ") or die "Error, couldn't open gunzip -c $bam_file";
    }
    elsif ($bam_file =~ /\.bam$/) {
        open($fh, "samtools view $bam_file |") or die "Error, cannot samtools view $bam_file";
    }
    else {
        die "Error, dont recognize file type: $bam_file";
    }
    
    my %read_to_target;  # read -> gene
    my %gene_to_cell;    # gene -> cell
    
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
            my $read = $x[0];    # ex.   S1-bc1816-CAGTTTGC-ACACTAAG:GGATTTT:0
            my $hit_acc = $x[2]; # ex.   c1072462_g3_i5|c1072462_g3_i5
            if ($hit_acc eq '*') { next; }

            
            # capture read to target gene mapping
            $read_to_target{$read}->{$hit_acc} = 1;
        

            $read =~ /^(\w+-\w+)-/ or die "Error, cannot extract cell info from $read";
            my $cell = $1;  # ex. S1-bc1816
            $gene_to_cell{$hit_acc}->{$cell} = 1;
            
            
        }
        close $fh;
    }

    
    print STDERR "\n\nBuilding count matrix\n";

    # note, reads encode cell and UMI info.
    
    my %gene_sample_counter;
    my %cells;
    
    {

        my %hit_to_umi; # avoid multi-counting UMIs / gene.
        
        foreach my $read (keys %read_to_target) {
            my @hits = keys %{$read_to_target{$read}};

            # exclude those reads that map to too many target genes
            if (scalar @hits > $max_read_mappings) { next; }
            
            my $num_hits = scalar(@hits);
            
            my ($sample, $barcode, @rest) = split(/-/, $read);  # ex. S1-bc1816-CAGTTTGC-ACACTAAG:GGATTTT:0
            $sample =~ s/_//; # seurat uses first _ for sample to cell delineation
            my $sample_barcode = "${sample}_${barcode}";
            
            # get read UMI
            my @pts = split(/:/, $read);
            my $umi = $pts[1];
            
            foreach my $hit (@hits) {

                # read umi is only counted once / cell
                if (! $hit_to_umi{$hit}->{$umi}) {
                    
                    # split read among the multiple hits.
                    $gene_sample_counter{$hit}->{$sample_barcode} += 1/$num_hits;
                    $hit_to_umi{$hit}->{$umi} = 1;
                }
            }

            
            # increment read/cell counter
            $cells{$sample_barcode} += 1;
        }
    }
    
    
    # remove cells w/ fewer than min_reads
    {

        my @all_cells = reverse sort { $cells{$a} <=> $cells{$b} } keys %cells;
        my $total_cells = scalar @all_cells;
        my @remove_cells;
        if ($max_top_cells > 0 && scalar(@all_cells) > $max_top_cells) {
            @remove_cells = @all_cells[$max_top_cells..$#all_cells];
            @all_cells = @all_cells[0..($max_top_cells-1)];
            print STDERR "-capturing just the top $max_top_cells / $total_cells total cells\n";
        }

        my $num_removed = 0;
        foreach my $cell (@all_cells) {
            my $num_reads = $cells{$cell};
            if ($num_reads < $min_reads) {
                push (@remove_cells, $cell);
                $num_removed++;
            }
        }
        print STDERR "-removed an additional $num_removed cells due to not meeting $min_reads min reads cutoff.\n" if $num_removed;
        
        foreach my $cell (@remove_cells) {
            delete $cells{$cell};
        }
        
        
    }
    
    
    print STDERR "\n\nPrinting count matrix\n";
    # output matrix
    {
        my @cells = sort keys %cells;

        print "\t" . join("\t", @cells) . "\n";
        foreach my $gene (sort keys %gene_sample_counter) {
            my @vals = ($gene);

            my $num_cells_expr = 0;
            foreach my $cell (@cells) {
                my $count = $gene_sample_counter{$gene}->{$cell} || 0;
                if ($count > 0) {
                    $count = sprintf("%.1f", $count);
                }
                push (@vals, $count);
                if ($count) { $num_cells_expr++; }
            }

            print join("\t", @vals) . "\n" if ($num_cells_expr >= $min_cells_per_gene);
        }
    }
    

    exit(0);
}

