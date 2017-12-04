#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fastq_reader;

my $usage = "\n\n\tusage: $0 file.fq.gz\n\n";

my $fq = $ARGV[0] or die $usage;

main: {

    my $fq_reader = new Fastq_reader($fq);

    my $counter = 0;
    
    while (my $fq_record = $fq_reader->next()) {

        my $fq_txt = $fq_record->get_fastq_record();

        my @lines = split(/\n/, $fq_txt);

        $counter++;
        $lines[0] .= ":$counter";

        print join("\n", @lines) . "\n";
    }


    exit(0);
}
