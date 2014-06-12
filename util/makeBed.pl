#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $size = $ARGV[1];

open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;

    my @F = split("\t", $_);

    print $F[1] . "\t" . ($F[2] - $size - 1) . "\t" . ($F[3] + $size) . "\t" . join(";", @F) . "\n";

}
close(IN);

