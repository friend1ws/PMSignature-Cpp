#! /usr/local/bin/perl

use strict;


my $input = $ARGV[0];
my $num = $ARGV[1];

open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;

    my @F = split("\t", $_);

    my $seq = substr($F[2], 5 - $num, 2 * $num + 1);
    $F[1] =~ s/ /_/g;
    
    next if ($seq =~ /N/);
    print $F[0] . "\t" . $F[1] . "\t" . $seq . "\t" . $F[3] . "\n";
}
close(IN);



