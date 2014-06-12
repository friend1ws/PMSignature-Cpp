#! /usr/local/bi/nperl

use strict;

my $input = $ARGV[0];
my $exon = $ARGV[1];
my $type = $ARGV[2];
my $size = $ARGV[3];

my %key2strand = ();
open(IN, $exon) || die "cannot open $!";
while(<IN>) {
    s/[\r\n]//g;
    my @F = split("\t", $_);
    if (not exists $key2strand{$F[3]}) {
        $key2strand{$F[3]} = $F[9];
    } elsif ($key2strand{$F[3]} ne $F[9]) {
        $key2strand{$F[3]} = "0";
    }
}
close(IN);


my $num = 0;
open(IN, $input) || die "cannot open $!";
while(<IN>) {

    s/[\r\n]//g;

    my @F = split("\t", $_);
    my @FF = split(";", $F[0]);

    # my $sample = join(";", @FF[0 .. ($#FF - 2)]);
    my $sample = $FF[0];
    my $chr = $FF[1];
    my $pos = $FF[2];

    my $ref = $FF[$#FF - 1];
    my $alt = $FF[$#FF];

    my $strand = exists $key2strand{$F[0]} ? $key2strand{$F[0]} : "0";

    if ($ref ne substr($F[1], $size, length($ref))) {
        die $num . " th inconsistent\n";
    }

    print $type . "_" . $num . "\t" . $sample . "\t" . $F[1] . "\t" . $alt . "\t" . $chr . "\t" . $pos . "\t" . $strand . "\n";

    $num = $num + 1;
}
close(IN);

