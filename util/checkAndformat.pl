#! /usr/local/bi/nperl

use strict;

my $input = $ARGV[0];
my $type = $ARGV[1];
my $size = $ARGV[2];

my $num = 0;
open(IN, $input) || die "cannot open $!";
while(<IN>) {

    s/[\r\n]//g;

    my @F = split("\t", $_);
    my @FF = split(";", $F[0]);

    my $sample = join(";", @FF[0 .. ($#FF - 2)]);
    my $ref = $FF[$#FF - 1];
    my $alt = $FF[$#FF];

    if ($ref ne substr($F[1], $size, 1)) {
        die $num . " th inconsistent\n";
    }

    print $type . "_" . $num . "\t" . $sample . "\t" . $F[1] . "\t" . $alt . "\n";

    $num = $num + 1;
}
close(IN);

