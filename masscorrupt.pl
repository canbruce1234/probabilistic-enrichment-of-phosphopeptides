#!/usr/bin/perl -w
#
################################################################
#
# Script to determine how many peptides are incorrectly placed
# above or below the p=0.9 line, according to their phosphorylation
# status, if their masses have a gaussian error of 1, 5, 20 or
# 100 ppm.
#
# Counting peptides according to whether they
# are truly above (TA), truly below (TB) falsely above (FA) and
# falsely below (FB), because their mass is in error.  Also,
# enumerated peptides as phosphorylated & below line (TP=true-positive),
# nophos & below line (FP= False Positive), phosphosrylate & above line
# (FN = False Negative), and nophos & above line (TN = True negative),
# according to their phosphorylation status and error level.
#
################################################################

use POSIX;
$| = 1;     # Do not buffer output

$infile  = '/home/phospho/corrupt_test.in';
$outfile = '/home/phospho/corrupt_test.out';

open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

open( IN, $infile ) || die "can't open $infile: $!";

$time = `date`;
print "started at $time \n";

while (<IN>) {
	push @AoA, [split];
}
close IN;
$time = `date`;
print "finished reading data at $time \n";

for ( $i = 0 ; $i < 100 ; $i++ ) {
	
	print".";
	if ( ($z=$i+1) % 10 == 0 ) {
		print"$z "; 
		$time = `date`;
		chomp $time;
		print " ($time) \n";}

	$sum_mass1   = 0;
	$sum_mass5   = 0;
	$sum_mass20  = 0;
	$sum_mass100 = 0;
	$sum_dev1    = 0;
	$sum_dev5    = 0;
	$sum_dev20   = 0;
	$sum_dev100  = 0;
	$count       = 0;

	$TB_1_phos0   = 0;
	$TA_1_phos0   = 0;
	$FB_1_phos0   = 0;
	$FA_1_phos0   = 0;
	$TB_5_phos0   = 0;
	$TA_5_phos0   = 0;
	$FB_5_phos0   = 0;
	$FA_5_phos0   = 0;
	$TB_20_phos0  = 0;
	$TA_20_phos0  = 0;
	$FB_20_phos0  = 0;
	$FA_20_phos0  = 0;
	$TB_100_phos0 = 0;
	$TA_100_phos0 = 0;
	$FB_100_phos0 = 0;
	$FA_100_phos0 = 0;
	$TB_1_phos1   = 0;
	$TA_1_phos1   = 0;
	$FB_1_phos1   = 0;
	$FA_1_phos1   = 0;
	$TB_5_phos1   = 0;
	$TA_5_phos1   = 0;
	$FB_5_phos1   = 0;
	$FA_5_phos1   = 0;
	$TB_20_phos1  = 0;
	$TA_20_phos1  = 0;
	$FB_20_phos1  = 0;
	$FA_20_phos1  = 0;
	$TB_100_phos1 = 0;
	$TA_100_phos1 = 0;
	$FB_100_phos1 = 0;
	$FA_100_phos1 = 0;
	$TB_1_phos2   = 0;
	$TA_1_phos2   = 0;
	$FB_1_phos2   = 0;
	$FA_1_phos2   = 0;
	$TB_5_phos2   = 0;
	$TA_5_phos2   = 0;
	$FB_5_phos2   = 0;
	$FA_5_phos2   = 0;
	$TB_20_phos2  = 0;
	$TA_20_phos2  = 0;
	$FB_20_phos2  = 0;
	$FA_20_phos2  = 0;
	$TB_100_phos2 = 0;
	$TA_100_phos2 = 0;
	$FB_100_phos2 = 0;
	$FA_100_phos2 = 0;
	$TB_1_phos3   = 0;
	$TA_1_phos3   = 0;
	$FB_1_phos3   = 0;
	$FA_1_phos3   = 0;
	$TB_5_phos3   = 0;
	$TA_5_phos3   = 0;
	$FB_5_phos3   = 0;
	$FA_5_phos3   = 0;
	$TB_20_phos3  = 0;
	$TA_20_phos3  = 0;
	$FB_20_phos3  = 0;
	$FA_20_phos3  = 0;
	$TB_100_phos3 = 0;
	$TA_100_phos3 = 0;
	$FB_100_phos3 = 0;
	$FA_100_phos3 = 0;
	$FP_1_phos0   = 0;
	$TN_1_phos0   = 0;
	$FP_5_phos0   = 0;
	$TN_5_phos0   = 0;
	$FP_20_phos0  = 0;
	$TN_20_phos0  = 0;
	$FP_100_phos0 = 0;
	$TN_100_phos0 = 0;
	$TP_1_phos1   = 0;
	$FN_1_phos1   = 0;
	$TP_5_phos1   = 0;
	$FN_5_phos1   = 0;
	$TP_20_phos1  = 0;
	$FN_20_phos1  = 0;
	$TP_100_phos1 = 0;
	$FN_100_phos1 = 0;
	$TP_1_phos2   = 0;
	$FN_1_phos2   = 0;
	$TP_5_phos2   = 0;
	$FN_5_phos2   = 0;
	$TP_20_phos2  = 0;
	$FN_20_phos2  = 0;
	$TP_100_phos2 = 0;
	$FN_100_phos2 = 0;
	$TP_1_phos3   = 0;
	$FN_1_phos3   = 0;
	$TP_5_phos3   = 0;
	$FN_5_phos3   = 0;
	$TP_20_phos3  = 0;
	$FN_20_phos3  = 0;
	$TP_100_phos3 = 0;
	$FN_100_phos3 = 0;
	$sum_mass1    = 0;
	$sum_mass5    = 0;
	$sum_mass20   = 0;
	$sum_mass100  = 0;
	$sum_dev1     = 0;
	$sum_dev5     = 0;
	$sum_dev20    = 0;
	$sum_dev100   = 0;
	$count        = 0;

	for $row (@AoA) {
		$phos_num = @$row[0];
		$mass     = @$row[1];

		$corrupt_mass_1   = $mass * ( 1 + 0.000001 * &ltqnorm(rand) );
		$corrupt_mass_5   = $mass * ( 1 + 0.000005 * &ltqnorm(rand) );
		$corrupt_mass_20  = $mass * ( 1 + 0.000020 * &ltqnorm(rand) );
		$corrupt_mass_100 = $mass * ( 1 + 0.000100 * &ltqnorm(rand) );

		$border0   = 0.000457 * $mass - 0.0448;
		$border1   = 0.000457 * $corrupt_mass_1 - 0.0448;
		$border5   = 0.000457 * $corrupt_mass_5 - 0.0448;
		$border20  = 0.000457 * $corrupt_mass_20 - 0.0448;
		$border100 = 0.000457 * $corrupt_mass_100 - 0.0448;

		$massDeff_0 =
		  floor(
			0.0004978 * $mass - 0.588105 - floor( 0.0004978 * $mass - 0.588105 )
			  - ( $mass - int $mass ) ) +
		  floor( 0.0004978 * $mass - 0.588105 + 1 ) + ( $mass - int $mass );

		$massDeff_1 =
		  floor( 0.0004978 * $corrupt_mass_1 - 0.588105 -
			  floor( 0.0004978 * $corrupt_mass_1 - 0.588105 ) -
			  ( $corrupt_mass_1 - int $corrupt_mass_1 ) ) +
		  floor( 0.0004978 * $corrupt_mass_1 - 0.588105 + 1 ) +
		  ( $corrupt_mass_1 - int $corrupt_mass_1 );

		$massDeff_5 =
		  floor( 0.0004978 * $corrupt_mass_5 - 0.588105 -
			  floor( 0.0004978 * $corrupt_mass_5 - 0.588105 ) -
			  ( $corrupt_mass_5 - int $corrupt_mass_5 ) ) +
		  floor( 0.0004978 * $corrupt_mass_5 - 0.588105 + 1 ) +
		  ( $corrupt_mass_5 - int $corrupt_mass_5 );

		$massDeff_20 =
		  floor( 0.0004978 * $corrupt_mass_20 - 0.588105 -
			  floor( 0.0004978 * $corrupt_mass_20 - 0.588105 ) -
			  ( $corrupt_mass_20 - int $corrupt_mass_20 ) ) +
		  floor( 0.0004978 * $corrupt_mass_20 - 0.588105 + 1 ) +
		  ( $corrupt_mass_20 - int $corrupt_mass_20 );

		$massDeff_100 =
		  floor( 0.0004978 * $corrupt_mass_100 - 0.588105 -
			  floor( 0.0004978 * $corrupt_mass_100 - 0.588105 ) -
			  ( $corrupt_mass_100 - int $corrupt_mass_100 ) ) +
		  floor( 0.0004978 * $corrupt_mass_100 - 0.588105 + 1 ) +
		  ( $corrupt_mass_100 - int $corrupt_mass_100 );

	  SWITCH11: {
			if (   ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$FP_1_phos0;
				++$TB_1_phos0;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$FB_1_phos0;
				++$FP_1_phos0;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$FA_1_phos0;
				++$TN_1_phos0;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$TA_1_phos0;
				++$TN_1_phos0;
				last SWITCH11;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$TB_1_phos1;
				++$TP_1_phos1;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$FB_1_phos1;
				++$TP_1_phos1;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$FA_1_phos1;
				++$FN_1_phos1;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$TA_1_phos1;
				++$FN_1_phos1;
				last SWITCH11;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$TB_1_phos2;
				++$TP_1_phos2;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$FB_1_phos2;
				++$TP_1_phos2;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$FA_1_phos2;
				++$FN_1_phos2;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$TA_1_phos2;
				++$FN_1_phos2;
				last SWITCH11;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$TB_1_phos3;
				++$TP_1_phos3;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 <= $border1 ) )
			{
				++$FB_1_phos3;
				++$TP_1_phos3;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$FA_1_phos3;
				++$FN_1_phos3;
				last SWITCH11;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_1 > $border1 ) )
			{
				++$TA_1_phos3;
				++$FN_1_phos3;
				last SWITCH11;
			} else {
				$nothing = 1;
			}

		}    # end SWITCH11

	  SWITCH12: {
			if (   ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$FP_5_phos0;
				++$TB_5_phos0;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$FB_5_phos0;
				++$FP_5_phos0;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$FA_5_phos0;
				++$TN_5_phos0;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$TA_5_phos0;
				++$TN_5_phos0;
				last SWITCH12;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$TB_5_phos1;
				++$TP_5_phos1;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$FB_5_phos1;
				++$TP_5_phos1;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$FA_5_phos1;
				++$FN_5_phos1;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$TA_5_phos1;
				++$FN_5_phos1;
				last SWITCH12;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$TB_5_phos2;
				++$TP_5_phos2;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$FB_5_phos2;
				++$TP_5_phos2;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$FA_5_phos2;
				++$FN_5_phos2;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$TA_5_phos2;
				++$FN_5_phos2;
				last SWITCH12;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$TB_5_phos3;
				++$TP_5_phos3;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 <= $border5 ) )
			{
				++$FB_5_phos3;
				++$TP_5_phos3;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$FA_5_phos3;
				++$FN_5_phos3;
				last SWITCH12;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_5 > $border5 ) )
			{
				++$TA_5_phos3;
				++$FN_5_phos3;
				last SWITCH12;
			} else {
				$nothing = 1;
			}

		}    # end SWITCH12

	  SWITCH13: {
			if (   ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$FP_20_phos0;
				++$TB_20_phos0;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$FB_20_phos0;
				++$FP_20_phos0;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$FA_20_phos0;
				++$TN_20_phos0;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$TA_20_phos0;
				++$TN_20_phos0;
				last SWITCH13;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$TB_20_phos1;
				++$TP_20_phos1;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$FB_20_phos1;
				++$TP_20_phos1;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$FA_20_phos1;
				++$FN_20_phos1;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$TA_20_phos1;
				++$FN_20_phos1;
				last SWITCH13;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$TB_20_phos2;
				++$TP_20_phos2;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$FB_20_phos2;
				++$TP_20_phos2;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$FA_20_phos2;
				++$FN_20_phos2;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$TA_20_phos2;
				++$FN_20_phos2;
				last SWITCH13;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$TB_20_phos3;
				++$TP_20_phos3;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 <= $border20 ) )
			{
				++$FB_20_phos3;
				++$TP_20_phos3;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$FA_20_phos3;
				++$FN_20_phos3;
				last SWITCH13;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_20 > $border20 ) )
			{
				++$TA_20_phos3;
				++$FN_20_phos3;
				last SWITCH13;
			} else {
				$nothing = 1;
			}

		}    # end SWITCH13

	  SWITCH14: {
			if (   ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$FP_100_phos0;
				++$TB_100_phos0;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$FB_100_phos0;
				++$FP_100_phos0;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$FA_100_phos0;
				++$TN_100_phos0;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 0 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$TA_100_phos0;
				++$TN_100_phos0;
				last SWITCH14;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$TB_100_phos1;
				++$TP_100_phos1;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$FB_100_phos1;
				++$TP_100_phos1;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$FA_100_phos1;
				++$FN_100_phos1;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 1 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$TA_100_phos1;
				++$FN_100_phos1;
				last SWITCH14;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$TB_100_phos2;
				++$TP_100_phos2;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$FB_100_phos2;
				++$TP_100_phos2;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$FA_100_phos2;
				++$FN_100_phos2;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 2 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$TA_100_phos2;
				++$FN_100_phos2;
				last SWITCH14;
			} else {
				$nothing = 1;
			}

			if (   ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$TB_100_phos3;
				++$TP_100_phos3;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 <= $border100 ) )
			{
				++$FB_100_phos3;
				++$TP_100_phos3;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 <= $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$FA_100_phos3;
				++$FN_100_phos3;
				last SWITCH14;
			} elsif ( ( $mass < 4000 )
				&& ( $phos_num == 3 )
				&& ( $massDeff_0 > $border0 )
				&& ( $massDeff_100 > $border100 ) )
			{
				++$TA_100_phos3;
				++$FN_100_phos3;
				last SWITCH14;
			} else {
				$nothing = 1;
			}
		}    # end SWITCH14

		$sum_mass1   = $sum_mass1 + $corrupt_mass_1 / $mass;
		$sum_mass5   = $sum_mass5 + $corrupt_mass_5 / $mass;
		$sum_mass20  = $sum_mass20 + $corrupt_mass_20 / $mass;
		$sum_mass100 = $sum_mass100 + $corrupt_mass_100 / $mass;
		$sum_dev1    = $sum_dev1 + ( ( $mass - $corrupt_mass_1 ) / $mass )**2;
		$sum_dev5    = $sum_dev5 + ( ( $mass - $corrupt_mass_5 ) / $mass )**2;
		$sum_dev20   = $sum_dev20 + ( ( $mass - $corrupt_mass_20 ) / $mass )**2;
		$sum_dev100  =
		  $sum_dev100 + ( ( $mass - $corrupt_mass_100 ) / $mass )**2;
		++$count;

		#	}    # end while
	}   #end for loop to read @AoA

		$mean1   = $sum_mass1 / $count;
		$mean5   = $sum_mass5 / $count;
		$mean20  = $sum_mass20 / $count;
		$mean100 = $sum_mass100 / $count;
		$sd1     = sqrt( $sum_dev1 / $count );
		$sd5     = sqrt( $sum_dev5 / $count );
		$sd20    = sqrt( $sum_dev20 / $count );    
		$sd100   = sqrt( $sum_dev100 / $count );

		print OUT
		  "$TB_1_phos0 $TA_1_phos0 $FB_1_phos0 $FA_1_phos0 ",      # fields 1-4
		  "$TB_5_phos0 $TA_5_phos0 $FB_5_phos0 $FA_5_phos0 ",      # fields 5-8
		  "$TB_20_phos0 $TA_20_phos0 $FB_20_phos0 $FA_20_phos0 ",  # fields 9-12
		  "$TB_100_phos0 $TA_100_phos0 $FB_100_phos0 $FA_100_phos0 "
		  ,                                                       # fields 13-16
		  "$TB_1_phos1 $TA_1_phos1 $FB_1_phos1 $FA_1_phos1 ",     # fields 17-20
		  "$TB_5_phos1 $TA_5_phos1 $FB_5_phos1 $FA_5_phos1 ",     # fields 21-24
		  "$TB_20_phos1 $TA_20_phos1 $FB_20_phos1 $FA_20_phos1 ", # fields 25-28
		  "$TB_100_phos1 $TA_100_phos1 $FB_100_phos1 $FA_100_phos1 "
		  ,                                                       # fields 29-32
		  "$TB_1_phos2 $TA_1_phos2 $FB_1_phos2 $FA_1_phos2 ",     # fields 33-36
		  "$TB_5_phos2 $TA_5_phos2 $FB_5_phos2 $FA_5_phos2 ",     # fields 37-40
		  "$TB_20_phos2 $TA_20_phos2 $FB_20_phos2 $FA_20_phos2 ", # fields 41-44
		  "$TB_100_phos2 $TA_100_phos2 $FB_100_phos2 $FA_100_phos2 "
		  ,                                                       # fields 45-48
		  "$TB_1_phos3 $TA_1_phos3 $FB_1_phos3 $FA_1_phos3 ",     # fields 49-52
		  "$TB_5_phos3 $TA_5_phos3 $FB_5_phos3 $FA_5_phos3 ",     # fields 53-56
		  "$TB_20_phos3 $TA_20_phos3 $FB_20_phos3 $FA_20_phos3 ", # fields 57-60
		  "$TB_100_phos3 $TA_100_phos3 $FB_100_phos3 $FA_100_phos3 "
		  ,                                     # fields 61-64
		  "$FP_1_phos0 $TN_1_phos0 ",           # fields 65-66
		  "$FP_5_phos0 $TN_5_phos0 ",           # fields 67-68
		  "$FP_20_phos0 $TN_20_phos0 ",         # fields 69-70
		  "$FP_100_phos0 $TN_100_phos0 ",       # fields 71-72
		  "$TP_1_phos1 $FN_1_phos1 ",           # fields 73-744
		  "$TP_5_phos1 $FN_5_phos1 ",           # fields 75-76
		  "$TP_20_phos1 $FN_20_phos1 ",         # fields 77-78
		  "$TP_100_phos1 $FN_100_phos1 ",       # fields 79-80
		  "$TP_1_phos2 $FN_1_phos2 ",           # fields 81-82
		  "$TP_5_phos2 $FN_5_phos2 ",           # fields 83-84
		  "$TP_20_phos2 $FN_20_phos2 ",         # fields 85-86
		  "$TP_100_phos2 $FN_100_phos2 ",       # fields 87-88
		  "$TP_1_phos3 $FN_1_phos3 ",           # fields 89-90
		  "$TP_5_phos3 $FN_5_phos3 ",           # fields 91-92
		  "$TP_20_phos3 $FN_20_phos3 ",         # fields 93-94
		  "$TP_100_phos3 $FN_100_phos3 ",       # fields 95-96
		  "$mean1 $mean5 $mean20 $mean100 ",    # fields 97-100
		  "$sd1 $sd5 $sd20 $sd100 ",            # fields 101-104
		  "$count ",                            # field 105
		  "\n";

		# end for loop
	}

	close OUT;
$time = `date`;
print "finished at $time \n";

	sub ltqnorm ($) {

		#
		# Lower tail quantile for standard normal distribution function.
		#
		# This function returns an approximation of the inverse cumulative
		# standard normal distribution function.  I.e., given P, it returns
		# an approximation to the X satisfying P = Pr{Z <= X} where Z is a
		# random variable from the standard normal distribution.
		#
		# The algorithm uses a minimax approximation by rational functions
		# and the result has a relative error whose absolute value is less
		# than 1.15e-9.
		#
		# Author:      Peter J. Acklam
		# Time-stamp:  2000-07-19 18:26:14
		# E-mail:      pjacklam@online.no
		# WWW URL:     http://home.online.no/~pjacklam

		my $p = shift;
		die "input argument must be in (0,1)\n" unless 0 < $p && $p < 1;

		# Coefficients in rational approximations.
		my @a = (
			-3.969683028665376e+01, 2.209460984245205e+02,
			-2.759285104469687e+02, 1.383577518672690e+02,
			-3.066479806614716e+01, 2.506628277459239e+00
		);
		my @b = (
			-5.447609879822406e+01, 1.615858368580409e+02,
			-1.556989798598866e+02, 6.680131188771972e+01,
			-1.328068155288572e+01
		);
		my @c = (
			-7.784894002430293e-03, -3.223964580411365e-01,
			-2.400758277161838e+00, -2.549732539343734e+00,
			4.374664141464968e+00,  2.938163982698783e+00
		);
		my @d = (
			7.784695709041462e-03, 3.224671290700398e-01,
			2.445134137142996e+00, 3.754408661907416e+00
		);

		# Define break-points.
		my $plow  = 0.02425;
		my $phigh = 1 - $plow;

		# Rational approximation for lower region:
		if ( $p < $plow ) {
			my $q = sqrt( -2 * log($p) );
			return (
				(
					( ( ( $c[0] * $q + $c[1] ) * $q + $c[2] ) * $q + $c[3] ) *
					  $q + $c[4]
				) * $q + $c[5]
			  ) /
			  ( ( ( ( $d[0] * $q + $d[1] ) * $q + $d[2] ) * $q + $d[3] ) * $q +
				  1 );
		}

		# Rational approximation for upper region:
		if ( $phigh < $p ) {
			my $q = sqrt( -2 * log( 1 - $p ) );
			return -(
				(
					( ( ( $c[0] * $q + $c[1] ) * $q + $c[2] ) * $q + $c[3] ) *
					  $q + $c[4]
				) * $q + $c[5]
			  ) /
			  ( ( ( ( $d[0] * $q + $d[1] ) * $q + $d[2] ) * $q + $d[3] ) * $q +
				  1 );
		}

		# Rational approximation for central region:
		my $q = $p - 0.5;
		my $r = $q * $q;
		return (
			(
				( ( ( $a[0] * $r + $a[1] ) * $r + $a[2] ) * $r + $a[3] ) * $r +
				  $a[4]
			) * $r + $a[5]
		  ) * $q / (
			(
				( ( ( $b[0] * $r + $b[1] ) * $r + $b[2] ) * $r + $b[3] ) * $r +
				  $b[4]
			) * $r + 1
		  );
	}
