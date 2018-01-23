#!/usr/bin/perl -w
####################################################
# Trypsin digestion human proteins from IPI, and calculate their masses.
# Generates all possible phosphorylations (up to 3 sites per peptide),
# and up to 2 missed trypsin cuts per peptide.
# Peptides with methioinine are calculated with all possible variants that are not/partially/fully
# oxidized. Cysteines are assumed to be carbamidomethylated. 
# In addition to the exact mass for each peptide, 4 corrupted mass values are generated, 
# corresponding to 1, 5 ,20, 100 ppm measurement error. 
#
####################################################
#use Math::BigInt;
#use strict; 
use Memoize; 
memoize('myfactorial');

$infile  = '/home/cbruce/ipi.HUMAN.v3.07.fasta';
#$infile  = '/home/cbruce/test1163.fasta';
$outfile   = '/home/cbruce/human_peptides2.out';
$outfile2  = '/home/cbruce/human_proteins3.out';
$badfile   = '/home/cbruce/human_badproteins2.out';

$maxMissed = 2;
$maxPhos   = 3;

open( IN,   $infile )    || die "can't open $infile: $!";
open( OUT,  ">$outfile" )  || die "can't ope $outfile: $!";
open( OUT2, ">$outfile2" ) || die "can't open $outfile2: $!";
open( BAD,  ">$badfile" )  || die "can't open $badfile: $!";
$proteinSeq = '';
$IPI        = '';

# $countProtein = 0;
while ( defined( $line = <IN> ) ) {
	if ( $line =~ m/^>/ ) {
		# print "IPI#= $IPI \n";
		processProtein( $IPI, $proteinSeq, $maxMissed, $maxPhos );
#		$line =~ /^>IPI:(.*?)\|(.*\|)*(.*:.*?)*(.*\=.*? )*(.*)$/; # $2  holds other ids
#		$line =~ /^>IPI:(.*?)\|.*? Tax_Id=\d* (.*)$/;
#		$description = $2;
		$all_swissprot ='';
		$all_trembl='';
 		$line =~ /^>IPI:(.*?)\|(SWISS-PROT:(.*?)\|){0,1}(TREMBL:(.*?)\|){0,1}.*? Tax_Id=\d* (.*)$/;
		$IPI         = $1;
		$swissprot_string   = $3 if defined $3;
		$trembl_string = $5 if defined $5;
		# $description = $6 if defined $6;

		if ($trembl_string =~ /\;/) {
			@all_trembl = split /;/, $trembl_string;
			} else {
			@all_trembl = ($trembl_string);
        }
		if (scalar(@all_trembl) > 0) { 
			foreach $trembl (@all_trembl) {		
			print 	OUT2	"$IPI\t\t$trembl\t\n";
		    }
		}
		if ($swissprot_string =~ /\;/) {
			@all_swissprot = split /;/, $swissprot_string;
			} else {
			@all_swissprot = ($swissprot_string);
        }
		 if (scalar(@all_swissprot) > 0) { 
			foreach $swissprot (@all_swissprot) {
			print 	OUT2	"$IPI\t$swissprot\t\t\n";
			}
		 }
		next;
	}
}

close(IN);
close(OUT2);
close(OUT);
close(BAD);

sub processProtein {
	my ( $IPI, $sequence, $maxMissed, $maxPhos ) = @_;
	my $begin   = 0;
	my $end     = 0;
	my $peptide = '';

	if ( $sequence eq '' ) { return }

	@sites = getsites($sequence);    #cutting sites in protein for trypsin
	$peptideCount = 0;

	# loop to get peptides with 0, 1, 2, ...missed cutting sites:
	for ( $missed = 0 ; $missed < $maxMissed + 1 ; $missed++ ) {
		@peptides = trypsinDigest( $sequence, $missed, @sites );
		
		for $row (@peptides) {
			$countSTY = 0;   #set counter of phosphorylatable sites to zero
			( $begin, $end, $peptide ) = @$row;
	 	  if ( $peptide =~ /[^ARNDCEQGHILKMFPSTWYV]/ ) {
		 	print BAD "$IPI $missed $begin $end $peptide \n";
 	 	  next;
	 	  }

			$_ = $peptide;
			$countA = tr/A//; # return value of tr/// operator is number of characters matched in old string
			$countR = tr/R//;
			$countN = tr/N//;
			$countD = tr/D//;
			$countC = tr/C//;
			$countE = tr/E//;
			$countQ = tr/Q//;
			$countG = tr/G//;
			$countH = tr/H//;
			$countI = tr/I//;
			$countL = tr/L//;
			$countK = tr/K//;
			$countM = tr/M//;
			$countF = tr/F//;
			$countP = tr/P//;
			$countS = tr/S//;
			$countT = tr/T//;
			$countW = tr/W//;
			$countY = tr/Y//;
			$countV = tr/V//;
			$countSTY = $countS + $countT + $countY ;
			
			
			# my $countM = 0;
			# my @allMsites = ();
			# while ( $peptide =~ m/M/g ) {
			#	push @allMsites, pos $peptide;
			#	$countM = scalar @allMsites;
			# }
			
			for ( $oxy = 0 ; $oxy < $countM + 1 ; $oxy++ ) {
				$num_oxyMet = $oxy;
				
				# record whether peptide is not/partially/fully (N/P/F) oxidized on methionines
				if($num_oxyMet == 0) {$stat_oxyMet = 'N'}
				elsif($num_oxyMet == $countM) {$stat_oxyMet = 'F'}
				else {$stat_oxyMet = 'P'}
				
				# calculate number of isomers of this peptide with this number of Met oxidation
				# nCk = n!/(k!(n-k!))
				# $oxyMet_isomers = myfactorial($countM)/(myfactorial($num_oxyMet)*myfactorial($countM-$num_oxyMet));

				#phosphorylation block:

				# my @allPSites = ();
				# while ( $peptide =~ m/[STY]/g ) {
				#	push @allPSites, pos $peptide;
				#	$countSTY = scalar @allPSites;
				# }

                
                # Max number of phoshates on this peptide to be the least of $maxPhos or $countSTY
				if ( $maxPhos < $countSTY ) {
					$thisMaxPhos = $maxPhos;
				} else {
					$thisMaxPhos = $countSTY;
				}

				for ( $phos = 0 ; $phos < $thisMaxPhos + 1 ; $phos++ ) {

					$mass = massOfPeptide($peptide) 
					        + $phos * 79.96633041 
					        + $num_oxyMet * 15.994914622 ;
					        
					$massDef = massDeficitOfPeptide($peptide)
							+ $phos * (-0.03366959)
							+ $num_oxyMet * (-0.0050853780);
					        
					$corrupt_mass_1 = $mass * ( 1 + 0.000001 * &ltqnorm(rand));
					$corrupt_mass_5 = $mass * ( 1 + 0.000005 * &ltqnorm(rand));
					$corrupt_mass_20 = $mass * ( 1 + 0.000020 * &ltqnorm(rand));
					$corrupt_mass_100 = $mass * ( 1 + 0.000100 * &ltqnorm(rand));
					        
					# calculate number of isomers of this peptide with this number of phosphates
					# nCk = n!/(k!(n-k!))
					# $phos_isomers = myfactorial($thisMaxPhos)/(myfactorial($phos)*myfactorial($thisMaxPhos-$phos));

					++$peptideCount;
					
					print  OUT
 "$IPI $peptideCount $mass $massDef $corrupt_mass_1 $corrupt_mass_5 $corrupt_mass_20 $corrupt_mass_100 ",
 "$missed $begin $end $stat_oxyMet $num_oxyMet ",
 "$countA $countR $countN $countD $countC $countE $countQ $countG $countH $countI ",
 "$countL $countK $countM $countF $countP $countS $countT $countW $countY $countV ",
 "$phos ",
 "$peptide \n";

					# end phosphorylation block

				}
			}
		}
	}
}

sub getsites {
	my ($sequence) = @_;
	@sites = ();
	while ( $sequence =~ /([KR])/g ) {
		if ( substr( $', 0, 1 ) eq 'P' ) { next; }
		push( @sites, pos($sequence) );
	}
	return @sites;
}

sub trypsinDigest {
	my ( $sequence, $missedSites, @sites ) = @_;
	my $begin = 1;
	my $end   = 0;
	my $peptide;
	my @line = ();
	my $begin_idx = 0;
	my $end_idx = 1;
	#	my %peptides = ();
	my @peptides = ();

	for ( $i = 0 ; $i < scalar(@sites)-$missedSites +1 ; $i++ ) {
			 if (scalar(@sites) == 0){
			 	$begin = 1; 
			 	$end = length($sequence); 
			 	$peptide = $sequence; 
			 	@line = ( $begin + 1, $end, $peptide );
			 	push @peptides, [@line];
			 	last;}
		 $end_idx = $i + $missedSites;
		 if ($end_idx > scalar(@sites)-1 ) {
		 	$end = length($sequence);
	     	$begin = $sites[$begin_idx]+1;
		  } 
		 else {
			$end = $sites[$end_idx];
		}
		$peptide = substr( $sequence, $begin-1, $end - $begin+1);

		@line = ( $begin, $end, $peptide );
		push @peptides, [@line];

		 $begin_idx = $i;
		 if (defined ($sites[$begin_idx] )) {$begin = $sites[$begin_idx]+1 } 
	}

	return @peptides;
}

sub phosphorylate {
	my ( $sequence, $numPhos ) = @_;
	my $mypeptide = '';
	my @peptides  = ();
	my @allPSites = ();
	use Math::Combinatorics;
	@indexes = ();

	while ( $sequence =~ m/[STY]/g ) {
		push @allPSites, pos $sequence;
	}

	for ( $i = 0 ; $i < scalar(@allPSites) ; $i++ ) {
		push @indexes, $i;
	}

	@allCombos = combine( $numPhos, @indexes );
	for $row (@allCombos) {

		@combo = @$row;
		@sites = ();
		for ( $i = 0 ; $i < scalar(@combo) ; $i++ ) {
			;
			$index = $combo[$i];
			push( @sites, $allPSites[$index] );
		}

		$mypeptide = $sequence;
		for ( $i = 0 ; $i < scalar(@sites) ; $i++ ) {
			$residueNum = $sites[$i];

			substr( $mypeptide, $residueNum - 1, 1 ) =~ tr/STY/sty/;

		}

		push @peptides, $mypeptide;
	}

	return @peptides;

}

sub massOfPeptide {
	my ($sequence) = @_;
	%mass = (
		A => 71.037113787,
		R => 156.101111026,
		N => 114.042927446,
		D => 115.026943031,
		C => 103.009184477 + 57.021463723 ,    # carbamidomethylation
		E => 129.042593095,
		Q => 128.058577510,
		G => 57.021463723,
		H => 137.058911861,
		I => 113.084063979,
		L => 113.084063979,
		K => 128.094963016,
		M => 131.040484605,     # not oxidized
		F => 147.068413915,
		P => 97.052763851,
		S => 87.032028409,
		T => 101.047678473,
		W => 186.079312952,
		Y => 163.063328537,
		V => 99.068413915,
		s => 87.032028409 + 79.96633041,       # Ser + HPO3
		t => 101.047678473 + 79.96633041,      # Thr + HPO3
		y => 163.063328537 + 79.96633041       # Tyr + HPO3
	);
	$mass_H2O = 18.010564690;
	$mass_Hyd = 1.007276470;                   # H+, not H

	@residues = split( //, $sequence );
	$total_mass = 0;
	foreach $residue (@residues) {

		# print "res=$residue $mass{$residue}\n";
		$total_mass = $total_mass + $mass{$residue};
	}
	$total_mass = $total_mass + $mass_H2O + $mass_Hyd;
	return $total_mass;
}

sub massDeficitOfPeptide {
	my ($sequence) = @_;
	%mass = (
		A => 0.037113787,
		R => 0.101111026,
		N => 0.042927446,
		D => 0.026943031,
		C => 0.009184477 + 0.021463723 ,    # carbamidomethylation
		E => 0.042593095,
		Q => 0.05857751,
		G => 0.021463723,
		H => 0.058911861,
		I => 0.084063979,
		L => 0.084063979,
		K => 0.094963016,
		M => 0.040484605,     # not oxidized
		F => 0.068413915,
		P => 0.052763851,
		S => 0.032028409,
		T => 0.047678473,
		W => 0.079312952,
		Y => 0.063328537,
		V => 0.068413915
	);
	$mass_H2O = 0.01056469;
	$mass_Hyd = 0.00727647;                   # H+, not H

	@residues = split( //, $sequence );
	$total_mass = 0;
	foreach $residue (@residues) {

		$total_mass = $total_mass + $mass{$residue};
	}
	$total_mass = $total_mass + $mass_H2O + $mass_Hyd;
	return $total_mass;
}

sub oxymet {
	my $sequence = @_;
	use Math::Combinatorics;

	# get position of methionines
	@Msites = ();
	while ( $sequence =~ /M/g ) {
		push( @Msites, pos($sequence) );
	}

	# create array of (1,2,..,N) where N is number of Methionines in sequence
	@indexes = ();
	for ( $i = 0 ; $i < scalar(@Msites) ; $i++ ) {
		push @indexes, $i;
	}

	for ( $oxy = 0 ; $oxy < scalar(@Msites) + 1 ; $oxy++ ) {

		# get combination of all arrays of $oxy elements for array @indexes
		@allCombos = combine( $oxy, @indexes );
		for $row (@allCombos) {

			#   $row is one element of the array of arrays @allCombos
			@combo = @$row;
			@sites = ();

# convert array of indexes to array of locations in peptide; e.g., (0,1,2) -> (3, 12, 23)
			for ( $i = 0 ; $i < scalar(@combo) ; $i++ ) {
				;
				$index = $combo[$i];
				push( @sites, $Msites[$index] );
			}

			$mypeptide = $sequence;

			# convert M to m in sequence at positions corresponding to @sites
			for ( $i = 0 ; $i < scalar(@sites) ; $i++ ) {
				$residueNum = $sites[$i];

				#    $oldstr = substr($mypeptide, $residueNum, 1, 'x');
				substr( $mypeptide, $residueNum - 1, 1 ) =~ tr/M/m/;

			}

			push @peptides, $mypeptide;
		}
	}

	return @peptides;
}

sub myfactorial {
    my $n = shift;
    my $s = 1;
    $s *= $n-- while $n > 0;
   return $s;
}

sub ltqnorm ($) {
    #
    # (code copied from below URL) 
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
    my @a = (-3.969683028665376e+01,  2.209460984245205e+02,
             -2.759285104469687e+02,  1.383577518672690e+02,
             -3.066479806614716e+01,  2.506628277459239e+00);
    my @b = (-5.447609879822406e+01,  1.615858368580409e+02,
             -1.556989798598866e+02,  6.680131188771972e+01,
             -1.328068155288572e+01 );
    my @c = (-7.784894002430293e-03, -3.223964580411365e-01,
             -2.400758277161838e+00, -2.549732539343734e+00,
              4.374664141464968e+00,  2.938163982698783e+00);
    my @d = ( 7.784695709041462e-03,  3.224671290700398e-01,
              2.445134137142996e+00,  3.754408661907416e+00);

    # Define break-points.
    my $plow  = 0.02425;
    my $phigh = 1 - $plow;

    # Rational approximation for lower region:
    if ( $p < $plow ) {
       my $q  = sqrt(-2*log($p));
       return ((((($c[0]*$q+$c[1])*$q+$c[2])*$q+$c[3])*$q+$c[4])*$q+$c[5]) /
               (((($d[0]*$q+$d[1])*$q+$d[2])*$q+$d[3])*$q+1);
    }

    # Rational approximation for upper region:
    if ( $phigh < $p ) {
       my $q  = sqrt(-2*log(1-$p));
       return -((((($c[0]*$q+$c[1])*$q+$c[2])*$q+$c[3])*$q+$c[4])*$q+$c[5]) /
                (((($d[0]*$q+$d[1])*$q+$d[2])*$q+$d[3])*$q+1);
    }

    # Rational approximation for central region:
    my $q = $p - 0.5;
    my $r = $q*$q;
    return ((((($a[0]*$r+$a[1])*$r+$a[2])*$r+$a[3])*$r+$a[4])*$r+$a[5])*$q /
           ((((($b[0]*$r+$b[1])*$r+$b[2])*$r+$b[3])*$r+$b[4])*$r+1);
}
