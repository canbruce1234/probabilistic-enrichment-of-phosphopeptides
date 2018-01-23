#!/usr/bin/perl -w
# For each peptide, calculates number of peptides having the same mass within a given mass error (1, 5, 20, 100 ppm)

($in_filename, $out_filename, $accuracy) = @ARGV;

#$in_filename  = '/home/mysql/unq_pept_test.out';
#$out_filename  = '/home/mysql/neighb_test_50ppm.out';
#$in_filename = '/home/mysql/unq_pept_masses2.out';
#$out_filename  = '/home/mysql/neighb_5ppm.out';
#$accuracy = 0.000005;

print "infile: $in_filename out: $out_filename accuracy: $accuracy \n";

open( IN,   "$in_filename" )    || die "can't open input file $in_filename: $!";
open( OUT,  ">$out_filename" )    || die "can't open output file $out_filename: $!";

$line = <IN>;
($new_id,  $new_mass, $new_phos) = unpack("A8 x1 A18 a2 x1", $line);

@neighborhood = ({	id 		=> $new_id,	
					mass  	=> $new_mass,
					phos	=> $new_phos,
					all_cnt => 0,	
					np_cnt 	=> 0, 
					},);
					
while ( defined( $line = <IN> ) ) {
	($new_id, $new_mass, $new_phos) = unpack("A8  x1 A18 a2 x1", $line);

	#create new neighborhood, which gets id and mass of newly read peptide
	push @neighborhood, ( {	id 		=> $new_id, 
							mass 	=> $new_mass, 
							phos 	=> $new_phos, 
							all_cnt => 0, 
							np_cnt 	=> 0, 
						},);

	$range = sprintf("%.6f", $accuracy * $new_mass);
		
	$sum_done = 0;	
	for $nghbr_idx (0 .. $#neighborhood-1) {
		
			if ( abs($new_mass - $neighborhood[$nghbr_idx]{mass}) <= $range ) {
					if ($neighborhood[$nghbr_idx]{phos} != 0 ) {
						if ($neighborhood[$#neighborhood]{phos} != 0 ) {
							++$neighborhood[$nghbr_idx]{all_cnt} ;
			 				++$neighborhood[$#neighborhood]{all_cnt};
			 			 }
			 			 elsif ($neighborhood[$#neighborhood]{phos} == 0 ) {
							++$neighborhood[$nghbr_idx]{all_cnt} ;
			 				++$neighborhood[$#neighborhood]{all_cnt};
							++$neighborhood[$nghbr_idx]{np_cnt};
			 			 }
					}			
			 		elsif ($neighborhood[$nghbr_idx]{phos} == 0 ) {
						if ($neighborhood[$#neighborhood]{phos} != 0 ) {
							++$neighborhood[$nghbr_idx]{all_cnt} ;
			 				++$neighborhood[$#neighborhood]{all_cnt};
			 				++$neighborhood[$#neighborhood]{np_cnt} ;
			 			 }
			 			 elsif ($neighborhood[$#neighborhood]{phos} == 0 ) {
							++$neighborhood[$nghbr_idx]{all_cnt} ;
			 				++$neighborhood[$#neighborhood]{all_cnt};
							++$neighborhood[$nghbr_idx]{np_cnt};
							++$neighborhood[$#neighborhood]{np_cnt};
			 			 }			 			 	
					}
			}
			elsif ($new_mass - $neighborhood[$nghbr_idx]{mass} > $range) {
				print OUT "$neighborhood[$nghbr_idx]{id} $neighborhood[$nghbr_idx]{mass} $range $neighborhood[$#neighborhood]{phos} $neighborhood[$nghbr_idx]{all_cnt} $neighborhood[$nghbr_idx]{np_cnt}\n";
				++$sum_done;
			}			
	}
	while (0<$sum_done) {
			shift @neighborhood;
			--$sum_done;
	}

}

# All items have been read into the @neighborhood array, 
# No more neighborhoods can be defined by adding new values outside the $range value
# Now, print remaining items in the @neighborhood array.
 
for $nghbr_idx ( 0 .. $#neighborhood ) {
	   print OUT "$neighborhood[$nghbr_idx]{id} $neighborhood[$nghbr_idx]{mass} $range $neighborhood[$#neighborhood]{phos} $neighborhood[$nghbr_idx]{all_cnt} $neighborhood[$nghbr_idx]{np_cnt}\n";
}


