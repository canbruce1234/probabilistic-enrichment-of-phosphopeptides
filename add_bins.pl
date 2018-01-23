#!/usr/bin/perl  -w
# Populate table BINS with rows corresponding to the pixels of Figure 2

$" = "\t";
use DBI;
     my $dbh = DBI->connect( "dbi:mysql:phospho", "root", $password )
         or die "Can't connect to mysql database: $DBI::errstr\n";
        $mass=0;
        while ($mass < 10000){
        	$fract = 0.00;
        	$mass =$mass+10;
        	while ($fract < 1.00) {
        		my $stmt= $dbh->prepare( 
                	"insert into bins(mass_bin, fract_bin) values ($mass, $fract);"
                	) or die "can't prepare stmt: $DBI::errstr";
         		$stmt->execute()  or die "Can't execute SQL statement: $DBI::errstr\n";
                $fract = $fract +0.01; 
        	}
        print "$mass "
        }

        
       warn "Data fetching terminated early by error: $DBI::errstr\n"
           if $DBI::err;


### Now, disconnect from the database
     $dbh->disconnect
         or warn "Disconnection failed: $DBI::errstr\n";
