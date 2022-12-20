#!/usr/bin/perl # -w

use strict;

my (@eachRead, $S1, $S2, $S3, $S4, $S5, $S6, $S7, $S8, $Total);


open(FIN, "genes.read_group_tracking") or die "can not open file\n";

open(FOUT, ">FPKM+4x2.txt") or die "can not open this file\n";


<FIN>;


while (<FIN>)

{ 
	@eachRead = split(/\t/, $_);	# space	

	$S1  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S2  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S3  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S4  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S5  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S6  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S7  = $eachRead[6]; 

	@eachRead = split(/\t/, <FIN>);	# space	

	$S8  = $eachRead[6]; 

	



	print FOUT qq($eachRead[0]\t$S1\t$S2\t$S3\t$S4\t$S5\t$S6\t$S7\t$S8\n);
} 


close(FIN);

close(FOUT);

exit(0);