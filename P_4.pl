#!/usr/bin/perl -w
use strict;
use Getopt::Long();

sub usage {
	print 
		"Usage:","\n",
		"\t","perl task4.pl -i data.txt -o output.txt -mr-min 1000 -mr-max 1550  -bin 50 -error 1.2Da -ion-with P M -ion-start M","\n",
		"\t","-man/-help","\t\t","read complete manual page","\n",
		"Mandatory Parameters:","\n",
		"\t","-i","\t\t\t","input file","\n",
		"\t","-o","\t\t\t","output file","\n",
		"Options:","\n",
		"\t","-mr-max","\t\t\t","upper limit of mass range","\n",
		"\t","-mr-min","\t\t\t","lower limit of mass range","\n",
		"\t","-bin","\t\t\t","bin size","\n",
		"\t","-error","\t\t\t","give value of mass error (ppm/Da), default error=0. ex. 0.5Da","\n",
		"\t","-ion-start","\t\t","only calculate sequences start with specific amino acids","\n",
		"\t","-ion-with","\t\t","only calculate sequences include specific amino acids","\n";
		
}

usage();
### read opts
my ($filename_input, $filename_output,$mass_range_max,$mass_range_min,$bin_size,$error,@ion_with,$ion_with_each);
my $ion_start="";
Getopt::Long::GetOptions(
	'i=s' => \$filename_input,
	'o=s' => \$filename_output, 
	'mr-max=s' => \$mass_range_max, 
	'mr-min=s' => \$mass_range_min, 
	'bin=s' => \$bin_size,
	'error=s' => \$error, 
	'ion-start=s' => \$ion_start,
	'ion-with=s{,}' => \@ion_with
);
#or usage("file name must be specified.") 
#unless defined $filename_input;
#print $filename_input,"\n",$filename_output;



### read file 

open(IFILE,$filename_input) or die "unable to open the file\n";
#create a list to store data from the input
my (@mass_data,@line);
my ($line,$line_array,$mass_data);


while ($line_array = <IFILE>) {
	chomp($line_array);
	@line = split ("\t",$line_array);
	#push (@_,$line);
	$line++;
    for ( my $i=0; $i< @line; $i++ ) {
	    $mass_data[$line-1][$i] = $line[$i];
	}
}
close IFILE;
#print $mass_data[0][0],"\t",$mass_data[0][5],"\n";

###################################################### 

### -mr-max -mr-min Given mass range and report number of peptides

my $total_pep;
for ( my $i=0;$i< @mass_data;$i++ ) {
	if ($mass_data[$i][2] >= $mass_range_min and $mass_data[$i][2] <= $mass_range_max) {
		$total_pep += $mass_data[$i][1];
	}
}

### match ion pattern
sub get_mass_start_with {
	my (@data) = @_;
	my (@newdata);
	my ($amino) = $ion_start;
	my ($num) = 0 ;
	for (my $i =0; $i<@data;$i++) {
		if ($data[$i][5] =~ /^$amino/) {
			for (my $j=0; $j<6; $j++) {
				$newdata[$num][$j] = $data[$i][$j];
			}
			$num++;
			#print $flag,"\t",scalar(@newdata),"\n";
		}
	} 
	return @newdata;
	print $num;
}

sub get_mass_include {
	my (@data) = @_;
	my (@newdata);
	my ($amino) = $ion_with_each;
	my ($num) = 0;
	for (my $i =0; $i<@data;$i++) {
		if ($data[$i][5] =~ /$amino/) {			
			for (my $j=0; $j<6; $j++) {
				$newdata[$num][$j] = $data[$i][$j];
				
			}
			$num++;
			#print $num,"\t",scalar(@newdata),"\n";
		}
	} 
	return @newdata;
}

if ($ion_start ne "") {
	@mass_data = get_mass_start_with(@mass_data);
	#print $mass_data[1][0],"\t",$mass_data[1][1],"\n";
	#print scalar(@mass_data);
}

for (my $i=0;$i<@ion_with;$i++) {
	$ion_with_each = $ion_with[$i];
	#print $ion_with_each,"\t";
	#print scalar(@ion_with),"\n";
	@mass_data = get_mass_include(@mass_data);
	#print scalar(@mass_data),"\n";
}
#print $mass_data[0][0];


### convert theoretic mass value to exact mass value by using Da/ppm mass accuracy
my ($mass_error,$type);
if (substr($error,length($error)-3,length($error)-1) eq "ppm") {
	$mass_error  = substr($error,0,length($error)-4);
	$type = "ppm";
} else {
	$mass_error  = substr($error,0,length($error)-3);
	$type = "Da";
}

sub Da_error_to_mass {
	my ($mass_cal) = @_;
	my ($mass_exact) = $mass_cal + $mass_error;
	return $mass_exact;
} 

sub ppm_error_to_mass {
	my ($mass_cal) = @_;
	my ($mass_exact) = (10^6 * $mass_cal) / (10^6 - $mass_error);
	return $mass_exact;
}

# process getopt of mass error

if ($type eq "ppm") {
	for (my $i=0;$i< @mass_data;$i++) {
		$mass_data[$i][2] = ppm_error_to_mass($mass_data[$i][2]);
	}
} elsif ($type eq "Da") {
	for (my $i=0;$i< @mass_data;$i++) {
		$mass_data[$i][2] = Da_error_to_mass($mass_data[$i][2]);
	}
}


### report peptides counts in m/z ranges by bins
#calculate m/z values
for (my $i=0;$i< @mass_data;$i++) {
	$mass_data[$i][6] = $mass_data[$i][2]/$mass_data[$i][3];
}
#sort mass data by m/z values
#print scalar(@mass_data);
#for (my $i=0;$i< 5;$i++) {
#	for (my $j=0;$j<7;$j++) {
#		print $mass_data[$i][$j],"\t";
#	}
#	print "\n";
#}



my @new_mass_data = sort {$a->[6] <=> $b->[6]} @mass_data;
@mass_data = @new_mass_data;



###### sliding windows

my @pepcounts_mz; # array for saving results
my $first_half_pep_counts = 0;
my $second_half_pep_counts = 0;
my $bin_counts = 0;
my $min = $mass_range_min - ($bin_size/2);
my $max = $mass_range_min;


for (my $i=1;$i<@mass_data-1;$i++) {
	# see if the bin is in the user-defined range
	if ($max <= $mass_range_max + $bin_size/2) {
		# sum up the peptides counts in the bin
		if ($mass_data[$i][6] >= $min and $mass_data[$i][6] < $max ) {
			$second_half_pep_counts += $mass_data[$i][1];
		} elsif ($mass_data[$i][6] >= $max ) {
			# see if it is the first bin in the loop
			if ($bin_counts == 0) {
				$first_half_pep_counts = $second_half_pep_counts;
				$bin_counts ++;
				$min = $max;
				$max = $min + ($bin_size/2);
				$i--;
			} else {
				$pepcounts_mz[$bin_counts-1][0] = $min;
				$pepcounts_mz[$bin_counts-1][1] = $first_half_pep_counts + $second_half_pep_counts;
				$bin_counts ++;
				$first_half_pep_counts = $second_half_pep_counts;
				$second_half_pep_counts = 0;
				$min = $max;
				$max = $min + ($bin_size/2);				
				$i--;
			}
		}
	}
}

# generate the remained bins in case that the user defined range beyonds the m/z value
while ($max <= $mass_range_max + $bin_size/2) {
	$pepcounts_mz[$bin_counts-1][0] = $min;
	$pepcounts_mz[$bin_counts-1][1] = 0;
	$bin_counts++;
	$min = $max;
	$max = $min + ($bin_size/2);
}

print scalar(@mass_data);

######################################################
#output
open (OFILE,'>',$filename_output);
print OFILE "m/z values","\t","peptides number","\t","(bin size = ${bin_size})","\n";
	for (my $i=0;$i<@pepcounts_mz;$i++) {
		print OFILE $pepcounts_mz[$i][0],"\t",$pepcounts_mz[$i][1],"\n";
	}
close OFILE;












