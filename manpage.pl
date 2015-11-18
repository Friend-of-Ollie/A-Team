#!/usr/bin/perl -w
use strict;
use Getopt::Long;
sub manpage {
	print 
	"\n","Task 4: the ion statistic calculator","\n\n",
	"This program is working on the relationship between m/z value and peptide counts.","\n",
	"The program take in the mass file which should contains protein name, peptide counts, mass-to-charge value, ","\n",
	"the charge of ion,the number of cleavages and sequence. Finally, the program will output a file which contains","\n",
	"two columns represents m/z values and relevant peptide counts. The output file can be easily used for plotting","\n",
	"a histogram in R or any other tools you would like to use. ", "\n\n",
	"There are three main analysis in this program:","\n\n",
	"\t","1) ion pattern match","\n",
	"\t","Program can identify sequences with specified amino acid in the start of the sequence or sequence contains","\n",
	"\t","one or several specified amino acids. User should use shortcut name of amino acids in the command line by","\n",
	"\t","-ion-start and -ion-with options. If -ion-start and -ion-with option are not specified, the program won't","\n",
	"\t","do this ion-pattern-match part.","\n\n",
	"\t","2) mass accuracy modelling","\n",
	"\t","Users should indicate the mass error value if the mass error is known. -error option is used to define the","\n",
	"\t","mass error value. The program can accecpt two units of mass error: Da or ppm. When user defines the mass error,","\n",
	"\t","the units should be indicated too or the program will take it as an error.","\n",
	"\t","When the unit is 'Da', the program will use the following formula to modify the mass-to-charge value:","\n",
	"\t\t","[ exact mass = accurate mass + mass error ].","\n",
	"\t","When the unit is 'ppm', program will using the following formula:","\n",
	"\t\t",,"[ mass error = (exact mass - accurate mass)* 10^6 / exact mass ]","\n",
	"\t","The modified exact mass will still in unit Da.","\n\n",
	"\t","3) calculate peptide counts in bins by sliding window","\n",
	"\t","Users can define the mass range by -mr-max and -mr-min and define the bin size by -bin. The program will use","\n",
	"\t","a sliding window which is as big as the defined bin size, and go over all the data in the mass range. Every","\n",
	"\t","bin represents the center position of the bin. For example, if the bin size is 50, the gap between two center","\n",
	"\t","position will be 25 which mean there is an half-bin-size overlap between every neighbouring bins. The program","\n",
	"\t","will count up all peptides number in every bin and output them.","\n",
	"\t","Options not in mandatory parameter can be remained undefined.","\n",
	"\n\n",
	"SYNOPSIS:","\n",
	"\t","perl P_4.pl [ -i input_file ] [ -o output_file ]","\n",
	"\t","perl P_4.pl [ -i input_file ] [ -o output_file ] [ -mr-min number ] [ -mr-max number ] [ -bin number]","\n",
	"\t\t","    [ -error number[ppm/Da] ] [ -ion-start amino_acids ] [ -ion-with amino_acids ]","\n",
	"\t","perl P_4.pl [ -man ] [ -h ] [ -help ]","\n",
	"\t","Example:","\n",
	"\t","perl P_4.pl -i data.txt -o output.txt -mr-min 1000 -mr-max 1550  -bin 50 -error 1.2Da -ion-with P M -ion-start M","\n\n",
	"MANUAL:","\n",
	"\t","-man","\t\t\t","Show complete manual page","\n",
	"\t","-h/help","\t\t\t","Show brief usage","\n\n",
	"MANDATORY PARAMETERS:","\n",
	"\t","-i","\t\t\t","Input file","\n",
	"\t","-o","\t\t\t","Output file","\n\n",
	"OPTIONS:","\n",
	"\t","-mr-max","\t\t\t","Upper limit of mass range (Default: the minimum m/z value of the input data)","\n",
	"\t","-mr-min","\t\t\t","Lower limit of mass range (Default: the maximum m/z value of the input data)","\n",
	"\t","-bin","\t\t\t","Bin size (Default: bin size = 50)","\n",
	"\t","-error","\t\t\t","Give value of mass error (ppm/Da) ex.0.5Da/0.5ppm (Default: mass error = 0 Da)","\n",
	"\t","-ion-start","\t\t","Only calculate sequences start with specific amino acids","\n",
	"\t","-ion-with","\t\t","Only calculate sequences include specific amino acids. Accepts single or several amino acids.","\n\n";
		
}

1;


