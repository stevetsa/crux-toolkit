#!/net/gs/vol3/software/bin/perl

# convert an input file from ms2 format to mgf format
# ms2 format is described at http://noble.gs.washington.edu/proj/crux/ms2-format.html
# mgf format is described at http://noble.gs.washington.edu/~andrones/internal/proj/ms2denoising/doc/html/PSMAnnotator.html

# If the ms2 file contains 2 Z rows per spectrum, the mgf file will have 2 spectra 

use strict;
use warnings;

if ($#ARGV != 1)
{
    print "Usage: $0 <input_file.ms2> <output_file.mgf>\n";
    exit;
}

my $input = $ARGV[0];
my $output = $ARGV[1];

open (INPUT, $input) || die ("Cannot open $input");
open (OUTPUT, ">$output") || die ("Cannot open $output");
my $line;


# All the data for a spectrum. Write all the data in variables first,
# to be able to write it several times, in case there are several Z lines
my $data = "";  # the data string. 
my $title = "";
my $pepmass = "";
my %mass = ();

my $charge;
 
while ($line = <INPUT>)
{
    # ignore the header lines and the charge (in)dependent lines
    next if ($line =~ /^[HID]/);
    if ($line =~ /^S\s/)  # new scan, means new data string
    {
        foreach $charge (sort(keys(%mass)))
        {
            save_spectrum ($title, $charge, $mass{$charge}, $data);
        }
        
        $data = "";
        $title = "";
        %mass = ();
        # get the scan number and the precursor mass/charge
        my @cols = split(/\s+/,$line);
        # remove the heading zeros
        my $scan = $cols[1];
        $scan =~ s/^0+//;
        $title = "TITLE=S_${scan}_";
    }
    elsif ($line =~ /^Z\s/)   # new supposed charge
    {
        my @cols = split(/\s+/,$line);
        $mass{$cols[1]} = $cols[2]/$cols[1];
    }
    elsif ($line =~ /\d+(\.\d+)*\s\d+(\.\d+)*/)   # data line
    {
        chomp($line);
        #print "<$line>";
        $data .= "$line\n";
    }
    else
    {
        die ("Wrong format of line!!!\n$line");
    }
}

# now make sure I write the last ones
foreach $charge (sort(keys(%mass)))
{
    save_spectrum ($title, $charge, $mass{$charge}, $data);
}

close (OUTPUT);
close (INPUT);


######################
sub save_spectrum
{
    my ($title, $charge, $pepmass, $data) = @_;
    print OUTPUT "BEGIN IONS\n";
    print OUTPUT "$title$charge\n";
    print OUTPUT "SEQ=\n";
    print OUTPUT "CHARGE=+$charge\n";
    print OUTPUT "PEPMASS=$pepmass\n";
    print OUTPUT "$data";   # has \n
    print OUTPUT "END IONS\n\n";
}



