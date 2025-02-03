#!/usr/bin/env perl

use strict;
use warnings "all";
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

my($Help, $Usage, $me);
my(@opt_table, %opt, @args, $tmpdir);

$me = &basename($0);
%opt = (
   'verbose'   => 1,
   'clobber'   => 0,
   'fake'      => 0,
   'text'      => 'use -title',
);



$Help = <<HELP;
| $me builds a jpg verification image for the hippocampus
| using a stereotaxic MRI minc file and a correponding label file
| Problems or comments should be sent to: louis.collins\@mcgill.ca
HELP

$Usage = "Usage: $me [options] minc labels jpg\n" . "where:\n" .
" inputs:\n".
"   minc   - N3 corrected MRI in stx space\n".
"   labels - label file containing HC labels\n".
" output:\n".
"   jpg - the name of the jpg image\n\n";


@opt_table = (
   ["-verbose", "boolean", 0, \$opt{verbose},
      "be verbose" ],
   ["-clobber", "boolean", 0, \$opt{clobber},
      "clobber existing check files" ],
   ["-fake", "boolean", 0, \$opt{fake},
      "do a dry run, (echo cmds only)" ],
   ["-title", "string", 1, \$opt{text},
      "do a dry run, (echo cmds only)" ],
);

# Check arguments
&Getopt::Tabular::SetHelp($Help, $Usage);
&GetOptions (\@opt_table, \@ARGV) || exit 1;


die $Usage if(! ($#ARGV == 2));

# Get the arguments.

my $mri = shift;
my $label = shift;
my $outfile = shift;

die "$me: Couldn't find input mri file: $mri\n\n" if (!-e $mri);
die "$label: Couldn't find input label file: $label\n\n" if (!-e $label);
die "$me: $outfile exists, use -clobber to overwrite\n\n" if (-e $outfile && !$opt{clobber});


# make tmp dir
$tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

# resample mri and labels

my $tmpmri = $tmpdir . "/resampled_mri.mnc";
my $tmplab2 = $tmpdir . "/resampled_lab2.mnc";
my $tmplab = $tmpdir . "/resampled_lab.mnc";

&do_cmd('mincresample', '-clobber', $mri, $tmpmri,
	'-step', 1.0, 1.0, 1.0,
	'-start', -65.000000, -64.000000, -55.000000,
	'-nelements', 131, 98, 84);

&do_cmd('mincresample', '-nearest','-clobber', $label, $tmplab2,
	'-step', 1.0, 1.0, 1.0,
	'-start', -65.000000, -64.000000, -55.000000,
	'-nelements', 131, 98, 84);

&do_cmd('minclookup', '-discrete',
  #'-lut', '11 1; 12 7; 21 3; 22 5',
  '-lut', '11 4; 12 6; 21 8; 22 13',
	$tmplab2,$tmplab);

# normalize MRI to 0..100

my $normmri = $tmpdir . "/normed_mri.mnc";

&do_cmd('mincnorm','-cutoff',0.5,$tmpmri, $normmri);


# make jpg

my $grey     = $tmpdir . "/grey_mri.mnc";
my $labelvol = $tmpdir . "/coloured_labels.mnc";
my $result   = $tmpdir . "/result.mnc";

&do_cmd('minclookup','-clobber','-grey','-range',20,90,$normmri, $grey);
&do_cmd('minclookup','-clobber','-discrete','-lookup_table','/home/bic/louis/lib/luts/labels.map', $tmplab, $labelvol);
&do_cmd('mincmath',  '-clobber','-nocheck_dimensions','-max', $grey, $labelvol, $result);



my @cmd;
my @mont_cmd;
my @ext_cmd;

#z
foreach (26, 30, 34, 38, 42, 44, 48 ) {
    @cmd = ('mincpik', '-scale','1','-transverse','-slice',$_,'-clobber'); #linux:add -clobber
    push(@cmd, @ext_cmd) if @ext_cmd;
    push(@cmd, $result, "$tmpdir/T$_.miff");
    &do_cmd(@cmd);

    push(@mont_cmd, "$tmpdir/T$_.miff");
}

#x
foreach  (28, 31, 34, 37, 40, 43, 47) {
    @cmd = ('mincpik','-scale','1','-sagittal','-slice',$_,'-clobber');
    push(@cmd, @ext_cmd) if @ext_cmd;
    push(@cmd, $result, "$tmpdir/S$_.miff");
    &do_cmd(@cmd);

    push(@mont_cmd, "$tmpdir/S$_.miff");
}

foreach  (102, 99, 96, 93, 90, 87, 84) {
    @cmd = ('mincpik','-scale','1','-sagittal','-slice',$_,'-clobber');
    push(@cmd, @ext_cmd) if @ext_cmd;
    push(@cmd, $result, "$tmpdir/S$_.miff");
    &do_cmd(@cmd);

    push(@mont_cmd, "$tmpdir/S$_.miff");
}

#y
foreach  (55, 50, 45, 40, 35, 30, 25 ) {
    @cmd = ('mincpik','-scale','1','-coronal','-slice',$_,'-clobber');
    push(@cmd, @ext_cmd) if @ext_cmd;
    push(@cmd, $result, "$tmpdir/C$_.miff");
    &do_cmd(@cmd);

    push(@mont_cmd, "$tmpdir/C$_.miff");
}

my $smalltilesize = 200;

&do_cmd('montage',
	'-tile', '7x4',
	'-background', 'grey10',
	'-geometry', $smalltilesize . 'x' . $smalltilesize . '+1+1',
	@mont_cmd,
	"$tmpdir/mont.miff");

    # Add the title
@cmd = ('convert', '-box', 'white');
#'-font', '7x13',
#'-fill', 'white',
#'-draw', "text 2,15 \"@opt{text}\"");

push(@cmd,'-draw', "text 2,15 \"$opt{text}\"") if $opt{text};

#push(@cmd, @more_cmd) if @more_cmd;

&do_cmd(@cmd,"$tmpdir/mont.miff", $outfile);



#  &do_cmd('mincpik','-clobber','-triplanar','-sagittal_offset',20,$result, $outfile);


exit 0;


sub do_cmd {
   print STDOUT "@_\n" if $opt{verbose};
   if(!$opt{fake}){
      system(@_) == 0 or die;
   }
}
