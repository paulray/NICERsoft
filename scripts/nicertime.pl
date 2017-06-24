#!/usr/bin/perl -s
#
# nicertime - NICER time conversion
#
# This task performs time conversion between various
# useful time systems in use by the NICER project.
#
# All time inputs and outputs are in UTC time.  This task converts
# between the standard ISO calendar notation (YYYY-MM-DDThh:mm:ss),
# day of year notation (YYYY:DDD:hh:mm:ss) and plain NICER seconds
# since epoch (MET).
# 
# Type nicertime alone on the command line for more information.
#

use Time::Local;
use POSIX;

# YYYY:DDD:hh:mm:ss  where YYYY can be two-digit year; seconds are optional
my $re_YYYY_DDD_hh_mm_ss   = '(\d\d\d\d|\d\d)[-:/](\d\d\d)[-T:\/](\d\d):(\d\d)(:(\d\d))?';
# YYYY-MM-ddT:hh:mm:ss
my $re_YYYY_MM_DD_hh_mm_ss = '(\d\d\d\d|\d\d)[-:/](\d\d)[-:/](\d\d)T(\d\d):(\d\d)(:(\d\d))?';
# Pure decimal MET
my $re_MET = '\d+(\.\d*)?';

my $nicer_time = nicer_set_epoch(2014);
nicer_set_epoch($epoch) if ($epoch);

if (scalar(@ARGV) == 0 && $stdin == 0 && $now == 0) {
  print "USAGE: nicertime [-iso] [-doy] [-met] [-stdin] [time value]\n".
        "  [time value] can be any of the following forms:\n".
	"    YYYY-MM-DDThh:mm:ss - ISO-standard time\n".
	"      YY-MM-DDThh:mm:ss - ISO-standard time with two digit year\n".
	"    YYYY:DDD:hh:mm:ss   - Day-of-year notation\n".
	"      YY:DDD:hh:mm:ss   - Day-of-year notation with two digit year\n".
	"    numerical value     - NICER mission elapsed time in seconds\n".
	"    -now                - use the current time\n".
	"  NOTES: all times are in UTC\n".
	"         specification of seconds is optional\n".
	"         separator character between years days and months can be\n".
	"         a colon(:), dash(-) or slash(/)\n".
	"  -iso - output the time as ISO standard YYYY-MM-DDThh:mm:ss\n".
	"  -doy - output the time as day-of-year YYYY-MM-DDThh:mm:ss\n".
	"  -met - output the time as Mission elapsed time in seconds\n".
	"  If none of -iso -doy and -met are specified then all forms\n".
	"  are printed in a more verbose form.\n".
        "  -stdin - if set, then instead of using a value from the command line\n".
	"     read standard input line by line and convert values\n";
  exit 0;
}

if ($stdin) {
  while ($line=<STDIN>) {
    chomp($line);
    if ($line =~ /^(\s*)(?<val>$re_YYYY_DDD_hh_mm_ss|$re_YYYY_MM_DD_hh_mm_ss|$re_MET)(?<rest>.*)$/) {
      my $val = $+{val};
      my $metval = parse_timespec($val);
      do_output($metval,0);
      print $+{rest}."\n";
    } else {
      print $line."\n";
    }
  }
} elsif ($now) {
  my $metval = parse_timespec("__NOW__");
  do_output($metval,1);
} else {
  my $metval = parse_timespec($ARGV[0]);
  do_output($metval,1);
}

sub do_output ($$) {
  my ($metval,$cr) = (@_);
  if ($cr) { $cr = "\n"; } else { $cr = ""; }
  if    ($iso) { print iso_output($metval).$cr; }
  elsif ($doy) { print doy_output($metval).$cr; }
  elsif ($met) { print met_output($metval).$cr; }
  else         { detail_output($metval); }
}

sub detail_output ($) {
  my ($met) = (@_);
  print "ISO = ".iso_output($met)."\n";
  print "DOY = ".doy_output($met)."\n";
  print "MET = ".met_output($met)."\n";
}
sub iso_output ($) {
  my ($met) = (@_);
  my ($ss,$mm,$hh,$day,$mon,$year,$wday,$yday) =
    gmtime($met - $nicer_time->{gmtime_offset});
  $year += 1900;
  $mon += 1;
  return sprintf("%04d-%02d-%02dT%02d:%02d:%02d",
		 $year, $mon, $day, $hh, $mm, floor($ss));
}
sub doy_output ($) {
  my ($met) = (@_);
  my ($ss,$mm,$hh,$day,$mon,$year,$wday,$yday) =
    gmtime($met - $nicer_time->{gmtime_offset});
  $year += 1900;
  $yday += 1;
  return sprintf("%04d:%03d:%02d:%02d:%02d",
		 $year, $yday, $hh, $mm, floor($ss));
}
sub met_output ($) {
  my ($met) = (@_);
  return $met;
}

sub parse_gti ($$$$) {
  my ($gtistr,$gref,$tminr,$tmaxr) = (@_);
  my (@gtis) = split(/,/,$gtistr);
  my (@gti);
  my ($tmin, $tmax);

  foreach my $gti (@gtis) {
    if ($gti =~ m/^(?<tstart>$re_YYYY_DDD_hh_mm_ss)-(?<tstop>$re_YYYY_DDD_hh_mm_ss)/ ||
	$gti =~ m/^(?<tstart>$re_YYYY_MM_DD_hh_mm_ss)-(?<tstop>$re_YYYY_MM_DD_hh_mm_ss)/ ||
	$gti =~ m/^(?<tstart>$re_MET)-(?<tstop>$re_MET)/) {
      my ($tstartv,$tstopv);
      $tstartv = parse_timespec($+{"tstart"});
      $tstopv  = parse_timespec($+{"tstop"});
      push @gti, $tstartv."-".$tstopv;
      push @$gref, [$tstartv, $tstopv];
      $tmin = $tstartv if (!defined($tmin) || $tstartv<$tmin);
      $tmax = $tstopv  if (!defined($tmax) || $tstopv >$tmax);
    } else {
      die "ERROR: unknown GTI '$gti'";
    }
  }

  $_[2] = $tmin;
  $_[3] = $tmax;
  return join(",",@gti);
}

sub parse_timespec ($) {
  my ($timespec) = (@_);
  my ($t);
  if ($timespec =~ m/^$re_YYYY_DDD_hh_mm_ss/) { # YYYY:DDD:hh:mm:ss
    my ($year,$doy,$hh,$mm,$ss);
    $year = $1; $doy = $2; $hh = $3; $mm = $4; $ss = $6+0;
    $year += 2000 if ($year < 50);
    $year += 1900 if ($year > 50 && $year < 100);
    $t = timegm($ss,$mm,$hh,1,0,$year-1900) + ($doy-1)*86400;
    $t += $nicer_time->{gmtime_offset};
  } elsif ($timespec =~ m/^$re_YYYY_MM_DD_hh_mm_ss/) { # YYYY-MM-ddT:hh:mm:ss
    my ($year,$mon,$day,$hh,$mm,$ss);
    $year = $1; $mon = $2; $day = $3; $hh = $4; $mm = $5; $ss = $7+0;
    $year += 2000 if ($year < 50);
    $year += 1900 if ($year > 50 && $year < 100);
    $t = timegm($ss,$mm,$hh,$day,$mon-1,$year-1900);
    $t += $nicer_time->{gmtime_offset};
  } elsif ($timespec =~ m/^__NOW__$/) {
    $t = timegm(gmtime());
    $t += $nicer_time->{gmtime_offset};
  } elsif ($timespec =~ m/^$re_MET/) {  # Pure decimal MET number
    $t = $timespec + 0.0;
  } else {
    die "ERROR: unknown time specification $timespec";
  }
  return $t;
}

sub nicer_set_epoch ($) {
  my ($epoch) = (@_);
  my %nicer_time;

  my $MJD_1970 = 40587;
  my $MJD_GPS =  44244;
  my $MJD_2014 = 56658;
  my $MJD_2017 = 57754;
  my $LEAP_1970 = 8;
  my $LEAP_GPS = 19;
  my $LEAP_2014 = 35;
  my $LEAP_2017 = 37;

  if ($epoch == 2017) {
    $nicer_time{leapsec} = $LEAP_2017;
    $nicer_time{mjdrefi} = $MJD_2017;
    $nicer_time{epoch_str} = "2017-01-01T00:00:00";
  } elsif ($epoch == 2014) {
    $nicer_time{leapsec} = $LEAP_2014;
    $nicer_time{mjdrefi} = $MJD_2014;
    $nicer_time{epoch_str} = "2014-01-01T00:00:00";
  } else {
    die "ERROR: unknown epoch $epoch";
  }

  # Number of leap seconds elapsed since the NICER epoch
  # WARNING: this needs to be updated if there are more leap seconds!
  $nicer_time{nicer_leapsec} = $LEAP_2017 - $nicer_time{leapsec};

  $nicer_time{gps_offset_day} = $nicer_time{mjdrefi} - $MJD_GPS;
  $nicer_time{"1970_offset_day"} = $nicer_time{mjdrefi} - $MJD_1970;
  $nicer_time{gps_offset_leapsec} = $nicer_time{leapsec} - $LEAP_GPS + $nicer_time{nicer_leapsec};
  $nicer_time{"1970_offset_leapsec"} = $nicer_time{leapsec} - $LEAP_1970 + $nicer_time{nicer_leapsec};

  $nicer_time{gps_offset} = (-$nicer_time{gps_offset_day}*86400 - $nicer_time{gps_offset_leapsec});
  $nicer_time{"1970_offset"} = (-$nicer_time{"1970_offset_day"}*86400 - $nicer_time{"1970_offset_leapsec"});
  $nicer_time{"gmtime_offset"} = (-$nicer_time{"1970_offset_day"}*86400 + $nicer_time{nicer_leapsec});

  $nicer_time{mjdreff} = (($nicer_time{leapsec} + 32.184)/86400.0);

  return \%nicer_time;
}

