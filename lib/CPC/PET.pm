package CPC::PET;

use 5.006;
use strict;
use warnings;
use Carp qw(carp croak cluck confess);
use List::Util qw(all);
use Scalar::Util qw(looks_like_number reftype);

=head1 NAME

CPC::PET - Calculate potential evapotranspiration (PET)

=head1 VERSION

Version 0.50

=cut

our $VERSION = '0.50';


=head1 SYNOPSIS

This module provides exportable functions to calculate potential evapotranspiration using different methods.

    # Export functions to calculate the PET using the Thornthwaite Equation
    use CPC::PET qw(get_thornthwaite_pet get_thornthwaite_tei);

=head1 EXPORT

=over 4

=item * get_thornthwaite_pet

Calculate the PET using the Thornthwaite Equation

=item * get_thornthwaite_tei

Calculate the TEI parameter used in the Thornthwaite Equation from monthly average temperatures

=back

=head1 FUNCTIONS

=head2 get_pet_thornthwaite

Calculate the potential evapotranspiration using the Thornthwaite Equation. Requires five arguments:

=over 4

=item 1

YDAY - the day of the year (1 - 366). When calculating for a range of days, use the midpoint's date

=item 2

NDAYS - number of days in the period

=item 3

LAT - the latitude(s) of the target location(s) (in degrees)

=item 4

TEMP - temperature(s) (deg C)

=item 5

TEI - the Thornthwaite climatological "temperature efficiency index" for the target location(s) based on monthly average temperatures

=back

The LAT, TEMP, and TEI arguments can be either numeric scalars to calculate the PET for a single location, or array refs of matching dimensions to calculate the PET for multiple locations. A separate function called get_thornthwaite_tei is provided in this package to calculate the TEI from monthly average temperature data.

Usage:

    my $pet_scalar = get_thornthwaite_pet($yday,$ndays,$lat,$temp,$tei);    # For a single location
    my $pet_ref    = get_thornthwaite_pet($yday$ndays,\@lat,\@temp,\@tei);  # For multiple locations (all for the same time period)
    my @pet_array  = @{$pet_ref};

Limitations or known bugs:

=over 4

=item 1

This function currently only supports northern hemisphere locations

=back

=cut

sub get_thornthwaite_pet {
    my $function = "CPC::PET::get_thornthwaite_pet";

    # --- Get args ---

    unless(@_ >= 4) { croak "$function: Five arguments are required"; }
    my $yday        = shift;
    my $ndays       = shift;
    my $lat         = shift;
    my $temperature = shift;
    my $tei         = shift;

    # --- Validate args ---

    unless(
        looks_like_number($yday)        and
        int($yday)  == $yday            and
        $ndays > 0                      and
        $ndays < 367
    ) { croak "$function: YDAY must be an integer between 1 and 366"; }

    unless(
        looks_like_number($ndays)       and
        int($ndays) == $ndays           and
        $ndays > 0
    ) { croak "$function: NDAYS arg must be a positive integer"; }

    my $nvals;

    if(defined reftype($lat)) {  # LAT, TEMP, and TEI are assumed to be array refs of equal dimensions - validate

        if(reftype($lat) eq 'ARRAY' and reftype($temperature) eq 'ARRAY' and reftype($tei) eq 'ARRAY') {
            $nvals = scalar(@{$lat});
            unless(defined reftype($temperature) and defined reftype($tei))         { croak "$function: Invalid argument - expected a ref arg but the arg is not a ref";              }
            unless(reftype($temperature) eq 'ARRAY' and reftype($tei) eq 'ARRAY')   { croak "$function: Invalid argument - expected an ARRAY ref but the arg is not an ARRAY ref";    }
            unless(scalar(@{$temperature}) == $nvals and scalar(@{$tei}) == $nvals) { croak "$function: Invalid argument - array dimension mismatch between LAT, TEMP, and TEI args"; }
        }
        else { croak "$function: Invalid LAT, TEMP, or TEI argument(s) - not a scalar or ARRAY ref";  }

    }
    else {  # LAT, TEMP, and TEI args are assumed to be scalars - convert to array refs with one element each
        $nvals = 1;
        if(defined reftype($temperature) or defined reftype($tei)) { croak "$function: Invalid argument - expected a scalar and found a ref"; }
        my @lat         = ($lat);         $lat         = \@lat;
        my @temperature = ($temperature); $temperature = \@temperature;
        my @tei         = ($tei);         $tei         = \@tei;
    }

    for(my $i=0; $i<$nvals; $i++) {
        unless(looks_like_number($$lat[$i]) and $$lat[$i] >= 0 and $$lat[$i] <= 90) { croak "$function: Invalid LAT value found: ".$$lat[$i]." - must be numeric and between 0 and 90"; }
        unless(looks_like_number($$temperature[$i]))                                { $$temperature[$i] = 'NaN'; }
        unless(looks_like_number($$tei[$i]))                                        { $$tei[$i]         = 'NaN'; }
    }

    # --- Calculate potential evapotranspiration ---

    my $pet = [];

    LOC: for(my $i=0; $i<$nvals; $i++) {

        unless(looks_like_number($$lat[$i]) and $$lat[$i] >= 0 and $$lat[$i] <= 90) {
            carp "$function: Invalid LAT value: ".$$lat[$i]." - setting PET to NaN for location $i";
            $$pet[$i] = 'NaN';
            next LOC;
        }

        unless(looks_like_number($$temperature[$i])) {
            carp "$function: Non-numeric TEMP value: ".$$temperature[$i]." - setting PET to NaN for location $i";
            $$pet[$i] = 'NaN';
            next LOC;
        }
        unless(looks_like_number($$tei[$i]))         {
            carp "$function: Non-numeric TEI value: ".$$tei[$i]." - setting PET to NaN for location $i";
            $$pet[$i] = 'NaN';
            next LOC;
        }

        if($$temperature[$i] <= 0) {

            # --- No PET is possible if the topsoil is frozen ---

            $$pet[$i] = 0;
        }
        elsif($$temperature[$i] > 0 and $$temperature[$i] < 26.5) {

            # --- Calculate the hours of sunlight ---

            my $lat_radians       = (pi/180)*$$lat[$i];
            my $solar_declination = 0.409*sin((2*pi/365)*$yday - 1.39);
            my $sunset_angle      = acos(-tan($lat_radians)*tan($solar_declination));
            my $daylength         = (24/pi)*$sunset_angle;

            # --- Calculate alpha term ---

            my $alpha             = (6.75e-7)*($$tei[$i]**3) - (7.71e-5)*($$tei[$i]**2) + 0.01792*$$tei[$i] + 0.49239;

            # --- Apply the Thornthwaite formula ---

            my $pet_gross         = 16*((10*$$temperature[$i])/($$tei[$i]))**$alpha;
            $$pet[$i]             = ($pet_gross*($daylength/12)*($ndays/30))/25.4;  # Converted to inches
        }
        else {

            # --- Apply the formula in Huang et al. 1995 ---

            $$pet[$i]             = ((-415.85 + 32.25*$$temperature[$i] - 0.43*($$temperature[$i]**2))*($ndays/30))/25.4;  # Converted to inches
        }

    }  # :LOC

    return $pet;
}

=head2 get_thornthwaite_tei

Calculate the climatological "temperature efficiency index" parameter used in the Thornthwaite Equation for potential evapotranspiration. Requires 12 arguments: the monthly average temperatures from January through December (month order does not matter). These can either be 12 scalars to compute the TEI for a single location, or 12 array refs of equal dimensions to compute the TEI for a list of locations. Non-numeric values will result in a TEI value set to NaN.

Usage:

    my $tei_scalar = get_thornthwaite_tei($jan,$feb,$mar,$apr,$may,$jun,$jul,$aug,$sep,$oct,$nov,$dec);             # For a single location
    my $tei_ref    = get_thornthwaite_tei(\@jan,\@feb,\@mar,\@apr,\@may,\@jun,\@jul,\@aug,\@sep,\@oct,\@nov,\@dec); # For multiple locations
    my @tei_array  = @{$tei_ref};

=cut

sub get_thornthwaite_tei {
    my $function = "CPC::PET::get_thornthwaite_tei";

    # --- Get args ---

    my @args;
    unless(@_ >= 12) { croak "$function: 12 arguments are required"; }
    for(my $i=0; $i<12; $i++) { $args[$i] = shift; }

    # --- Validate args ---

    my @monthly_temperatures;
    my $testmon = $args[0]; # Use first arg to determine whether we are working with scalars or array refs
    my $nvals;

    if(defined reftype($testmon)) {  # Args assumed to be all array refs of equal dimensions - validate this

        if(reftype($testmon) eq 'ARRAY') {
            $monthly_temperatures[0] = $testmon;
            $nvals                   = scalar(@{$testmon});

            for(my $mon=1; $mon<12; $mon++) {
                unless(defined reftype($args[$mon]))     { croak "$function: Invalid argument - expected a ref arg but the arg is not a ref"; }
                unless(reftype($args[$mon]) eq 'ARRAY')  { croak "$function: Invalid argument - expected an ARRAY ref but the arg is not an ARRAY ref"; }
                unless(scalar(@{$args[$mon]}) == $nvals) { croak "$function: Invalid argument - array dimension mismatch with the other args"; }
                $monthly_temperatures[$mon] = $args[$mon];
            }

        }
        else { croak "$function: Invalid argument - not a scalar or ARRAY ref"; }

    }
    else {  # Args assumed to be all numeric scalars - validate this
        $nvals = 1;

        for(my $mon=0; $mon<12; $mon++) {
            unless(looks_like_number($args[$mon])) { croak "$function: Invalid argument - expected a numeric scalar but the argument is not numeric"; }
            # Convert the arguments to array refs of dim 1
            $monthly_temperatures[$mon] = \($args[$mon]);
        }

    }

    # --- Calculate the TEI ---

    my $tei = [];  # Array ref

    for(my $i=0; $i<$nvals; $i++) {
        $$tei[$i] = 0;

        for(my $mon=0; $mon<12; $mon++) {
            my $temp   = ${$monthly_temperatures[$mon]}[$i];
            if(looks_like_number($temp) { $temp = $temp > 0 ? $temp : 0; }
            else                        { $temp = 'NaN'; }
            $$tei[$i] += ($temp/5)**1.514;
        }

    }

    # --- Return the TEI ---

    if($nvals == 1) { return $$tei[0]; }
    else            { return $tei;     }
}

=head1 AUTHOR

Adam Allgood, C<< <adam.allgood at noaa.gov> >>

=over 4

Meteorologist

L<Climate Prediction Center (CPC)|https://www.cpc.ncep.noaa.gov>

L<National Weather Service (NWS)|https://www.weather.gov>

=back

=head1 BUGS

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CPC::PET

=head1 REFERENCES

Thornthwaite, C. W., 1948: An approach toward a rational classification of climate. I<Geographical Review>, B<38> 55-94.

=head1 LICENSE

As a work of the United States Government, this project is in the
public domain within the United States.

Additionally, we waive copyright and related rights in the work
worldwide through the CC0 1.0 Universal public domain dedication.

=head2 CC0 1.0 Universal Summary

This is a human-readable summary of the L<Legal Code (read the full text)|https://creativecommons.org/publicdomain/zero/1.0/legalcode>.

=head3 No Copyright

The person who associated a work with this deed has dedicated the work to
the public domain by waiving all of his or her rights to the work worldwide
under copyright law, including all related and neighboring rights, to the
extent allowed by law.

You can copy, modify, distribute and perform the work, even for commercial
purposes, all without asking permission.

=head3 Other Information

In no way are the patent or trademark rights of any person affected by CC0,
nor are the rights that other persons may have in the work or in how the
work is used, such as publicity or privacy rights.

Unless expressly stated otherwise, the person who associated a work with
this deed makes no warranties about the work, and disclaims liability for
all uses of the work, to the fullest extent permitted by applicable law.
When using or citing the work, you should not imply endorsement by the
author or the affirmer.

=cut

1; # End of CPC::PET
