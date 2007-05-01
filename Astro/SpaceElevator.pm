package Astro::SpaceElevator;
use utf8; # because I'm crazy

# This module handles all of the details of the ECI coordinate system,
# which is a cartesian system with the origin at the earth's
# center. Also handles conversions from latitude/longitude into ECI
# (based on a reference geoid, not a sphere) and finds the location of
# the sun.
use Astro::Coord::ECI;
use Astro::Coord::ECI::Sun;
use Astro::Coord::ECI::Utils qw{PI rad2deg deg2rad};

# this module lets me do vector and matrix math at a high level, which
# makes the code easier to read.
use Math::MatrixReal;

=head1 NAME

Astro::SpaceElevator - Model a Space Elevator

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

=head1 SYNOPSIS

    use Astro::SpaceElevator;

    my $elevator = Astro::SpaceElevator->new(0, 120, 100_000, time());
    print "The elevator leaves the Earth's shadow at " . ($elevator->shadows())[1] . "km above the base.\n";

=head1 METHODS

=head2 new

=over

Creates a new elevator object. Takes four arguments: the latitude, longitude and height (in km) of the elevator, and a time in seconds since the epoch, in the GMT timezone.

=back

=cut

sub new
{
    my ($class, $lat, $lon, $height, $time) = @_;
    $lat = deg2rad $lat;
    $lon = deg2rad $lon;

    my $self = bless({lat    => $lat,
                      lon    => $lon,
                      height => $height,
                     },
                     $class);
    $self->time($time);

    return $self;
}

=head2 time

=over

Gets the time associated with the model. If you supply an argument, it uses that as a new time to update all of the time-dependant aspects of the model.

=back

=cut

sub time
{
    my ($self, $time) = @_;

    if ($time)
    {
        $self->{time} = $time;
        $self->{sun}  = Astro::Coord::ECI::Sun->universal($time);
        $self->{base} = _geodetic($self->{lat}, $self->{lon}, 0, $time);
    }

    return $self->{time};
}

sub _geodetic
{
    my ($lat, $lon, $elev, $time) = @_;
    return Astro::Coord::ECI->geodetic($lat, $lon, $elev)->universal($time);
}

sub _eci
{
    if (ref $_[0])
    {
        my ($vector, $time) = @_;
        return Astro::Coord::ECI->eci($vector->element(1,1),
                                      $vector->element(2,1),
                                      $vector->element(3,1))->universal($time);
    }
    
    my ($x, $y, $z, $time) = @_;
    return Astro::Coord::ECI->eci($x, $y, $z)->universal($time);
}

sub _vector
{
    my ($x, $y, $z) = @_;
    return Math::MatrixReal->new_from_cols([[$x, $y, $z]]);
}

sub _between
{
    my ($a, $x, $b) = @_;
    return $a if ($x < $a);
    return $x if ($a <= $x and $x <= $b);
    return $b if ($b < $x);
}

my $I = Math::MatrixReal->new_diag([1, 1, 1]);

=head2 shadows

=over

Returns a two element list containing the height in kilometers of the
Earth's umbral and penumbral shadows on the elevator.

=back

=cut

sub shadows
{
    my $self = shift;

    my $earth_diameter = 6372.797;
    my $earth_radius = $earth_diameter / 2;

    my $time = $self->{time};
    my $sun = $self->{sun};
    my $base = $self->{base};
    my $height = $self->{height};

    my $direction = _vector(_geodetic($self->{lat}, $self->{lon}, $self->{height}, $time)->eci);
    $direction *= 1 / $direction->length(); # the matrix module has a slight problem with straight division.

    # first we check to see if the sun has risen over the base station.
    my ($azimuth, $elevation, $range) = $base->azel($sun, 1);
    if ($elevation <= 0)
    {
        # Figure out the dimensions of the umbra, which is a function of
        # the Earth's distance from the sun. The penumbra is congruent to
        # the umbra, but reflected around the plane of the terminator
        my $umbra_temp = _vector($sun->eci);
        my $umbraV = -$umbra_temp * ($earth_radius / ($sun->get('diameter') / 2));
        my $umbraA = $umbra_temp * (1 / $umbra_temp->length());
        my $umbraΘ = PI/2 + _eci($umbraV, $time)->dip();

        my $baseV = _vector($base->eci);

        my @umbra = _intersect($umbraV, $umbraA, $umbraΘ, $baseV, $direction);
        my @penumbra = _intersect(-$umbraV, -$umbraA, $umbraΘ, $baseV, $direction);

        return (_between(0, $umbra[-1], $height),
                _between(0, $penumbra[-1], $height));
    }
    else
    {
        return (0, 0);
    }
}

sub _intersect
{
    # I've taken this algorithm for intersecting lines and cones from
    # http://www.geometrictools.com/Documentation/IntersectionLineCone.pdf,
    # though I'm not rejecting intersections with the anti-cone, as
    # they indicate a transition into an annular eclipse.

    # the cone is defined using a Vertex, an Axis and the angle (Θ)
    # between the axis and the edge of the cone. The axis is a
    # normalized vector that indicates the direction the cone is
    # facing.

    my ($V, $A, $Θ, $P, $D) = @_;
    my $Δ = $P - $V;
    my $M = ($A * ~$A) - ((cos $Θ)**2 * $I);

    my $c0 = (~$Δ * $M * $Δ)->element(1, 1);
    my $c1 = (~$D * $M * $Δ)->element(1, 1);
    my $c2 = (~$D * $M * $D)->element(1, 1);

    if ($c2 != 0)
    {
        my $δ = $c1**2 - $c0*$c2;
        if ($δ > 0)
        {
            return ('segment',
                    (-$c1 + sqrt($δ)) / $c2,
                    (-$c1 - sqrt($δ)) / $c2);
        }
        elsif ($δ = 0)
        {
            return ('ray',
                    -$c1 / $c2);
        }
    }

    # TODO: there are a few other cases that need to be included, but
    # as I've modelled the elevator here, they can't occur. If a more
    # detailed geometric model of the elevator is used to account for
    # oscillation or other effects then they need to be filled in.

    return ('none', 0);
}

=head1 AUTHOR

Daniel Brooks, C<< <db48x at yahoo.com> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-astro-spaceelevator at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Astro-SpaceElevator>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 COPYRIGHT

Copyright © 2007 by Daniel Brooks. All rights reserved.

=head1 LICENSE

This program is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut

1;
