#! /usr/bin/perl
use Math::Trig;
$particle_data = $ARGV[0];
$force_factor = 0.01;
if ($#ARGV >= 1){
	$force_factor = $ARGV[1];
}
$y_section = 0;
printf "$#ARGV";
if ($#ARGV == 2){
	$y_section = $ARGV[2];
	printf "section $y_section\n";
	exit;
}
# Read the header
#
# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$interaction_data = "int_${name}.dat";
printf "$interaction_data\n";
$output = "y_$name.yap";
open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
$first=1;
while (1){
	&InParticles;
	&InInteractions;
	last unless defined $line;
	&OutYaplotData;
	printf "$shear_rate\n";
}
close (OUT);

sub readHeader{
	$line = <IN_particle>; ($buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $Lz) = split(/\s+/, $line);
	printf "$np, $VF, $Lx, $Ly, $Lz\n";
}

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
    ($buf, $shear_rate) = split(/\s+/, $line);
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_particle> ;
        ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz ) = split(/\s+/, $line);
		$radius[$i] = $a;
        $posx[$i] = $x;
        $posy[$i] = $y;
        $posz[$i] = $z;
		if ($radius_max < $a){
			$radius_max = $a
			
		}
    }
}

sub InInteractions {
	$line = <IN_interaction>;
	($buf, $shear_rate, $num_interaction) = split(/\s+/, $line);
	if ($buf != "#"){
		exit;
	}
	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction> ;
		($i, $j, $f_lub, $fc_n, $fc_t, $gap, $fric_st) = split(/\s+/, $line);
		$int0[$k] = $i;
		$int1[$k] = $j;
		$F_lub[$k] = $f_lub;
		$Fc_n[$k] = $fc_n;
		$Ft_t[$k] = $fc_t;
		$Gap[$k] = $gap;
	}
}

sub OutYaplotData{
	if ($first == 0){
		printf OUT "\n";
	} else {
		$first = 0;
	}
	
	printf OUT "y 1\n";
    printf OUT "@ 2\n";
	$r = $radius[0];
	printf OUT "r $r\n";
    for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $radius[$i];
			printf OUT "r $r\n";
		}
		if ($y_section == 0 ||
			abs($posy[$i]) < $y_section ){
				printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
			}
    }
	
    printf OUT "y 2\n";
    printf OUT "@ 6\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        if ($Fc_n[$k] > 0){
            $string_with = $force_factor*($Fc_n[$k]);
			if ( $y_section == 0
				|| abs($yi) < $y_section
				|| abs($yj) < $y_section){
					printf OUT "r ${string_with}\n";
					OutString($int0[$k],  $int1[$k]);
				}
        }
    }
    printf OUT "y 3\n";
    printf OUT "@ 3\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        $force = $F_lub[$k] + $Fc_n[$k];
        if ($force < 0){
            $string_with = -${force_factor}*${force};
			if ( $y_section == 0
				|| abs($yi) < $y_section
				|| abs($yj) < $y_section){
					
					printf OUT "r ${string_with}\n";
					&OutString($int0[$k],  $int1[$k]);
				
				}
        }
    }
    printf OUT "@ 4\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        $force = $F_lub[$k] + $Fc_n[$k];
        if ($force > 0){
            $string_with = ${force_factor}*${force};
			if ( $y_section == 0
				|| abs($yi) < $y_section
				|| abs($yj) < $y_section){
					printf OUT "r ${string_with}\n";
					&OutString($int0[$k],  $int1[$k]);
				}
			
        }
    }
   
}

sub OutString {
    ($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i];
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j];
    $zj = $posz[$j];
	
	if (abs($xi-$xj) < $radius_max*5
		&&  abs($yi-$yj) < $radius_max*5
		&&  abs($zi-$zj) < $radius_max*5){
			printf OUT "s $xi $yi $zi $xj $yj $zj\n";
		}
}

