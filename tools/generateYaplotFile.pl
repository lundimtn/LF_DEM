#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [output_interval] [xz_shift]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;
use Getopt::Long;

my $particle_data = $ARGV[0];
my $yap_radius = 1;
my $force_factor = 0.0001;
my $output_interval = 1;
my $xz_shift = 0;
my $axis = 0;
my $reversibility_test = 0;
my $monodisperse = 0;

GetOptions(
'forcefactor=f' => \$force_factor,
'interval=i' => \$output_interval,
'shift=f' => \$xz_shift,
'axis' => \$axis,
'reversibility' => \$reversibility_test,
'monodisperse' => \$monodisperse);

printf "force_factor = $force_factor\n";
printf "output_interval = $output_interval\n";
printf "xz_shift = $xz_shift\n";
printf "axis = $axis\n";
printf "reversibility = $reversibility_test\n";
printf "monodisperse = $monodisperse\n";

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$interaction_data = "int_${name}.dat";
$output = "y_$name.yap";

open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");

&readHeader;
&yaplotColor;

$cnt_interval = 0;
$first = 1;
$first_int = 1;
$checkpoint = 1;
$shear_strain_previous = 0;
$shearrate_positive = 1;

while (1) {
	if ($cnt_interval == 0 ||
		$cnt_interval % $output_interval == 0) {
			$output = 1;
		} else {
			$output = 0;
		}
	&InParticles;
	last unless defined $line;
	if ($shearrate_positive > 0) {
		if ($shear_strain < $shear_strain_previous) {
			$shearrate_positive = -1;
		}
		$checkpoint = 0;
	} else {
		if ($shear_strain > $shear_strain_previous) {
			$shearrate_positive = 1;
			$checkpoint = 1;
		} else {
			$checkpoint = 0;
		}
	}
	$shear_strain_previous = $shear_strain;
	&InInteractions;
	if ($reversibility_test) {
		if ($first || $checkpoint == 1) {
			&keepInitialConfig;
		}
	}
	if ($output == 1) {
		&OutYaplotData;
	}
	$cnt_interval ++;
}

close (OUT);

close (IN_particle);
close (IN_interaction);

##################################################################
sub keepInitialConfig {
	for ($i = 0; $i < $np; $i ++){
		$posx_init[$i] = $posx[$i];
		$posy_init[$i] = $posy[$i];
		$posz_init[$i] = $posz[$i];
		$ang_init[$i] = $ang[$i];
		$radius_init[$i] = $radius[$i];
	}
}

sub readHeader {
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $flwtyp) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $dataunit) = split(/\s+/, $line);
	
	if ($Ly==0) {
		$number_of_header = 8;
	} else {
		$number_of_header = 7;
	}
	for ($i = 0; $i<$number_of_header; $i++) {
		$line = <IN_particle>;
		printf "$line";
	}
	if ($Ly == 0) {
		$number_of_header_int = 20;
	} else {
		$number_of_header_int = 20;
	}
	for ($i = 0; $i<$number_of_header_int; $i++) {
		$line = <IN_interaction>;
	}
	#	$xo = $Lx/2;
	#	$yo = $Ly/2;
	#	$zo = $Lz/2;
	$xo = 0;
	$yo = 0;
	$zo = 0;
}

sub yaplotColor {
	printf OUT "\@0 0 0 0 \n";
	#printf OUT "\@1 50 100 205 \n";
	printf OUT "\@1 25 50 102 \n";
	#printf OUT "\@1 255 255 255 \n";
	printf OUT "\@2 200 200 200 \n";
	printf OUT "\@3 50 150 255 \n";
	printf OUT "\@4 50 200 50 \n";
	printf OUT "\@5 255 100 100 \n";
	printf OUT "\@6 50 200 50 \n";
	printf OUT "\@7 255 255 0 \n";
	printf OUT "\@8 255 255 255\n";
	printf OUT "\@9 150 150 150\n";
	#printf OUT "\@8 224 143 0 \n";
	#printf OUT "\@9 67 163 230 \n";
	#printf OUT "\@8 253 105 6 \n";
	#printf OUT "\@9 109 109 109 \n";
	printf OUT "\@10 250 250 250 \n";
	printf OUT "\@11 240 240 240 \n";
	printf OUT "\@12 230 230 230 \n";
	printf OUT "\@13 220 220 220 \n";
	printf OUT "\@14 210 210 210 \n";
	printf OUT "\@15 200 200 200 \n";
	printf OUT "\@16 190 190 190 \n";
	printf OUT "\@17 180 180 180 \n";
	printf OUT "\@18 170 170 170 \n";
	printf OUT "\@19 160 160 160 \n";
	printf OUT "\@20 150 150 150 \n";
	printf OUT "\@21 140 140 140 \n";
	printf OUT "\@22 130 130 130 \n";
	printf OUT "\@23 120 120 120 \n";
	printf OUT "\@24 110 110 110 \n";
	printf OUT "\@25 100 100 100 \n";
	printf OUT "\@26 90 90 90 \n";
	printf OUT "\@27 80 90 90 \n";
	printf OUT "\@28 70 70 70\n";
	printf OUT "\@29 60 60 60 \n";
	printf OUT "\@30 50 50 50 \n";
	printf OUT "\@31 40 40 40 \n";
	printf OUT "\@32 30 30 30 \n";
	printf OUT "\@33 20 20 20 \n";
	printf OUT "\@34 10 10 10 \n";
	printf OUT "\@35 0 0 0 \n";
}

sub InParticles {
	$radius_max = 0;
	##  Snapshot Header
	$j = 0;
	while (1) {
		$line = <IN_particle>;
		($buf, $val) = split(" : ", $line);
		($buf1) = split(/\s+/, $buf);
		printf "val $buf1 \n";
		if ($buf1 eq '0') {
			last;
		} else {
			$ssHeader[$j++] = $val;
		}
		last unless defined $line;
	}
	$shear_strain = $ssHeader[0];
	$shear_disp = $ssHeader[1];
	$shear_rate = $ssHeader[2];
	$target_stress = $ssHeader[3];
	printf "=== $shear_strain $shear_disp  $shear_rate $target_stress \n";
	for ($i = 0; $i < $np; $i ++){
		if ($i > 0) {
			$line = <IN_particle>;
		}
		if ($output == 1) {
			# LF_DEM version v3.0-53-g9f8a2585-dirty
			# np 500
			# VF 0.78
			# Lx 51.6458
			# Ly 0
			# Lz 51.6458
			# flow_type shear
			# data in r units.
			#1: particle index
			#2: radius
			#3: position x
			#4: position z
			#5-7: velocity (x, y, z)
			#8-10: angular velocity (x, y, z)
			#11: angle
			#				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			#	$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			($ip, $a, $x, $z, $vx, $vz, $vy, $ox, $oz, $oy, $angle) = split(/\s+/, $line);
			#
			$ang[$i] = $angle;
			$radius[$i] = $a;
			if ($xz_shift) {
				$x += $xz_shift;
				$z += $xz_shift;
				if ($x > $Lx/2) {
					$x -= $Lx;
				}
				if ($z > $Lz/2) {
					$z -= $Lz;
				}
			}
			$posx[$i] = $x-$xo;
			$posy[$i] = $y-$yo;
			$posz[$i] = $z-$zo;
			$velx[$i] = $vx;
			$vely[$i] = $vy;
			$velz[$i] = $vz;
			$omegax[$i] = $ox;
			$omegay[$i] = $oy;
			$omegaz[$i] = $oz;
			if ($radius_max < $a) {
				$radius_max = $a;
			}
		}
	}
}

sub InInteractions{
	#$line = <IN_interaction>;
	#($buf, $shear_strain_i, $num_interaction) = split(/\s+/, $line);
	#printf "int $buf $shear_strain_i $num_interaction\n";
	while (1) {
		$line = <IN_interaction>;
		($buf, $val) = split(" : ", $line);
		($buf1) = split(/\s+/, $buf);
		printf "int $buf1 \n";
		if ($buf1 eq '0') {
			last;
		} else {
			$ssHeader[$j++] = $val;
		}
		last unless defined $line;
	}

	# 1, 2: numbers of the interacting particles
	# 3: 1=contact, 0=apart
	# 4, 5, 6: normal vector
	# 7: dimensionless gap = s - 2, s = 2r/(a1+a2)
	# 8: lubrication force
	# 9: Normal part of contact force
	# 10: Tangential part of contact force
	# 11: Colloidal force
	# 12: Viscosity contribution of contact xF
	# 13: N1 contribution of contact xF
	# 14: N2 contribution of contact xF
	$k = 0;
	while (true) {
		if ($k > 0) {
			$line = <IN_interaction>;
		}
		#1: particle 1 label
		#2: particle 2 label
		#3: contact state (0 = no contact, 1 = frictionless contact, 2 = non-sliding frictional, 3 = sliding frictional)
		#4-6: normal vector, oriented from particle 1 to particle 2
		#7: dimensionless gap = s-2, s = 2r/(a1+a2)
		#8: normal part of the lubrication force (positive for compression)
		#9-11: tangential part of the lubrication force
		#12: norm of the normal part of the contact force
		#13-15: tangential part of the contact force
		#16: norm of the normal repulsive force
		#17: Viscosity contribution of contact xF
		
		($i, $j, $contact, $nx, $nz, $ny, #1---6
		$gap, $f_lub_norm, # 7, 8
		$f_lub_tan_x, $f_lub_tan_z, $f_lub_tan_y, # 9, 10, 11
		$fc_norm, # 12
		$fc_tan_x, $fc_tan_z, $fc_tan_y, # 13, 14, 15
		$fr_norm, $s_xF) = split(/\s+/, $line); # 16, 17
	
		if ($i eq '#' || $i eq NU) {
			($buf, $shear_strain) = split(" : ", $line);
			last;
		}
		if (! defined $i) {
			last;
		}
		
		if ($output == 1) {
			$int0[$k] = $i;
			$int1[$k] = $j;
			$contactstate[$k] = $contact;
			if ($contact > 0) {
				#$force[$k] = $fc_norm + $f_lub_norm + $fr_norm;
				$force[$k] = $fc_norm + $f_lub_norm + $fr_norm;
				#$force[$k] = $fc_norm;
			} else {
				$force[$k] = $f_lub_norm + $fr_norm;
				#$force[$k] = 0;
			}
			$F_lub[$k] = $f_lub_norm;
			$Fc_n[$k] = $fc_norm;
			$Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
			$F_lub_t[$k] = sqrt($f_lub_tan_x**2+ $f_lub_tan_z**2 + $f_lub_tan_y**2);
			$S_bf[$k] =  $s_xF;
			$force[$k] += $Fc_t[$k] + $F_lub_t[$k];
			$nrvec_x[$k] = $nx;
			$nrvec_y[$k] = $ny;
			$nrvec_z[$k] = $nz;
			$Gap[$k] = $gap;
			$k++;
		}
	}
	printf "---- $k\n";
	$num_interaction = $k;
}

sub OutYaplotData{
	if ($first == 0) {
		printf OUT "\n";
	} else {
		$first = 0;
	}
	printf OUT "y 1\n";
	printf OUT "@ 8\n";
	## visualize particles
	if ($monodisperse) {
		printf OUT "r $radius[0]\n";
		for ($i = 0; $i < $np; $i++) {
			printf OUT "c $posx[$i] $posy[$i] $posz[$i]  \n";
		}
	} else {
		for ($i = 0; $i < $np; $i++) {
			$rr = $yap_radius*$radius[$i];
			printf OUT "r $rr\n";
			printf OUT "c $posx[$i] $posy[$i] $posz[$i]  \n";
		}
	}
	
	## visualize contact network
	printf OUT "y 2\n";
	printf OUT "r 0.2\n";
	printf OUT "@ 2\n"; # static
	for ($k = 0; $k < $num_interaction; $k ++) {
		if ($contactstate[$k] >= 2) {
			&OutString2($int0[$k], $int1[$k]);
		}
	}
	## visualize force chain network
	printf OUT "y 4\n";
	printf OUT "@ 7\n";
	printf "num $num_interaction \n";
	for ($k = 0; $k < $num_interaction; $k ++) {
		#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
		#if (1 || $force[$k] >= 0) {
		#	&OutString_width($int0[$k], $int1[$k], $force_factor*$force[$k], 0.01);
		#}
		if ($force[$k] > 0) {
			&OutString_width($int0[$k], $int1[$k], $force_factor*$force[$k], 0.01);
		} else {
			#	&OutString_width($int0[$k], $int1[$k], -$force_factor*$force[$k], 0.01);
		}
	}
	#	printf OUT "y 3\n";
	#	printf OUT "@ 5\n";
	#	for ($k = 0; $k < $num_interaction; $k ++) {
	#		#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
	#		if ($force[$k] < 0) {
	#			&OutString_width($int0[$k], $int1[$k], -$force_factor*$force[$k], 0.02);
	#		}
	#	}
	
	## visualize rotation in 2D
	if ($Ly == 0) {
		if (0) {
			printf OUT "y 6\n";
			printf OUT "@ 1\n";
			for ($i = 0; $i < $np; $i++) {
				OutCross($i);
			}
		}
	}
	if ($reversibility_test) {
		printf OUT "y 5\n";
		printf OUT "@ 7\n";
		$r = $radius_init[0];
		printf OUT "r $r\n";
		for ($i = 0; $i < $np; $i++) {
			if ($i >= 1 && $radius_init[$i] != $radius_init[$i-1]) {
				$r = $yap_radius*$radius_init[$i];
				printf OUT "r $r\n";
			}
			if ($Ly == 0) {
				printf OUT "c $posx_init[$i] 0.01 $posz_init[$i] \n";
			} else {
				printf OUT "c $posx_init[$i] $posy_init[$i] $posz_init[$i] \n";
			}
		}
	}
	#	&OutBoundaryBox;
}

sub OutBoundaryBox {
	$x0 = -$Lx/2;
	$x1 = -$Lx/2 + $shear_disp / 2;
	$z1 = $Lz/2;
	$x2 = $Lx/2;
	$z2 = 0;
	$x3 = $Lx/2 - $shear_disp / 2;
	$z3 = -$Lz/2;
	
	printf OUT "y 7\n";
	printf OUT "@ 6\n";
	
	if ($Ly == 0) {
		$lx2 = $Lx/2;
		$ly2 = $Ly/2;
		$lz2 = $Lz/2;
		#printf OUT "p 4 -$lx2 $yb $lz2 $lx2 $yb $lz2 $lx2 $yb -$lz2 -$lx2 $yb -$lz2\n";
		printf OUT "l -$lx2 0 $lz2    $lx2 0 $lz2\n";
		printf OUT "l -$lx2 0 -$lz2   $lx2 0 -$lz2\n";
		printf OUT "l -$lx2 0 -$lz2  -$lx2 0 $lz2\n";
		printf OUT "l $lx2 0 $lz2   $lx2 0 -$lz2\n";
	} else {
		$lx2 = $Lx/2;
		$ly2 = $Ly/2;
		$lz2 = $Lz/2;
		printf OUT "l -$lx2 -$ly2 $lz2    $lx2 -$ly2 $lz2\n";
		printf OUT "l -$lx2 -$ly2 -$lz2   $lx2 -$ly2 -$lz2\n";
		printf OUT "l -$lx2 -$ly2 -$lz2  -$lx2 -$ly2 $lz2\n";
		printf OUT "l  $lx2 -$ly2 $lz2    $lx2 -$ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 $lz2    $lx2 $ly2 $lz2\n";
		printf OUT "l -$lx2 $ly2 -$lz2   $lx2 $ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 -$lz2  -$lx2 $ly2 $lz2\n";
		printf OUT "l  $lx2 $ly2 $lz2    $lx2 $ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 $lz2    -$lx2 -$ly2 $lz2\n";
		printf OUT "l $lx2 $ly2 $lz2    $lx2 -$ly2 $lz2\n";
		printf OUT "l $lx2 $ly2 -$lz2    $lx2 -$ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 -$lz2    -$lx2 -$ly2 -$lz2\n";
	}
	
	if ($axis) {
		printf OUT "l -$lx2 0 0 $lx2 0 0\n";
		printf OUT "l 0 0 -$lz2 0 0 $lz2\n";
	}
	
}

sub OutString_width {
	($i, $j, $w, $delta) = @_;
	$xi = $posx[$i];
	$yi = $posy[$i] - $delta;
	$zi = $posz[$i];
	$xj = $posx[$j];
	$yj = $posy[$j] - $delta;
	$zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2 ;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1) {
		printf OUT "r $w\n";
		printf OUT "s $xi $yi $zi $xj $yj $zj \n";
	}
}

sub OutString2{
	$offset2 = -0.012;
	($i, $j) = @_;
	$xi = $posx[$i];
	$yi = $posy[$i]+$offset2;
	$zi = $posz[$i];
	$xj = $posx[$j];
	$yj = $posy[$j]+$offset2;
	$zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1) {
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
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
		&&  abs($zi-$zj) < $radius_max*5) {
			printf OUT "s $xi $yi $zi $xj $yj $zj\n";
		}
}

sub OutCross {
	$offset1 = -0.01;
	($i) = @_;
	$a = $radius[$i];
	$xi = $posx[$i];
	$yi = $posy[$i]+$offset1;
	$zi = $posz[$i];
	$angle = $ang[$i];
	$ux =  $a*cos($angle);
	$uz = -$a*sin($angle);
	$xa = $xi-$ux;
	$ya = $yi+$offset1;
	$za = $zi-$uz;
	$xb = $xi+$ux;
	$yb = $yi+$offset1;
	$zb = $zi+$uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
	$ux = $a*sin($angle);
	$uz = $a*cos($angle);
	$xa = $xi-$ux;
	$ya = $yi+$offset1;
	$za = $zi-$uz;
	$xb = $xi+$ux;
	$yb = $yi+$offset1;
	$zb = $zi+$uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
}
