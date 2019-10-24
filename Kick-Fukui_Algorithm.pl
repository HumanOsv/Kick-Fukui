#!/usr/bin/env perl


#
# Write code: 
#            Diego Inostrosa & Osvaldo Yañez Osses
#            contact: osvyanezosses@gmail.com


# Modules
use strict;
use Benchmark; #entrega cuando demora el proceso, cuanto de CPU utilizó, etc
use Cwd qw(abs_path);
use File::Basename;
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
# install: sudo cpan Parallel::ForkManager
use Parallel::ForkManager;
# install: sudo cpan Math::Matrix
use Math::Matrix;

#
#
# Number of process 
my $ncpus = 100;
# Weight spherical restriction algorithm 
my $weight_spherical = 2.0;
#
#
###########################
# Configure file variables
my $num_atoms_global;
my $num_basin_global;
my $Num_of_geometries;
my $XYZSize;
my $Box_x;
my $Box_y;
my $Box_z;
my $decision;
my $sort_option;
my $duplicate_mol;
my $output_file;
my $Final_File;
#
my @ValueJ  = ();
#
my $Tam1 = 0;
my $Tam2 = 0;
#
my %HashCoordsAttra = ();
my %HashDist        = ();
my %HashMult        = ();

# Default is 1 unless factor is given.
my $Box_Multiplication_Factor = 1;
my $BoxAutomatic              = 1;
#
my %Atomic_radii = ( 'H'  => '1.00', 'Li' => '1.45', 'Be' => '1.05', 'B'  => '0.85',
					 'C'  => '0.70', 'N'  => '0.65', 'O'  => '1.35', 'F'  => '0.50', 
					 'Na' => '1.80', 'Mg' => '1.50', 'Al' => '1.25', 'Si' => '1.10',
					 'P'  => '1.00', 'S'  => '1.00', 'Cl' => '1.00', 'K'  => '2.20',
					 'Ca' => '1.80', 'Sc' => '1.60', 'Ti' => '1.40', 'V'  => '1.35',
					 'Cr' => '1.40', 'Mn' => '1.40', 'Fe' => '1.40', 'Co' => '1.35',
					 'Ni' => '1.35', 'Cu' => '1.35', 'Zn' => '1.35', 'Ga' => '1.30',
					 'Ge' => '1.25', 'As' => '1.15', 'Se' => '1.15', 'Br' => '1.15',
					 'Rb' => '2.35', 'Sr' => '2.00', 'Y'  => '1.80', 'Zr' => '1.55',
					 'Nb' => '1.45', 'Mo' => '1.45', 'Tc' => '1.35', 'Ru' => '1.30',
					 'Rh' => '1.35', 'Pd' => '1.40', 'Ag' => '1.60', 'Cd' => '1.55',
					 'In' => '1.55', 'Sn' => '1.45', 'Sb' => '1.45', 'Te' => '1.40',
					 'I'  => '1.40', 'Cs' => '2.60', 'Ba' => '2.15', 'La' => '1.95',
					 'Ce' => '1.85', 'Pr' => '1.85', 'Nd' => '1.85', 'Pm' => '1.85',
					 'Sm' => '1.85', 'Eu' => '1.85', 'Gd' => '1.80', 'Tb' => '1.75',
					 'Dy' => '1.75', 'Ho' => '1.75', 'Er' => '1.75', 'Tu' => '1.75',
					 'Yb' => '1.75', 'Lu' => '1.75', 'Hf' => '1.55', 'Ta' => '1.45',
					 'W'  => '1.35', 'Re' => '1.35', 'Os' => '1.30', 'Ir' => '1.35',
					 'Pt' => '1.35', 'Au' => '1.35', 'Hg' => '1.50', 'Tl' => '1.90',
					 'Bi' => '1.60', 'Po' => '1.90', 'Ra' => '2.15', 'Ac' => '1.95',
					 'Th' => '1.80', 'Pa' => '1.80', 'U'  => '1.75', 'Np' => '1.75',
					 'Pu' => '1.75', 'Am' => '1.75', default=> '0.8' );

###################################
# Regular expression
sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
###################################
# Read .frag file 
sub FragFileReader{
	my ($taffFile)  = @_;
	#open tf file
	my @information = ();
	if(open(FRAG,$taffFile)){
		print "MESSAGE Reading File $taffFile\n";
		@information=<FRAG>;
	}else{
		print "ERROR Can't open $taffFile\n";
		exit(1);
	}
	close(FRAG);
	# Variable creations
	my $aux_na         = 0;
	my $aux_ba         = 0;
	my @array_integral = ();
 	my @array_X        = ();
 	my @array_Y        = ();
 	my @array_Z        = ();
 	my @array_atoms    = ();
	my @array_atoms_x  = ();
	my @array_atoms_y  = ();
	my @array_atoms_z  = ();
	#set data in varbiales
	foreach my $line (@information){
		chomp($line);
		my @data=split(" ",$line);
		if($data[0] eq "X"){
			#if($data[4])
			push @array_integral, $data[4];
	 		push @array_X, $data[1];
	 		push @array_Y, $data[2];
	 		push @array_Z, $data[3];
	 		$aux_ba++;

		}else{
			push @array_atoms, $data[0];
			push @array_atoms_x, $data[1];
			push @array_atoms_y, $data[2];
			push @array_atoms_z, $data[3];
			$aux_na++;
		}
	}
	my $MaxIntegral = max(@array_integral);
	my $CutUmbral   = $MaxIntegral/1000.0;
	#
	# Counter of atoms and basins
	my @array_set_atoms  = ([@array_atoms],
                           [@array_atoms_x],
                           [@array_atoms_y],
                           [@array_atoms_z]);
	my @array_set_basin  = (0,
                           [@array_integral],
                            0,
                           [@array_X],
                           [@array_Y], 
                           [@array_Z],
                            0,
                            0);
    #
	# This code block and eliminate negative and low integals values
	# dump (@array_set_basin);
	for (my $i = 0; $i < $aux_ba; $i++) {
		#print "Viendo $array_set_basin[1][$i] en $i\n";
		if($array_set_basin[1][$i] < $CutUmbral){
		#	print "Eliminando.. $array_set_basin[1][$i]\n";
			splice @{$array_set_basin[1]}, $i, 1;
			splice @{$array_set_basin[3]}, $i, 1;
			splice @{$array_set_basin[4]}, $i, 1;
			splice @{$array_set_basin[5]}, $i, 1;
			$aux_ba--;
			$i--;
		}
	}
	# print "AFTER\n";
	# dump (@array_set_basin);
	$num_atoms_global=$aux_na;
	$num_basin_global=$aux_ba;
	#
	my @matrix = ([@array_set_atoms],
                  [@array_set_basin]);
	return @matrix;
}
###################################
# Box size manipulation
sub BoxCoordinates{
	my ($f1,$f2) = @_;
	#
	my @frag1    = @{$f1};
	my @frag2    = @{$f2};
	#llama al otro metodo y luego suma los tamaños
	my $dist1 = BoxSizeCalculation(@frag1);
	my $dist2 = BoxSizeCalculation(@frag2);
	#
    $Tam1 = $dist1;
    $Tam2 = $dist2;
    #
	my $sum = $dist2+$dist1;
	# Box size = (frag1+frag2) * 0.8 (Default).
	my $Box = $sum * $Box_Multiplication_Factor;
	# This means, maximum size of each fragments, the we choose only X of the space
	$Box_x = $Box;
	$Box_y = $Box;
	$Box_z = $Box;
}
###################################
# Get optimal size
sub BoxSizeCalculation{
	my @frag = @_;
	my @x    = ();
	my @y    = ();
	my @z    = ();
	my @at   = ();
	foreach my $line(@frag){
		($at[++$#at],$x[++$#x],$y[++$#y],$z[++$#z]) = split(' ',$line);
	}
	my @dist = ();
	push(@dist, MinMax(@x));
	push(@dist, MinMax(@y));
	push(@dist, MinMax(@z));
	my @sorted = sort { $b <=> $a } @dist;
	return $sorted[0];
}
###################################
# Get maximun or minimun from an array, Box use
sub MinMax{
	# We got the biggest distance (X Y Z)
	my @array = @_;
	my($min, $max);
	for(@array){
		$min = $_ if !$min || $_ < $min;
    	$max = $_ if !$max || $_ > $max
	};
	my $dist = $max-$min;
	return $dist;
}
###################################
# Read files
sub clear_data{
	my (@data)=@_;
	my @datos=();
	foreach my $info(@data){
		my @tmp=split(' ', $info);
		$datos[++$#datos]=[@tmp];
	}
	# Print $datos[0][0]
	# First [] indicates the position, or line in the .frag file.
	# Second [] indicates data in specific:
	#	0 = Atom type.
	#	1 = X
	#	2 = Y
	#	3 = Z
	return @datos;
}
###################################
# Verification
sub verification{
	my ($a1, $a2, $dist)=@_;
	# hash values	
	my $v1  = $Atomic_radii{$a1} || $Atomic_radii{default}; 
	my $v2  = $Atomic_radii{$a2} || $Atomic_radii{default};
	my $sum = $v1+$v2;  
	my $resultado;
	# Steric effects if radii_1 + radii_2 < distance
	if($dist <= $sum){
		# Steric problem	
		$resultado = 1; 
	}else{
		$resultado = 0;
	}
	return $resultado;
}
###################################
# Steric Impediment
sub steric_impediment{
	# array are send by reference
	my ($cords1,$cords2)=@_;
	# arrays	
	my @coords1=@{$cords1}; 
	my @coords2=@{$cords2};
	# get size
	my $distance;
	my ($final_trial, $resultado)=(0, 0);
	# all the atoms from the first structure
	for (my $i=0; $i<scalar(@coords1);$i++){
		# second structure
		for(my $j=0; $j<scalar(@coords2);$j++){
			# only atoms, not basins
			$distance=distance_point($coords1[$i][1],$coords1[$i][2],$coords1[$i][3],$coords2[$j][1],$coords2[$j][2],$coords2[$j][3]);
			$final_trial=$final_trial+verification($coords1[$i][0],$coords2[$j][0],$distance);
			if($final_trial==1){
				$resultado=1;
				last;
			}
		}
	}
	return $resultado;
}
###################################
# Distance between one point to all other in a coord xyz
sub medium_point {
	# array coords basin 1 and basin 2    
	my ($coords_1,$coords_2) = @_;
    my @arr1 = @{$coords_1};
    my @arr2 = @{$coords_2};
    #    
    my $med_x = ( ($arr1[0] + $arr2[0]) / 2 );
    my $med_y = ( ($arr1[1] + $arr2[1]) / 2 );
    my $med_z = ( ($arr1[2] + $arr2[2]) / 2 );
    #
    my @medium = ($med_x, $med_y, $med_z);
    return @medium;
}
###################################
# Distance between one point to all other in a coord xyz
sub distance_point {
	# array coords basin 1 and basin 2
	my ($p1,$p2,$p3, $x1, $y1, $z1) = @_;
	# measure distance between two point
	my $dist = sqrt(
				($x1-$p1)**2 +
				($y1-$p2)**2 +
				($z1-$p3)**2
				); 
	return $dist;
}
###################################
# Parsing name atom 
sub parsing_name_atom {
	my ($atom_set) = @_;
	my @array = ();
	@array = split('_', $atom_set);	
	# return array	
	return $array[0];
}
###################################
# This function opens a file, and returns the info in an array of lines
sub read_file {
	# filename
	my ($input_file) = @_;
	my @lines = ();
	my @array = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	@lines = <FILE>;
	close (FILE);
	# loop
	foreach my $a (@lines){
		chomp ($a);
		$array[++$#array] = $a;
	}
	# return array	
	return @array;
}
###################################
# compute the center of mass
sub measure_center {
	my ($coord_x,$coord_y,$coord_z) = @_;
	my $num_data = scalar (@{$coord_x});
	my @array  = ();
	my $weight = 1;
	# variable sum
	my $sum_weight = 0;
	my $sum_x = 0;
	my $sum_y = 0;
	my $sum_z = 0;
	for ( my $j = 0 ; $j < $num_data ; $j = $j + 1 ){
		$sum_weight+= $weight;
		$sum_x+= $weight * @$coord_x[$j];
		$sum_y+= $weight * @$coord_y[$j];
		$sum_z+= $weight * @$coord_z[$j];		
	}
	my $com_x = $sum_x / $sum_weight;
	my $com_y = $sum_y / $sum_weight;
	my $com_z = $sum_z / $sum_weight;
	# array
	@array = ($com_x,$com_y,$com_z);
	# return array	
	return @array;
}
###################################
# Returns the additive inverse of v(-v)
sub vecinvert {
	my ($center_mass) = @_;
	my @array         = ();
	foreach my $i (@$center_mass) {
		my $invert        = $i * -1;
		$array[++$#array] = $invert; 
	}	
	# return array	
	return @array;
}
###################################
# Returns the vector sum of all the terms
sub vecadd {
	my ($coord_x,$coord_y,$coord_z,$vecinvert_cm ) = @_;
	my $num_data = scalar (@{$coord_x});
	my @array   = ();
	my $sum_coord_x;
	my $sum_coord_y;
	my $sum_coord_z;
	# array 
	my @array_x = ();
	my @array_y = ();
	my @array_z = ();
	for ( my $i = 0 ; $i < $num_data ; $i = $i + 1 ){	
		$sum_coord_x = @$coord_x[$i]+@$vecinvert_cm[0] ; 
		$sum_coord_y = @$coord_y[$i]+@$vecinvert_cm[1] ;
		$sum_coord_z = @$coord_z[$i]+@$vecinvert_cm[2] ;
		# save array
		$array_x[++$#array_x] = $sum_coord_x;
		$array_y[++$#array_y] = $sum_coord_y;
		$array_z[++$#array_z] = $sum_coord_z;
	}
	@array = ( [@array_x], 
              [@array_y], 
              [@array_z] ); 
	# return array	
	return @array;
}
###################################
# Focus on the origin of molecular systems
sub CenterMolecule {
	my ($input, $energy, $total_num)=@_;
	my @problematic_molecule=@{$input};
	my @coord_x=@{$problematic_molecule[1]};
	my @total_array = ();
	#
	my @array_center_mass = measure_center(\@{$problematic_molecule[1]},\@{$problematic_molecule[2]},\@{$problematic_molecule[3]});
	my @array_vecinvert   = vecinvert(\@array_center_mass);
	my @array_catersian   = vecadd (\@{$problematic_molecule[1]},\@{$problematic_molecule[2]},\@{$problematic_molecule[3]},\@array_vecinvert);
	#
	my @at    = @{$problematic_molecule[0]};
	my @cx    = @{$array_catersian[0]};
	my @cy    = @{$array_catersian[1]};
	my @cz    = @{$array_catersian[2]};
	#
	open(FILE, ">>Mol_Center.xyz");
	print FILE "$total_num\n";
	print FILE "$energy\n";
	for (my $i=0; $i < scalar(@at) ; $i++) {
		my $x_coord = sprintf '%.6f', $cx[$i];
		my $y_coord = sprintf '%.6f', $cy[$i];
		my $z_coord = sprintf '%.6f', $cz[$i];
		print FILE "$at[$i]\t$x_coord\t$y_coord\t$z_coord\n";
	}
	close (FILE);
}
###################################
# Generation of random co-ordinates
sub gen_xyz {
	my $x = rand($Box_x);
	my $y = rand($Box_y);
	my $z = rand($Box_z);
	my $x_coord = sprintf '%.6f', $x;
	my $y_coord = sprintf '%.6f', $y;
	my $z_coord = sprintf '%.6f', $z;
	my @coords =($x_coord, $y_coord, $z_coord);
	return @coords;
}
###################################
# phi, theta, psi
sub gen_ptp {
	my $pi     = 3.14159265;
	my $phi    = sprintf '%.6f', rand()*2*$pi;
	my $theta  = sprintf '%.6f', rand()*2*$pi;
	my $psi    = sprintf '%.6f', rand()*2*$pi;
	my @angles =($phi, $theta, $psi);
	return @angles;
}
###################################
# Generating .cart files
sub cart_file_write {
	my ($name_without_extension,$num_atom, $num_basin, @matrix)=@_;
	open(FILE, ">$name_without_extension.cart");
	for(my $i=0; $i<$num_atom;$i++){
		print FILE "$matrix[0][0][$i]\t$matrix[0][1][$i]\t$matrix[0][2][$i]\t$matrix[0][3][$i]\n";
	}
	for(my $j=0; $j<$num_basin;$j++){
		print FILE "X\t$matrix[1][3][$j]\t$matrix[1][4][$j]\t$matrix[1][5][$j]\n"
	}
	close (FILE);
}
###################################
# multiplication of attractors 
sub basin_integral_multi{
	my ($m1, $m2,$basin1,$basin2)=@_;
	my @matrix1=@{$m1};
	my @matrix2=@{$m2};
	my @array_multi=();
	for(my $i=0; $i<$basin1; $i++){
		for(my $j=0; $j<$basin2; $j++){
			my $tmp=$matrix1[1][1][$i]*$matrix2[1][1][$j];
			push(@array_multi,$tmp);
		}
	}
	return @array_multi;
}
###################################
# Cartesian distance measuring
sub distance_medition{
	my ($ba1, $ba2)=@_;
	my @big_array1=@{$ba1};
	my @big_array2=@{$ba2};
	my @dist_total=();
    #
    my $count = 0;
    #
	for (my $i = 0; $i < scalar (@{$big_array1[2]}); $i++) {
		for (my $j = 0; $j < scalar (@{$big_array2[2]}); $j++) {	
			my $tmp_dist = distance_point ( $big_array1[1][$i],$big_array1[2][$i],$big_array1[3][$i], $big_array2[1][$j],$big_array2[2][$j],$big_array2[3][$j]);
			push (@dist_total,$tmp_dist);
            #
            my $key = $count;
            $HashCoordsAttra{$key} ="$big_array1[1][$i],$big_array1[2][$i],$big_array1[3][$i] = $big_array2[1][$j],$big_array2[2][$j],$big_array2[3][$j]" ;
            $HashDist{$key}        = $tmp_dist;
            $count++;
		}
	}
	return @dist_total;
}
###################################
# Here, magics happends
# Essencially the fragments coordinates are randomly rotated and traslatated
sub position_change_fragment{
	my(@FragLines)  = @_;
	# Rotate
	my @base_angles = gen_ptp();
	my $phi         = $base_angles[0];
	my $theta       = $base_angles[1];
	my $psi         = $base_angles[2];
	# Translate
	my @base_xyz    = gen_xyz();
	my $base_x      = $base_xyz[0];
	my $base_y      = $base_xyz[1];
	my $base_z      = $base_xyz[2];
	# Do the trig
	my $cos_phi     = sprintf '%.6f', cos($phi);
	my $cos_theta   = sprintf '%.6f', cos($theta);
	my $cos_psi     = sprintf '%.6f', cos($psi);
	my $sin_phi     = sprintf '%.6f', sin($phi);
	my $sin_theta   = sprintf '%.6f', sin($theta);
	my $sin_psi     = sprintf '%.6f', sin($psi);
	# Make the rotation matrix
	my $D = new Math::Matrix ([$cos_phi,$sin_phi,0],[-$sin_phi,$cos_phi,0],[0,0,1]);
	my $C = new Math::Matrix ([1,0,0],[0,$cos_theta,$sin_theta],[0,-$sin_theta,$cos_theta]);
	my $B = new Math::Matrix ([$cos_psi,$sin_psi,0],[-$sin_psi,$cos_psi,0],[0,0,1]);
	my $A = $B->multiply($C)->multiply($D);
	#
	my @ar_x   = ();
	my @ar_y   = ();
	my @ar_z   = ();
	my @coords = ();
	#
	while (my $Fline = shift (@FragLines)) {
		my @Cartesians              = split '\s+', $Fline;
		my ($Atom_label, @orig_xyz) = @Cartesians;
		my $matrix_xyz              = new Math::Matrix ([$orig_xyz[0],$orig_xyz[1],$orig_xyz[2]]);
		my $trans_xyz               = ($matrix_xyz->transpose);
		my $rotated_xyz             = $A->multiply($trans_xyz);
		my @new_xyz                 = split '\n+',$rotated_xyz;
		# Save rotated fragment
		push(@new_xyz,$Atom_label);
		#
		my $new_x = $base_x + $new_xyz[0];
		my $new_y = $base_y + $new_xyz[1];
		my $new_z = $base_z + $new_xyz[2];
		if($new_xyz[3] eq "X"){
			push (@ar_x,$new_x);
			push (@ar_y,$new_y);
			push (@ar_z,$new_z);
		}else{
			push (@coords,"$Atom_label\t$new_x\t$new_y\t$new_z");	
		}	
	}
    #
	my @bigarray = ([@coords],[@ar_x],[@ar_y],[@ar_z]);
	return @bigarray;
}
###################################
# Equation
sub equation{
	my ($multi, $dist) = @_;
	my @b_multi        = @{$multi};
	my @b_dist         = @{$dist};
	my $sum=0;
	for (my $i = 0; $i < scalar (@b_multi); $i++) {
		#print "f+ * f- = $b_multi[$i] with distance of $b_dist[$i]\t\n";
		my $division = ( $b_multi[$i] / $b_dist[$i] );
		#print "f+ * f- = $b_multi[$i] with distance of $b_dist[$i]\t Division: $division\n";
		$sum+=$division;
	}
	return $sum;
}
###################################
# write xyz file
sub write_xyz_file_center{
	my ($pid,$arr_center, $atom) = @_;
    my @arr = @{$arr_center};
    open (XYZ,">>tmp$pid.tmp");
    print XYZ "$atom\t$arr[0]\t$arr[1]\t$arr[2]\n";
    close (XYZ);
}
###################################
# write xyz file data
sub write_xyz_file_header{
	#write only the header for each tmp file
	my ($pid, $atoms_total, $sum)=@_;
	open(XYZ,">>tmp$pid.tmp");
	print XYZ "$atoms_total\n";
	print XYZ "E = $sum\n";
	close (XYZ);
}
###################################
# write xyz file
sub write_xyz_file{
	#write the coordinate information for each tmp file
	my ($pid,  $Attract, @big_array )=@_;
    #
    my @tmpr = @{$Attract};
    #
	open(XYZ,">>tmp$pid.tmp");
	my @coords	=clear_data(@{$big_array[0]});
	for(my $i=0; $i<scalar(@coords); $i++){
		#
		#my $c_x = sprintf '%.6f', $coords[$i][1];
		#my $c_y = sprintf '%.6f', $coords[$i][2];
		#my $c_z = sprintf '%.6f', $coords[$i][3];
		#
		print XYZ "$coords[$i][0]\t$coords[$i][1]\t$coords[$i][2]\t$coords[$i][3]\n";
	}
	close (XYZ);
}
###################################
# sort in descending order energies
sub sort_energies{
	# Simple sort, uses hash for sorting
	my @files     = glob "*.tmp";
	my ($na, $me) = @_;
	my @energies  = ();
	my @sorted_numbers;
	my %super;
	my $natom;
	foreach my $file (@files){
		# Open files .tmp
		open (TMP, "<$file") || die "cannot open $file in sort subroutine\n";
		my @lines=<TMP>;
		close(TMP);
		#
		$natom     = $lines[0];
		my $energy = $lines[1];
		$energy=~s/(E = )?//g;
		chomp($energy);
		push(@energies, $energy);
		my $coords="";
		#
        foreach $a  (2..($#lines)){
			$coords=$coords.$lines[$a];
		}
		$super{$energies[$#energies]}=$coords;	
		unlink "$file";
	}
	#
	@sorted_numbers = sort { $b <=> $a } @energies;
	@ValueJ = @sorted_numbers;
	# write inputfile .XYZ
	open(NEW, ">tmp_kick-fukui.xyz");
	my $aux = 0;
	#
	foreach my $item (@sorted_numbers){
		if($aux < $XYZSize){
			print NEW "$natom", "E = $item\n",$super{$item};
		}else{
			last;
		}
		$aux++;
	}
	close(NEW);
}
#####################################
# Normal order
sub unify_tmp_files_noSort{
	my ($na, $me)=@_;
	my @files =glob "*.tmp";
	open(NEW, ">tmp_kick-fukui.xyz");
	#
	my @energies = ();
	#
	foreach my $file (@files){
		# Open files .tmp
		open (TMP, "<$file") || die "cannot open $file in sort subroutine\n";
		my @lines=<TMP>;
		close(TMP);
		#
		my $energy = $lines[1];
		$energy=~s/(E = )?//g;
		chomp($energy);
		push(@energies, $energy);
	}
	#
	@ValueJ = sort { $b <=> $a } @energies;
	#
	my $aux = 0;
	foreach my $file (@files){
		open (TMP, "<$file") || die "cannot open $file in sort subroutine\n";
		my @lines=<TMP>;
		close(TMP);
		#
		if($aux < $XYZSize){
			print NEW @lines;
		}else{
			last;
		}
		$aux++;
	}
	close(NEW);
	#
	foreach my $file (@files){
		unlink "$file";
	}
}
#####################################
# Logo Kick-Fukui
sub logo {
	print "\n";	
	print "      ____  __.__        __            ___________     __         .__        \n";
	print "     |    |/ _|__| ____ |  | __        \\_   _____/_ __|  | ____ __|__|      \n";
	print "     |      < |  |/ ___\\|  |/ /  ______ |    __)|  |  \\  |/ /  |  \\  |    \n";
	print "     |    |  \\|  \\  \\___|    <  /_____/ |     \\ |  |  /    <|  |  /  |   \n";
	print "     |____|__ \\__|\\___  >__|_ \\         \\___  / |____/|__|_ \\____/|__|  \n";
	print "             \\/       \\/     \\/             \\/             \\/         	\n";
	print "\n\n";
	print " Search global minimum structures of atomic clusters and molecules in the gas phase,\n"; 
	print "       using as a Coulomb interaction between Fukui functions-based method \n";
	print "\n";	
} 
###################################
# Configuration input file 
sub ConfigFileReader{
	#
	my $configfile = $ARGV[0];
	#
	if (not defined $configfile) {
		die "\n Kick-Fukui_Algorithm must be run with:\n\n Usage :\n\tperl $0 [configure-file]\n
		\n Help :\n\tperl $0 -help \n\n\n";
		exit(1);
	}
	#
	if($configfile=~/-help/i) {
		show_help();
		exit (1);
	}
	# Open configure file
	open(CONFIG,$configfile) or die "ERROR Configuration file doesn't exist\n\n";
	my @lines=<CONFIG>;
	close(CONFIG);
	#
	my $aux = 0;
	my ($frag1, $frag2);
	foreach my $line(@lines){
		chomp($line);
		my $letter = substr($line, 0, 1);
		my ($key, $value)=split("=",$line);
		# Check repository 
		# the char # marks a commentary
		if($letter ne "#" && $line ne ""){
			$value = trim($value);
			#
			if($key=~/initial_species/i){
				$Num_of_geometries=$value;
				$aux++;
			}
			#
			if ($key=~/final_species/i) {
				$XYZSize=$value;
				$aux++;
			}
			#
			if($key=~/mol_1/i){
				$frag1 = $value;
				$aux++;
			}
			#
			if($key=~/mol_2/i){
				$frag2 = $value;
				$aux++;
			}
			#
			if($key=~/box_size_factor/i){
				my @box=split(",",$value);
				if(scalar(@box)==1){
					$Box_Multiplication_Factor=$value;
				}elsif(scalar(@box)==3){
					if(looks_like_number($box[0]) && looks_like_number($box[1]) && looks_like_number($box[2])){
						print "MESSAGE Box Size given by user\n";
						$BoxAutomatic=2;
						$Box_x = $box[0];
						$Box_y = $box[1];
						$Box_z = $box[2];
					}else{
						print "ERROR Box coordinates are not numeric\n";
						exit(1);
					}
				}else{
					print "ERROR Box Size error, usage box_size_factor = X,Y,Z \n";
					exit(1);
				}
			}
			#
			if($key=~/energy_order/i){
				$sort_option = $value;
				if (($sort_option=~/yes/i )) {
					print "MESSAGE Choose Sort by Descending Order  \n";
				} elsif (($sort_option=~/no/i )) {
					print "MESSAGE Select Stochastic Order \n";
				} else {
					print "ERROR Set the Correct Option energy_order = YES/NO \n";
					exit (1);
				}
			}
			#
			if($key=~/duplicate_species/i){
				$duplicate_mol = $value;
				if (($duplicate_mol=~/yes/i )) {
					print "MESSAGE Grigoryan-Springborg's Algorithm \n";
				} elsif (($duplicate_mol=~/no/i )) {
					print "MESSAGE Don't Select Duplicates \n";
				} else {
					print "ERROR Set The Correct Option duplicate_species YES/NO \n";
					exit (1);
				}
			}
			#
			if($key=~/output_file_name/i){
				$output_file=$value;
				$aux++;
			}
			
		}
	}
	if($aux!=5){
		print "ERROR Wrong number of params in $ARGV[0]\naux = $aux\n";
		exit(1);
	}
	return ($frag1, $frag2);
}
####################################
# Summary Kick-Fukui Algorithm
sub resume_log_file {
	my ($na,$me,$timeS,$timeE,$timeT)=@_;
	#
	print "\n";
	print "  <> Summary Kick-Fukui Algorithm <> \n\n";
	# Sample Size
	print "    * Sample Size: $Num_of_geometries\n";	
	# Molecules f+ and f- 
	print "    * Molecules \n";
	print "        f+ : $na.frag \n";
	print "        f- : $me.frag \n";
	# Box size
	my $c_x = sprintf '%.2f', $Box_x;
	my $c_y = sprintf '%.2f', $Box_y;
	my $c_z = sprintf '%.2f', $Box_z;
	print "    * Box size : $c_x Å x $c_y Å x $c_z Å\n";
	# Output file xyz
	print "    * Output File: $Final_File\n";
	#
	my $J_value = sprintf '%.6f', $ValueJ[0];
	print "    * Maximun Coulombic-Interaction Value (J) : $J_value\n";
	#
	print "\n  <> Performance Data <> \n\n";
	# Total number of process
	print "    * N° Process : $ncpus\n";	
	# Start of the program
	print "    * Start : $timeS\n";
	# Final execution of the programme 
	print "    * End : $timeE\n";
	# Execution time  	
	print "    * Execution time :$timeT\n\n";
	#
	print "MESSAGE Normal Termination Kick-Fukui Algorithm\n\n\n";
	#
}
###################################
# Help 
sub show_help {
	print " <> To use Kick-Fukui Algorithm <> \n\n";
	#
	print "\n";
	print " * The Initial Population size (1000N, N = Atoms number) \n";
	print "\t initial_species = 1000 \n";
	print "\n";
	print " * Number of final geometries \n";
	print "\t final_species = 100 \n";
	print "\n";
	print " * * Name of molecular species in the simulation \n";
	print " * File with the nucleophilic Fukui function \n";
	print "\t Mol_1 = f+FileCoords.frag \n";
	print " * File with the electrophilic Fukui function \n";
	print "\t Mol_2 = f-FileCoords.frag \n";
	print "\n";
	print " * Box size the sum of the sides from both molecules multiply for a factor (Default 1) \n";
	print " * The factor size or size of the box (in Angstroms) length, width, and height \n";
	print "\t box_size_factor = 0.8 or box_size_factor = 10, 10, 10\n";
	print "\n";
	print " * Order Coulombic-Interaction Value (J) in a descending (YES) or scholastic (NO) manner \n";
	print "\t energy_order = YES \n";
	print "\n";
	print " * Search for duplicate molecular species (YES/NO) \n";
	print "\t duplicate_species = NO \n";
	print "\n";
	print " * Output xyz file name (Default name Kick-Fukui-Output.xyz) \n";
	print "\t output_file_name = NameOutputFile.xyz \n";
	print "\n";
	
}



###################################
# MAIN
my @matrix_1;
my $atoms_num_1;
my $basin_num_1;
# extension 1
my $without_extension_1;
#
my @matrix_2;
my $atoms_num_2;
my $basin_num_2;
# extension 2
my $without_extension_2;
#
logo ();
my $tiempo_inicial  = new Benchmark; #funcion para el tiempo de ejecucion del programa
my $datestringStart = localtime();
# delete .tmp files
unlink glob "*.tmp";
unlink ("$output_file.xyz");
# 
# Read configure file
my ($file_name_1, $file_name_2) = ConfigFileReader();
#
# Read .frag files
# 
@matrix_1    = FragFileReader($file_name_1);
$atoms_num_1 = $num_atoms_global;
$basin_num_1 = $num_basin_global;
# extension
$file_name_1=(split("/",$file_name_1))[-1];
( $without_extension_1 = $file_name_1) =~ s/\.[^.]+$//;
#
@matrix_2    = FragFileReader($file_name_2);
$atoms_num_2 = $num_atoms_global;
$basin_num_2 = $num_basin_global;
# extension
$file_name_2=(split("/",$file_name_2))[-1];
($without_extension_2 = $file_name_2) =~ s/\.[^.]+$//;
#
#
# Creation of .cart files
cart_file_write($without_extension_1,$atoms_num_1,$basin_num_1,@matrix_1);
cart_file_write($without_extension_2,$atoms_num_2,$basin_num_2,@matrix_2);
# Reading .cart file 
open(FRAG1,"$without_extension_1.cart" ) or die "ERROR Unable to open fragment file: $without_extension_1.cart";
open(FRAG2,"$without_extension_2.cart" ) or die "ERROR Unable to open fragment file: $without_extension_2.cart";
my @f1=<FRAG1>;
my @f2=<FRAG2>;
close(FRAG1);
close(FRAG2);
unlink glob "*.cart";
#
# Parallel processing
my $pm = new Parallel::ForkManager($ncpus);
if ($BoxAutomatic == 1){
	BoxCoordinates(\@f1,\@f2);
}
#
# The integral value for each attractor doesn't change moving the fragments.
# So the multiplication between atractor (integral) is only done once.
my @multi = basin_integral_multi(\@matrix_1, \@matrix_2, $basin_num_1, $basin_num_2);
#
foreach my $iteration(1 .. $Num_of_geometries){
	$pm->start($iteration) and next;
	# All children process have their own random
	srand();
	# To verify steric impediment
	my $decision = 1;		
	#
	my @resume1;
	my @resume2;
	#
	while($decision!=0){
		# Randomize atoms and basin positions
		@resume1 = position_change_fragment(@f1);
		@resume2 = position_change_fragment(@f2);
		# Get atoms coordinates
		my @coords_1 = clear_data(@{$resume1[0]});
		my @coords_2 = clear_data(@{$resume2[0]});
		# Verify for steric impediment, 1 yes, 0 no;             
		$decision = steric_impediment(\@coords_1, \@coords_2);
	}
	###################################
	# Spherical restriction algorithm 
	# 
	# Attractors X=@{$resume1[1]}, Y=@{$resume1[2]},Z=@{$resume1[3]}
	# measure center of mass 
	my @array_center_mass_attrac1 = measure_center(\@{$resume1[1]},\@{$resume1[2]},\@{$resume1[3]});
	my @array_center_mass_attrac2 = measure_center(\@{$resume2[1]},\@{$resume2[2]},\@{$resume2[3]});
	#
	my @medium = medium_point (\@array_center_mass_attrac1,\@array_center_mass_attrac2);
	#
	# Sphere radius
	my @ValuesTam   = sort {$a cmp $b} ($Tam1,$Tam2);
	#
	my $Raddi_Value = ($ValuesTam[0]/2) * $weight_spherical;
	#
	my $Radii_Xmin  = ($medium[0] + (-$Raddi_Value));
	my $Radii_Ymin  = ($medium[1] + (-$Raddi_Value));
	my $Radii_Zmin  = ($medium[2] + (-$Raddi_Value));
	my @Ra_min      = ($Radii_Xmin,$Radii_Ymin,$Radii_Zmin);
	#
	my $Radii_Xmax  = ($medium[0] + $Raddi_Value);
	my $Radii_Ymax  = ($medium[1] + $Raddi_Value);
	my $Radii_Zmax  = ($medium[2] + $Raddi_Value);
	my @Ra_max      = ($Radii_Xmax,$Radii_Ymax,$Radii_Zmax);
	#    
	# The measurement of the distance of the attractors is carried out
	my @dist = distance_medition(\@resume1,\@resume2);
	# Hash for multiplication
	for ( my $i =0; $i < scalar (@multi); $i=$i+1 ){
		my $key = $i;
		$HashMult{$key}=$multi[$i];
	}
	#
	my @Mult_1    = ();
	my @Dist_1    = ();
	my @Attracs_1 = ();
	my @Attracs_2 = ();
	# A point (x,y,z) is within the sphere, with center (cx,cy,cz) and radius r if
	# ( x-cx )^2 + (y-cy)^2 + (z-cz)^2 < r^2     
	foreach my $key (keys %HashMult){
		my ($Arr1,$Arr2)                   = split '\=', $HashCoordsAttra{$key};
		my ($axis_x_1,$axis_y_1,$axis_z_1) = split ',', $Arr1;
		my ($axis_x_2,$axis_y_2,$axis_z_2) = split ',', $Arr2;
		#
		my $diameter = ($Raddi_Value)**2;
		#
		my $value_1 = (($axis_x_1 - $medium[0])**2) + (($axis_y_1 - $medium[1])**2) + (($axis_z_1 - $medium[2] )**2);
		my $value_2 = (($axis_x_2 - $medium[0])**2) + (($axis_y_2 - $medium[1])**2) + (($axis_z_2 - $medium[2] )**2);
		#
		if ( ($value_1 <= $diameter ) && ( $value_2 <= $diameter ) ) {
			push (@Mult_1, $HashMult{$key});
			push (@Dist_1, $HashDist{$key});
			# Attractors
			push (@Attracs_1,"$axis_x_1\t$axis_y_1\t$axis_z_1");
			push (@Attracs_2,"$axis_x_2\t$axis_y_2\t$axis_z_2");
		}
	}
	#
	# End of restricted spherical algorithm
	# print Dumper(\%HashMult);
	# print Dumper(\%HashDist);
	# print Dumper(\%HashCoordsAttra);
	#
	# Maximum Matching value
	#
	my $sum = equation(\@Mult_1,\@Dist_1);
	# Write tmp file with Maximum Matching value and coordinates
	write_xyz_file_header($iteration,($atoms_num_1+$atoms_num_2),$sum);
	write_xyz_file($iteration,\@Attracs_1,@resume1);
	write_xyz_file($iteration,\@Attracs_2,@resume2);
	#
	# write_xyz_file_center ($iteration,\@array_center_mass_attrac1, "He");
	# write_xyz_file_center ($iteration,\@array_center_mass_attrac2, "He");
	#
	# write_xyz_file_center ($iteration,\@medium, "Rn");
	#
	# write_xyz_file_center ($iteration,\@Ra_min, "Xe");
	# write_xyz_file_center ($iteration,\@Ra_max, "Xe");
	#
	# Delete data hash
	%HashCoordsAttra = ();
	%HashDist        = ();
	%HashMult        = ();
	#
	$pm->finish;
}
#
$pm->wait_all_children;
#
#
# Sort by Coulombic integral value
#
if ($sort_option=~/yes/i ) {
	sort_energies($without_extension_1, $without_extension_2);
	print "MESSAGE Delete files .tmp\n";
	$Final_File = "$output_file\_sort.xyz";
} else {
	unify_tmp_files_noSort($without_extension_1, $without_extension_2);
	print "MESSAGE Delete files .tmp\n";
	$Final_File = "$output_file\_stochastic.xyz";
}
#
#
# Center molecule at origin
my @FilesCenter = read_file("tmp_kick-fukui.xyz");
my $count_cen   = 0;
my $string_type = "E =";
my @array_cen   = ();
my @value_j     = ();
foreach my $i (@FilesCenter) {
	if ($i =~ m/$string_type/) {
		push (@array_cen, $count_cen);	
		push (@value_j, $i);	
	}
	$count_cen++;
}
for (my $i=0; $i < scalar(@array_cen) ; $i++) {
	my (@at,@cx,@cy,@cz) = ();
	#
	my $sum_cen_1 = $array_cen[$i] + 1;
	my $sum_cen_2 = $array_cen[$i] + ($atoms_num_1 + $atoms_num_2) ;
	my $Coords_cen; 
	for my $x ($sum_cen_1 .. $sum_cen_2) {
		my ($name_atom, $axis_x, $axis_y, $axis_z) = split (" " , $FilesCenter[$x]);
		push (@at,$name_atom);
		push (@cx,$axis_x);
		push (@cy,$axis_y);
		push (@cz,$axis_z);
	}
	my @mol_center = ([@at],[@cx],[@cy],[@cz]);
	CenterMolecule (\@mol_center, $value_j[$i], ($atoms_num_1 + $atoms_num_2));
}
#
#
# Grigoryan and Springborg's algorithm 
my $root                = dirname(abs_path($0));
my $DupliScript         = "$root/DupGrigoryanSpringborg.pl";
my $threshold_duplicate = 0.009;
my @final_coords = ();
if (($duplicate_mol=~/yes/i )) {
	#
	print "MESSAGE Find Structures Similar\n";
	#
	system ("$DupliScript $threshold_duplicate Mol_Center.xyz");
	#
	my @filesdupli = read_file("FinalCoords_GS.xyz");
	#
	my $count  = 0;
	my $word_1 = "E =";
	my @linearray_1 = ();
	foreach my $i (@filesdupli) {
		if ( ($i =~ m/$word_1/) ) {
			push (@linearray_1,$count);		
		}
		$count++;
	}
	my @array_keys  = ();
	for (my $i=0; $i < scalar(@linearray_1) ; $i++) { 
		my $tmpsum_1     = $linearray_1[$i] + 1;
		my $tmpsum_2     = $linearray_1[$i] + ($atoms_num_1 + $atoms_num_2) ;
		my $concatCoords; 
		for my $x ($tmpsum_1 .. $tmpsum_2) {			
			$concatCoords.="$filesdupli[$x]\n";
		}
		push (@final_coords,$concatCoords);
	}
	print "MESSAGE Structures similar end\n";
}
#
#
# Renombrar archivos
if (($duplicate_mol=~/yes/i )) {
	rename ("$root/FinalCoords_GS.xyz", "$root/$Final_File") || die ( "ERROR Error in renaming" ); 
} else {
	rename ("$root/Mol_Center.xyz", "$root/$Final_File") || die ( "ERROR Error in renaming" );	
}
#
# Delete tmp files
unlink ("$root/FinalCoords_GS.xyz");
unlink ("$root/tmp_kick-fukui.xyz");
unlink ("$root/Mol_Center.xyz");
#
# Time in console is printed
my $tiempo_final  = new Benchmark;
my $datestringEnd = localtime();
my $tiempo_total  = timediff($tiempo_final, $tiempo_inicial);
#
# Outputfile is done
resume_log_file($without_extension_1,$without_extension_2,$datestringStart,$datestringEnd,timestr($tiempo_total));
