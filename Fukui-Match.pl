#!/usr/bin/env perl

# Notas de version:
# Caja Automatica (mol1 + mol2)*0.8
# taf reader 	
# Lector de numeros de geometrias como ultimo argumento
# calculo automatico en el numeor de procesos en relacion 1 CPU= 32 procesos
# Control de argumetnos de entrada.
# Muestra help
# Verificacion de Archivos se produce en el programa Osiris.pl
# Commented
# El factor de multiplicacion de la caja es definida por el usuario, DEFAULT=1
# Agregado la opcion de ordenar de peor a mejor (-s)
# Lector de parametros de entrada por medio de un archivo de configuracion
# Caja manual o automatica
# Elimnar los valores muy pequeños apra hacer la integral

#Modules.
use strict;
use Parallel::ForkManager;
use Math::Matrix;
use Benchmark; #entrega cuando demora el proceso, cuanto de CPU utilizó, etc
use Cwd qw(abs_path);
use File::Basename;
#use Data::Dump qw(dump ddx);
use Data::Dumper;
#
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
###########################
# SIMILITUD
my $Sim = "DeleteSimil.pl";
#
#
my $Tam1 =0;
my $Tam2 =0;
#
my %HashCoordsAttra = ();
my %HashDist        = ();
my %HashMult        = ();
###########################
#variables.
my $num_atoms_global;
my $num_basin_global;
my $Num_of_geometries;
my $XYZSize;
my $ncpus;
my $Box_x;
my $Box_y;
my $Box_z;
my $decision;
my $sort_option=1;
#Default is 1 unless factor is given.
my $Box_Multiplication_Factor=1;
my $BoxAutomatic=1;
my $program_name="Stochastic Fukui Function";
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
###########################

sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
sub convert_bohr_armstrong{my $s=shift; my $multi=$s*0.529177249; return $multi}
############################
#TAFF FIle TF reader
sub TaffFileReader{
	my ($taffFile)=@_;
	#open tf file
	my @information=();
	if(open(TAFF,$taffFile)){
		print "Opening: $taffFile\n";
		@information=<TAFF>;
	}else{
		print "Can't open $taffFile\n";
		exit(1);
	}
	close(TAFF);
	#variable creations
	my $aux_na=0;
	my $aux_ba=0;
	my @array_integral=();
 	my @array_X=();
 	my @array_Y=();
 	my @array_Z=();
 	my @array_atoms=();
	my @array_atoms_x=();
	my @array_atoms_y=();
	my @array_atoms_z=();
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
	my $MaxIntegral=max(@array_integral);
	my $CutUmbral=$MaxIntegral/1000.0;
	#counter of atoms and basins
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
	# This code block is to eliminate negative and low integals values
	#dump (@array_set_basin);
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
	#print "AFTER\n";
	#dump (@array_set_basin);
	$num_atoms_global=$aux_na;
	$num_basin_global=$aux_ba;

	my @matrix = ([@array_set_atoms],
              [@array_set_basin]);
	return @matrix;	#return variable
	
}
#############################
#Box size manipulation:
sub BoxCoordinates{
	my ($f1,$f2)=@_;		#Read fragmenst coordinates, and use them to reate the optimal box size
	my @frag1=@{$f1};
	my @frag2=@{$f2};
	#llama al otro metodo y luego suma los tamaños
	my $dist1=BoxSizeCalculation(@frag1);
	my $dist2=BoxSizeCalculation(@frag2);
    # OSVALDO
    $Tam1 = $dist1;
    $Tam2 = $dist2;
    #
	my $sum=$dist2+$dist1;
	my $Box=$sum*$Box_Multiplication_Factor;		#Box size= (frag1+frag2) * 0.8.
	$Box_x=$Box;			#this means , maxiimum size of each fragkments, the we choose only te X% of the space
	$Box_y=$Box;
	$Box_z=$Box;
	print "Box Size (Å):$Box\n";
	
}
###########################
# Get optimal size
sub BoxSizeCalculation{
	my @frag=@_;
	my @x=();
	my @y=();
	my @z=();
	my @at=();
	foreach my $line(@frag){
		($at[++$#at],$x[++$#x],$y[++$#y],$z[++$#z])=split(' ',$line);
	}
	my @dist=();
	push(@dist, MinMax(@x));
	push(@dist, MinMax(@y));
	push(@dist, MinMax(@z));
	my @sorted = sort { $b <=> $a } @dist;
	return $sorted[0];
}#Get max or min from an array, Box use
sub MinMax{
	#we got the biggest distance (X Y Z)
	my @array=@_;
	my($min, $max);
	for(@array){
		$min = $_ if !$min || $_ < $min;
    	$max = $_ if !$max || $_ > $max
	};
	my $dist=$max-$min;
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
	#print $datos[0][0]
	#first [] indicates the position, or line in the .cart file.
	#second [] indicates data in specific:
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
	my $v1=$Atomic_radii{$a1} || $Atomic_radii{default}; 
	my $v2=$Atomic_radii{$a2} || $Atomic_radii{default};
	my $sum= $v1+$v2;  
	my $resultado;
	# steric effects if radio1+radio2 < distance
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
	#  arrays	
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
    #
#    my $radii_1 = ($Tam1/2);
#    my $radii_2 = ($Tam2/2);
    #
#    my $value1_x = ( $arr1[0]/$radii_1 );
#    my $value1_y = ( $arr1[1]/$radii_1 );
#    my $value1_z = ( $arr1[2]/$radii_1 );
#    my $value2_x = ( $arr2[0]/$radii_2 );
#    my $value2_y = ( $arr2[1]/$radii_2 );
#    my $value2_z = ( $arr2[2]/$radii_2 );
    #
#    my $med_x = ( ($value1_x + $value2_x) / 2 );
#    my $med_y = ( ($value1_y + $value2_y) / 2 );
#    my $med_z = ( ($value1_z + $value2_z) / 2 );
    #
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
##################
sub read_file { #this function opens a file, and returns the info in an array of lines.
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
sub parsing_data {
	my ($num_line_min,$num_line_max,@array_data) = @_;	
	my @array = ();
	for ( my $i = $num_line_min ; $i < $num_line_max ; $i = $i + 1 ){
		$array[++$#array] = $array_data[$i];
	}
	# return array	
	return @array;
}
sub data_set{		#AKA the bas file reader
	my ($filename)=@_;
	
	my $matrix=();
	my $count_total=0;
	my $num_spaces=2;
	my $line_1;
	my $line_2;
	my $num_basin=();  #posicion inicio de info de basin en array
	my @array_data = read_file($filename);
	foreach my $value (@array_data){
		if ( $value =~/NUMBER/ && 
	        $value =~/OF/    && 
	        $value =~/BASINS/    && 
	        $value =~ /\:/ ){
				# array num basin
				$num_basin= (split(':', $value))[-1];
		}
		if ( $value =~/ATOM/ && 
	        $value =~/x/    && 
	        $value =~/y/    && 
	        $value =~ /z/ ){
				$line_1 = $count_total; 
		}      
		if ( $value =~/BASIN/    && 
			  $value =~/VOLUME/   && 
	        $value =~/INTEGRAL/ && 
	        $value =~ /MAXIMUM/ && 
	        $value =~/ATOM/ ){
				$line_2 = $count_total;
		}
		$count_total++;
	}
	
	my $expr_basin=($num_basin+$line_2+$num_spaces);
	my @basin_set=parsing_data(($line_2+$num_spaces),$expr_basin, @array_data);
	#all the info is separate in an specific array.
	my @array_VOLUME   = ();
	my @array_INTEGRAL = ();
	my @array_ATOMS    = ();	
	my @array_DIST     = ();
	my @array_MAXIMUM  = ();	
	my @array_axes_X   = ();	
	my @array_axes_Y   = ();	
	my @array_axes_Z   = ();	
	#
	my $i=-1;
	foreach my $basin (@basin_set){
		my $dataset=substr($basin,8,length($basin));
		my @array_tmp=split(' ',$dataset);
		my $set_count = scalar(grep {defined $_} @array_tmp); 

		if($set_count==9){
			$array_VOLUME[++$#array_VOLUME]     = $array_tmp[0];
			$array_INTEGRAL[++$#array_INTEGRAL] = $array_tmp[1];
			$array_MAXIMUM[++$#array_MAXIMUM]   = $array_tmp[2];
			$array_axes_X[++$#array_axes_X]     = convert_bohr_armstrong ($array_tmp[3]);
			$array_axes_Y[++$#array_axes_Y]     = convert_bohr_armstrong ($array_tmp[4]);
			$array_axes_Z[++$#array_axes_Z]     = convert_bohr_armstrong ($array_tmp[5]);
			$array_ATOMS[++$#array_ATOMS]       = $array_tmp[6];
			$array_DIST[++$#array_DIST]         = $array_tmp[7];
		}
	}
	#
	my @atom_set=parsing_data(($line_1+$num_spaces),$line_2,@array_data);
	#info to put
	my @array_atoms   = ();	
	my @array_coord_x = ();	
	my @array_coord_y = ();	
	my @array_coord_z = ();
	foreach my $atom(@atom_set){
		my @array_tmp2=split(' ',$atom);
		my $elements_count = scalar(grep {defined $_} @array_tmp2); 
		if($elements_count==4){
			$array_atoms[++$#array_atoms]     = parsing_name_atom($array_tmp2[0]);
			$array_coord_x[++$#array_coord_x] = convert_bohr_armstrong ($array_tmp2[1]);
			$array_coord_y[++$#array_coord_y] = convert_bohr_armstrong ($array_tmp2[2]);
			$array_coord_z[++$#array_coord_z] = convert_bohr_armstrong ($array_tmp2[3]);
		}
	}
	$num_atoms_global=scalar(@array_atoms);
	$num_basin_global=ltrim($num_basin);
	#array of arrays contruction.
	my @array_set_atoms  = ([@array_atoms],
                           [@array_coord_x],
                           [@array_coord_y],
                           [@array_coord_z]);
	my @array_set_basin  = ([@array_VOLUME],
                           [@array_INTEGRAL],
                           [@array_MAXIMUM],
                           [@array_axes_X],
                           [@array_axes_Y], 
                           [@array_axes_Z],
                           [@array_ATOMS],
                           [@array_DIST]);
	my @matrix = ([@array_set_atoms],
              [@array_set_basin]);
	# return array	
	return @matrix;
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
###########################################
sub cart_file_write{
	my ($name_without_extension,$num_atom, $num_basin, @matrix)=@_;
	open(FUK, ">$name_without_extension.cart");
	for(my $i=0; $i<$num_atom;$i++){
		print FUK "$matrix[0][0][$i]\t$matrix[0][1][$i]\t$matrix[0][2][$i]\t$matrix[0][3][$i]\n";
	}
	for(my $j=0; $j<$num_basin;$j++){
		print FUK "X\t$matrix[1][3][$j]\t$matrix[1][4][$j]\t$matrix[1][5][$j]\n"
	}
	close (FUK);
}
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
sub position_change_fragment{	#Here, magics happends. Essencially the fragments coordintaes are randomly rotated and traslatated. Then
	my(@FragLines)=@_;
	#rotate
	my @base_angles=gen_ptp();
	my $phi         = $base_angles[0];
	my $theta       = $base_angles[1];
	my $psi         = $base_angles[2];
	#translate
	my @base_xyz  = gen_xyz();
	my $base_x    = $base_xyz[0];
	my $base_y    = $base_xyz[1];
	my $base_z    = $base_xyz[2];
	# do the trig
	my $cos_phi     = sprintf '%.6f', cos($phi);
	my $cos_theta   = sprintf '%.6f', cos($theta);
	my $cos_psi     = sprintf '%.6f', cos($psi);
	my $sin_phi     = sprintf '%.6f', sin($phi);
	my $sin_theta   = sprintf '%.6f', sin($theta);
	my $sin_psi     = sprintf '%.6f', sin($psi);
	# make the rotation matrix
	my $D = new Math::Matrix ([$cos_phi,$sin_phi,0],[-$sin_phi,$cos_phi,0],[0,0,1]);
	my $C = new Math::Matrix ([1,0,0],[0,$cos_theta,$sin_theta],[0,-$sin_theta,$cos_theta]);
	my $B = new Math::Matrix ([$cos_psi,$sin_psi,0],[-$sin_psi,$cos_psi,0],[0,0,1]);
	my $A = $B->multiply($C)->multiply($D);
	#
	my @ar_x = ();
	my @ar_y = ();
	my @ar_z = ();
	my @coords=();
	#
	while (my $Fline = shift (@FragLines)) {
		my @Cartesians = split '\s+', $Fline;
		my ($Atom_label, @orig_xyz) = @Cartesians;
		my $matrix_xyz = new Math::Matrix ([$orig_xyz[0],$orig_xyz[1],$orig_xyz[2]]);
		my $trans_xyz = ($matrix_xyz->transpose);
		my $rotated_xyz = $A->multiply($trans_xyz);
		my @new_xyz = split '\n+',$rotated_xyz;
		push(@new_xyz,$Atom_label); #rotated fragment.
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
	my @bigarray=	([@coords],
								[@ar_x],
								[@ar_y],
								[@ar_z]);
	return @bigarray;
}
###################################
#equation.
sub equation{
	my ($multi, $dist)=@_;
	my @b_multi=@{$multi};
	my @b_dist=@{$dist};
	my $sum=0;
	for (my $i = 0; $i < scalar (@b_multi); $i++) {
			#print "f+ * f- = $b_multi[$i] con distancia de $b_dist[$i]\t\n";
			my $division = ( $b_multi[$i] / $b_dist[$i] );
			#print "f+ * f- = $b_multi[$i] con distancia de $b_dist[$i]\t Division: $division\n";
			$sum+=$division;
		}
	return $sum;
}

###################################
#write xyz file
sub write_xyz_file_center{
	my ($pid,$arr_center, $atom)=@_;
    my @arr = @{$arr_center};
    open(XYZ,">>tmp$pid.tmp");
    print XYZ "$atom\t$arr[0]\t$arr[1]\t$arr[2]\n";
    close (XYZ);
}
###################################
#write xyz file
sub write_xyz_file{
	#write the coordinate information for each tmp file
	my ($pid,  $Attract, @big_array )=@_;
    #
    my @tmpr = @{$Attract};
    #
	open(XYZ,">>tmp$pid.tmp");
	my @coords	=clear_data(@{$big_array[0]});
	for(my $i=0; $i<scalar(@coords); $i++){
		print XYZ "$coords[$i][0]\t$coords[$i][1]\t$coords[$i][2]\t$coords[$i][3]\n";
	}
#	for (my $i = 0; $i < scalar (@{$big_array[1]}); $i++) {
#	for (my $i = 0; $i < scalar (@tmpr); $i++) {
#		print XYZ "X\t$tmpr[$i]\n";
#	}
	close (XYZ);
}
###################################
#write xyz file data.
sub write_xyz_file_header{
	#write only the header for each tmp file
	my ($pid, $atoms_total, $sum)=@_;
	open(XYZ,">>tmp$pid.tmp");
	print XYZ "$atoms_total\n";
	print XYZ "E = $sum\n";
	close (XYZ);
}
###################################
sub sort_energies{
	#simple sort, uses hash for sorting
	my @files =glob "*.tmp";
	my ($na, $me)=@_;
	my @energies=();
	my @sorted_numbers;
	my %super;
	my  $natom;
	foreach my $file (@files){
		open (TMP, "<$file") || die "cannot open $file in sort subroutine\n";
		my @lines=<TMP>;
		close(TMP);
		$natom=$lines[0];
		my $energy=$lines[1];
		$energy=~s/(E = )?//g;
		chomp($energy);
		push(@energies, $energy);
		my $coords="";
		#
# Te comias el ultimo atractor por el -1        #foreach $a  (2..($#lines-1)){
        foreach $a  (2..($#lines)){
			$coords=$coords.$lines[$a];
		}
		$super{$energies[$#energies]}=$coords;	
		unlink "$file";
	}
	if($sort_option==1){
		@sorted_numbers = sort { $b <=> $a } @energies;
	}else{
		@sorted_numbers = sort { $a <=> $b } @energies;
	}
	open(NEW, ">$na$me.xyz");
	my $aux=0;
	#print "$aux and $XYZSize \n";
	foreach my $item (@sorted_numbers){
		if($aux<$XYZSize){
			print NEW "$natom", "E = $item\n",$super{$item};
		}else{
			last;
		}
		$aux++;
	}
	close(NEW);
#Osvaldo	my $root=dirname(abs_path($0));
	#print "ROOT: $root\n";
#Osvaldo	$Sim="$root/".$Sim;
	#`perl $Sim $na$me.xyz`;
#	unlink "$na$me.xyz";
}
#####################################
sub unify_tmp_files_noSort{
	my ($na, $me)=@_;
	my @files =glob "*.tmp";
	open(NEW, ">$na$me\_noSort.xyz");
	foreach my $file (@files){
		open (TMP, "<$file") || die "cannot open $file in sort subroutine\n";
		my @lines=<TMP>;
		close(TMP);
		print NEW @lines;
	}
	close(NEW);
}
####################################
sub resume_log_file{
	my ($na,$me, $timeS,$timeE,$timeT)=@_;
	
	open(LOG,">$na$me.txt") || die "cannot write $na$me.txt";
	print LOG "\t\t\t\t$program_name\n";
	#nombre del archivo
	print LOG "Output File:\t$na$me.xyz\n";
	#fragmntos
	print LOG "Fragments:\n\tf+:\t$na.bas\n\tf-:\t$me.bas\n";
	#tamaño de caja
	print LOG "Box size:\t $Box_x, $Box_y, $Box_z\n";
	#numero demuestra
	print LOG "\nProgram data:\nSample size:\t$Num_of_geometries\t";
	#procesadores
	print LOG "N° processes:\t$ncpus\n";
	#inicio, fin, tiempo total
	print LOG "\nStart: $timeS\nEnd: $timeE\nTotal Time: $timeT\n";
	
	close(LOG);	
}
sub show_help{
	print "\nTo use Rubik:\n\n";
	print "\tUse: \$ perl $0 [File Format] [f+ File] [f- File] [Number of Geometries] {OPTIONAL: Box size factor} {Sort mode}";
	print "\nWhere:\n[File Format]: -taff -> reads .tf files from TAFF program\n\t\t-bas -> reads .bas files\n";
	print "[f+ File]: File with the nucleophilic Fukui funcion information.\n";
	print "[f- File]: File with the electrophilic Fukui funcion information.\n" ;
	print "[Number of Geometries]: Number of geometries to assemble\n";
	print "{Box size Factor}: the box size in Rubik is the sum of the sides from both fragments multiply for a factor (Default 1)\n";
	print "{Sort mode}: By default the data is given sorted from best \"Maximum Matching\" to worst, using \"-s\" the data will be sorted less to more\n";
}
###################################
sub ConfigFileReader{
	open(CONFIG,"$ARGV[0]") or die "Configuration file doesn't exist\n";
	my @lines=<CONFIG>;
	close(CONFIG);
	my $aux=0;
	my ($type, $frag1, $frag2) ;
	foreach my $line(@lines){
		chomp($line);
		my $letter = substr($line, 0, 1);
		my ($key, $value)=split("=",$line);
		# Check repository 
		#1
		# the char # marks a commentary in Miconf.txt
		if($letter ne "#" && $line ne ""){
			$value=trim($value);
			if($key=~/NumberOfAssemblies/i){
				$Num_of_geometries=$value;
				$aux++;
			}
			if ($key=~/numberoffinalstructures/i) {
				$XYZSize=$value;
				$aux++;
			}
			if ($key=~/FileType/i) {
				$type=$value;
				$aux++;
			}
			if($key=~/fragment1/i){
				$frag1=$value;
				$aux++;
			}
			if($key=~/fragment2/i){
				$frag2=$value;
				$aux++;
			}
			if($key=~/boxsizefactor/i){
				my @box=split(",",$value);
				if(scalar(@box)==1){
					$Box_Multiplication_Factor=$value;
				}elsif(scalar(@box)==3){
					print "Box Size given by user\n";
					if(looks_like_number($box[0]) && looks_like_number($box[1]) && looks_like_number($box[2])){
						$BoxAutomatic=2;
						$Box_x=$box[0];
						$Box_y=$box[1];
						$Box_z=$box[2];
					}else{
						print "Box coordinates are not numeric\n";
						exit(1);
					}
				}else{
					print "Box Size error\nUsage:\nBoxSizeFactor=X,Y,Z\n";
					exit(1);
				}
			}
			if($key=~/sort/i){
				$sort_option=$value;
			}
		}
	}
	if($aux!=5){
		print "Wrong number of params in $ARGV[0]\naux = $aux\n";
		exit(1);
	}
	return ($type, $frag1, $frag2);
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
# MAIN

my @matrix_1;
my $atoms_num_1;
my $basin_num_1;
# extension
my $without_extension_1;
#
my @matrix_2   ;
my $atoms_num_2;
my $basin_num_2;
	# extension
my $without_extension_2 ;

print "                          ______        _     _ _     
                         (_____ \\      | |   (_) |    
                          _____) )_   _| |__  _| |  _ 
                         |  __  /| | | |  _ \\| | |_/ )
                         | |  \\ \\| |_| | |_) ) |  _ ( 
                         |_|   |_|____/|____/|_|_| \\_)
                      _______________________________________
                       A Stochastic Minimum Candidate Finder                                     
                      _______________________________________
";

my $tiempo_inicial = new Benchmark; #funcion para el tiempo de ejecucion del programa
my $datestringStart = localtime();
unlink glob "*.tmp";

my ($type, $file_name_1, $file_name_2)=ConfigFileReader();

if ( ($type eq "-bas") ){ 
	#if option para validadr archivos y formato.
	#$file_name_1 = $ARGV[1];
	@matrix_1    = data_set($file_name_1);
	$atoms_num_1 = $num_atoms_global;
	$basin_num_1 = $num_basin_global;
	# extension
	($without_extension_1 = $file_name_1) =~ s/\.[^.]+$//;
	#
	#$file_name_2 = $ARGV[2];
	@matrix_2    = data_set($file_name_2);
	$atoms_num_2 = $num_atoms_global;
	$basin_num_2 = $num_basin_global;
	# extension
	($without_extension_2 = $file_name_2) =~ s/\.[^.]+$//;

}elsif( ($type eq "-taff") ){
	#if option para validadr archivos y formato.
	#$file_name_1 = $ARGV[1];
	#
	@matrix_1    = TaffFileReader($file_name_1);
	$atoms_num_1 = $num_atoms_global;
	$basin_num_1 = $num_basin_global;
	# extension
	$file_name_1=(split("/",$file_name_1))[-1];
	( $without_extension_1 = $file_name_1) =~ s/\.[^.]+$//;
	#
	#$file_name_2 = $ARGV[2];
	#
	@matrix_2    = TaffFileReader($file_name_2);
	$atoms_num_2 = $num_atoms_global;
	$basin_num_2 = $num_basin_global;
	# extension
	$file_name_2=(split("/",$file_name_2))[-1];
	($without_extension_2 = $file_name_2) =~ s/\.[^.]+$//;
}else{
	print "Error en el numero de argumntos\n";
	show_help();
	exit(1);
}
	#CART FILE CREATION
	cart_file_write($without_extension_1,$atoms_num_1,$basin_num_1,@matrix_1);
	cart_file_write($without_extension_2,$atoms_num_2,$basin_num_2,@matrix_2);
	#NUMBER OF PROCESSORS AND GEOMETRY
	#chomp(my $cpu_count = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
	$ncpus=32*1;			##Si bien no encuentro un patron apartente, este parecice ser la mejro combnacion.
									## con 32 es mejor para 10.000 combinaciones, pero no el mejor para 1000.
	#READING CART FILE (REDUNDANCY)
	open(FRAG1,"$without_extension_1.cart" ) or die "Unable to open fragment file: $without_extension_1.cart";
	open(FRAG2,"$without_extension_2.cart" ) or die "Unable to open fragment file: $without_extension_2.cart";
	my @f1=<FRAG1>;
	my @f2=<FRAG2>;
	close(FRAG1);
	close(FRAG2);
	unlink glob "*.cart";

	my $pm= new Parallel::ForkManager($ncpus);
	if($BoxAutomatic==1){
		BoxCoordinates(\@f1,\@f2);
	}
    #
    #
	# The integral value for each attractor doesn't change moving the fragments.
	# So the multiplication between atractor (integral) is only done once.
    my @multi=basin_integral_multi(\@matrix_1, \@matrix_2, $basin_num_1, $basin_num_2);
    #
	foreach my $iteration(1 .. $Num_of_geometries){

		$pm->start($iteration) and next;
		srand();			# all children process havee their own random.
		my $decision=1;		# to verify steric impediment
        #
		my @resume1;
		my @resume2;
        #
		while($decision!=0){
			#randomize atoms/basin positions
			@resume1=position_change_fragment(@f1);
			@resume2=position_change_fragment(@f2);
			#get atoms coordinates
			my @coords_1 = clear_data(@{$resume1[0]});
			my @coords_2 = clear_data(@{$resume2[0]});
			# Verify for steric impediment, 1 yes, 0 no;             
			$decision = steric_impediment(\@coords_1, \@coords_2);
		}
        # # # # # # # #
        # Inicio colocar burbuja 
        # Attractors X=@{$resume1[1]}, Y=@{$resume1[2]},Z=@{$resume1[3]}
        # measure center of mass 
        my @array_center_mass_attrac1 = measure_center(\@{$resume1[1]},\@{$resume1[2]},\@{$resume1[3]});
        my @array_center_mass_attrac2 = measure_center(\@{$resume2[1]},\@{$resume2[2]},\@{$resume2[3]});
        #
        my @medium = medium_point (\@array_center_mass_attrac1,\@array_center_mass_attrac2);
        # Radius sphere
        my @ValuesTam   = sort {$a cmp $b} ($Tam1,$Tam2);
        #
        my $weight      = 1.6;
        my $Raddi_Value = ($ValuesTam[0]/2) * $weight;
        #my $Raddi_Value = 2;
        #
        my $Radii_Xmin = ($medium[0] + (-$Raddi_Value));
        my $Radii_Ymin = ($medium[1] + (-$Raddi_Value));
        my $Radii_Zmin = ($medium[2] + (-$Raddi_Value));
        my @Ra_min = ($Radii_Xmin,$Radii_Ymin,$Radii_Zmin);
        #
        my $Radii_Xmax = ($medium[0] + $Raddi_Value);
        my $Radii_Ymax = ($medium[1] + $Raddi_Value);
        my $Radii_Zmax = ($medium[2] + $Raddi_Value);
        my @Ra_max = ($Radii_Xmax,$Radii_Ymax,$Radii_Zmax);
        #    
        # Distance medition between atractor is done.
        my @dist = distance_medition(\@resume1,\@resume2);
        # Hash for mult
        for ( my $i =0; $i < scalar (@multi); $i=$i+1 ){
            my $key = $i;
            $HashMult{$key}=$multi[$i];
        }
        #
        my @Mult_1 = ();
        my @Dist_1 = ();
        my @Attracs_1 = ();
        my @Attracs_2 = ();
        #A point (x,y,z) is inside the sphere with center (cx,cy,cz) and radius r if
        #( x-cx )^2 + (y-cy)^2 + (z-cz)^2 < r^2     
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
        # Fin colocar burbuja
#        print Dumper(\%HashMult);
#        print Dumper(\%HashDist);
#        print Dumper(\%HashCoordsAttra);
        # Maximum Matching value
        #
        my $sum = equation(\@Mult_1,\@Dist_1);        
        # Write tmp file with Maximum Matching value and coordinates
        write_xyz_file_header($iteration,($atoms_num_1+$atoms_num_2),$sum);
        write_xyz_file($iteration,\@Attracs_1,@resume1);
        write_xyz_file($iteration,\@Attracs_2,@resume2);
        #
        #
    #    write_xyz_file_center ($iteration,\@array_center_mass_attrac1, "He");
    #    write_xyz_file_center ($iteration,\@array_center_mass_attrac2, "He");
        #
    #    write_xyz_file_center ($iteration,\@medium, "Rn");
        #
    #    write_xyz_file_center ($iteration,\@Ra_min, "Xe");
    #    write_xyz_file_center ($iteration,\@Ra_max, "Xe");
        # Delete data hash
        %HashCoordsAttra = ();
        %HashDist        = ();
        %HashMult        = ();
        #
		$pm->finish;
	}
	#paralel
	$pm->wait_all_children;
	######################################################
	#sort by E.
	#unify_tmp_files_noSort($without_extension_1, $without_extension_2);
	sort_energies($without_extension_1, $without_extension_2);
	# Time in console is printed
	my $tiempo_final = new Benchmark;
	my $datestringEnd = localtime();
	my $tiempo_total = timediff($tiempo_final, $tiempo_inicial);
	print "\n\tTiempo de ejecucion: ",timestr($tiempo_total),"\n";
	print "\n";
	# Log File is written
	resume_log_file($without_extension_1,$without_extension_2,$datestringStart,$datestringEnd,timestr($tiempo_total));




