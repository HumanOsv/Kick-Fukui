#!/usr/bin/perl -s

sub read_cart_file{
	my ($input_file)=@_;
	my @data=();
	open(FILE, "<", $input_file) or die "Can't open";
	my @lines=<FILE>;
	close(FILE);
	foreach $i (@lines){
		chomp($i);
		$data[++$#data]=$i;
	}
	return @data;
}

my ($file,$energy) = @ARGV;
if (not defined $file) {
	die "\nObtener Basines:\n\nUsage:\n\tperl 01Basin_Orden.pl [file.txt]\n";
	exit(1);  
}
#read and parse
my @data   = read_cart_file($file); 
my $word_1 = "BASIN";
my $word_2 = "VOLUME";
my $word_3 = "INTEGRAL";
my $word_4 = "VALUE";
#
my $word_5 = "TOTAL";
#
my @array_energy = ();
my @linearray_1  = ();
my @linearray_2  = ();
my $count = 0;
#
foreach my $i (@data) {
	if ( ($i =~ m/$word_1/) && ($i =~ m/$word_2/) && ($i =~ m/$word_3/) && ($i =~ m/$word_4/) ) {
		my $lines = $count;
		push (@linearray_1,$lines);
	}
	if ( ($i =~ m/$word_5/) ) {
		my $lines = $count;
		push (@linearray_2,$lines);
	} 
	$count++;
}

#	my @array_tabs = ();
#	@array_tabs    = split ('\s+',$i);
#	push (@linearray_1,$lines);
#
#push (@array_energy,$array_tabs[0]);
#
my @new_data = ();
for my $x ( $linearray_1[$#linearray_1] .. $linearray_2[$#linearray_2] ) {
	push (@new_data,$data[$x]);
}
#
shift (@new_data);
shift (@new_data);
#
pop (@new_data);
pop (@new_data);
#
print "Numero de Basines\n";
foreach my $i (@new_data){
	my @array_tabs = ();
	@array_tabs    = split ('\s+',$i);
	print "$array_tabs[0]\n";
}
print "Coordenadas para Monte-Fukui\n\n";
foreach my $i (@new_data){
	my @array_tabs = ();
	@array_tabs    = split ('\s+',$i);
	#
	my $ma_x  = sprintf '%.6f',$array_tabs[3];
	my $ma_y  = sprintf '%.6f',$array_tabs[4];
	my $ma_z  = sprintf '%.6f',$array_tabs[5];
	my $basin = sprintf '%.6f',$array_tabs[2];
	print "X  $ma_x  $ma_y  $ma_z  $basin\n";
}
print "\n\n";
