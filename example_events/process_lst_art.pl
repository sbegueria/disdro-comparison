use Time::Local;
use POSIX 'strftime';

# https://www.researchgate.net/figure/233065459_fig2_Fig-3-Empirical-equation-to-estimate-the-raindrop-fall-velocity-in-still-air-using
sub terminal_velocity
{
  my ($d) = @_;
  return 9.5 * (1.0 - exp( - 0.54 * ( $d ** 1.13 ) ));
}

$argc = $#ARGV + 1;
if(!($argc == 2 || $argc == 3))
{
	print "Usage: perl process.pl folder outfile (middle|uniform|linear)\n";
	exit;
}

$folder  = $ARGV[0];
$outfile = $ARGV[1];
$prob = lc($ARGV[2]);

if(!($prob eq 'middle' || $prob eq 'uniform' || $prob eq 'linear'))
{
	$prob = 'uniform';
}

# Read directory tree
@fs = ();
@ys = sort <$folder/*>;
#@ys=("$folder/2015");
foreach $y (@ys){
	if(-d $y && $y =~ /\/\d+$/){
		@ms = sort <$y/*>;
		#@ms=("$folder/2015/05");
		foreach $m (@ms){
			if(-d $m && $m =~ /\/\d+$/){
				@ds = sort <$m/*>;
				#@ds=("$folder/2015/05/18");
				foreach $d (@ds){
					if(-d $d && $d =~ /\/\d+$/){
						@ps = sort <$d/*\.*>;
						foreach $p (@ps){
							if(-f $p && $p =~ /(\d\d\d\d)(\d\d)(\d\d)\d+\./){
								push(@fs,$p);
							}
						}
					}
				}
			}
		}
	}
}

# Test disdrometer model
$f = $fs[0];
print $f."\n";
if( $f =~ /\.gz$/ ){
	$x = `gunzip -c $f`;
	$x =~ s/[\r\n]+/\n/g;
                @ls = split(/(?<=\n)/, $x);
}
else{
	open F, $f;
	@ls = <F>;
	close F;
}
@x = split(/\;/, $ls[1]);

if($x[1] =~ /^(\d\d\d\d)$/)
{
	$iserial = $1;
	print "Thies disdrometer ($iserial)\n";
	$size{'0436'} = 46.65314E-04;
	$size{'0655'} = 49.04051E-04;
	$alpha = $size{$iserial};
	if($alpha eq "")
	{
		print "Unknown sensor area\n";
		exit;
	}
	$disdro = 0;
}
else
{
	if($x[12] =~ /^(\d\d\d\d\d\d)$/)
	{
		$iserial = $1;
		print "Parsivel disdrometer ($iserial)\n";
		$alpha = 54.0E-04; # superficie de medida (teórica), m2
		$disdro = 1;
	}
	else
	{
		print "Unknown disdrometer\n";
		exit;
	}
}

# Date interval
if($fs[0] =~ /(\d\d\d\d)(\d\d)(\d\d)\d+\./){
	$date_i = "$1/$2/$3 00:00";
	$secs_i = secs($date_i);
	print "Dates from $1/$2/$3";
}
if($fs[$#fs] =~ /(\d\d\d\d)(\d\d)(\d\d)\d+\./){
	$date_f = "$1/$2/$3 23:59";
	$secs_f = secs($date_f);
	print " to $1/$2/$3\n";
}

# Thies disdrometer
if($disdro == 0)
{
	@d1 = (0.125, 0.25, 0.375, 0.5, 0.750, 1, 1.250, 1.5, 1.75,
			2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8);
	@d2 = (0.25, 0.375, 0.5, 0.750, 1, 1.250, 1.5, 1.75, 2, 2.5, 
			3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5);
	for($i=0;$i<22;$i++){
		$dmed[$i]=0.5*($d1[$i]+$d2[$i]);
		$dmed3[$i]=$dmed[$i]*$dmed[$i]*$dmed[$i];
		$rmed3[$i]=$dmed3[$i]/8.0;
		$dmed6[$i]=$dmed3[$i]*$dmed3[$i];
	}
	@v1 = (0, 0.2, 0.4, 0.6, 0.8, 1, 1.4, 1.8, 2.2, 2.6, 3, 3.4,
			4.2, 5, 5.8, 6.6, 7.4, 8.2, 9, 10);
	@v2 = (0.2, 0.4, 0.6, 0.8, 1, 1.4, 1.8, 2.2, 2.6, 3, 3.4,
			4.2, 5, 5.8, 6.6, 7.4, 8.2, 9, 10, 11);
	for($i=0;$i<20;$i++){
		$vmed[$i]=0.5*($v1[$i]+$v2[$i]);
		$vmed2[$i]=$vmed[$i]*$vmed[$i];
	}
	for($i=0;$i<20;$i++){
		$w1med[$i] = 4.85*$dmed[$i]*exp(-0.195*$dmed[$i]);
		$w1med2[$i] = $w1med[$i]*$w1med[$i];
	}
	for($i=0;$i<20;$i++){
		$w2med[$i] = 0.051*$dmed[$i]*$dmed[$i]*$dmed[$i]-
			0.912*$dmed[$i]*$dmed[$i]+5.038*$dmed[$i]-0.254;
		$w2med2[$i] = $w2med[$i]*$w2med[$i];
	}
	for($i=0;$i<20;$i++){
		$w3med[$i] = 9.65-10.3*exp(-0.6*$dmed[$i]);
		$w3med2[$i] = $w3med[$i]*$w3med[$i];
	}

	# Constants
	$pi = 3.141592653589793;
	$tau = 60.0; # sampling interval in seconds
	$ta = $tau*$alpha;

	open F, 'filter_thies.csv';
	@ls = <F>;
	close F;
	shift @ls;
	@mask = ();
	foreach $l (@ls)
	{
		$l =~ s/[\r\n\"]+//g;
		$l =~ s/^.*?\,//g;
		$l =~ s/true/1/ig;
		$l =~ s/false/0/ig;
		@xs = split(/\s*\,\s*/, $l);
		foreach $x (@xs)
		{
			push(@mask, $x);
		}
	}
	if(int($#mask+1) != 22*20 )
	{
		print "Using default mask\n";
		@mask = ();
		for($i=0; $i<22*20; $i++)
		{
			push(@mask, 1);
		}
	}
	else
	{
		print "Mask read from filter_thies.csv\n";
	}

	open F, 'margin_thies.csv';
	@ls = <F>;
	close F;
	$ls[0] =~ s/[\r\n\"]+//g;
	@xs = split(/\s*\,\s*/, $ls[0]);
	@margin = ();
	foreach $x (@xs)
	{
		push(@margin, $x);
	}
	if(int($#margin+1) != 22 )
	{
		print "Using default margin factors\n";
		@margin = ();
		for($i=0; $i<22; $i++)
		{
			push(@margin, 1.0);
		}
	}
	else
	{
		print "Margin factors read from margin_thies.csv\n";
	}

	# Procesado de los archivos
	$chunk = "";
	foreach $f (@fs){
		if($f =~ /(\d\d\d\d)(\d\d)(\d\d)(\d\d)\./){
			$dy = $1;
			$dm = $2;
			$dd = $3;
			$dh = $4;
			if( $f =~ /\.gz$/ ){
				$x = `gunzip -c $f`;
				$x =~ s/[\r\n]+/\n/g;
		                    @ls = split(/(?<=\n)/, $x);
			}
			else{
				open F, $f;
				@ls = <F>;
				close F;
			}
			if( $chunk ne "" ){
				@bs = split(/\;/, $ls[0]);
				if($#bs != 524 && $#bs != 525)
				{
					$ls[0] = $chunk . $ls[0];
					@bs = split(/\;/, $ls[0]);
				}
			}
			@bs = split(/\;/, $ls[$#ls]);
			if( $#bs != 524 && $#bs != 525 ){
				$chunk = pop(@ls);
			}
			else{
				$chunk = "";
			}
			@d = split(/\;/, $ls[0]);
			if($d[4] =~ /\d+\:(\d+)\:\d+/){
				$imin = $1;
			}

			$line = 0;
			foreach $l (@ls){
				$l =~ s/[\r\n]+//g;
				if( $l ne '' ){

					@d = split(/\;/, $l);

					$id = $d[0];
					$serial = $d[1];

					if($d[4] =~ /\d+\:(\d+)\:\d+/){
						$amin = $1;
						if($amin>=$imin){
							$amin=$amin-$imin;
						}
						else{
							$amin=$amin+60-$imin;
						}
						if($amin<10){$amin='0'.$amin;}
					}
					$time = "$dy/$dm/$dd $dh:$amin:00";

					#print "$time $#d\n";

					# Fix time based on Artila ntp
					$shift = 0;
					if( $d[524] =~ /^\d\d\d\d\-\d\d\-\d\d \d\d\:\d\d\:\d\d$/ )
					{
						# Round to nearest minute
						$shift = secs2($d[524]);
						$time = date2(int($shift/60+0.5)*60);
						$shift = $shift - int($shift/60+0.5)*60;
						#print $d[524]."\t".$time."\n";
					}
					$dt = $tau;
					$qual = $d[18];
					$synop = $d[9]; # Synop code 4677
					$metar = $d[11];

					# $d[12]  R (rain intensity, mm h-1)
					# $d[16]  MOR visibility, m
					# $d[17]  Z (radar reflectivity, dB mm6 m-3)
					# $d[18]  Measurement quality (%)
					# $d[49]  M (number of particles detected, -)

					if($d[519] ne '')
					{
						# $d[519] T (air temperature, ºC)
						$d[519] =~ s/^99999$/0/g;
					}
				
					if($d[520] ne '')
					{
						# $d[520] RH (relative air humidity, %)
						$d[520] =~ s/^99999$/0/g;
					}

					if($d[521] ne '')
					{
						# $d[521] W (wind velocity, m/s)
						$d[521] =~ s/^9999$/0/g;
					}

					if($d[522] ne '')
					{
						# $d[522] WD (wind direction, º)
						$d[522] =~ s/^999$/0/g;
					}

					$error = 0;
					for($i=22; $i<=34; $i++)
					{
						if($d[$i] == 1)
						{
							$error++;
							$ierror = $i;
						}
					}
					if($error == 1)
					{
						$error = 0;
						#$error = $ierror + 2; # TODO Descomentar
					}
					elsif($error > 1)
					{
						$error = 0;
						#$error = 37; # TODO Descomentar
					}

					if($error != 0)
					{
						$data{$time} = "err$error,$#d";
					}
					elsif($synop !~ /^[\d\+\-\.]+$/){
						$data{$time} = "err3,$#d";
					}
					elsif($d[12] =~ /[^\d\+\-\.eE]/){
						$data{$time} = "err4,$#d";
					}
					elsif($d[12] =~ /9999\.999/){
						$data{$time} = "err5,$#d";
					}
					elsif($l =~ /(OK|Version)/){
						$data{$time} = "err6,$#d";
					}
					elsif($l !~ /^[\w\+\-\.\:\;\r\n\s\t ]+$/i || $l =~ /^\;/){
						$data{$time} = "err7,$#d";
					}
					#elsif( $synop != 0 ){
					elsif ($serial ne "" && $serial ne $iserial)
                                        {
						#print "Telegram with wrong ID ($serial) found in $f\n$l\n";
						$data{$time} = "err7,$#d";
                                        }
					elsif ($#d < 523 || $#d > 525)
					{
						$data{$time} = "err7,$#d";
					}
					else{
						$satura = 0;

						for($i=0;$i<22;$i++){$dsdd[$i]=0.0;}
						for($j=0;$j<20;$j++){$dsdv[$j]=0.0;}
						# Lectura matrix DSD (velocidades x diametros)
						$k=0;
						for($i=0;$i<22;$i++){
							for($j=0;$j<20;$j++){
								$dsd[$j][$i] = $d[79+$k] * $mask[$k] * $margin[$i];
								$dsdv[$j] = $dsdv[$j] + $d[79+$k] * $mask[$k] * $margin[$i];
								$dsdd[$i] = $dsdd[$i] + $d[79+$k] * $mask[$k] * $margin[$i];
		
								if($dsd[$j][$i] eq '999'){
									$satura=1;
								}
								$k++;
							}
						}

						# Cálculo del número de partículas N (-)
						$k = 0;
						for($j=0;$j<20;$j++){
							for($i=0;$i<22;$i++){
								$k = $k + $dsd[$j][$i];
							}
						}
						$np = int($k);

						# Cálculo R (rain intensity, mm h-1)
						$r=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i];
							}
							$r += $x*$dmed3[$i];
						}
						$r = 6E-4*$pi/$ta*$r;

						# Cálculo P (M precipitation amount, mm)
						$p = $r/60.0;

						# Cálculo vector ND (number density, m-3 mm-1)
						@nd = ();
						@nd2 = ();
						for($i=0;$i<22;$i++){
							$nd[$i] = 0.0;
							for($j=0;$j<20;$j++){
								$nd[$i] += $dsd[$j][$i]/$vmed[$j]/$ta;
							}
							if($p>0.0){
								$nd[$i] = $nd[$i] / $p;
							}
							$nd2[$i] = 'nd'.($i+1);
						}

						# Cálculo M (water content, g m3)
						$m=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]/$vmed[$j];
							}
							$m += $x*$dmed3[$i];
						}
						$m = 1E-3*$pi/6.0/$ta*$m;

						# Cálculo Z (radar reflectivity, dB mm6 m-3)
						$z=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]/$vmed[$j];
							}
							$z += $x*$dmed6[$i];
						}
						$z = 10*log10(1/$ta*$z);

						# Cálculo E (kinetic energy, J m-2 mm-1)
						$e=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]*$vmed2[$j];
							}
							$e += $x*$dmed3[$i];
						}
						if($p>0.0){
							$e = 1.0/$alpha/$p/12.0*$pi*1.E-6*$e;
						}
						else{
							$e = 0.0;
						}

						# Cálculo ET (kinetic energy ac. to Uplinger, J m-2 mm-1)
						$et1=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]*$w1med2[$i];
							}
							$et1 += $x*$dmed3[$i];
						}
						if($p>0.0){
							$et1 = 1.0/$alpha/$p/12.0*$pi*1.E-6*$et1;
						}else{
							$et1 = 0.0;
						}

						# Cálculo ET (kinetic energy ac. to Van Dijk, J m-2 mm-1)
						$et2=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]*$w2med2[$i];
							}
							$et2 += $x*$dmed3[$i];
						}
						if($p>0.0){
							$et2 = 1.0/$alpha/$p/12.0*$pi*1.E-6*$et2;
						}
						else{
							$et2 = 0.0;
						}

						# Cálculo ET (kinetic energy ac. to Atlas, J m-2 mm-1)
						$et3=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]*$w3med2[$i];
							}
							$et3 += $x*$dmed3[$i];
						}
						if($p>0.0){
							$et3 = 1.0/$alpha/$p/12.0*$pi*1.E-6*$et3;
						}
						else{
							$et3 = 0.0;
						}

						# Cálculo V (MOR visibility, m)
						$v=0.0;
						for($i=0;$i<22;$i++){
							$x=0.0;
							for($j=0;$j<20;$j++){
								$x += $dsd[$j][$i]/$vmed[$j];
							}
							$v += $x*($dmed[$i]/1000)*($dmed[$i]/1000);
						}
						if($v>0.0){
							$v = 3/(($pi/2)*(1/$ta)*$v);
						}
						else{
							$v = 0.0;
						}

						$d1 = ''; # D at 10%
						$d2 = ''; # D at 25%
						$d3 = ''; # D at 50%
						$d4 = ''; # D at 75%
						$d5 = ''; # D at 90%
						$d6 = 0.0; # Avergage D

						@list = ();

						if($prob eq 'middle')
						{
							$x = 0.0;
							for($i=0;$i<22;$i++){
								for($j=0;$j<20;$j++){
									$d6 += $dsd[$j][$i]*$dmed[$i];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $dmed[$i]);
									}
								}
							}
						}
						elsif($prob eq 'uniform')
						{
							$x = 0.0;
							for($i=0;$i<22;$i++){
								for($j=0;$j<20;$j++){
									$d6 += $dsd[$j][$i]*$dmed[$i];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $d1[$i] + rand() * ($d2[$i] - $d1[$i]));
									}
								}
							}
						}
						elsif($prob eq 'linear')
						{
							$x = 0.0;
							for($i=0;$i<22;$i++){

								if($i==0)
								{
									$m = $dsdd[$i+1] - $dsdd[$i];
								}
								elsif($i==21)
								{
									$m = $dsdd[$i] - $dsdd[$i-1];
								}
								else
								{
									$m = ($dsdd[$i+1] - $dsdd[$i-1]) * 0.5;
								}

								for($j=0;$j<20;$j++){
									$d6 += $dsd[$j][$i]*$dmed[$i];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										$done = 0;
										while($done == 0)
										{
											$x1 = rand();
											$x2 = rand();
											$p1 = $dsdd[$i] - $m * 0.5 + $x1 * $m;
											$p2 = $x2 * ($dsdd[$i] + $m * 0.5);
											if($p2 < $p1)
											{
												push(@list, $d1[$i] + $x1 * ($d2[$i] - $d1[$i]));
												$done = 1;
											}
										}
									}
								}
							}
						}

						@list = sort { $a <=> $b } @list;

						if($x>0.0)
						{
							$d6 = $d6 / $x;

							$d1 = $list[int(0.10*$x)];
							$d2 = $list[int(0.25*$x)];
							$d3 = $list[int(0.50*$x)];
							$d4 = $list[int(0.75*$x)];
							$d5 = $list[int(0.90*$x)];
						}
						else
						{
							$d1 = 0.0;
							$d2 = 0.0;
							$d3 = 0.0;
							$d4 = 0.0;
							$d5 = 0.0;
							$d6 = 0.0;
						}

						$v1 = ''; # V at 10%
						$v2 = ''; # V at 25%
						$v3 = ''; # V at 50%
						$v4 = ''; # V at 75%
						$v5 = ''; # V at 90%
						$v6 = 0.0; # Average V

						@list = ();

						if($prob eq 'middle')
						{
							$x = 0.0;
							for($j=0;$j<20;$j++){
								for($i=0;$i<22;$i++){
									$v6 += $dsd[$j][$i]*$vmed[$j];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $vmed[$j]);
									}
								}
							}
						}
						elsif($prob eq 'uniform')
						{
							$x = 0.0;
							for($j=0;$j<20;$j++){
								for($i=0;$i<22;$i++){
									$v6 += $dsd[$j][$i]*$vmed[$j];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $v1[$j] + rand() * ($v2[$j] - $v1[$j]));
									}
								}
							}
						}
						elsif($prob eq 'linear')
						{
							$x = 0.0;
							for($j=0;$j<20;$j++){

								if($j==0)
								{
									$m = $dsdv[$j+1] - $dsdv[$j];
								}
								elsif($j==19)
								{
									$m = $dsdv[$j] - $dsdv[$j-1];
								}
								else
								{
									$m = ($dsdv[$j+1] - $dsdv[$j-1]) * 0.5;
								}

								for($i=0;$i<22;$i++){
									$v6 += $dsd[$j][$i]*$vmed[$j];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										$done = 0;
										while($done == 0)
										{
											$x1 = rand();
											$x2 = rand();
											$p1 = $dsdv[$j] - $m * 0.5 + $x1 * $m;
											$p2 = $x2 * ($dsdv[$j] + $m * 0.5);
											if($p2 < $p1)
											{
												push(@list, $v1[$j] + $x1 * ($v2[$j] - $v1[$j]));
												$done = 1;
											}
										}
									}
								}
							}
						}

						@list = sort { $a <=> $b } @list;

						if($x>0.0){
							$v6 = $v6 / $x;

							$v1 = $list[int(0.10*$x)];
							$v2 = $list[int(0.25*$x)];
							$v3 = $list[int(0.50*$x)];
							$v4 = $list[int(0.75*$x)];
							$v5 = $list[int(0.90*$x)];
						}
						else
						{
							$v1 = 0.0;
							$v2 = 0.0;
							$v3 = 0.0;
							$v4 = 0.0;
							$v5 = 0.0;
							$v6 = 0.0;
						}

						#print "$time,";
						#print secs($time);
						$r = sprintf("%.3e", $r);
						$p = sprintf("%.3e", $p);
						$m = sprintf("%.3e", $m);
						$z = sprintf("%.3e", $z);
						$e = sprintf("%.3e", $e);
						$et1 = sprintf("%.3e", $et1);
						$et2 = sprintf("%.3e", $et2);
						$et3 = sprintf("%.3e", $et3);
						$v = sprintf("%.3e", $v);

						$d1 = sprintf("%.3f", $d1);
						$d2 = sprintf("%.3f", $d2);
						$d3 = sprintf("%.3f", $d3);
						$d4 = sprintf("%.3f", $d4);
						$d5 = sprintf("%.3f", $d5);
						$d6 = sprintf("%.3f", $d6);

						$v1 = sprintf("%.3f", $v1);
						$v2 = sprintf("%.3f", $v2);
						$v3 = sprintf("%.3f", $v3);
						$v4 = sprintf("%.3f", $v4);
						$v5 = sprintf("%.3f", $v5);
						$v6 = sprintf("%.3f", $v6);

						$data{$time} = "$synop,$r,$p,$m,$z,$e,$v,$d[12],$d[17],0,$d[16],$d[18],$d[519],$d[520],$d[521],$d[522],$d[49],$np,$d[38],$d[40],$d[41],$d[36],$d1,$d2,$d3,$d4,$d5,$d6,$v1,$v2,$v3,$v4,$v5,$v6,".join(',', @nd).",$shift,$line,$satura,$#d"; 

					}
					#else{
					#	$data{$time} = "$synop,0,0,0,0,0,0,0,0,0,0,0,0,0,$d[519],$d[520],$d[521],$d[522],0"; 
					#}

				}

				$line++;
			}
		}
	}

	open F, '>'.$outfile;

	print F "'type','serial','time','seconds','synop','r','p','m','z','e','mor','r_meas',";
	print F "'z_meas','e_meas','mor_meas','qual','tmp','rh','w','wd','np_meas','np','lcurrent','ocontrol','power',";
	print F "'tmp_int','d10','d25','d50','d75','d90','dmean','v10','v25','v50','v75','v90','vmean','".join('\',\'', @nd2)."','t_shift','now','err','ncol'\n";

	# Recorremos fechas
	for($i=$secs_i;$i<=$secs_f;$i+=60){
			$x = get_date($i);
			# Comprobamos si hay dato para esta fecha
			if($data{$x} eq ''){
				$err = 1;
				$num = 0;
				#print F "Thi,$iserial,$x,$i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,$err,$num\n";
				print F "Thi,$iserial,$x,$i,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,";
				print F "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,$err,$num\n";
			}
			elsif($data{$x} =~ /^err(\d+)\,([^\,]+)/)
			{
				$err = $1;
				$num = $2+1;
				#print F "Thi,$iserial,$x,$i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,$err,$num\n";
				print F "Thi,$iserial,$x,$i,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,";
				print F "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,$err,$num\n";
			}
			else{
				# Si todo ha ido bien ponemos un valor de 0 en el error status
				if($data{$x} =~ /\,0\,([^\,]+)$/){
					$num = $1+1;
					$data{$x} =~ s/\,[^\,]+$//g;
					$data{$x} =~ s/\,[^\,]+$//g;

					$data{$x} =~ s/\,/></g;
					$data{$x} =~ s/<[0\.e\+\-]+>/<0>/g;
					$data{$x} =~ s/^[0\.e\+\-]+>/0>/g;
					$data{$x} =~ s/<[0\.e\+\-]+$/<0/g;
					$data{$x} =~ s/></\,/g;

					print F "Thi,$iserial,$x,$i,".$data{$x}.",0,$num\n";
				}
				# Si se ha producido saturación en la matrix DSD ponemos un valor de 2 en el error status
				elsif($data{$x} =~ /\,1\,([^\,]+)$/){
					$num = $1+1;
					$data{$x} =~ s/\,[^\,]+$//g;
					$data{$x} =~ s/\,[^\,]+$//g;

					$data{$x} =~ s/\,/></g;
					$data{$x} =~ s/<[0\.e\+\-]+>/<0>/g;
					$data{$x} =~ s/^[0\.e\+\-]+>/0>/g;
					$data{$x} =~ s/<[0\.e\+\-]+$/<0/g;
					$data{$x} =~ s/></\,/g;

					print F "Thi,$iserial,$x,$i,".$data{$x}.",2,$num\n";
				}
			}
	}
	close F;
}

# Parsivel disdrometer
if($disdro == 1)
{
	# Diámetros y velocidades:
	@dmed = (0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,
		1.375,1.625,1.875,2.125,2.375,2.750,3.250,3.750,4.250,4.750,5.500,
		6.50,7.50,8.50,9.50,11.00,13.00,15.00,17.00,19.00,21.50,24.50);
	@vmed = (0.050,0.150,0.250,0.350,0.450,0.550,0.650,0.750,0.850,0.950,
		1.100,1.300,1.500,1.700,1.900,2.200,2.600,3.000,3.400,3.800,4.400,
		5.200,6.000,6.800,7.600,8.80,10.40,12.00,13.60,15.20,17.60,20.80);
	for($i=0;$i<32;$i++){
		$dmed3[$i]=$dmed[$i]*$dmed[$i]*$dmed[$i];
		$rmed3[$i]=$dmed3[$i]/8.0;
		$dmed6[$i]=$dmed3[$i]*$dmed3[$i];
	}
	for($i=0;$i<32;$i++){
		$vmed2[$i]=$vmed[$i]*$vmed[$i];
	}

	@d1 = (0.00, 0.12, 0.25, 0.37, 0.50, 0.62, 0.75, 
			0.87, 1.00, 1.12, 1.28, 1.50, 1.75, 2.00, 
			2.25, 2.56, 3.00, 3.50, 4.00, 4.50, 5.13, 
			6.00, 7.00, 8.00, 9.00, 10.25, 12.00, 
			14.00, 16.00, 18.00, 20.25, 23.00);
	@d2 = (0.12, 0.25, 0.37, 0.50, 0.62, 0.75, 0.87, 
			1.00, 1.12, 1.28, 1.50, 1.75, 2.00, 2.25, 
			2.56, 3.00, 3.50, 4.00, 4.50, 5.13, 6.00, 
			7.00, 8.00, 9.00, 10.25, 12.00, 14.00, 
			16.00, 18.00, 20.25, 23.00, 25.75);

	@v1 = (0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 
			0.70, 0.80, 0.90, 1.03, 1.20, 1.40, 1.60, 
			1.80, 2.05, 2.40, 2.80, 3.20, 3.60, 4.10, 
			4.80, 5.60, 6.40, 7.20, 8.20, 9.60, 11.20, 
			12.80, 14.40, 16.40, 19.20);
	@v2 = (0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 
			0.80, 0.90, 1.03, 1.20, 1.40, 1.60, 1.80, 
			2.05, 2.40, 2.80, 3.20, 3.60, 4.10, 4.80, 
			5.60, 6.40, 7.20, 8.20, 9.60, 11.20, 12.80, 
			14.40, 16.40, 19.20, 22.00);

	# Constantes
	$pi = 3.141592653589793;
	$tau = 60.0; # sampling time interval, seconds
	$ta = $tau*$alpha;

	open F, 'filter_parsivel.csv';
	@ls = <F>;
	close F;
	shift @ls;
	@mask = ();
	foreach $l (@ls)
	{
		$l =~ s/[\r\n\"]+//g;
		$l =~ s/^.*?\,//g;
		$l =~ s/true/1/ig;
		$l =~ s/false/0/ig;
		@xs = split(/\s*\,\s*/, $l);
		foreach $x (@xs)
		{
			push(@mask, $x);
		}
	}
	if(int($#mask+1) != 32*32 )
	{
		print "Using default mask\n";
		@mask = ();
		for($i=0; $i<32*32; $i++)
		{
			push(@mask, 1);
		}
	}
	else
	{
		print "Mask read from filter_parsivel.csv\n";
	}

	open F, 'margin_parsivel.csv';
	@ls = <F>;
	close F;
	$ls[0] =~ s/[\r\n\"]+//g;
	@xs = split(/\s*\,\s*/, $ls[0]);
	@margin = ();
	foreach $x (@xs)
	{
		push(@margin, $x);
	}
	if(int($#margin+1) != 32 )
	{
		print "Using default margin factors\n";
		@margin = ();
		for($i=0; $i<32; $i++)
		{
			push(@margin, 1.0);
		}
	}
	else
	{
		print "Margin factors read from margin_parsivel.csv\n";
	}

	# Procesado de los archivos
	$last = '';
	foreach $f (@fs){
		if($f =~ /(\d\d\d\d)(\d\d)(\d\d)(\d\d)\./){
			$dy = $1;
			$dm = $2;
			$dd = $3;
			$dh = $4;
			if( $f =~ /\.gz$/ ){
				$x = `gunzip -c $f`;
				$x =~ s/[\r\n]+/\n/g;
		                    @ls = split(/(?<=\n)/, $x);
			}
			else{
				open F, $f;
				@ls = <F>;
				close F;
			}

			if( $last ne '' )
			{
				$ls[0] = $last.$ls[0];
				$last = '';
			}
			if( ! ( $ls[$#ls] =~ /\n$/ ) )
			{
				$last = pop @ls;
			}

			@d = split(/\;/, $ls[0]);
			if($d[19] =~ /\d+\:(\d+)\:\d+/){
				$imin = $1;
			}

			$line = 0;
			foreach $l (@ls){
				$l =~ s/[\r\n]+//g;
				if( $l ne '' ){

					@d = split(/\;/, $l);

					$id = $d[22];
					$serial = $d[12];

					if($d[19] =~ /\d+\:(\d+)\:\d+/){
						$amin = $1;
						if($amin>=$imin){
							$amin=$amin-$imin;
						}
						else{
							$amin=$amin+60-$imin;
						}
						if($amin<10){$amin='0'.$amin;}
					}
					$time = "$dy/$dm/$dd $dh:$amin:00";


					# Fix time based on Artila ntp
					$shift = 0;
					if( $d[1118] =~ /^\d\d\d\d\-\d\d\-\d\d \d\d\:\d\d\:\d\d$/ )
					{
						# Round to nearest minute
						$shift = secs2($d[1118]);
						$time = date2(int($shift/60+0.5)*60);
						$shift = $shift - int($shift/60+0.5)*60;
						#print $d[1118]."\t".$time."\n";
					}
					$dt = $d[8];
					$synop = $d[3]; # Synop code 4677

		            #$r = $d[25]; # R (rain intensity, mm h-1)
		            #if($r eq "-9.999"){
		            #       $r = $d[26];
		            #}
		            $r = $d[0]; # R (rain intensity, mm h-1)

					$p = $r/60; # P (precipitation amount, mm)

					$m = $d[7]; #MOR visibility, m
					$z = $d[6]; #Z (radar reflectivity, dB mm6 m-3)
					if($z eq "-9.999"){
						$z = 0;
					}
				
					$e = $d[29]; # Kinetic energy (KJ)
					if($e>0.0 && $p>0.0){
						$e = 1000*$e/$alpha/$p/1.E7; # Kinetic energy (J m-2 mm-1)); el escalado /1.E7 es para obtener valores comparables con Thies y con la literatura, pero no debería estar
					}else{
						$e = 0;
					}

					$na = 0;

					# Strange matrix value
        	        #$k=0;
	                #for($i=0;$i<1023;$i++){
	                #	$k = $k + $d[94+$i];
	                #}
					#if($k == 0 && $d[94+1023] == 256)
					#{
					#	$d[94+1023] = 0;
					#}

					# Strange matrix value
					if($d[94+1023] == 256)
					{
						$d[94+1023] = 0;
					}

					$error = 0;
					if(int($d[24]) == 1)
					{
						$error = 21;
					}
					elsif(int($d[24]) == 2)
					{
						$error = 22;
					}
					elsif(int($d[24]) == 3)
					{
						$error = 23;
					}

					if($error != 0)
					{
						$data{$time} = "err$error,$#d";
					}
					elsif($synop !~ /^[\d\+\-\.]+$/){
						$data{$time} = "err3,$#d";
					}
					elsif($r =~ /[^\d\+\-\.eE]/){
						$data{$time} = "err4,$#d";
					}
					elsif($r =~ /9999\.999/){
						$data{$time} = "err5,$#d";
					}
					elsif($l =~ /(OK|Version)/){
						$data{$time} = "err6,$#d";
					}
					elsif($l !~ /^[\w\+\-\.\:\;\r\n\s\t ]+$/i || $l =~ /^\;/){
						$data{$time} = "err7,$#d";
					}
					elsif($serial ne "" && $serial ne $iserial)
					{
						#print "Telegram with wrong ID ($serial) found in $f\n";
						$data{$time} = "err7,$#d";
					}
                                        elsif ($#d < 1117 || $#d > 1118)
                                        {
                                                $data{$time} = "err7,$#d";
                                        }
					else{

		                for($j=0;$j<32;$j++){$dsdv[$j]=0.0;}
						for($i=0;$i<32;$i++){$dsdd[$i]=0.0;}

		                # Lectura matrix DSD (velocidades x diametros)
		                $k=0;
		                for($j=0;$j<32;$j++){
		                    for($i=0;$i<32;$i++){
								$dsd[$j][$i] = $d[94+$k] * $mask[$j+$i*32] * $margin[$i];
								$dsdv[$j] = $dsdv[$j] + $d[94+$k] * $mask[$j+$i*32] * $margin[$i];
								$dsdd[$i] = $dsdd[$i] + $d[94+$k] * $mask[$j+$i*32] * $margin[$i];

		                        $k++;
		                    }
		                }

						# Cálculo del número de partículas
						$k = 0;
		                for($i=0;$i<32;$i++){
		                    for($j=0;$j<32;$j++){
								$k = $k + $dsd[$j][$i];
							}
						}
						$np = int($k);

		                # Cálculo R (rain intensity, mm h-1)
		                $cr=0.0;
		                for($i=0;$i<32;$i++){
		                    $x=0.0;
		                    for($j=0;$j<32;$j++){
		                        $x += $dsd[$j][$i];
		                    }
		                    $cr += $x*$dmed3[$i];
		                }
		                $cr = 6E-4*$pi/$ta*$cr;

						$cp = $cr/60; # P (precipitation amount, mm)

						# Cálculo vector ND (number density, m-3 mm-1)
						@nd = ();
						@nd2 = ();
						for($i=0;$i<32;$i++){
							$nd[$i] = 0.0;
							for($j=0;$j<32;$j++){
								$nd[$i] += $dsd[$j][$i]/$vmed[$j]/$ta;
							}
							if($cp>0.0){
								$nd[$i] = $nd[$i] / $cp;
							}
							$nd2[$i] = 'nd'.($i+1);
						}

		         		# Cálculo M (water content, g m3)
		                $cm=0.0;
		                for($i=0;$i<32;$i++){
		                    $x=0.0;
		                    for($j=0;$j<32;$j++){
		                        $x += $dsd[$j][$i]/$vmed[$j];
		                    }
		                    $cm += $x*$dmed3[$i];
		                }
		                $cm = 1E-3*$pi/6.0/$ta*$cm;

		                # Cálculo Z (radar reflectivity, dB mm6 m-3)
		                $cz=0.0;
		                for($i=0;$i<32;$i++){
		                    $x=0.0;
		                    for($j=0;$j<32;$j++){
		                        $x += $dsd[$j][$i]/$vmed[$j];
		                    }
		                    $cz += $x*$dmed6[$i];
		                }
		                $cz = 10*log10(1/$ta*$cz);

						# Cálculo E (kinetic energy, J m-2 mm-1)
						$ce=0.0;
						for($i=0;$i<32;$i++){
							$x=0.0;
							for($j=0;$j<32;$j++){
								$x += $dsd[$j][$i]*$vmed2[$j];
							}
							$ce += $x*$dmed3[$i];
						}
						if($cp>0.0){
							$ce = 1.0/$alpha/$cp/12.0*$pi*1.E-6*$ce;
						}
						else{
							$ce = 0.0;
						}

						# Cálculo V (MOR visibility, m)
						$cv=0.0;
						for($i=0;$i<32;$i++){
							$x=0.0;
							for($j=0;$j<32;$j++){
								$x += $dsd[$j][$i]/$vmed[$j];
							}
							$cv += $x*($dmed[$i]/1000)*($dmed[$i]/1000);
						}
						if($cv>0.0){
							$cv = 3/(($pi/2)*(1/$ta)*$cv);
						}
						else{
							$cv = 0.0;
						}

						$d1 = ''; # D at 10%
						$d2 = ''; # D at 25%
						$d3 = ''; # D at 50%
						$d4 = ''; # D at 75%
						$d5 = ''; # D at 90%
						$d6 = 0.0; # Avergage D

						@list = ();

						if($prob eq 'middle')
						{
							$x = 0.0;
							for($i=0;$i<32;$i++){
								for($j=0;$j<32;$j++){
									$d6 += $dsd[$j][$i]*$dmed[$i];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $dmed[$i]);
									}
								}
							}
						}
						elsif($prob eq 'uniform')
						{
							$x = 0.0;
							for($i=0;$i<32;$i++){
								for($j=0;$j<32;$j++){
									$d6 += $dsd[$j][$i]*$dmed[$i];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $d1[$i] + rand() * ($d2[$i] - $d1[$i]));
									}
								}
							}
						}
						elsif($prob eq 'linear')
						{
							$x = 0.0;
							for($i=0;$i<32;$i++){

								if($i==0)
								{
									$m = $dsdd[$i+1] - $dsdd[$i];
								}
								elsif($i==31)
								{
									$m = $dsdd[$i] - $dsdd[$i-1];
								}
								else
								{
									$m = ($dsdd[$i+1] - $dsdd[$i-1]) * 0.5;
								}

								for($j=0;$j<32;$j++){
									$d6 += $dsd[$j][$i]*$dmed[$i];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										$done = 0;
										while($done == 0)
										{
											$x1 = rand();
											$x2 = rand();
											$p1 = $dsdd[$i] - $m * 0.5 + $x1 * $m;
											$p2 = $x2 * ($dsdd[$i] + $m * 0.5);
											if($p2 < $p1)
											{
												push(@list, $d1[$i] + $x1 * ($d2[$i] - $d1[$i]));
												$done = 1;
											}
										}
									}
								}
							}
						}

						@list = sort { $a <=> $b } @list;

						if($x>0.0){
							$d6 = $d6 / $x;

							$d1 = $list[int(0.10*$x)];
							$d2 = $list[int(0.25*$x)];
							$d3 = $list[int(0.50*$x)];
							$d4 = $list[int(0.75*$x)];
							$d5 = $list[int(0.90*$x)];
						}
						else
						{
							$d1 = 0.0;
							$d2 = 0.0;
							$d3 = 0.0;
							$d4 = 0.0;
							$d5 = 0.0;
							$d6 = 0.0;
						}

						$v1 = ''; # V at 10%
						$v2 = ''; # V at 25%
						$v3 = ''; # V at 50%
						$v4 = ''; # V at 75%
						$v5 = ''; # V at 90%
						$v6 = 0.0; # Average V

						@list = ();

						if($prob eq 'middle')
						{
							$x = 0.0;
							for($j=0;$j<32;$j++){
								for($i=0;$i<32;$i++){
									$v6 += $dsd[$j][$i]*$vmed[$j];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $vmed[$j]);
									}
								}
							}
						}
						elsif($prob eq 'uniform')
						{
							$x = 0.0;
							for($j=0;$j<32;$j++){
								for($i=0;$i<32;$i++){
									$v6 += $dsd[$j][$i]*$vmed[$j];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										push(@list, $v1[$j] + rand() * ($v2[$j] - $v1[$j]));
									}
								}
							}
						}
						elsif($prob eq 'linear')
						{
							$x = 0.0;
							for($j=0;$j<32;$j++){

								if($j==0)
								{
									$m = $dsdv[$j+1] - $dsdv[$j];
								}
								elsif($j==31)
								{
									$m = $dsdv[$j] - $dsdv[$j-1];
								}
								else
								{
									$m = ($dsdv[$j+1] - $dsdv[$j-1]) * 0.5;
								}

								for($i=0;$i<32;$i++){
									$v6 += $dsd[$j][$i]*$vmed[$j];
									$x += $dsd[$j][$i];

									for($k=0; $k<$dsd[$j][$i]; $k++)
									{
										$done = 0;
										while($done == 0)
										{
											$x1 = rand();
											$x2 = rand();
											$p1 = $dsdv[$j] - $m * 0.5 + $x1 * $m;
											$p2 = $x2 * ($dsdv[$j] + $m * 0.5);
											if($p2 < $p1)
											{
												push(@list, $v1[$j] + $x1 * ($v2[$j] - $v1[$j]));
												$done = 1;
											}
										}
									}
								}
							}
						}

						@list = sort { $a <=> $b } @list;

						if($x>0.0){
							$v6 = $v6 / $x;

							$v1 = $list[int(0.10*$x)];
							$v2 = $list[int(0.25*$x)];
							$v3 = $list[int(0.50*$x)];
							$v4 = $list[int(0.75*$x)];
							$v5 = $list[int(0.90*$x)];
						}
						else
						{
							$v1 = 0.0;
							$v2 = 0.0;
							$v3 = 0.0;
							$v4 = 0.0;
							$v5 = 0.0;
							$v6 = 0.0;
						}

						$r = sprintf("%.3e", $r);
						$p = sprintf("%.3e", $p);
						$e = sprintf("%.3e", $e);

						$cr = sprintf("%.3e", $cr);
						$cm = sprintf("%.3e", $cm);
						$cz = sprintf("%.3e", $cz);
						$cv = sprintf("%.3e", $cv);
						$ce = sprintf("%.3e", $ce);

						$d1 = sprintf("%.3f", $d1);
						$d2 = sprintf("%.3f", $d2);
						$d3 = sprintf("%.3f", $d3);
						$d4 = sprintf("%.3f", $d4);
						$d5 = sprintf("%.3f", $d5);
						$d6 = sprintf("%.3f", $d6);

						$v1 = sprintf("%.3f", $v1);
						$v2 = sprintf("%.3f", $v2);
						$v3 = sprintf("%.3f", $v3);
						$v4 = sprintf("%.3f", $v4);
						$v5 = sprintf("%.3f", $v5);
						$v6 = sprintf("%.3f", $v6);

						$data{$time} = "$synop,$cr,$cp,$cm,$cz,$ce,$cv,$r,$z,$e,$m,$na,$na,$na,$na,$na,$d[10],$np,$na,$d[9],$d[16],$d[11],$d1,$d2,$d3,$d4,$d5,$d6,$v1,$v2,$v3,$v4,$v5,$v6,".join(',', @nd).",$shift,$line,0,".int($#d+1);
					}

				}else{
					$data{$time}='';
				}

				$line++;
			}
		}
	}

	open F, '>'.$outfile;

	print F "'type','serial','time','seconds','synop','r','p','m','z','e','mor','r_meas',";
	print F "'z_meas','e_meas','mor_meas','qual','tmp','rh','w','wd','np_meas','np','lcurrent','ocontrol','power',";
	print F "'tmp_int','d10','d25','d50','d75','d90','dmean','v10','v25','v50','v75','v90','vmean','".join('\',\'', @nd2)."','t_shift','now','err','ncol'\n";

	# Recorremos fechas
	for($i=$secs_i;$i<=$secs_f;$i+=60){
			$x = get_date($i);
			# Comprobamos si hay dato para esta fecha
			if($data{$x} eq ''){
				#print F "Par,$iserial,$x,$i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0\n";
				print F "Par,$iserial,$x,$i,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,";
				print F "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,0\n";
			}
			elsif($data{$x} =~ /^err(\d+)\,([^\,]+)/)
			{
				$err = $1;
				$num = $2+1;
				#print F "Par,$iserial,$x,$i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,$err,$num\n";
				print F "Par,$iserial,$x,$i,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,";
				print F "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,$err,$num\n";
			}
			else{
				# Si todo ha ido bien ponemos un valor de 0 en el error status
				$data{$x} =~ s/\,/></g;
				$data{$x} =~ s/<[0\.e\+\-]+>/<0>/g;
				$data{$x} =~ s/^[0\.e\+\-]+>/0>/g;
				$data{$x} =~ s/<[0\.e\+\-]+$/<0/g;
				$data{$x} =~ s/></\,/g;

				print F "Par,$iserial,$x,$i,".$data{$x}."\n";
			}
	}
	close F;
}

exit;

# Segundos desde 1970
sub secs{
	my $x=$_[0];
	if($x =~ /^(\d+)\/(\d+)\/(\d+) (\d+)\:(\d+)/){
		my $y=$1;
		my $m=$2;
		my $d=$3;
		my $h=$4;
		my $s=$5;
		return timegm(0,$s,$h,$d,$m-1,$y-1900);
	}
}

sub secs2{
	my $x=$_[0];
	if($x =~ /^(\d+)[^\d](\d+)[^\d](\d+) (\d+)[^\d](\d+)[^\d](\d+)/){
		my $y=$1;
		my $m=$2;
		my $d=$3;
		my $h=$4;
		my $mm=$5;
		my $s=$6;
		return timegm($s,$mm,$h,$d,$m-1,$y-1900);
	}
}

sub date2{
	strftime '%Y/%m/%d %H:%M:%S', gmtime $_[0];
}

# Fecha a partir de segundos desde 1970
sub get_date{

        my ($sec,$min,$hour,$mday,$mon,
        $year,$wday,$yday,$isdst) = gmtime($_[0]);

		$year=$year+1900;
		$mon++;
		if($mon<10){ $mon='0'.$mon; }
		$mday=int($mday);
		if($mday<10){ $mday='0'.$mday; }
		$min=int($min);
		if($min<10){ $min='0'.$min; }
		$hour=int($hour);
		if($hour<10){ $hour='0'.$hour; }
		$sec=int($sec);
		if($sec<10){ $sec='0'.$sec; }
        return "$year/$mon/$mday $hour:$min:$sec";

}

sub log10 {
	my $n = shift;
	if($n<=0.0){ return 0.0; }
	return log($n)/log(10);
}
