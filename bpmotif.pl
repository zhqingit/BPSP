#!/usr/bin/perl
############################################################################
#Find the motif by using mixture model

use strict;
use warnings;
use Switch;
use List::Util qw/sum/;
use Getopt::Long;

sub get_PWMbp_file;
sub scores_BP;
sub initial_paras;
sub upZ;
sub upPwm;
sub upLamda;
sub EM;
sub E;
sub M;
sub dis_paras;
sub in_pwm;
sub loglik;
sub logg;
sub ifcon;

sub print_usage(){
	print "Progarm: bpmotif\n",
		"Version: 0.0.1 (05/02/2015)\n",
		"\n",
		"Usage: bpmotif.pl \n",
		"\n",
		"Required options:\n",
		"  --pmotif  FILE  prior motif file\n",
		"  --motif   FILE  inferred motif saved\n",
		"  --bpolyn  FILE  frequencies of heptanucleotides in branch point and background regions\n",
}
if ($#ARGV==-1){
	print_usage();
	exit;
}

my ($fbpwm,$fout,$fbpolyn);

GetOptions(
	'pmotif=s' => \$fbpwm,
	'motif=s' => \$fout,
	'bpolyn=s' => \$fbpolyn,
);
#my $fbpwm = $ARGV[0];
#my $fout = $ARGV[1];
#my $fbpolyn = $ARGV[2];

my %c_polyn = ();
#my $fbpolyn = "out/human/bp.polyn";
open IN, $fbpolyn || die "Can't open the file!";
while(<IN>){
	next if (/^back/);
	next if (/N/);
	chomp;
	my @line = split /\t/;
	$c_polyn{$line[0]}{"bp"} = $line[2];
	$c_polyn{$line[0]}{"back"} = $line[1];
	$c_polyn{$line[0]}{"ppt"} = $line[3];
}
close IN;



my %pwmB = ();
my %pwmBP = ();
my %pwmPPT = ();
my %beta = ();
my %Z = ();
my $lamdaBP = 0.5;
my $lamdaB = 0.5;
my $base = "TACTAAC";
my $ie = 10e-6;

print("Initing......\n");
init($fbpwm);
print_pwm();

print("EM......\n");
EM(100,\%c_polyn);
print_pwm();

open OUT, ">$fout" || die "Can't open the file!";
print OUT "\t",join("\t",qw(A C G T)),"\n";
foreach my $p (sort keys %pwmBP){
	print OUT $p;
	foreach my $n (sort keys %{$pwmBP{$p}}){
#		if ($pwmBP{$p}{$n}<0.001){
#			print OUT "\t",0.001;
#		}else{
			print OUT "\t",$pwmBP{$p}{$n};
#		}
	}
	print OUT "\n";
}

print "######################################\n";
foreach my $m (sort keys %c_polyn){
#	print $m,"\t",$Z{"BP"}{$m},"\t",$Z{"BACK"}{$m},"\n";
}


sub EM{
	my ($times,$c_seq) = @_;
#	for (my $i=1;$i<=$times;$i++){
		#print $i,"\ttimes\n";
	my $d = 100;
	my $i =0 ;
	while($d>$ie){
		print ++$i,"\n";
		my %pwmBold = %pwmB;
		E($c_seq);
		M($c_seq);
		$d = dis_e(\%pwmBold,\%pwmB);
#		print_pwm();
	}
}

sub M{
	my ($c_seq) = @_;
	$lamdaBP=0;
	$lamdaB=0;
	my $sumotif = 0;

	%pwmBP = ();
	%pwmB = ();

	##foreach my $m (@NSS){
	foreach my $m (sort keys %{$c_seq}){
		$lamdaBP += $Z{"BP"}{$m}*$$c_seq{$m}{"bp"};
		$lamdaB += $Z{"BACK"}{$m}*$$c_seq{$m}{"bp"};
		$sumotif += $$c_seq{$m}{"bp"};


		#my $sumppt = 0; #newppt	
		my $cm = $$c_seq{$m}{"bp"};
		my $zbp = $Z{"BP"}{$m};
		my $zb = $Z{"BACK"}{$m};
		my @n = split(//,$m);
		for(my $i=0;$i<=$#n;$i++){
			my $nn = $n[$i];
			$pwmB{$nn} += $zb*$cm;
			$pwmBP{$i+1}{$nn} += $zbp*$cm;
		}
		#$pwmPPT{$m} += $zppt*$cm; #newppt
		#$sumppt += $zppt*$cm; #newppt
        }

        $lamdaBP=$lamdaBP/$sumotif;
        $lamdaB=$lamdaB/$sumotif;


	foreach my $p(sort keys %pwmBP){
                my $sumBP=0;
                foreach my $nn(sort keys %{$pwmBP{$p}}){
                        $sumBP += $pwmBP{$p}{$nn};
                }
                foreach my $nn(sort keys %{$pwmBP{$p}}){
#$pwmBP{$p}{$nn} = $pwmBP{$p}{$nn}/$sumBP;
                        $pwmBP{$p}{$nn} = ($pwmBP{$p}{$nn}+$beta{$nn})/($sumBP+$c_polyn{$base}{"bp"});
                }
        }

	my $sumB =0;
        foreach my $nn (sort keys %pwmB){
		$sumB += $pwmB{$nn};
        }
        foreach my $nn (sort keys %pwmB){
                $pwmB{$nn} /= $sumB;
        }
}

sub E{
	my ($c_seq) = @_;
	#foreach my $m (@NSS){
	foreach my $m (sort keys %{$c_seq}){

		my @n = split(//,$m);
		my $scoreBP=1;
                my $scoreB=1;

		for(my $i=0;$i<=$#n;$i++){
			my $nn = $n[$i];
			$scoreBP *= $pwmBP{$i+1}{$nn};
			$scoreB *= $pwmB{$nn};
		}
		
		my $zbp=$lamdaBP*$scoreBP;
                my $zb=$lamdaB*$scoreB;
		my $sumZ=$zbp+$zb;
		$Z{"BP"}{$m} = $zbp/$sumZ;
		$Z{"BACK"}{$m} = $zb/$sumZ;
	}
	return(1);
}

sub init{
	my ($f_bp) = @_;
	my %cpwmB = ();
	my %cpwmPPT = ();

	my $betatot = 0;
	my $sumB = 0;
	#my $sumppt = 0; #forppt
	foreach my $m (sort keys %c_polyn){
		my @n = split(//,$m);
		for(my $i=0;$i<=$#n;$i++){
			my $nn = $n[$i];
			$cpwmB{$nn} += $c_polyn{$m}{"back"};
			$sumB += $c_polyn{$m}{"back"};
			#$cpwmPPT{$i}{$nn} += $$c_seq{"ppt"}{$m};
			$beta{$nn} += $c_polyn{$m}{"bp"};
			$betatot += $c_polyn{$m}{"bp"};
		}
	}
	foreach my $n (sort keys %beta){
		$beta{$n} = int $beta{$n}/$betatot*$c_polyn{$base}{"bp"}+0.5;
	}

	foreach my $n (sort keys %cpwmB){
		$pwmB{$n} = $cpwmB{$n}/$sumB;
	}
=forppt
	foreach my $p (sort keys %cpwmPPT){
		my $sumPPT = 0;
		foreach my $n (sort keys %{$cpwmPPT{$p}}){
			$sumPPT += $cpwmPPT{$p}{$n};
		}
		foreach my $n (sort keys %{$cpwmPPT{$p}}){
			$pwmPPT{$p}{$n} = $cpwmPPT{$p}{$n}/$sumPPT;
		}
	}
=cut
	open IN, $f_bp || die "Can't open the file!";
	my @n = ();
	while(<IN>){
		chomp;
		if (/^\t|#/){
			@n = split /\t/;
			next;
		}
		my @line = split /\t/;
		for (my $i=1;$i<=$#n;$i++){
			$pwmBP{$line[0]}{$n[$i]} = $line[$i];
		}
	}
}


sub print_pwm{
	print "LAMDA:\n";
	print "back:\t$lamdaB\tbranch point:\t$lamdaBP\n";
	print "BP:\n";
	foreach my $p (sort keys %pwmBP){
		print $p;
		foreach my $n (sort keys %{$pwmBP{$p}}){
#		if ($pwmBP{$p}{$n}<=0.001){
#				print "\t",0.001;
#			}else{
				print "\t",$pwmBP{$p}{$n};
#			}
		}
		print "\n";
	}

	print "BACK:\n";
	foreach my $n (sort keys %pwmB){
		print $n,"\t",$pwmB{$n},"\n";
	}
}

sub dis_e{
	my ($v1,$v2) = @_;
	my $d = 0;
	foreach my $n (sort keys %$v1){
		$d += ($$v1{$n}-$$v2{$n})*($$v1{$n}-$$v2{$n});
	}
	return(sqrt($d));
}
