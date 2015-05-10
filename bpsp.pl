#!/usr/bin/perl
#
#
#
#
#
#
use strict;
use warnings;
use Switch;
use List::Util qw/sum/;
use Getopt::Long;


sub ppt_max;
sub find_ppt;
sub score_ppt;
sub if_ppt;
######-------------get all the parameters from files-------- 

sub print_usage(){
        print "Progarm: bpsp\n",
                "Version: 0.0.1 (05/08/2015)\n",
                "\n",
                "Usage: bpsp.pl \n",
                "\n",
                "Required options:\n",
                "  --onlyM   INT   1 : only use motif; 0 : use both motif and PPT. [1]\n",
                "  --nBP     INT   Reported BPs. [3]\n",
                "  --motif   FILE  motif file\n",
                "  --intron  FILE  intron file\n",
                "  --PPT     FILE  PPT score\n",
                "  --out     FILE  output file\n",
}
if ($#ARGV==-1){
        print_usage();
        exit;
}

my ($fpwm,$fseq,$fppt,$fout);
my $onlyBP = 0;
my $nbp = 3;

GetOptions(
        'onlyM=i' => \$onlyBP,
        'nBP=i' => \$nbp,
        'motif=s' => \$fpwm,
        'intron=s' => \$fseq,
        'PPT=s' => \$fppt,
        'out=s' => \$fout,
);

#my $fpwm = $ARGV[0];
#my $onlyBP = $ARGV[1];
#my $onlyTNA = $ARGV[2];
#my $nbp = $ARGV[3];

my @NN= qw(A C G T);
#my $S_BACK = -300;
#my $E_BACK = -287;
#my $S_BP = -34;
#my $E_BP = -21;
#my $S_PPT = -16;
#my $E_PPT = -3;
my $L_BP_MOTIF = 7;
my $L_PPT_MOTIF = 8;
my $L_PPT = 25;
my $L_WINDOW = 40;
#my $L_WINDOW = 25;
my $OFFSET = 10;
my $BP_BASE="TACTAAC";
my $PPT_BASE="TTTTTTTT";
my $onlyTNA = 0;


my @NSS = gene_subseq($L_BP_MOTIF,"max");
my $pwmBP = get_pwmBP($fpwm);
my $bpS = pwmBP2Score(\@NSS,$pwmBP,$BP_BASE);
#my $pptS = get_pptS("out/human/ppt.score.8");
my $pptS = get_pptS($fppt);

#my $fseq = "src/human/BPS_exp1.txt";
open OUT, ">$fout" || die "Can't open the file!";
open IN, $fseq || die "Can't open the file!";
while (<IN>){
	next if (/^#/);
	chomp;
	my @line = split /\t/;
#my $seq = $line[5];
	my $seq = $line[1];
	my $ge = $line[0];
	my ($p_bp,$BP,$BP_sc,$PPT_sc,$SC) = get_BPPT_score($seq,$ge);
#print OUT join("\t",@line[0..4]),"\t",join(",",@$p_bp),"\t",join(",",@$BP),"\t",join(",",@$BP_sc),"\t",join(",",@$PPT_sc),"\t",join(",",@$SC),"\n";
	print OUT $ge,"\t",join(",",@$p_bp),"\t",join(",",@$BP),"\t",join(",",@$BP_sc),"\t",join(",",@$PPT_sc),"\t",join(",",@$SC),"\n";
}
close IN;
close OUT;

sub get_pwmBP{
	my ($f_pwmBP) = @_;
	my %pwmBP = ();
	open IN, $f_pwmBP || die "Can't open the file!";
        my $k = 0;
	my @n = ();
        while(<IN>){
                chomp;
		if (/^\t/){
			@n = split /\t/;
		}else{
			my @line = split /\t/;
			for (my $i=1;$i<=$#n;$i++){
				$pwmBP{$line[0]}{$n[$i]} = $line[$i];
			}
		}
        }
	close IN;
#	foreach my $p (sort keys %pwmBP){
#		print $p;
#		foreach my $n (sort keys %{$pwmBP{$p}}){
#			print "\t",$n,"\t",$pwmBP{$p}{$n};
#		}
#		print "\n";
#	}
	return (\%pwmBP);
}

sub get_pptS{
	my ($f_pptS) = @_;
	my %pptS = ();
	open IN, $f_pptS || die "Can't open the file!";
	while(<IN>){
		next if (/^#/);
		chomp;
		my @line = split /\t/;
		$pptS{$line[0]} = $line[4];
	}
	close IN;
	return(\%pptS);
}

sub pwmBP2Score{
	my ($NSS,$pwmBP,$bp_base) = @_;

#my $bp_base_score=1;
	my $bp_base_score=0;
	my @BP=split(//,$bp_base);
	for(my $i=0;$i<=$#BP;$i++){
		my $n=$BP[$i];
		switch($n){
#case 'A'{$bp_base_score*=$$pwmBP{$i+1}{"A"};}
#			case 'C'{$bp_base_score*=$$pwmBP{$i+1}{"C"};}
#			case 'G'{$bp_base_score*=$$pwmBP{$i+1}{"G"};}
#			case 'T'{$bp_base_score*=$$pwmBP{$i+1}{"T"};}
			case 'A'{$bp_base_score+=log($$pwmBP{$i+1}{"A"});}
			case 'C'{$bp_base_score+=log($$pwmBP{$i+1}{"C"});}
			case 'G'{$bp_base_score+=log($$pwmBP{$i+1}{"G"});}
			case 'T'{$bp_base_score+=log($$pwmBP{$i+1}{"T"});}
		}
	}

	my %bpS =();
	foreach my $m (sort @$NSS){
#my $score=1;
		my $score=0;
		my @BP=split(//,$m);
		for(my $i=0;$i<=$#BP;$i++){
			my $n=$BP[$i];
			switch($n){
#case 'A'{$score*=$$pwmBP{$i+1}{"A"};}
#				case 'C'{$score*=$$pwmBP{$i+1}{"C"};}
#				case 'G'{$score*=$$pwmBP{$i+1}{"G"};}
#				case 'T'{$score*=$$pwmBP{$i+1}{"T"};}
				case 'A'{$score+=log($$pwmBP{$i+1}{"A"});}
				case 'C'{$score+=log($$pwmBP{$i+1}{"C"});}
				case 'G'{$score+=log($$pwmBP{$i+1}{"G"});}
				case 'T'{$score+=log($$pwmBP{$i+1}{"T"});}
			}
		}
#$bpS{$m} = $score/$bp_base_score;
		$bpS{$m} = $score-$bp_base_score;
	}
	return (\%bpS);
}
sub gene_subseq{
        my ($l,$flag) = @_;
        my @NN= qw(A C G T);
        my @NSS=@NN;
        my @NS=@NN;
        my @NNS=();

        if ($l eq 1){
                @NSS=@NN;
        }else{
                for (my $i=2;$i<=$l;$i++){
                        foreach my $N (@NS){
                                foreach my $n (@NN){
                                        push(@NNS,$N.$n);       
                                }       
                        }
                        @NS=@NNS;
                        @NSS=(@NSS,@NNS);
                        @NNS=();
                }       
        }
        if ($flag eq "full"){
                return(@NSS);
        }elsif($flag eq "max"){
                return(@NS);
        }else{
                return(0);
        }
}


sub get_AGEZ{ #find the AGEZ
	my ($_) = @_;
	my $seq = $_;
	my $sL = length($_);
	my @posAG=();
	my $pAG = -1;
	push (@posAG, pos()-length($&)), "\n" , while(/AG/g);

	for (my $i=$#posAG-1;$i>=1;$i--){
		if ($posAG[$#posAG]-$posAG[$i] > 12) {$pAG=$posAG[$i]-14; last;}
	}
	##if ($pAG <0 ) {$pAG=$sL+$pstart};
	#$pAG=$sL if ($pAG < 0);
	$pAG=0 if ($pAG < 0);
	return ($pAG,substr($seq,$pAG,$sL-$pAG));
}


############################################################################
# calculate the scores of candidate BP and PPT, and the distance between them
# in consTNA sequences
sub get_BPPT_score{
	my ($seq,$ge) = @_;
	our $maxppt=();
	our $pptstart=();
	our $pptend=();

	#$seq="AGGCGTACCTTAACTACTAACAAGGGGTTTTTTTTTTTTTTTTTTTTGGGGGGAAAAAACAAAGGAC";
	my $sL=length($seq);
	#find the AGEZ
	my ($pAG,$s_agez) = get_AGEZ($seq);
#	print $pAG,"\t",$s_agez,"\n";
	my $l_s_agez = length($s_agez);

	my @BP = ("NA")x$nbp;
	my @BP_sc = ("NA")x$nbp;
	my @PPT_sc = ("NA")x$nbp;
	my @P_bp = (0)x$nbp;
	my $PP = 0;
	#my $SC = -1*$L_WINDOW;
	my @SC = (-1*$L_WINDOW)x$nbp;
#print $l_s_agez,"\t",$l_s_agez-$L_BP_MOTIF,"\t-1\n";
	for(my $p=0;$p<=$l_s_agez-2*$L_BP_MOTIF-$L_PPT_MOTIF;$p++){
		my $subseq = ();
		if ($l_s_agez-$p < $L_WINDOW){
#$subseq = $s_agez;
			$subseq = substr($s_agez,$p,$l_s_agez-$p);
		}else{
			$subseq = substr($s_agez,$p,$L_WINDOW);
		}
		my ($bp_can,$bp_sc,$ppt_sc,$sc) = find_bppt($subseq);
		if ($onlyTNA == 1){
			next if (!($bp_can =~ /\w{3}T\wA\w/));
		}

		for(my $ii=0;$ii<$nbp;$ii++){
			my $ssc = $SC[$ii];
			if ($sc>$ssc){
				splice(@SC,$ii,0,$sc);
#print $sc,"\t",$ssc,"\t",$ii,"\t",join(",",@SC),"---------\n";
				splice(@SC,$nbp,1);
#print $sc,"\t",$ssc,"\t",$ii,"\t",join(",",@SC),"++++++++++\n";
				splice(@BP,$ii,0,$bp_can);
				splice(@BP,$nbp,1);
				splice(@BP_sc,$ii,0,$bp_sc);
				splice(@BP_sc,$nbp,1);
				splice(@PPT_sc,$ii,0,$ppt_sc);
				splice(@PPT_sc,$nbp,1);
				my $p_bp = $sL-($p+$L_BP_MOTIF-1+$pAG)+1;
				splice(@P_bp,$ii,0,$p_bp);
				splice(@P_bp,$nbp,1);
				last;
			}
		}
	}
#	my $p_bp = $sL-($PP+$L_BP_MOTIF-1+$pAG)+1;
#	return($p_bp,$BP,$BP_sc,$PPT_sc,$SC);
	return(\@P_bp,\@BP,\@BP_sc,\@PPT_sc,\@SC);
}

sub find_bppt{
	my ($seq) = @_;
	my $l = length($seq);
	my $bp = substr($seq,0,$L_BP_MOTIF);
	my $scppt = 0;
	my $bps = 0;
	$bps = $$bpS{$bp} if (exists($$bpS{$bp}));
	if ($onlyBP){
		return($bp,$bps,$scppt,$bps);
	}
	my $n=0;
	for (my $p=$L_BP_MOTIF;$p<=$l-$L_PPT_MOTIF;$p++){
		++$n;
		my $s = substr($seq,$p,$L_PPT_MOTIF);
		if (exists($$pptS{$s})){
			$scppt += $$pptS{$s};
		}
	}
	$scppt = $scppt/$n;
#return($bp,$bps,$scppt,$bps*$scppt);
	return($bp,$bps,$scppt,$bps+$scppt);
}


