#!/usr/bin/perl;
use warnings;
use strict;
#use Math::GSL::CDF qw/:binomial/;
#use Math::GSL::Randist qw/:binomial/;
use binomtest;
use Getopt::Std;

## infor for revel
our $DEBUG=0;
our $Permutation_Times=1000;
#our $D_col;
our $D_0;
our $D_1;
our $D_step;
our $D_bin; #=($D_1-$D_0)/$D_step+1;
## infor for frequency
#our $F_col;
our $F_0; #=0;
our $F_1; #=0.001;
our $F_step; #=0.00005;
our $F_bin;#=($F_1-$F_0)/$F_step+1;
#print "Totally bins ".$F_bin."\n";


sub Btest {
	### do binomial test, input caseN, controlN, probability; return pvalue
	my $num=$_[0];
	my $denum=$_[1];
	my $mp=$_[2];
	my $r=-1;

	if(@_ <3){print "$num:$denum:$mp\tERROR\n";return ($r);}
	$r=binomtest::binomtest($num,$num+$denum,$mp);
	if($DEBUG==1){print "Btest: $num:$denum:$mp:$r\n";}
	return ($r);
}

sub shuffle {
  my $cnt=$_[0];
  my $mp=$_[1];
  my @arr=();
  for(my $i=0;$i<$cnt;$i++){
     my $r=rand(1);
	 if($r>$mp){push(@arr,0);}else{push(@arr,1)}
  }

  return (\@arr);
}

sub  inti_matrix {
		my ($d_bin,$f_bin)=@_;
		my @mat=();
		for(my $i=0;$i<$d_bin;$i++){
			for(my $j=0;$j<$f_bin;$j++){
					$mat[$i][$j]=0;
			}
		}
		return(\@mat);
}
sub rowsum {
	my @arr=@{$_[0]};
	my $sum=0;
	foreach my $c(@arr){
		if($c eq ""){print "error: at 131";exit;}
		$sum+=$c;
	}
	return ($sum);
}


sub permutation {
	my $phe=$_[0];
	my $matd=$_[1];
	my $adfreqs=$_[2];
	my $adrevels=$_[3];
	my $mp=$_[4];
	my $rawp=$_[5];
	my $permut=0;
	my $TIMES=$Permutation_Times;
	my $T=0;
	while($T<$TIMES){
		 my @newphe=@{shuffle(scalar(@$phe),$mp)};
		 if($DEBUG==1){
     		print "raw_phe: ".join(",",@$phe)."\n";
		 	print "new_phe: ".join(",",@newphe)."\n";
		# exit;
		}
		 if(rowsum(\@newphe)==0){$T=$T+1;next;}
		 if(rowsum(\@newphe)==@newphe){$permut+=1;$T=$T+1;next;}
		 my ($ad_cs_mat,$ad_ctr_mat)=@{summary($matd,\@newphe,$adfreqs,$adrevels)};
		 my ($pvalue,$PF,$PR)=@{choose_best_cell($ad_cs_mat,$ad_ctr_mat,$mp)};
		 if($pvalue <=$rawp){
			 $permut+=1;
			 if($DEBUG==1){
				print "$pvalue < = $rawp, $permut\n";
			}
		 }
		 $T=$T+1;
		 if($T==$TIMES-1){
			 if($permut< 10){
				 $TIMES=10*$TIMES;
				 if($DEBUG==1){
					 print "increase the permutation times to $TIMES\n";
				 }
			 }
		 }
	 }
	 my @rs=(($permut+1)/($TIMES+1),$TIMES,);
	 return \@rs;

}

sub choose_best_cell {
	## get cs_mat, ctr_mat, might compress cs_mat, ctr_mat, as many replicates cells.
	## compare per cell using bimo
	## it might be able to be compressed.
	## return a pvalue
	my ($add_cs,$add_ctr,$mp)=@_;
	my @cs_mat=@$add_cs;
	my @ctr_mat=@$add_ctr;
	my $lcs=0;
	my $lctr=0;
	my $bestp=1;
	my $bestF=0;
	my $bestR=0;
	my %tests=();
	for(my $i=0;$i<$D_bin;$i++){
		for(my $j=0;$j<$F_bin;$j++){
			my $Ncs=@{$cs_mat[$i]}[$j];
			my $Nctr=@{$ctr_mat[$i]}[$j];
			if($lcs==$Ncs && $lctr==$Nctr){next;}
			if(exists($tests{$Ncs.":".$Nctr})){next;}
			my $pvalue=Btest($Ncs,$Nctr,$mp);
			$tests{$Ncs.":".$Nctr}=1;
			if($pvalue <$bestp){
				$bestp=$pvalue;
				$bestR=$i;
				$bestF=$j;
			}
			if($DEBUG==1){
				print "$i,$j: ".$Ncs."\t".$Nctr."$pvalue\t$bestp\n";
			}
			$lcs=$Ncs;
			$lctr=$Nctr;
		}
	}
	my @rs=($bestp,$bestR,$bestF);
	if($DEBUG==1){
		print "final: "."$bestp\n";
	}
	return (\@rs);
}


sub summary {
	my ($matd,$phed,$adfreqs,$adrevels)=@_;
	my @mat=@$matd;
	my @phe=@$phed;
	my @freqs=@$adfreqs;
	my @revels=@$adrevels;
	my @cs_mat=@{inti_matrix($D_bin,$F_bin)};
	my @ctr_mat=@{inti_matrix($D_bin,$F_bin)};
	my @cs_ids=();
	my @ctr_ids=();
	for(my $id=0;$id<@phe;$id++){
		if($phe[$id]==1){
			push(@cs_ids,$id);
		}else{
			push(@ctr_ids,$id)
		}
	}
	if($DEBUG==1){
		print "newcs".join(",",@cs_ids)."\n";
		print "newctr".join(",",@ctr_ids)."\n";
	}
	## go through each row
	## check freq, check revel
	## fille cells in the matrix

	for(my $i=0;$i<@freqs;$i++){
		my $freq=$freqs[$i];
		my $revel=$revels[$i];
		my $FN=int(($freq-$F_0)/$F_step);
		my $RN=int(($revel-$D_0)/$D_step);
		my @cs_var=@{$mat[$i]}[@cs_ids];
		my @ctr_var=@{$mat[$i]}[@ctr_ids];
		if($cs_var[0] eq ""){
			print "Error: $i.\t.$freq\tfreq length:". scalar(@freqs)."\n";
			exit;
		}

	# print join(",",@cs_var)."\n";
		my $cs_sum=rowsum(\@cs_var);
		my $ctr_sum=rowsum(\@ctr_var);
		if($DEBUG==1){
		 print $i."\t".$freq."\t".$revel."\tALL:".join("\t",@{$mat[$i]})."\n";
	 	 print $i."\t".$freq."\t".$revel."\tcase\t".$cs_sum."\t:".join("\t",@cs_var)."\n";
	 	 print $i."\t".$freq."\t".$revel."\tcontrol\t".$ctr_sum."\t:".join("\t",@ctr_var)."\n";
	  }
		for(my $rid=0; $rid<$RN;$rid++){
			for(my $fid=$FN;$fid<$F_bin;$fid++){  ## fill in cells
				$cs_mat[$rid][$fid]+=$cs_sum;
				$ctr_mat[$rid][$fid]+=$ctr_sum;
			#	print "changed ".($D_0 + $rid * $D_step)."\t".($F_0 + $fid * $F_step)."\tcase:control: = ".($cs_mat[$rid][$fid]).":".($ctr_mat[$rid][$fid])."\n";
			}
		}
	}
	if($DEBUG==1){
		 print "Final summary table :\n";
		 print "Revel/Freq\t";
		  for(my $mm=0;$mm<$F_bin;$mm++){
			  print +($F_0 + $mm * $F_step )."\t";
	    }
			print "\n";
			for(my $nn=0;$nn<$D_bin;$nn++){
				print +($D_0 + $nn * $D_step).":\t";
				for(my $mm=0;$mm<$F_bin;$mm++){
					print $cs_mat[$nn][$mm].":".$ctr_mat[$nn][$mm]."\t";
			}
				print "\n";
			}
	}

	my @rt=(\@cs_mat,\@ctr_mat);
	return (\@rt);
}


sub gene_burden {
	 my ($matd,$phe,$adfreqs,$adrevels,$mp)=@_;
#	 my $T=0;
	 my @newphe=@$phe;
#	 print "NN ".@newphe."\n";
	 my $RR=0;
	 my $bestF=0;
	 my $bestR=0;

	 my $permut=0;
	 my ($ad_cs_mat,$ad_ctr_mat)=@{summary($matd,$phe,$adfreqs,$adrevels)};
	 my @cs_mat=@$ad_cs_mat;
	 my @ctr_mat=@$ad_ctr_mat;
	 my ($rawp,$PR,$PF)=@{choose_best_cell($ad_cs_mat,$ad_ctr_mat,$mp)};
	 if($DEBUG==1){
	 		print "R: ".$rawp."\t".$PR."\t".$PF."\n";
 	}
	# $rawp=$pvalue;
	 $RR=1/$mp*($cs_mat[$PR][$PF]/(0.05+$ctr_mat[$PR][$PF]));
	 $bestF=$F_0+$F_step*$PF;
	 $bestR=$D_0+$D_step*$PR;
	 my ($permut_p, $T)=@{permutation($phe,$matd,$adfreqs,$adrevels,$mp,$rawp)};
	 if($DEBUG==1){
	 	 print "permutation: ".$rawp."\t".$permut_p."\t".$T."\n";
	 }
	 my @rs=($permut_p, $T,$rawp,$cs_mat[$PR][$PF],$ctr_mat[$PR][$PF],$RR,$bestF,$bestR);  #$cutoff,$p,$permutp,$Te,$OR,$Nc,$Ns
	return \@rs;
}



sub process {
	### read variants
	## read ped
	###  load variants to variants matrix
	my $fin=$_[0]; ## variants file
	my $chr=$_[1]; ## seperate by chromosome
	my @phe=@{$_[2]}; ## phenotype for cases and controls
  my $mp=$_[3];  ## fraction of cases
  my $keycol=$_[4]; ## REVEL column
	my $fkeycol=$_[5]; ## frequency column
	my @includes=@{$_[6]}; ## samples included in this test
	my %col_filter=%{$_[7]};
#	my $fout=$_[7]; ## the output file
 # my $lastgene=""; ## gene pointer
	my @matrix=(); ## matrix used for storing revel for each person, it has to be changed to 2D, carriersXvariants, first two columns are frequency and revel, the
	my %pheno_col=();
	my @newphe=();
	my @freqs=();
	my @revels=();
	my $row=0;
  # print "score: ".$start."\t".$step."\t$NR\n";
  open IN, "tabix $fin $chr|"; ## could change to get region from the cannoical transcript. reading region from outside, then parallel each gene


  while(my $line=<IN>){
       chomp($line);
       my @sets=split(/\t+/,$line);
		 my @c=@sets[@includes]; ## genotype loading
		 my $revel=$sets[$keycol];  ## get revel
		 my $freq=$sets[$fkeycol]; ## get frequency
    #   my $gene=$sets[5];
		# if($gene =~ /;/){next}
		 foreach my $fkey (keys %col_filter){
				if(index($sets[$fkey], $col_filter{$fkey}."|\.")==-1){
					print $sets[$fkey]." not match with  ". $col_filter{$fkey}."\n";
					next;
				}
			}
			if($revel eq "."){$revel=1;}
      	if($revel< $D_0){next;}
			if($freq eq "."){$freq=0;}
		  if($freq > $F_1){next;}

	    my @arr= (0) x scalar(keys %pheno_col);
	    for(my $i=0;$i<@c;$i++){
			  if($c[$i]==0){next;}
			  my $colTo=0;
			  if(!exists($pheno_col{$i})){
			  	 	$colTo=scalar(keys %pheno_col);
					$pheno_col{$i}=$colTo;
					for(my $r=0;$r<$row;$r++){
					  $matrix[$r][$colTo]=0;
				   }
				   push(@newphe,$phe[$i]);
					push(@arr,$c[$i]);
			  }else{
			  		$colTo=$pheno_col{$i};
					$arr[$colTo]=$c[$i];
		  	  }
			  $matrix[$row]=\@arr;
			  if($DEBUG==1){
				  print $row.":\t";
				  for(my $ss=0;$ss<scalar(keys %pheno_col);$ss++){
					print $matrix[$row][$ss]."\t";
			  }
				print "\n";
			}
		  #$matrix[$row]=\@arr;
	}
   push @freqs, $freq;
   push @revels, $revel;
	#$lastgene=$gene;
	$row=$row+1;
   # last;
  }
  close IN;
	if(scalar(keys %pheno_col)>0){
		#print "process $lastgene\n";
		return(gene_burden(\@matrix,\@newphe,\@freqs,\@revels,$mp));
	#	my ($permut_p, $T,$rawp,$cspf,$ctpf,$RR,$bestF,$bestR)= @{gene_burden(\@matrix,\@newphe,\@freqs,\@revels,$mp)}; #
	#	print OUT $lastgene."\t".$permut_p."\t".$T."\t".$rawp."\t".$cspf."\t".$ctpf."\t".$RR."\t".$bestF."\t".$bestR."\n";
 	}
	#close OUT;
}


sub load_info {
	print "input format: inputFile pedfile(case+control)  region KeyColumn1 keyColum2 OutPut_prefix\n ";
	print "input format pedfile: sampleName	Affected_state(1 or 2)\n ";
	print "input format region : chromosome (1..22) X Y\n";
   print "input format KeyColumn1: the 1st column name for variable threshold\n";
	print "input format KeyColumn2: the 2nd column name for variable threshold\n";
	print "example: perl VT_multiD.pl  DNMT3A_TET2.CCHMC.vcf.gz  IPAH.APAH.ped  2 REVEL gnomAD_exome_ALL test\n";
	if(@_<1){print "check the input parameters please!\n";exit;}
	my %mapInfo=%{$_[0]};
	if(!exists($mapInfo{"variantFile"})||!exists($mapInfo{"pedFile"}) || !exists($mapInfo{"dmis_query"}) ||!exists($mapInfo{"freq_query"})){
		 print join("\t", keys %mapInfo)."\n"; 
		 print "please check the profile, it must have \"variantFile\", pedFile and 2 query column setting \n";
		 exit;
	 }
	my $fin=$mapInfo{"variantFile"};
	my $fped=$mapInfo{"pedFile"};
	my @str_q1=split(/\t+/,$mapInfo{"dmis_query"});
	my @str_q2=split(/\t+/,$mapInfo{"freq_query"});
	if(@str_q1 <4 || @str_q2 <4){print "please check the profile, it must have 3 scores in each query key \n";exit;}
	my $key1=$str_q1[0];
	$D_0=$str_q1[1];
	$D_1=$str_q1[2];
	$D_step=$str_q1[3];
	$D_bin=int(($D_1-$D_0)/$D_step)+1;
	my $key2=$str_q2[0];
	$F_0=$str_q2[1];
	$F_1=$str_q2[2];
	$F_step=$str_q2[3];
	$F_bin=int(($F_1-$F_0)/$F_step)+1;
	print $key1.":\tfrom ".$D_0."\t to ".$D_1."\t by $D_step\t #bin".$D_bin."\n";
	print $key2.":\tfrom ".$F_0."\t to ".$F_1."\t by $F_step\t #bin".$F_bin."\n";
	my %filter=();
	foreach my $key(("ExonicFunCol","FunCol")){
		if(exists($mapInfo{$key})){
			my @fun_info=split(/\t+/,$mapInfo{$key});
			$filter{$fun_info[0]}=$fun_info[1];
		}
	}
	#$chr=
	#my ($fin,$fped,$key1,$key2)=@_;

	my %ped=();
	my @phe=();
	open PED, "$fped";
	while(my $line=<PED>){
	   chomp($line);
	   my @sets=split(/\t+/,$line);
	   $ped{$sets[0]}=$sets[1]-1;
	 }
	close PED;
	my %col_filter=();
	open IN,"tabix $fin  -H|tail -n 1 |" or die "$fin cannot find!\n";
#	print "tabix $fin $chr -H\n";#open OUT,">$fout";
	my $line=<IN>;
	close IN;
	chomp($line);
	my @sets=split(/\t+/,$line);
	#my $outkey=$sets[$D_key];
	my $mp=0;
	my @includes=();
	my $D_col=-1;
	my $F_col=-1;

	for(my $i=0;$i<@sets;$i++){
		if($sets[$i] eq $key1){
			$D_col=$i;
			print $key1.": at column $i\n";
		}
		if($sets[$i] eq $key2){
			$F_col=$i;
			print $key2.": at column $i\n";
		}
		if(exists($filter{$sets[$i]})){
				$col_filter{$i}=$filter{$sets[$i]};
		}

	  if(exists($ped{$sets[$i]})){
	     push(@phe,$ped{$sets[$i]});
		   push(@includes,$i);
	    }else{
	      next;
	    }
	   #print $sets[$i]."\t".$ped{$sets[$i]}."\n";
	    $mp+=$ped{$sets[$i]};
	  }
 if($D_col==-1||$D_col==-1){
	 print $key1." or $key2 "." were not found in the header\n";
 }
	print "Total samples:".@phe.", cases: $mp\n";
#	print "samples:". join(",",@includes).", cases: $mp\n";
	$mp=$mp/(@phe);
   my @rs=(\@phe,$mp,$D_col,$F_col,\@includes,\%col_filter);
	return (\@rs);
}



sub transcript_test {

	my ($fin,$ftrans,$phe,$mp,$D_col,$F_col,$includes,$col_filter,$fout)=@_;

	open OUT, " > $fout";
	print "Result outputs to $fout\n";
  	print OUT  "GeneName"."\tpermut.p\tTimes\traw.P\tNCase\tNControl\tRR\tBestFreq\tBestDmis\n" ;
	open RIN, "less $ftrans|";
	while (my $line=<RIN>){
		chomp($line);
		my @sets=split(/\t+/,$line); ## gene	chr:start-end chr_start-end
		my $region=$sets[1]; #.":".$sets[1]."-".$sets[2];
		my $gene=$sets[0];
		my ($permut_p, $T,$rawp,$cspf,$ctpf,$RR,$bestF,$bestR)=@{process($fin,$region,$phe,$mp,$D_col,$F_col,$includes,$col_filter)};
   	print OUT $gene."\t".$permut_p."\t".$T."\t".$rawp."\t".$cspf."\t".$ctpf."\t".$RR."\t".$bestF."\t".$bestR."\n";
	}
	close OUT;
}


if(@ARGV<1){
	print "test profile required!\n";
	exit;
}

## could change to get region from the cannoical transcript. reading region from outside, then parallel each gene
## change to read a setting profile
## change to case and control from several files. readling from differnt file and store to the same genotype matrix
##

print "Input: ".join(" ",@ARGV)."\n";
#my %options=();
#getopts("iprdfo",\%options);


my $profile=$ARGV[0];
open PIN,"$profile";
my %mapInfo=();
while(my $line=<PIN>){
	chomp($line);
	if($line eq ""){next}
	my @info=split(/\#/,$line);
	my @sets=split(/\t+/,$info[0]);
		print @sets."\t".$sets[0]."\n";
	$mapInfo{$sets[0]}=join("\t",@sets[1..(@sets-1)]);

}
close PIN;


#our $R = Statistics::R->new();
#my $fout="$prefix.$chr.$key1.$key2.txt";

my ($phe,$mp,$D_col,$F_col,$includes,$col_filter)=@{load_info(\%mapInfo)};
my $fout=$mapInfo{"OutName"};
transcript_test($mapInfo{"variantFile"},$mapInfo{"regionFile"},$phe,$mp,$D_col,$F_col,$includes,$col_filter,$fout);
