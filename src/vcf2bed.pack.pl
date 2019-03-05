
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use INFO2map;


print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: vcf_file\n";
my $file=$ARGV[0];
my @outf=split(/\//,$file);
my %samples=();
open IN,"zcat $file|" or die "Error: $file cannot find!\n";
open OUT, ">$outf[@outf-1].2.bed";
print "Please find out in $file.2.bed\n";
#my @header=();
my @csq=();
my $wi=0;
my $lc=0;
my @header=("FILTER","AC","AF","AN","BaseQRankSum","ClippingRankSum","DB","DP","FS","InbreedingCoeff","MQ","MQRankSum","QD","ReadPosRankSum","SOR","VQSLOD","VQSR_culprit","VQSR_NEGATIVE_TRAIN_SITE","VQSR_POSITIVE_TRAIN_SITE","GQ_HIST_ALT","DP_HIST_ALT","AB_HIST_ALT","GQ_HIST_ALL","DP_HIST_ALL","AB_HIST_ALL","AC_Male","AC_Female","AN_Male","AN_Female","AF_Male","AF_Female","GC_Male","GC_Female","GC_raw","AC_raw","AN_raw","GC","AF_raw","Hom_AFR","Hom_AMR","Hom_ASJ","Hom_EAS","Hom_FIN","Hom_NFE","Hom_OTH","Hom","Hom_raw","AC_AFR","AC_AMR","AC_ASJ","AC_EAS","AC_FIN","AC_NFE","AC_OTH","AN_AFR","AN_AMR","AN_ASJ","AN_EAS","AN_FIN","AN_NFE","AN_OTH","AF_AFR","AF_AMR","AF_ASJ","AF_EAS","AF_FIN","AF_NFE","AF_OTH","STAR_AC","STAR_AC_raw","STAR_Hom","POPMAX","AC_POPMAX","AN_POPMAX","AF_POPMAX","DP_MEDIAN","DREF_MEDIAN","GQ_MEDIAN","AB_MEDIAN","AS_RF","AS_FilterStatus","AS_RF_POSITIVE_TRAIN","AS_RF_NEGATIVE_TRAIN","CSQ","GC_AFR","GC_AMR","GC_ASJ","GC_EAS","GC_FIN","GC_NFE","GC_OTH","Hom_Male","Hom_Female","segdup","lcr");

while(my $line=<IN>){
	chomp($line);
	#print $line."\n";
	else{
	
		$wi+=1;
		my @sets=split(/\s+/,$line);
		my $ref=$sets[3];
		my @alts=$sets[4];
		if($sets[4] =~/,/){
			@alts=split(",",$sets[4]);
			
		}
		#print "$chr\t$start\t$end\t$ref\t$sets[4]\n";
		my %names=%{INFO2map::loadINFO($sets[7])};
		my @indexs=(); ## cannot work for mutiple allele.
		my %sub_names=%{INFO2map::extractINFO(\%names,scalar(@alts),$index)};
#			if($sub_names{"AC_NFE"}==0){next;}
		#	if($sub_names{"AF"} ne "."&&$sub_names{"AF"}>0.001){next;}
			my $nalt=$alt;
                        my $nref=$ref;
                        if(length($alt)>1 && length($ref)>1){
                                if(length($alt)==length($ref)){
                                        $nref=substr($ref,0,1);
                                        $nalt=substr($alt,0,1);
                                }elsif(length($alt)>length($ref)){
                                        $nref=substr($ref,0,1);
                                        $nalt=substr($alt,0,length($alt)-length($ref)+1)
                                }else{

                                        $nalt=substr($alt,0,1);
                                        $nref=substr($ref,0,length($ref)-length($alt)+1)

                                }
                        }
                	$end=$start+length($nref)-1;
		        my $outstr="$chr\t$start\t$end\t$nref\t$nalt";
#			my $outstr="$chr\t$start\t$end\t$ref\t$alt";
		#	print "outsearch\t".$outstr."\t".$index."\ttotal\t".@alts."\n";
			foreach my $item(@header){  ## simplify the output [(0..8,30..56,92,98..@header)]
			#	print $item."\n";
				if(exists($sub_names{$item})){
					if($sub_names{$item} eq ""){$sub_names{$item}="."}
					$outstr=$outstr."\t".$sub_names{$item};
				#	if($item eq "CSQ"){
				#		print $sub_names{$item}."\n";
					
				#			my @AAVs=split(/\|/,$sub_names{$item});
						#	print @AAVs."\t".$item."\t".$sub_names{$item}."\n";	exit;
						#if(@AAVs<1){print $outstr."\n";exit;}
				#			for(my $km=0;$km<@AAVs;$km++){
				#				if(ex)
				#				$sub_names{$csq[$km]}=$AAVs[$km];
				#			}				
			#		}
				}else{
					$outstr=$outstr."\t.";
				}
			}
			my $global_accept=1;
			my $flag=0;
			#if($sub_names{"Consequence"} =~ /splice_donor|splice_acceptor_/){$flag=1;}
			#if($flag==0 && $sub_names{"Consequence"} =~ /UTR|intron|RNA|downstream|intergenic|non_coding|regulatory|TF|upstream/){next;}
			AF:{
				#print $sub_names{"Consequence"}."\n" ;
				foreach my $skey("ExAC_MAF","ExAC_AFR_MAF","ExAC_AMR_MAF","ExAC_EAS_MAF","ExAC_FIN_MAF","ExAC_NFE_MAF","ExAC_OTH_MAF","ExAC_SAS_MAF"){
					if($sub_names{$skey}) {
						my @af_sets=split(/:|\&/,$sub_names{$skey});
						if(@af_sets<2){next;}
						my $accept=1;
						for(my $jk=0;$jk<@af_sets;$jk=$jk+2){
							if($af_sets[$jk] eq $alt) {
								if($af_sets[$jk+1]>0.01){$accept=0;$global_accept=0;last AF;}
							}
							
						}
						if($accept==0){$global_accept=0;last AF;}
					}
			 }
	    }
		 if($global_accept==0){next}
		#	print "outsearch\t".$outstr."\t".$index."\ttotal\t".@alts."\n";
		  print OUT $outstr."\n";
			
		}
	}
	
}


close IN;
close OUT;
