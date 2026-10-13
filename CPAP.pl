#!usr/bin/perl
###################################################
#		CPAP - Poly(A)site Prediction Tool
###################################################

###################################################
#	Declaring file names
###################################################
$filename="./Datasets/seq_positiveDataset.fa";
$filename1="./Datasets/oligo_scores.txt";
$filename2="./Datasets/meme_test_matrices.txt";
$filename3="./Datasets/2dproperty.txt";
$filename4="./Datasets/logit_output.txt";
$filename5="./Datasets/CSmatrix.txt";

open SEQFILE,"$filename"
	or die "Could not open $filename. Aborting!\n";
open OLIGOFILE,"$filename1"
	or die "Could not open $filename1. Aborting!\n";
open MEMEFILE,"$filename2"
	or die "Could not open $filename2. Aborting!\n";
open DINUCLEOFILE,"$filename3"
	or die "Could not open $filename3. Aborting.\n";
open LOGITFILE,"$filename4"
	or die "Could not open $filename4. Aborting.\n";
open CSFILE,"$filename5"
	or die "Could not open $filename5. Aborting.\n";
open  (WRITEFILE, ">test_introns_3_Negative_output_18thDec2015.txt")
	or die "Could not open output file. Aborting!\n";
###################################################
#	Declaring Variables
###################################################
@sequences=();
@seq_headers=();
@FUE=();
@NUE=();
@CE=();
@CS=();
@DSE=();
@di_ids=();
$p_value=0;
@param_names=();
@param_values=();
@dinucleo_scores=();
@nucleotides=('A','T','G','C');
my @CS_Matrix= ([(0) x 10], [(0) x 10],[(0) x 10], [(0) x 10]);
my @tri_count= ([(0) x 64],[(0) x 64], [(0) x 64]);
my @tetra_count= ([(0) x 256],[(0) x 256], [(0) x 256]);
@meme_CE1=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_CE2=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_CE3=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_CE4=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_CE5=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_DSE1=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_DSE2=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_DSE3=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_DSE4=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_DSE5=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_FUE1=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_FUE2=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_FUE3=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_FUE4=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_FUE5=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);
@meme_NUE1=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);	
@meme_NUE5=([(0)x4],[(0)x4],[(0)x4],[(0)x4],[(0)x4]);	
my @array=([(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16],[(0)x16]);
@scores_test=();
%tri_nucleo=();
%tetra_nucleo=();
%pro_values=(AA=>0,AC=>0,AG=>0,AT=>0,CA=>0,CC=>0,CG=>0,CT=>0,GA=>0,GC=>0,GG=>0,GT=>0,TA=>0,TC=>0,TG=>0,TT=>0);
$N=0;
$FUE_length=100;
$NUE_length=40;
$CE_length=40;
$CS_length=10;
$DSE_length=50;
$di_row=0;
$di_column=0;
$string="";
###################################################
#	Reading Input Files
###################################################
#==>Reading input files and storing sequencs in respective 	s

while(<SEQFILE>)
{
	chomp($_);
	if($_=~/>/)
	{
		push @seq_headers,$_;
	}
	else
	{
		push @sequences,$_;
	}
}
$N=$#seq_headers+1;
while(<DINUCLEOFILE>)
{
	if($_=~ /ID/)
	{
	}
	else
	{
		@values=split("\t",$_);
		$string=shift @values;
		$string.=shift @values;
		$string=~ s/ /_/;
		push @di_ids,$string;
		for($u=0;$u<16;$u++)
		{
			chomp($values[$u]);
			$array[$di_row][$di_column]=$values[$u];
			$di_column++;
		}
		$di_row++;
		$di_column=0;
	}
}
while(<LOGITFILE>)
{
	if($_=~/Estimate/)
	{
	}
	else
	{
		@values1=split(/\t/,$_);
		push @param_names,$values1[0];
		push @param_values,$values1[1];
	}
}
local $/=undef;
chomp(my $whole_CS_file=<CSFILE>);
local $/=undef;
chomp(my $whole_meme_file = <MEMEFILE>);
local $/=undef;
chomp(my $whole_oligo_file = <OLIGOFILE>);
###################################################
#	Subroutines
###################################################
$temp="";
$count_tri=0;
$count_tetra=0;
for($b=0;$b<4;$b++)
{
	$temp.=$nucleotides[$b];
	for($c=0;$c<4;$c++)
	{
		$temp.=$nucleotides[$c];
		for($d=0;$d<4;$d++)
		{
			$temp.=$nucleotides[$d];
			$tri{$temp}=0;
			$tri_nucleo{$temp}=$count_tri;
			$count_tri++;
			$temp=$nucleotides[$b].$nucleotides[$c];
		}
		$temp=$nucleotides[$b];
	}
	$temp="";
}
$temp="";
for($e=0;$e<4;$e++)
{
	$temp.=$nucleotides[$e];
	for($f=0;$f<4;$f++)
	{
		$temp.=$nucleotides[$f];
		for($g=0;$g<4;$g++)
		{
			$temp.=$nucleotides[$g];
			for($h=0;$h<4;$h++)
			{
				$temp.=$nucleotides[$h];
				$tetra{$temp}=0;
				$tetra_nucleo{$temp}=$count_tetra;
				$p++;
				$temp=$nucleotides[$e].$nucleotides[$f].$nucleotides[$g];
			}
			$temp=$nucleotides[$e].$nucleotides[$f];
		}
		$temp=$nucleotides[$e];
	}
	$temp="";
}
sub oligo
{
	@values2=split("TETRA",$whole_oligo_file);
	@lines1=split(/\n/,$values2[0]);
	$string=shift @lines1;		
	for($b=0;$b<3;$b++)
	{
		@counts=split(/\t/,$lines1[$b]);
		for($d=0;$d<64;$d++)
		{
			$tri_count[$b][$d]=$counts[$d];
		}
	}
	@lines2=split(/\n/,$values2[1]);
	$string=shift @lines2;
	$string=shift @lines2;
	for($c=0;$c<3;$c++)
	{
		@counts1=split(/\t/,$lines2[$c]);
		for($e=0;$e<256;$e++)
		{
				$tetra_count[$c][$e]=$counts1[$e];
		}
	}
}
sub meme_matrices
{
	@matrices=split(/\#/,$whole_meme_file);
	for($a=1;$a<=2;$a++)
	{
		@lines=split("\n",$matrices[$a]);
		for($b=1;$b<=$#lines;$b++)
		{
			chomp($lines[$b]);
			@columns=split(/\s+/,$lines[$b]);
			$row1=$b-1;
			for($k=0;$k<=$#columns;$k++)
			{
				if($a==1)
				{
					$meme_NUE1[$row1][$k]=$columns[$k];
				}
				elsif($a==2)
				{
					$meme_NUE5[$row1][$k]=$columns[$k];
				}		
			}
		}
	}
	for($a=3;$a<=7;$a++)
	{
		@lines=split("\n",$matrices[$a]);
		for($b=1;$b<=$#lines;$b++)
		{
			chomp($lines[$b]);
			@columns=split(/\s+/,$lines[$b]);
			$row1=$b-1;
			for($k=0;$k<=$#columns;$k++)
			{
				if($a==3)
				{
					$meme_FUE1[$row1][$k]=$columns[$k];
				}
				elsif($a==4)
				{
					$meme_FUE2[$row1][$k]=$columns[$k];
				}
				elsif($a==5)
				{
					$meme_FUE3[$row1][$k]=$columns[$k];
				}
				elsif($a==6)
				{
					$meme_FUE4[$row1][$k]=$columns[$k];
				}
				elsif($a==7)
				{
					$meme_FUE5[$row1][$k]=$columns[$k];
				}			
			}
		}
	}
	for($a=8;$a<=12;$a++)
	{
		@lines=split("\n",$matrices[$a]);
		for($b=1;$b<=$#lines;$b++)
		{
			chomp($lines[$b]);
			@columns=split(/\s+/,$lines[$b]);
			$row1=$b-1;
			for($k=0;$k<=$#columns;$k++)
			{
				if($a==8)
				{
					$meme_DSE1[$row1][$k]=$columns[$k];
					
				}
				elsif($a==9)
				{
					$meme_DSE2[$row1][$k]=$columns[$k];
				}
				elsif($a==10)
				{
					$meme_DSE3[$row1][$k]=$columns[$k];
				}
				elsif($a==11)
				{
					$meme_DSE4[$row1][$k]=$columns[$k];
				}
				elsif($a==12)
				{
					$meme_DSE5[$row1][$k]=$columns[$k];
				}			
			}
		}
	}
	for($a=13;$a<=17;$a++)
	{
		@lines=split("\n",$matrices[$a]);
		for($b=1;$b<=$#lines;$b++)
		{
			chomp($lines[$b]);
			@columns=split(/\s+/,$lines[$b]);
			$row1=$b-1;
			for($k=0;$k<=$#columns;$k++)
			{
				if($a==13)
				{
					$meme_CE1[$row1][$k]=$columns[$k];
				}
				elsif($a==14)
				{
					$meme_CE2[$row1][$k]=$columns[$k];
				}
				elsif($a==15)
				{
					$meme_CE3[$row1][$k]=$columns[$k];
				}
				elsif($a==16)
				{
					$meme_CE4[$row1][$k]=$columns[$k];
				}
				elsif($a==17)
				{
					$meme_CE5[$row1][$k]=$columns[$k];
				}			
			}
		}
	}
}
sub CSmatrix
{
	@rows_CS=split(/\n/,$whole_CS_file);
	for($row=0;$row<4;$row++)
	{
		@cs_values=split(/\t/,$rows_CS[$row]);
		for($column=0;$column<10;$column++)
		{
			$CS_Matrix[$row][$column]=$cs_values[$column];
		}
	}
}

###################################################
#	generating all matrices for scoring from 
#	respective input files
###################################################
oligo();
meme_matrices();
CSmatrix();
###################################################
#	Subroutines for Scoring the sequences
###################################################
# "@meme_DSE2\n";
sub dinucleo_score
{
	$seq=shift();
	$start=shift();
	$row_dinu=shift();
	$length=shift();
	$index=0;
	$sum=0;
	$score=0;
	$subSeq="";
	$dinucl="";
	@counts=();
	foreach $value(sort keys %pro_values)
	{
		$pro_values{$value}=$array[$row_dinu][$index];
		$index++;
	}
	for ($h=$start;$h<=($length-30);$h++)
	{
		#==>getting 50 nucleotide seq at a time
		$subSeq=substr($seq,$h,30);
		my @nucleo= split //, $subSeq;
		#==>counting occurence of each possible dinucleo combination of 
		#=>A,T,G,C with eachother and storing the count in a hash
		for($i=0;$i<$#nucleo;$i++)
		{
			$dinucl=$nucleo[$i]."".$nucleo[$i+1];
			$counter{$dinucl}=$counter{$dinucl}+1;
		}
		#==>finding multiplying each count by its repective value of 
		#==>dinucleo property and taking their sum in the end and dividing by 49		
		foreach $value (sort keys %counter)
		{
			$sum=$sum+($counter{$value}*$pro_values{$value});
		}
		$count=$sum/29;
		#==>summing values at each location of all sequences in to a separate 	
		#==>to get average value of property value in the end 
		push @counts,$count;
		#==> setting values to zero for next iteration
		$sum=0;
		foreach $value1 (keys %counter)
		{
			$counter{$value1}=0;
		}
	}
	return @counts;
}
sub tri_score
{
	$seq1=shift();
	$length2=shift();
	$row2=shift();
	$tri="";
	$sum3=0;
	for($j=0;$j<($length2-3);$j++)
	{
		$tri=substr($seq1,$j,3);
		$column1=$tri_nucleo{$tri};
		$sum3+=$tri_count[$row2][$column1];		
	}
	return $sum3;
}
sub tetra_score
{
	$seq=shift();
	$length=shift();
	$row3=shift();
	$tetra="";
	$sum1=0;
	for($k=0;$k<($length-4);$k++)
	{
		$tetra=substr($seq,$k,4);
		$column2=$tetra_nucleo{$tetra};
		$sum1+=$tetra_count[$row3][$column2];		
	}
	return $sum1;
}
sub cs_score
{
	$seq=shift();
	@csMatrix=@{shift()};
	$score1=0;
	for($l=0;$l<10;$l++)
	{
		$nucleo=substr($seq,$l,1);
		
		if($nucleo eq 'A')
		{
			$score1+=$csMatrix[0][$l];
		}
		elsif($nucleo eq 'C')
		{
			$score1+=$csMatrix[3][$l];
		}
		elsif($nucleo eq 'G')
		{
			$score1+=$csMatrix[2][$l];
		}
		elsif($nucleo eq 'T')
		{
			$score1+=$csMatrix[1][$l];
		}
	}
	return $score1;
}
sub meme_score
{
	$seq=shift();
	$length=shift();
	@meme_matrix=@{shift()};
	$size=$#meme_matrix+1;
	$max=0;
	$temp=0;
	$seq_temp="";
	my @arr = split //,$seq;
	for($p=0;$p<($length-$size);$p++)
	{
		for($q=$p;$q<($p+$size);$q++)
		{
			$pentamer.=$arr[$q];
		}
		for($index=0;$index<$size;$index++)
		{
			$nucleo=substr($pentamer,$index,1);
			if($nucleo eq 'A')
			{
				$temp+=$meme_matrix[$index][0];
			}
			elsif($nucleo eq 'C')
			{
				$temp+=$meme_matrix[$index][1];
			}
			elsif($nucleo eq 'G')
			{
				$temp+=$meme_matrix[$index][2];
			}
			elsif($nucleo eq 'T')
			{
				$temp+=$meme_matrix[$index][3];
			}
		}
		if($temp>$max)
		{
			$max=$temp;
			$seq_temp=$pentamer;
		}
		$temp=0;
		$pentamer="";
	}
	return $max;
}
###################################################
#	Scoring the sequences
###################################################
$FUE_score=0;
$NUE_score=0;
$CE_score=0;
$NUE_CE_score=0;
for($m=0;$m<$N;$m++)
{
	print WRITEFILE "$seq_headers[$m]\n";
	print WRITEFILE "$sequences[$m]\n";
	for($zeroes=0;$zeroes<140;$zeroes++)
	{
		print WRITEFILE "0,";
	}
	#==>Extracting NUE,FUE,CE,DSE,CS from the dataset foreach sequence
	for($site=140;$site<351;$site++)
	{
		push @FUE,substr($sequences[$m],$site-140,100);
		push @NUE,substr($sequences[$m],$site-40,40);
		push @CE,substr($sequences[$m],$site-20,40);
		push @DSE,substr($sequences[$m],$site,50);
		push @CS,substr($sequences[$m],$site-5,10);
	}
	#==>Computing Cleavage Site Matrix for Positive Sequences
	for($n=0;$n<$#FUE;$n++)
	{
		push @scores_test,meme_score($NUE[$n],$NUE_length,\@meme_NUE1);
		push @scores_test,meme_score($NUE[$n],$NUE_length,\@meme_NUE5);
		push @scores_test,meme_score($FUE[$n],$FUE_length,\@meme_FUE1);
		push @scores_test,meme_score($FUE[$n],$FUE_length,\@meme_FUE2);
		push @scores_test,meme_score($FUE[$n],$FUE_length,\@meme_FUE3);
		push @scores_test,meme_score($FUE[$n],$FUE_length,\@meme_FUE4);
		push @scores_test,meme_score($FUE[$n],$FUE_length,\@meme_FUE5);
		push @scores_test,meme_score($DSE[$n],$DSE_length,\@meme_DSE1);
		push @scores_test,meme_score($DSE[$n],$DSE_length,\@meme_DSE2);
		push @scores_test,meme_score($DSE[$n],$DSE_length,\@meme_DSE3);
		push @scores_test,meme_score($DSE[$n],$DSE_length,\@meme_DSE4);
		push @scores_test,meme_score($DSE[$n],$DSE_length,\@meme_DSE5);
		push @scores_test,meme_score($CE[$n],$CE_length,\@meme_CE1);
		push @scores_test,meme_score($CE[$n],$CE_length,\@meme_CE2);
		push @scores_test,meme_score($CE[$n],$CE_length,\@meme_CE3);
		push @scores_test,meme_score($CE[$n],$CE_length,\@meme_CE4);
		push @scores_test,meme_score($CE[$n],$CE_length,\@meme_CE5);
		push @scores_test,cs_score($CS[$n],\@CS_Matrix);
		push @scores_test,tri_score($NUE[$n],$NUE_length,0);
		push @scores_test,tri_score($FUE[$n],$FUE_length,1);
		push @scores_test,tri_score($DSE[$n],$DSE_length,2);
		push @scores_test,tetra_score($NUE[$n],$NUE_length,0);
		push @scores_test,tetra_score($FUE[$n],$FUE_length,1);
		push @scores_test,tetra_score($CE[$n],$CE_length,2);
		for($property=0;$property<=$#di_ids;$property++)
		{
			@dinucleo_scores=dinucleo_score($sequences[$m],$n,$property,190);
			for($rowFUE=0;$rowFUE<90;$rowFUE++)
			{
				$FUE_score+=$dinucleo_scores[$rowFUE];
			}
			for($rowNUE=90;$rowNUE<135;$rowNUE++)
			{
				$NUE_score+=$dinucleo_scores[$rowNUE];
			}
			for($rowCE=135;$rowCE<161;$rowCE++)
			{
				$CE_score+=$dinucleo_scores[$rowCE];
			}
			$NUE_CE_score=($NUE_score+$CE_score)/2;
			push @scores_test,$FUE_score;
			push @scores_test,$NUE_CE_score;			
			$FUE_score=0;
			$NUE_score=0;
			$CE_score=0;
			$NUE_CE_score=0;
			@dinucleo_scores=();
		}
		# print WRITEFILE $n+141,"\n";
		# foreach $value (@scores_test)
		# {
			# print WRITEFILE "$value\t";
		# }
		# print WRITEFILE "\n";
		@dinucleo_scores=();
		for($ind=1;$ind<=$#param_values;$ind++)
		{
			$p_value+=$param_values[$ind]*$scores_test[$ind-1];
		}
		$probability=1/(1+exp(-($p_value+$param_values[0])));
		#print  WRITEFILE $n+141,"\t";
		print WRITEFILE $probability.",";
		$p_value=0;
		@scores_test=();
	}	
	$p_value=0;
	@dinucleo_scores=();
	@scores_test=();
	@FUE=();
	@NUE=();
	@CE=();
	@DSE=();
	@CS=();
	$FUE_score=0;
	$NUE_score=0;
	$CE_score=0;
	$NUE_CE_score=0;
	for($zeroes=0;$zeroes<49;$zeroes++)
	{
		print WRITEFILE "0,";
	}
	print WRITEFILE "0";
	print WRITEFILE "\n";
}

###################################################
#	Closing all files
###################################################
close SEQFILE;
close MEMEILE;
close LOGITFILE;
close DINUCLEOFILE;
close OLIGOFILE;
close WRITEFILE;