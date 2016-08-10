package parseMpileup;
use warnings;
use strict;
use Data::Dumper;
use File::stat;
use List::Util qw[min max];

sub parse()
{
  my ($processName,$pileup,$filt,$param,$istumor)=@_;
  #print "istumor=$istumor\n";
  my @pileupLines=split(/\n/,$pileup);
  my $pileupData=[];
  my @NT=('A', 'C', 'G', 'T');
  my @currSeqIdx;
  my $nextSeqIdx;
  my $mismatchCount={};
  my $seqLength={};
  my $seqIdx={};
  my $first=1;
  my $i=0;
  #print Dumper($pileup);
  foreach my $line (@pileupLines)
  {
    #print $line . "\n";
    my ($chr,$pos,$ref,$readCount,$bases,$baseQual,$mapQual,$readPos)=split(/\t/,$line);
    my $currBases=$bases;
    my $currBaseQual=$baseQual;
    $currBaseQual=~s/\\\\/\t\t/g;
    $currBaseQual=~s/\t/\\/g;
    $currBaseQual=~s/\\'/'/g;
    #print "baseQual\n$currBaseQual\n";
    my $currMapQual=$mapQual;
    my @readPosArray=split(/\,/,$readPos);
    if($first)
    {
      @currSeqIdx=(1..$readCount);
      $nextSeqIdx=$readCount+1;
      $currBases=~s/\^.//g;
      $first=0;
    }
    elsif($pileupData->[$i-1]->{'Pos'}!=$pos-1)
    {
      @currSeqIdx=($currSeqIdx[-1]...$currSeqIdx[-1]+$readCount);
    }

    my $j=0;
    #print "pos: $pos depth: $readCount parsing: $bases\n";
    my %ntCount=();
    my %ntCountStrand=();
    my %ntBQsum=();
    my %ntMQsum=();
    my %ntReadPosSum=();
    my $lowQCBaseCount=0;
    foreach my $nt (@NT)
    {
      #print "$nt at pos: ";
      $ntCount{$nt}=0;
      $ntCountStrand{$nt}=0;
      $ntCountStrand{lc($nt)}=0;
      $ntBQsum{$nt}=0;
      $ntMQsum{$nt}=0;
      $ntReadPosSum{$nt}=0;
    }

    while($currBases)
    {
     #print "at $j: $currBases\n";
    #print "curr base qual has length" . length($currBaseQual) . "\n$currBaseQual\n";
    #print "curr map qual has length" . length($currMapQual) . "\n$currMapQual\n";
    #print "curr bases has length" . length($currBases) . "\n$currBases\n";
     #print "BQ: $currBaseQual\n";
     #print "MQ: $currMapQual\n";
      if($currBases=~/^\^/)
      {
        splice(@currSeqIdx,$j,0,$nextSeqIdx);
        #print "added new seq $nextSeqIdx\n";
        $nextSeqIdx++;
        $currBases=substr($currBases,2);
      }
      my $currBQ=min(ord(substr($currBaseQual,0,1))-$param->{'BQoffset'},$param->{'maxBQ'});
      my $currMQ=ord(substr($currMapQual,0,1))-$param->{'MQoffset'};
     #print "BQ is $currBQ and MQ is $currMQ with $currBases\n";
      if($currBases=~/^(.)([\+\-])(\d+)/gc)
      {
        #print "1: $1 2: $2 3: $3\n";
        my $indelStr=$2 . $3 . substr($currBases,pos($currBases),$3);
        my $indelLen=$3;
        my $b=$1;
        if($2 =~/\+/)
        {
          $seqLength->{$currSeqIdx[$j]}+=$indelLen;
        }
        #else
        #{
        #  $seqLength->{$currSeqIdx[$j]}-=$indelLen;
        #}
        if($currBQ>=$param->{'minBQ'} && $currMQ>=$param->{'minMQ'} && $b=~/([AGCTagct])/)
        {
          $mismatchCount->{$currSeqIdx[$j]}++;
          $seqIdx->{$pos}->{uc($b)}.=$currSeqIdx[$j] . ';';
          $ntCount{uc($b)}++;
          $ntBQsum{uc($b)}+=$currBQ;
          $ntMQsum{uc($b)}+=$currMQ;
          $ntReadPosSum{uc($b)}+=$readPosArray[0];
          $ntCountStrand{$b}++;
        }
        if($currMQ>=$param->{'minMQ'})
        {
          $mismatchCount->{$currSeqIdx[$j]}+=$indelLen;
          $seqIdx->{$pos}->{uc($indelStr)}.=$currSeqIdx[$j] . ';';
          if(exists($ntCount{uc($indelStr)}))
          {
            $ntCount{uc($indelStr)}++;
            $ntBQsum{uc($indelStr)}+=$param->{'defaultBQ'};
            $ntMQsum{uc($indelStr)}+=$currMQ;
            $ntReadPosSum{uc($indelStr)}+=$readPosArray[0];
            #print "currMQ is $currMQ and currBQ is $currBQ for $indelStr at $pos\n";
          }
          else
          {
            $ntCount{uc($indelStr)}=1;
            $ntBQsum{uc($indelStr)}=$param->{'defaultBQ'};
            $ntMQsum{uc($indelStr)}=$currMQ;
            $ntReadPosSum{uc($indelStr)}=$readPosArray[0];
            #print "currMQ is $currMQ and currBQ is $currBQ for $indelStr at $pos\n";
          }
          if(exists($ntCountStrand{$indelStr}))
          {
            $ntCountStrand{$indelStr}++
          }
          else
          {
            $ntCountStrand{$indelStr}=1;
          }
        }
        else
        {
          $lowQCBaseCount++;
        }
        $currBases=substr($currBases,length($indelStr)+1);
      }
      elsif($currBases=~/^([\.\,])/)
      {
        #print "found ref $ref at $j seqIdx=$currSeqIdx[$j]\n";
        #print "currMQ is $currMQ and currBQ is $currBQ for $1 at $pos\n";
        $seqLength->{$currSeqIdx[$j]}++;
        if($currBQ>=$param->{'minBQ'} && $currMQ>=$param->{'minMQ'})
        {
          $seqIdx->{$pos}->{uc($ref)}.=$currSeqIdx[$j] . ';';
          $ntCount{$ref}++;
          $ntBQsum{$ref}+=$currBQ;
          $ntMQsum{$ref}+=$currMQ;
          $ntReadPosSum{$ref}+=$readPosArray[0];
          #print "currMQ is $currMQ and currBQ is $currBQ for $1 at $pos\n";
          if($1=~/\./)
          {
            $ntCountStrand{$ref}++;
          }
          else
          {
            $ntCountStrand{lc($ref)}++;
          }
        }
        else
        {
          $lowQCBaseCount++;
        }
        $currBases=substr($currBases,1);
      }
      elsif($currBases=~/^([AGCTagct])/)
      {
        $seqLength->{$currSeqIdx[$j]}++;
        if($currBQ>=$param->{'minBQ'} && $currMQ>=$param->{'minMQ'})
        {
          $mismatchCount->{$currSeqIdx[$j]}++;
          $seqIdx->{$pos}->{uc($1)}.=$currSeqIdx[$j] . ';';
          $ntCount{uc($1)}++;
          $ntBQsum{uc($1)}+=$currBQ;
          $ntMQsum{uc($1)}+=$currMQ;
          $ntReadPosSum{uc($1)}+=$readPosArray[0];
          $ntCountStrand{$1}++;
          #print "currMQ is $currMQ and currBQ is $currBQ for $1 at $pos\n";
        }
        else
        {
          $lowQCBaseCount++;
        }
        $currBases=substr($currBases,1);
      }
      else
      {
        $seqLength->{$currSeqIdx[$j]}++;
        if($currMQ>=$param->{'minMQ'})
        {
          $ntCount{'*'}++;
          $ntMQsum{'*'}+=$currMQ;
          #$ntBQsum{'*'}+=$currBQ;
          #print "currMQ is $currMQ and currBQ is $currBQ for $1 at $pos\n";
        }
        else
        {
          $lowQCBaseCount++;
        }
        $currBases=substr($currBases,1);
      }
      if (length($currBaseQual)==0)
      {
        print "at $pos empty base qual\n";
      }
      $currBaseQual=substr($currBaseQual,1);
      $currMapQual=substr($currMapQual,1);
      shift(@readPosArray);
      if($currBases=~/^\$/)
      {
        #print "removing $currSeqIdx[$j]\n";
        $currBases=substr($currBases,1);
        splice(@currSeqIdx,$j,1);
        #print Dumper(@currSeqIdx);
      }
      else
      {
        $j++;
      }
    }
    #print "currentSeqIdx: " . Dumper(@currSeqIdx);
    #print "seqIdx: " . Dumper($seqIdx);
    #print "seqLength: " . Dumper($seqLength);
    #print "mismatchCount: " . Dumper($mismatchCount);
    #unless (length($baseQual)==$readCount)
    #{
    #print "found " . length($baseQual) . " BQ " . length ($mapQual) . " MQ for RD" . $readCount . " in:\n $line\n";
    #print Dumper(ord($baseQual));
    #}

    $pileupData->[$i]->{'Chr'}=$chr;
    $pileupData->[$i]->{'Pos'}=$pos;
    $pileupData->[$i]->{'readCount'}=$readCount;
    $pileupData->[$i]->{'readCountPass'}=$readCount-$lowQCBaseCount;
    #$pileupData->[$i]->{'RefCountF'}=($bases =~ tr/\.//);
    #$pileupData->[$i]->{'RefCountR'}=($bases =~ tr/\,//);
    $pileupData->[$i]->{'RefCountFpass'}=$ntCountStrand{$ref};
    $pileupData->[$i]->{'RefCountRpass'}=$ntCountStrand{lc($ref)};
    if(exists($ntCount{'*'}))
    {
      $pileupData->[$i]->{'delCount'}=$ntCount{'*'};
      $pileupData->[$i]->{'delMQ'}=$ntMQsum{'*'}/$ntCount{'*'};
      #$pileupData->[$i]->{'delBQ'}=$ntBQsum{'*'}/$ntCount{'*'};
      delete($ntCount{'*'});
    }
    $pileupData->[$i]->{'MismatchCountF'}=0;
    $pileupData->[$i]->{'MismatchCountR'}=0;
    #print "ntCountStrand for $pos:\n";
    #print Dumper(%ntCountStrand);
    foreach $b (keys %ntCountStrand)
    {
      if(uc($b) eq $ref)
      {
        next;
      }
      elsif($b eq uc($b))
      {
        $pileupData->[$i]->{'MismatchCountF'}+=$ntCountStrand{$b};
        #print "found f mismatch $b\n";
      }
      else
      {
        $pileupData->[$i]->{'MismatchCountR'}+=$ntCountStrand{$b};
        #print "found r mismatch $b\n"
      }
    }

    my @ntSort = (sort { $ntCount{$b} <=> $ntCount{$a} } keys %ntCount);

    #print "ntCountF: " . Dumper(%ntCountF);
    #print "ntCountF: " . Dumper(%ntCountR);
    #print Dumper($pileupData->[$i]);
    if(!exists($ntCountStrand{$ntSort[1]}))
    {
      $ntCountStrand{$ntSort[1]}=0;
    }
    if(!exists($ntCountStrand{lc($ntSort[1])}))
    {
      $ntCountStrand{lc($ntSort[1])}=0;
    }
    if(!exists($ntCountStrand{$ntSort[0]}))
    {
      $ntCountStrand{$ntSort[0]}=0;
    }
    if(!exists($ntCountStrand{lc($ntSort[0])}))
    {
      $ntCountStrand{lc($ntSort[0])}=0;
    }
    #print "ntSort: " . Dumper(@ntSort);
    #print "ntCountStrand: " . Dumper(%ntCountStrand);
    #if ($filt>0 && ($ntCountStrand{$ntSort[1]}/($readCount-$lowQCBaseCount)<$param->{'minPercentStrand'}  || $ntCountStrand{lc($ntSort[1])}/($readCount-$lowQCBaseCount)<$param->{'minPercentStrand'} || $ntCountStrand{$ntSort[0]}/($readCount-$lowQCBaseCount)<$param->{'minPercentStrand'}  || $ntCountStrand{lc($ntSort[0])}/($readCount-$lowQCBaseCount)<$param->{'minPercentStrand'} || ))
    #{
    # #print "strand biase filter $bases\n";
    # #print Dumper(%ntCountStrand);
    # #print Dumper($param);
    # #print Dumper(@ntSort);
    # #print "readCount $readCount, lowQCBaseCount $lowQCBaseCount\n";
    #  $i++;
    #  next;
    #}
    #print "ntBQsum: " . Dumper(%ntBQsum);
    #print "ntSort: " . Dumper(@ntSort);
    my $indelCount=0;
    my @tempNTSort;
    foreach my $j (0 .. $#ntSort)
    {
      if($ntSort[$j]=~/^[\-\+]/)
      {
        $indelCount++;
        #print "indelCount: $indelCount old " . Dumper(@ntSort);
        unless($indelCount>1)
        {
          push(@tempNTSort,$ntSort[$j]);
        }
      }
      else
      {
        push(@tempNTSort,$ntSort[$j]);
      }
    }

    #print "ntSort: " . Dumper(@ntSort);
    #print "tempNTSort: " . Dumper(@tempNTSort);
    @ntSort=@tempNTSort;
    if(($ntSort[0] eq $ref && $ntCount{$ntSort[1]}<$param->{'minBCount'} || $ntCount{$ntSort[0]}<$param->{'minBCount'}) && $filt>0)
    {
      $i++;
      next;
    }

    $pileupData->[$i]->{'Astring'}=$ntSort[0];
    if($ntCount{$ntSort[1]}>0)
    {
      $pileupData->[$i]->{'Bstring'}=$ntSort[1];
    }
    elsif($ntSort[0] eq $ref)
    {
      $pileupData->[$i]->{'Bstring'}='N';
    }
    else
    {
      $pileupData->[$i]->{'Bstring'}=$ref;
    }

    $pileupData->[$i]->{'RefString'}=$ref;
    if($ntSort[0]=~/^\-(\d+)/)
    {

     #print "found deletion: " . Dumper(@ntSort);
      $pileupData->[$i]->{'Ref'}=$ref . substr($ntSort[0],length($ntSort[0])-$1);
      $pileupData->[$i]->{'A'}=$ref;
      $pileupData->[$i]->{'B'}=$ntSort[1] . substr($ntSort[0],length($ntSort[0])-$1);
    }
    elsif($ntSort[1]=~/^\-(\d+)/)
    {
       #print "found deletion: " . Dumper(@ntSort);
      $pileupData->[$i]->{'Ref'}=$ref . substr($ntSort[1],length($ntSort[1])-$1);
      $pileupData->[$i]->{'A'}=$ntSort[0] . substr($ntSort[1],length($ntSort[1])-$1);
      $pileupData->[$i]->{'B'}=$ref;
    }
    elsif($ntSort[0]=~/^\+(\d+)/)
    {
       #print "found insertion: " . Dumper(@ntSort);
      $pileupData->[$i]->{'Ref'}=$ref;
      $pileupData->[$i]->{'A'}=$ref . substr($ntSort[0],length($ntSort[0])-$1);
      $pileupData->[$i]->{'B'}=$ntSort[1];
    }
    elsif($ntSort[1]=~/^\+(\d+)/)
    {
       #print "found insertion: " . Dumper(@ntSort);
      $pileupData->[$i]->{'Ref'}=$ref;
      $pileupData->[$i]->{'A'}=$ntSort[0];
      $pileupData->[$i]->{'B'}=$ref . substr($ntSort[1],length($ntSort[1])-$1);
    }
    elsif($ntCount{$ntSort[1]}>0)
    {
      $pileupData->[$i]->{'Ref'}=$ref;
      $pileupData->[$i]->{'A'}=$ntSort[0];
      $pileupData->[$i]->{'B'}=$ntSort[1];
    }
    elsif($ntSort[0] eq $ref)
    {
      $pileupData->[$i]->{'Ref'}=$ref;
      $pileupData->[$i]->{'A'}=$ntSort[0];
      $pileupData->[$i]->{'B'}='N';
    }
    else
    {
      $pileupData->[$i]->{'Ref'}=$ref;
      $pileupData->[$i]->{'A'}=$ntSort[0];
      $pileupData->[$i]->{'B'}=$ref;
    }

    $pileupData->[$i]->{'AcountF'}=$ntCountStrand{$ntSort[0]};
    $pileupData->[$i]->{'AcountR'}=$ntCountStrand{lc($ntSort[0])};
    if($ntCount{$ntSort[0]}>0)
    {
      #print "Calculating A BQ for $ntSort[0] sum= $ntBQsum{$ntSort[0]} count= $ntCount{$ntSort[0]}\n";
      $pileupData->[$i]->{'AmeanBQ'}=$ntBQsum{$ntSort[0]}/$ntCount{$ntSort[0]};
      $pileupData->[$i]->{'AmeanMQ'}=$ntMQsum{$ntSort[0]}/$ntCount{$ntSort[0]};
      $pileupData->[$i]->{'AmeanReadPos'}=$ntReadPosSum{$ntSort[0]}/$ntCount{$ntSort[0]};
    }
    else
    {
      $pileupData->[$i]->{'AmeanBQ'}=$param->{'defaultBQ'};
      $pileupData->[$i]->{'AmeanMQ'}=$param->{'defaultMQ'};
      $pileupData->[$i]->{'AmeanReadPos'}=0;
    }
    $pileupData->[$i]->{'BcountF'}=$ntCountStrand{$ntSort[1]};
    $pileupData->[$i]->{'BcountR'}=$ntCountStrand{lc($ntSort[1])};
    if($ntCount{$ntSort[1]}>0)
    {
      #print "Calculating B BQ for $ntSort[1] sum= $ntBQsum{$ntSort[1]} count= $ntCount{$ntSort[1]}\n";
      $pileupData->[$i]->{'BmeanBQ'}=$ntBQsum{$ntSort[1]}/$ntCount{$ntSort[1]};
      $pileupData->[$i]->{'BmeanMQ'}=$ntMQsum{$ntSort[1]}/$ntCount{$ntSort[1]};
      $pileupData->[$i]->{'BmeanReadPos'}=$ntReadPosSum{$ntSort[1]}/$ntCount{$ntSort[1]};
    }
    else
    {
      $pileupData->[$i]->{'BmeanBQ'}=$param->{'defaultBQ'};
      $pileupData->[$i]->{'BmeanMQ'}=$param->{'defaultMQ'};
      $pileupData->[$i]->{'BmeanReadPos'}=0;
    }
    $i++;
  }
  #print "SeqIdx: " . Dumper($seqIdx);
  #print "pileupData: " . Dumper($pileupData);
  if($istumor)
  {
    for (my $i = 0; $i <= $#$pileupData; $i++)
    {
      unless(exists($pileupData->[$i]->{'Ref'}) && exists($seqIdx->{$pileupData->[$i]->{'Pos'}}->{$pileupData->[$i]->{'Astring'}}))
      {
        next;
      }
      #print "analyzing pileup of $i " . Dumper($pileupData->[$i]);
      #print Dumper($seqIdx->{$pileupData->[$i]->{'Pos'}});
      my @AseqIdx=split(/\;/,$seqIdx->{$pileupData->[$i]->{'Pos'}}->{$pileupData->[$i]->{'Astring'}});
      my $mismatchSum=0;
      my $lengthSum=0;
      ###if doesn't exist 0 mismatches
      #print "mismatch counts: " . Dumper($mismatchCount);
      #print "seqIdx: " . Dumper($seqIdx);
      foreach my $idx (@AseqIdx)
      {
        ####need to subtract variant from mismatch count
        if(exists($mismatchCount->{$idx}))
        {
          $mismatchSum+=$mismatchCount->{$idx};
          #print "$idx has $mismatchCount->{$idx} with length $seqLength->{$idx}\n";
        }
	#print "$idx has length " . $seqLength->{$idx} . "\n";
        $lengthSum+=$seqLength->{$idx};
      }
      if($pileupData->[$i]->{'RefString'} eq $pileupData->[$i]->{'Astring'})
      {
        #print "A variant is Ref with $mismatchSum mismatches out of $lengthSum\n";
        #if($lengthSum>0)
        #{
          $pileupData->[$i]->{'AmeanPMM'}=$mismatchSum/$lengthSum;
        #}
        #else
        #{
         # print "lengthSum: $lengthSum\n";
        #  print "pileup: " . Dumper($pileupData->[$i]);
        #  print "AseqIdx: " . Dumper(@AseqIdx);
        #  print "seqlength: " . Dumper($seqLength);
        #}
      }
      elsif($pileupData->[$i]->{'Astring'}=~/^([\-\+])(\d+)/)
      {
        #print "1: $1 2: $2\n";
        my $sign=$1;
        my $indelLen=$2;
        #print "A variant is indel, old BcountF" . $pileupData->[$i]->{'BcountF'} . ", old BcountR" . $pileupData->[$i]->{'BcountR'} . "\n";
        #if($lengthSum>0)
        #{
          $pileupData->[$i]->{'AmeanPMM'}=($mismatchSum-$2*($pileupData->[$i]->{'AcountF'}+$pileupData->[$i]->{'AcountR'}))/$lengthSum;
        #}
        #else
        #{
        #  print "pileup: " . Dumper($pileupData->[$i]);
        #  print "AseqIdx: " . Dumper(@AseqIdx);
        #  print "seqlength: " . Dumper($seqLength);
        #}

        if($sign=~/\-/ && $pileupData->[$i]->{'Bstring'} eq $pileupData->[$i]->{'RefString'})
        {
          #my $readCountSum=0;
          #my $readCountPassSum=0;
          #my $refSumFpass=0;
          #my $refSumRpass=0;
          #my $delCountSum=0;
          #my $delMQSum=0;
          #my $delBQSum=0;
          my $misMatchSumF=0;
          my $misMatchSumR=0;
          for(my $j=$i+1;$j<=$i+$indelLen;$j++)
          {
            if(exists($pileupData->[$j]))
            {
              $misMatchSumF+=$pileupData->[$j]->{'MismatchCountF'};
              $misMatchSumR+=$pileupData->[$j]->{'MismatchCountR'};
              #print Dumper($pileupData->[$j]);
              #print "mismatch sumF $misMatchSumF mismatch sumR $misMatchSumR\n";
              #$refSumFpass+=$pileupData->[$j]->{'RefCountFpass'};
              #$refSumRpass+=$pileupData->[$j]->{'RefCountRpass'};
              #$readCountSum+=$pileupData->[$j]->{'readCount'};
              #$readCountPassSum+=$pileupData->[$j]->{'readCountPass'};
              #$delCountSum+=$pileupData->[$j]->{'delCount'};
              #$delMQSum+=$pileupData->[$j]->{'delMQ'};
              #$delBQSum=+$pileupData->[$j]->{'delBQ'};
            }
          }
          $pileupData->[$i]->{'BcountF'}=max($pileupData->[$i]->{'BcountF'}-$misMatchSumF/$indelLen,0);
          $pileupData->[$i]->{'BcountR'}=max($pileupData->[$i]->{'BcountR'}-$misMatchSumR/$indelLen,0);
          #if($pileupData->[$i]->{'AcountF'}+$pileupData->[$i]->{'AcountR'}!=$delCountSum/$indelLen)
          #{
          #  my $delCountDiff=$pileupData->[$i]->{'AcountF'}+$pileupData->[$i]->{'AcountR'}-$delCountSum/$indelLen;
          #  $pileupData->[$i]->{'AcountF'}-=$delCountDiff/2;
          #  $pileupData->[$i]->{'AcountR'}-=$delCountDiff/2;
          #}
          #$pileupData->[$i]->{'readCount'}=$readCountSum/$indelLen;
          #$pileupData->[$i]->{'readCountPass'}=$readCountPassSum/$indelLen;
          #$pileupData->[$i]->{'AmeanMQ'}=$delMQSum/$indelLen;          #print "new BCountF" . $pileupData->[$i]->{'BcountF'} . "\n";
          #$pileupData->[$i]->{'AmeanBQ'}=$delBQSum/$indelLen;
        }
      }
      else
      {
        #print "A variant is SNV " . $pileupData->[$i]->{'Astring'} . "and ref is" . $pileupData->[$i]->{'Ref'} . "\n";
        #if($lengthSum>0)
        #{
          $pileupData->[$i]->{'AmeanPMM'}=($mismatchSum-$pileupData->[$i]->{'AcountF'}-$pileupData->[$i]->{'AcountR'})/$lengthSum;
        #}
        #else
        #{
        #  print "pileup: " . Dumper($pileupData->[$i]);
        #  print "AseqIdx: " . Dumper(@AseqIdx);
        #  print "seqlength: " . Dumper($seqLength);
        #}
      }
      #print Dumper($pileupData->[$i]);
      my @BseqIdx=();
      if (defined($seqIdx->{$pileupData->[$i]->{'Pos'}}->{$pileupData->[$i]->{'Bstring'}}))
      {
        @BseqIdx=split(/\;/,$seqIdx->{$pileupData->[$i]->{'Pos'}}->{$pileupData->[$i]->{'Bstring'}});
      }
      $mismatchSum=0;
      $lengthSum=0;
      foreach my $idx (@BseqIdx)
      {
        if(exists($mismatchCount->{$idx}))
        {
          $mismatchSum+=$mismatchCount->{$idx};
          #print "$idx has $mismatchCount->{$idx} with length $seqLength->{$idx}\n";
        }
        $lengthSum+=$seqLength->{$idx};
      }
      if($lengthSum==0)
      {
        $pileupData->[$i]->{'BmeanPMM'}=0;
      }
      elsif($pileupData->[$i]->{'RefString'} eq $pileupData->[$i]->{'Bstring'})
      {
        #print "B variant is Ref\n";
        #if($lengthSum>0)
        #{
          $pileupData->[$i]->{'BmeanPMM'}=$mismatchSum/$lengthSum;
        #}
        #else
        #{
        #  print "pileup: " . Dumper($pileupData->[$i]);
        #  print "AseqIdx: " . Dumper(@AseqIdx);
        #  print "seqlength: " . Dumper($seqLength);
        #}
      }
      elsif($pileupData->[$i]->{'Bstring'}=~/^([\-\+])(\d+)/)
      {
        #print "B variant is indel, old AcountF" . $pileupData->[$i]->{'AcountF'} . "\n";
        #print "B 1: $1 2: $2\n";
        my $sign=$1;
        my $indelLen=$2;
        #print "sign: $sign indelLen: $indelLen\n";
        $pileupData->[$i]->{'BmeanPMM'}=($mismatchSum-$2*($pileupData->[$i]->{'BcountF'}+$pileupData->[$i]->{'BcountR'}))/$lengthSum;
        if($sign=~/\-/ && $pileupData->[$i]->{'Astring'} eq $pileupData->[$i]->{'RefString'})
        {
          #my $readCountSum=0;
          #my $readCountPassSum=0;
          #my $refSumFpass=0;
          #my $refSumRpass=0;
          #my $delCountSum=0;
          #my $delMQSum=0;
          #my $delBQSum=0;
          my $misMatchSumF=0;
          my $misMatchSumR=0;
          for(my $j=$i+1;$j<=$i+$indelLen;$j++)
          {
            if(exists($pileupData->[$j]))
            {
              $misMatchSumF+=$pileupData->[$j]->{'MismatchCountF'};
              $misMatchSumR+=$pileupData->[$j]->{'MismatchCountR'};
             #print Dumper($pileupData->[$j]);
              #$refSumFpass+=$pileupData->[$j]->{'RefCountFpass'};
              #$refSumRpass+=$pileupData->[$j]->{'RefCountRpass'};
              #$readCountSum+=$pileupData->[$j]->{'readCount'};
              #$readCountPassSum+=$pileupData->[$j]->{'readCountPass'};
              #$delCountSum+=$pileupData->[$j]->{'delCount'};
              #$delMQSum+=$pileupData->[$j]->{'delMQ'};
              #$delBQSum+=$pileupData->[$j]->{'delBQ'};
            }
          }
          $pileupData->[$i]->{'AcountF'}=max($pileupData->[$i]->{'AcountF'}-$misMatchSumF/$indelLen,0);
          $pileupData->[$i]->{'AcountR'}=max($pileupData->[$i]->{'AcountR'}-$misMatchSumR/$indelLen,0);
          #if($pileupData->[$i]->{'BcountF'}+$pileupData->[$i]->{'BcountR'}!=$delCountSum/$indelLen)
          #{
          #  my $delCountDiff=$pileupData->[$i]->{'BcountF'}+$pileupData->[$i]->{'BcountR'}-$delCountSum/$indelLen;
          #  $pileupData->[$i]->{'BcountF'}-=$delCountDiff/2;
          #  $pileupData->[$i]->{'BcountR'}-=$delCountDiff/2;
          #}
          #$pileupData->[$i]->{'readCount'}=$readCountSum/$indelLen;
          #$pileupData->[$i]->{'readCountPass'}=$readCountPassSum/$indelLen;
          #$pileupData->[$i]->{'BmeanMQ'}=$delMQSum/$indelLen;
          #$pileupData->[$i]->{'BmeanBQ'}=$delBQSum/$indelLen;
          #print "new ACountF" . $pileupData->[$i]->{'AcountF'} . "\n";
        }
      }
      else
      {
        #print "B variant is SNV\n";
        $pileupData->[$i]->{'BmeanPMM'}=($mismatchSum-$pileupData->[$i]->{'BcountF'}-$pileupData->[$i]->{'BcountR'})/$lengthSum;
      }
    }
  }

  #print Dumper($pileupData);
  #print Dumper($param);
  return $pileupData;
}

sub parseMulti()
{
  my ($processName,$pileup,$filt,$param,$istumor)=@_;
  my @pileupLines=split(/\n/,$pileup);
  my $splitPileup=[];
  my $pileupData=[];
  my $line=shift(@pileupLines);
  my @array=split(/\t/,$line);
  #print Dumper(@array);
  #print "parsing pileupline with $#array elements\n";
  for(my $i=0;$i<($#array-2)/5;$i++)
  {
    my $startPos=$i*5+3;
    my $endPos=$i*5+7;
    #print "joining element $i from $startPos  to $endPos:" .join("\t",@array[$startPos..$endPos]) . "\n";
    $splitPileup->[$i]=join("\t",@array[0..2]) . "\t" . join("\t",@array[$startPos..$endPos]) . "\n";
  }
  foreach $line (@pileupLines)
  {
    #my ($chr,$pos,$ref,$readCount,$bases,$baseQual,$mapQual,$readPos)=split(/\t/,$line);
    my @array=split(/\t/,$line);
    #print Dumper(@array);
    #print "parsing pileupline with $#array elements\n";
    for(my $i=0;$i<($#array-2)/5;$i++)
    {
      my $startPos=$i*5+3;
      my $endPos=$i*5+7;
      #print "joining element $i from $startPos  to $endPos:" .join("\t",@array[$startPos..$endPos]) . "\n";
      $splitPileup->[$i] .= join("\t",@array[0..2]) . "\t" . join("\t",@array[$startPos..$endPos]) . "\n";
    }
  }
  my $i=0;
  foreach my $currPileup (@$splitPileup)
  {
    #print "pileup for $i\n";
    $pileupData->[$i]=&parse('NA',$currPileup,$filt,$param,$istumor);
    $i++;
  }
  return $pileupData;
}

return 1;
