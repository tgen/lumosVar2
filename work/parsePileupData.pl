#!/usr/bin/perl
#use warnings;
use strict;
use Data::Dumper;
use File::stat;
use Math::BaseCalc;
use lib '/Users/rhalperin/perl5/lib/perl5/';
use YAML::Tiny;
use List::Util qw[min max];
use Cwd;
my $dir = getcwd;
use TGen::parseSNPvcf;
use TGen::parseMpileup;

my $configFile=$ARGV[0];
my $yaml = YAML::Tiny->read($configFile);
chomp($yaml->[0]);
my $param=$yaml->[0];

my $block=$ARGV[1];
my $tumorCount=$ARGV[2];

if ($tumorCount>1)
{
  my $calc = new Math::BaseCalc(digits => ['N','A','C','G','T']);
  my $mPileup  = `cd $param->{'samPath'}\n ./samtools mpileup -B -R -r $block -l $param->{'regionsFile'} --rf 0x2 --ff 0x400 -O -s -Q 0 -C $param->{'mpileupC'} -f $param->{'refGenome'} -b $param->{'bamList'}`;
  #print Dumper($mPileup);
  my $output=parseMpileup->parseMulti($mPileup,0,$param,1);
  #print Dumper($output);
  my ($chr,$range)=split(/\:/,$block);
  my $snpVCF=$param->{'snpVCFpath'} . $chr . $param->{'snpVCFname'};
  my $vcfText= `cd $param->{'tabixPath'} \n./tabix $snpVCF $block`;
  my $snpData = parseSNPvcf->getAFbyNT($vcfText,$param->{'pvFreq'});
  my $cosmicVCF=$param->{'cosmicVCF'};
  my $cosmicText= `cd $param->{'tabixPath'} \n./tabix $cosmicVCF $block`;
  my $cosmicData = parseSNPvcf->getCNT($cosmicText);
  my $i=0;
  my $dataMat='';
  my $rdMat='';
  my %posList;
  foreach my $sample (@$output)
  {
    #print "minBcount " . $param->{'minBCount'};
    foreach my $data (@$sample)
    {
      if(exists($data->{'Ref'}) && $data->{'Ref'} ne $data->{'A'} && $data->{'AcountF'}+$data->{'AcountR'}>=$param->{'minBCount'})
      {
        $posList{$data->{'Pos'}}=1;
        #print "found non ref pos " . $data->{'Pos'};
      }
      elsif($data->{'BcountF'}+$data->{'BcountR'} >= $param->{'minBCount'})
      {
        $posList{$data->{'Pos'}}=1;
        #print "found var pos " . $data->{'Pos'};
      }
      else
      {
        #print "skipping " . $data->{'Pos'} . "\n";
        #print "bcount sum " . $data->{'BcountF'}+$data->{'BcountR'} . "\n";
      }
    }
  }
  #print Dumper(%posList);
  foreach my $sample (@$output)
  {
    $i++;
    foreach my $data (@$sample)
    {
      #print Dumper($data);
      #$dataMat .= "\n" . "$i";
      $rdMat .= "\n\@$i\t" . $data->{'Chr'};
      $rdMat .= "\t" . $data->{'Pos'};
      $rdMat .= "\t" . $data->{'readCountPass'};
      #print Dumper($data->{'indels'});
      unless(exists($posList{$data->{'Pos'}}))
      {
        next;
      }
      $dataMat .= "\n$i\t" . $data->{'Chr'};
      $dataMat .= "\t" . $data->{'Pos'};
      $dataMat .= "\t" . $data->{'readCount'};
      $dataMat .= "\t" . $data->{'readCountPass'} . "\t";
      $dataMat .= $calc->from_base($data->{'Ref'});
      #for my $c (split //, $data->{'Ref'})
      #{
      #$dataMat .= $ntInt{$c};
      #}
      $dataMat .= "\t";
      $dataMat .= $calc->from_base($data->{'A'});
      #for my $c (split //, $data->{'A'})
      #{
      #  $dataMat .= $ntInt{$c};
      #}
      $dataMat .= "\t" . $data->{'AcountF'};
      $dataMat .= "\t" . $data->{'AcountR'};
      $dataMat .= "\t" . $data->{'AmeanBQ'};
      $dataMat .= "\t" . $data->{'AmeanMQ'};
      $dataMat .= "\t" . $data->{'AmeanPMM'};
      $dataMat .= "\t" . $data->{'AmeanReadPos'} . "\t";
      $dataMat .= $calc->from_base($data->{'B'});
      #for my $c (split //, $data->{'B'})
      #{
      #  $dataMat .= $ntInt{$c};
      #}
      $dataMat .= "\t" . $data->{'BcountF'};
      $dataMat .= "\t" . $data->{'BcountR'};
      $dataMat .= "\t" . $data->{'BmeanBQ'};
      $dataMat .= "\t" . $data->{'BmeanMQ'};
      $dataMat .= "\t" . $data->{'BmeanPMM'};
      $dataMat .= "\t" . $data->{'BmeanReadPos'};
      my $chrPos = $data->{'Chr'} . ':' . $data->{'Pos'};
      if(exists($snpData->{$chrPos}->{$data->{'A'}}))
      {
        $dataMat .= "\t" . $snpData->{$chrPos}->{$data->{'A'}};
      }
      elsif($data->{'Ref'} eq $data->{'A'})
      {
        my $AF = 1-3*$param->{'pvFreq'}-$param->{'pvFreqIndel'};
        $dataMat .= "\t$AF";
      }
      elsif(max(length($data->{'A'}),length($data->{'B'}))==1)
      {
        $dataMat .= "\t$param->{'pvFreq'}";
      }
      else
      {
        $dataMat .= "\t$param->{'pvFreqIndel'}";
      }
      if(exists($snpData->{$chrPos}->{$data->{'B'}}))
      {
        $dataMat .= "\t" . $snpData->{$chrPos}->{$data->{'B'}};
      }
      elsif($data->{'Ref'} eq $data->{'B'})
      {
        my $AF = 1-3*$param->{'pvFreq'}-$param->{'pvFreqIndel'};
        $dataMat .= "\t$AF";
      }
      elsif(max(length($data->{'A'}),length($data->{'B'}))==1)
      {
        $dataMat .= "\t$param->{'pvFreq'}";
      }
      else
      {
        $dataMat .= "\t$param->{'pvFreqIndel'}";
      }
      if(exists($cosmicData->{$chrPos}->{'CNT'}))
      {
        $dataMat .= "\t$cosmicData->{$chrPos}->{'CNT'}";
      }
      else
      {
        $dataMat .= "\t0";
      }
    }
  }
  #}
  print $dataMat;
  print $rdMat;
}

elsif ($tumorCount==1)
{
  my $calc = new Math::BaseCalc(digits => ['N','A','C','G','T']);
  my $mPileup  = `cd $param->{'samPath'}\n ./samtools mpileup -B -R -r $block -l $param->{'regionsFile'} --rf 0x2 --ff 0x400 -O -s -Q 0 -C $param->{'mpileupC'} -f $param->{'refGenome'} $param->{'bamFile'}`;
  my $output=parseMpileup->parse($mPileup,1,$param);
  my ($chr,$range)=split(/\:/,$block);
  my $snpVCF=$param->{'snpVCFpath'} . $chr . $param->{'snpVCFname'};
  my $vcfText= `cd $param->{'tabixPath'} \n./tabix $snpVCF $block`;
  my $snpData = parseSNPvcf->getAFbyNT($vcfText,$param->{'pvFreq'});
  my $cosmicVCF=$param->{'cosmicVCF'};
  my $cosmicText= `cd $param->{'tabixPath'} \n./tabix $cosmicVCF $block`;
  my $cosmicData = parseSNPvcf->getCNT($cosmicText);
  my $i=0;
  my $dataMat='';
  my $rdMat='';

  foreach my $data (@$output)
  {
    #print Dumper($data);
    #$dataMat .= "\n" . "$i";
    $rdMat .= "\n\@" . $data->{'Chr'};
    $rdMat .= "\t" . $data->{'Pos'};
    $rdMat .= "\t" . $data->{'readCountPass'};
    #print Dumper($data->{'indels'});
    unless(exists($data->{'Ref'}))
    {
      next;
    }
    $dataMat .= "\n" . $data->{'Chr'};
    $dataMat .= "\t" . $data->{'Pos'};
    $dataMat .= "\t" . $data->{'readCount'};
    $dataMat .= "\t" . $data->{'readCountPass'} . "\t";
    $dataMat .= $calc->from_base($data->{'Ref'});
    #for my $c (split //, $data->{'Ref'})
    #{
      #$dataMat .= $ntInt{$c};
    #}
    $dataMat .= "\t";
    $dataMat .= $calc->from_base($data->{'A'});
    #for my $c (split //, $data->{'A'})
    #{
    #  $dataMat .= $ntInt{$c};
    #}
    $dataMat .= "\t" . $data->{'AcountF'};
    $dataMat .= "\t" . $data->{'AcountR'};
    $dataMat .= "\t" . $data->{'AmeanBQ'};
    $dataMat .= "\t" . $data->{'AmeanMQ'};
    $dataMat .= "\t" . $data->{'AmeanPMM'};
    $dataMat .= "\t" . $data->{'AmeanReadPos'} . "\t";
    $dataMat .= $calc->from_base($data->{'B'});
    #for my $c (split //, $data->{'B'})
    #{
    #  $dataMat .= $ntInt{$c};
    #}
    $dataMat .= "\t" . $data->{'BcountF'};
    $dataMat .= "\t" . $data->{'BcountR'};
    $dataMat .= "\t" . $data->{'BmeanBQ'};
    $dataMat .= "\t" . $data->{'BmeanMQ'};
    $dataMat .= "\t" . $data->{'BmeanPMM'};
    $dataMat .= "\t" . $data->{'BmeanReadPos'};
    my $chrPos = $data->{'Chr'} . ':' . $data->{'Pos'};
    if(exists($snpData->{$chrPos}->{$data->{'A'}}))
    {
      $dataMat .= "\t" . $snpData->{$chrPos}->{$data->{'A'}};
    }
    elsif($data->{'Ref'} eq $data->{'A'})
    {
      my $AF = 1-3*$param->{'pvFreq'}-$param->{'pvFreqIndel'};
      $dataMat .= "\t$AF";
    }
    elsif(max(length($data->{'A'}),length($data->{'B'}))==1)
    {
      $dataMat .= "\t$param->{'pvFreq'}";
    }
    else
    {
      $dataMat .= "\t$param->{'pvFreqIndel'}";
    }
    if(exists($snpData->{$chrPos}->{$data->{'B'}}))
    {
      $dataMat .= "\t" . $snpData->{$chrPos}->{$data->{'B'}};
    }
    elsif($data->{'Ref'} eq $data->{'B'})
    {
      my $AF = 1-3*$param->{'pvFreq'}-$param->{'pvFreqIndel'};
      $dataMat .= "\t$AF";
    }
    elsif(max(length($data->{'A'}),length($data->{'B'}))==1)
    {
      $dataMat .= "\t$param->{'pvFreq'}";
    }
    else
    {
      $dataMat .= "\t$param->{'pvFreqIndel'}";
    }
    if(exists($cosmicData->{$chrPos}->{'CNT'}))
    {
      $dataMat .= "\t$cosmicData->{$chrPos}->{'CNT'}";
    }
    else
    {
      $dataMat .= "\t0";
    }
  }
  #}
  print $dataMat;
  print $rdMat;
}
else
{
  my %ntInt=('A' => 1, 'C' => 2, 'G' =>3, 'T' =>4, 'R' => 5, 'Y' => 6, 'K' => 7, 'M' => 8, 'S' => 9,
  'W' => 10, 'B' => 11, 'D' => 12, 'H' => 13, 'V' => 14, 'N' =>15, '-' => 16, '*' => 17);
  my $mPileup  = `cd $param->{'samPath'}\n ./samtools mpileup -r $block -l $param->{'regionsFile'} --rf 0x2 --ff 0x400 -O -s -Q 0 -C $param->{'mpileupC'} -f $param->{'refGenome'} -b $param->{'bamList'}`;
  my $output=parseMpileup->parseMulti($mPileup,0,$param);
  my $pvFreq=$param->{'pvFreq'};
  my ($chr,$range)=split(/\:/,$block);
  my $snpVCF=$param->{'snpVCFpath'} . $chr . $param->{'snpVCFname'};
  my $vcfText= `cd $param->{'tabixPath'}\n ./tabix $snpVCF $block`;
  my $snpData = parseSNPvcf->getAFbyNT($vcfText,$pvFreq);
  my $i=0;
  my $dataMat='';
  foreach my $sample (@$output)
  {
    $i++;
    foreach my $data (@$sample)
    {
      #print Dumper $data;
      $dataMat .= "\n" . "$i";
      $dataMat .= "\t" . $data->{'Chr'};
      $dataMat .= "\t" . $data->{'Pos'};
      $dataMat .= "\t" . $data->{'readCount'};
      $dataMat .= "\t" . $data->{'readCountPass'} . "\t";
      for my $c (split //, $data->{'Ref'})
      {
        $dataMat .= $ntInt{$c};
      }
      $dataMat .= "\t";
      for my $c (split //, $data->{'A'})
      {
        $dataMat .= $ntInt{$c};
      }
      $dataMat .= "\t" . $data->{'AcountF'};
      $dataMat .= "\t" . $data->{'AcountR'};
      $dataMat .= "\t" . $data->{'AmeanBQ'};
      $dataMat .= "\t" . $data->{'AmeanMQ'} . "\t";
      for my $c (split //, $data->{'B'})
      {
        $dataMat .= $ntInt{$c};
      }
      $dataMat .= "\t" . max($data->{'BcountF'},0);
      $dataMat .= "\t" . max($data->{'BcountR'},0);
      $dataMat .= "\t" . $data->{'BmeanBQ'};
      $dataMat .= "\t" . $data->{'BmeanMQ'};
      my $chrPos = $data->{'Chr'} . ':' . $data->{'Pos'};
      if(exists($snpData->{$chrPos}->{$data->{'A'}}))
      {
        $dataMat .= "\t" . $snpData->{$chrPos}->{$data->{'A'}};
      }
      elsif($data->{'Ref'}=~$data->{'A'})
      {
        my $AF = 1-3*$pvFreq;
        $dataMat .= "\t$AF";
      }
      else
      {
        $dataMat .= "\t$pvFreq";
      }
      if(exists($snpData->{$chrPos}->{$data->{'B'}}))
      {
        $dataMat .= "\t" . $snpData->{$chrPos}->{$data->{'B'}};
      }
      elsif($data->{'Ref'}=~$data->{'B'})
      {
        my $AF = 1-3*$pvFreq;
        $dataMat .= "\t$AF";
      }
      else
      {
        $dataMat .= "\t$pvFreq";
      }
    }
  }
  print $dataMat;
}
