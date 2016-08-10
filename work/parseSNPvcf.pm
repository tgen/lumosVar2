#!/usr/bin/perl
package parseSNPvcf;
use warnings;
use strict;
use Data::Dumper;
use File::stat;
use List::Util qw[min max];

sub getAF
{
  my ($processName,$vcfText)=@_;
  my @vcfLines=split(/\n/,$vcfText);
  my $vcfData={};
  my @NT=('A', 'C', 'G', 'T');
  foreach my $line (@vcfLines)
  {
    chomp $line;
    my @fields=split(/\t/,$line);
    my $chrPos=$fields[0] . ':' . $fields[1];
    $vcfData->{$chrPos}->{'ref'}=$fields[3];
    $vcfData->{$chrPos}->{'alt'}=$fields[4];
    $vcfData->{$chrPos}->{'qual'}=$fields[5];
    $vcfData->{$chrPos}->{'filter'}=$fields[6];
    #my @infoFields=split(/\;/,$info);
    $fields[7] =~m/AF\=(0\.[0-9]*)/;
    $vcfData->{$chrPos}->{'AF'}=$1;
  }
  return $vcfData;
}

sub getAFbyNT
{
  my ($processName,$vcfText,$pvFreq)=@_;
  my @vcfLines=split(/\n/,$vcfText);
  my $vcfData={};
  foreach my $line (@vcfLines)
  {
    chomp $line;
    my @fields=split(/\t/,$line);
    my $chrPos=$fields[0] . ':' . $fields[1];
    $vcfData->{$chrPos}->{'qual'}=$fields[5];
    $vcfData->{$chrPos}->{'filter'}=$fields[6];
    #my @infoFields=split(/\;/,$info);
    my $AF=$pvFreq;
    if ($fields[7] =~m/AF\=(0\.[0-9]*)/)
    {
      $AF=$1;
    }
    #$vcfData->{$chrPos}->{$fields[3]}=1-$AF-2*$pvFreq;
    $vcfData->{$chrPos}->{$fields[4]}=$AF;
    if(exists($vcfData->{$chrPos}->{$fields[3]}))
    {
      $vcfData->{$chrPos}->{$fields[3]}-=$AF;
    }
    else
    {
      $vcfData->{$chrPos}->{$fields[3]}=1-$AF-2*$pvFreq;
    }
  }
  return $vcfData;
}

sub getCNT
{
  my ($processName,$vcfText)=@_;
  my @vcfLines=split(/\n/,$vcfText);
  my $vcfData={};
  my @NT=('A', 'C', 'G', 'T');
  foreach my $line (@vcfLines)
  {
    chomp $line;
    my @fields=split(/\t/,$line);
    my $chrPos=$fields[0] . ':' . $fields[1];
    $vcfData->{$chrPos}->{'ref'}=$fields[3];
    $vcfData->{$chrPos}->{'alt'}=$fields[4];
    $vcfData->{$chrPos}->{'qual'}=$fields[5];
    $vcfData->{$chrPos}->{'filter'}=$fields[6];
    #my @infoFields=split(/\;/,$info);
    $fields[7] =~m/CNT\=([0-9]*)/;
    if(exists($vcfData->{$chrPos}->{'CNT'}))
    {
      $vcfData->{$chrPos}->{'CNT'}=max($1,$vcfData->{$chrPos}->{'CNT'});
    }
    else
    {
      $vcfData->{$chrPos}->{'CNT'}=$1;
    }
  }
  return $vcfData;
}


return 1;
