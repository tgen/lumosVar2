import sys, yaml
import os, subprocess

def die(msg):
  sys.stderr.write("%s: " % sys.argv[0])
  sys.stderr.write(msg)
  sys.stderr.write("\n")
  exit(1)

def chr_lists(autosomes, sexs):
  a_start, a_end = map(int, autosomes.split(":"))

  sex_list = sexs.split(",")

  autosomes_list = list(range(a_start, a_end + 1))

  return list(map(str, autosomes_list)), sex_list

def pipe_commands(command_list):
  print("process " + str(0) + " is: " +  command_list[0])
  print(command_list[0].split())
  plist = [subprocess.Popen(command_list[0].split(),stdout=subprocess.PIPE)]
  #print("first process is " + str(plist[0].pid))
  for i in range(1,len(command_list)):
    print("process " + str(i) + " is: " +  command_list[i])
    print(command_list[i].split())
    #print("input is :")
    #print(plist[i-1].stdout)
    plist.append(subprocess.Popen(command_list[i].split(), stdin=plist[i-1].stdout,stdout=subprocess.PIPE))
    #print("prev process is " + str(plist[i-1].pid))
    #print("curr process is " + str(plist[i].pid))
    #plist[i-1].stdout.close()
  #my_env=os.environ.copy()
  #print(my_env)
  #if i==0:
  #  p = subprocess.Popen(command_list[i].split(),subprocess.PIPE)
  #else:
  #  p =  subprocess.Popen(command_list[i].split(), stdin=plist[i-1].stdout.read())
  #plist.append(p)
  #if i<len(plist):
  #  pipe_commands(command_list,plist,i+1)
  #print("last process is " + str(plist[-1].pid))
  #(out,err)=plist[-1].communicate()
  #print(out)
  #print(err)
  return plist

def check_plist(plist,name):
  for i in range(0,len(plist)):
    p=plist[i]
    p.wait()
    print(name + "process " + str(i) + " finished")
    if not p.returncode == 0:
      die(name + "process " + str(i) + " failed:\n")
  return 0

if __name__ == "__main__":
  if len(sys.argv) < 2:
    die("missing positional argument: config")

  cfile=sys.argv[1]
  with open(cfile) as f:
    config = yaml.load(f)

  autosomes, sexs = chr_lists(config["autosomes"], config["sexChr"])
  chrs = autosomes + sexs

  mlist = list(map(int, config["M"].split(",")))
  flist = list(map(int, config["F"].split(",")))
  
  #print(flist)
  try:
    index = flist.index(2)
    diploid = "F"
  except ValueError:
    try:
      index = mlist.index(2)
      diploid = "M"
    except ValueError:
      die("diploid sex chr not found")
  
  testChrs=[autosomes[-1],sexs[index]]
  plist=[]
  for chrom in testChrs:
    print("finding common var pos for " + chrom)
    command_list=["bcftools view -q 0.2:nonmajor -v snps -R " + config["regionsFile"] + " " + config["snpVCFpath"] + chrom + config["snpVCFname"]]
    command_list.append("grep -v ^#")
    p=pipe_commands(command_list)
    plist.append(p)
    bedfile = open("commonVar.chr" + chrom + ".bed","w")
    while True:
      line = p[-1].stdout.readline().decode("utf8")
      if line != '':
        chrPos=line.split('\t',2)
        bedfile.write(chrPos[0] + "\t" + str(int(chrPos[1])-1) + "\t" + chrPos[1] + "\n")
      else:
        break
    bedfile.close()
  
  #print(plist) 
  for i in range(0,len(plist)):
    check_plist(plist[i],testChrs[i])
  
  plist=[]
  stats=[]
  for chrom in testChrs:
    command_list=["bcftools mpileup -Ou -b " + config["bamList"] +  " -B -f " + config["refGenome"] + " -R commonVar.chr" + chrom  + ".bed -I -Ou"]
    command_list.append("bcftools call -Ou -m")
    command_list.append("bcftools stats -s -")
    command_list.append("grep ^PSC")
    command_list.append("cut -f3,4,5,6,14")
    p=pipe_commands(command_list)
    plist.append(p)
    stats.append(p[-1].stdout)
  
  for i in range(0,len(plist)):
    check_plist(plist[i],testChrs[i])
    #print(stats[i].readline().decode("utf8"))
  out=open(config["bamList"] + ".guessSex.txt","w")
  header=["sample",testChrs[0] + "-nHomRef",testChrs[0] + "-nHomAlt",testChrs[0] + "-nHet",testChrs[0]+"-none"]
  header.extend([testChrs[1] + "-nHomRef",testChrs[1] + "-nHomAlt",testChrs[1] + "-nHet",testChrs[1]+"-none"]) 
  header.extend([testChrs[0] + "-Het/Hom",testChrs[1] + "-Het/Hom","sex"])
  out.write("\t".join(header) + "\n")
  sex_list=[]
  while True:
    lineA = stats[0].readline().decode("utf8").rstrip()
    lineS = stats[1].readline().decode("utf8").rstrip()
    print("autosome stats:" + lineA)
    print("sexChr stats:" + lineS)
    if lineA == '':
      break
    dataA = lineA.split("\t")
    dataS = lineS.split("\t")
    print("numHomA: " + str(int(dataA[1]+int(dataA[2])))
    print("numHomS: " + str(int(dataS[1]+int(dataS[2])))
    print("numHetA: " + str(int(dataA[3])))
    print("numHetS: " + str(int(dataS[3])))
    hetA=int(dataA[3])/(int(dataA[1])+int(dataA[2]))
    hetS=int(dataS[3])/(int(dataS[1])+int(dataS[2]))
    if hetS/hetA > 0.5:
      sex=diploid
    elif (hetS/hetA < 0.2) & (diploid=="F"):
      sex="M"
    elif (hetS/hetA < 0.2) & (diploid=="M"):
      sex="F"
    else:
      sex="U"
    out.write("\t".join(dataA) + "\t" + "\t".join(dataS[1:len(dataA)]) + "\t" + str(hetA) + "\t" + str(hetS) + "\t" + sex + "\n")
    sex_list.append(sex)
  out.write("\nsexList: " + ','.join(sex_list) + "\n")
  out.close()  

#module load samtools/1.7
#PATH=$PATH:/home/rhalperin/bin/samtools-1.5/bin/

#bcftools view -q 0.2:nonmajor -v snps -R $BED ${VCFPATH}${DCHR}${VCFEXT} | grep -v ^# | awk '{ print $1 "\t" ($2-1) "\t" $2}' >commonVar.chr${DCHR}.bed
#bcftools mpileup -Ou -b $BAMLIST -B -f $REF -R commonVar.chr${DCHR}.bed -I -Ou | bcftools call -Ou -m | bcftools stats -i '%QUAL>20' -s - |  grep -B 1 ^PSC | cut -f3,4,5,6,14 >${BAMLIST}.callCounts.chr${DCHR}.txt &

#bcftools view -q 0.2:nonmajor -v snps -R $BED ${VCFPATH}${TCHR}${VCFEXT} | grep -v ^# | awk '{ print $1 "\t" ($2-1) "\t" $2}' >commonVar.chr${TCHR}.bed
#bcftools mpileup -Ou -b $BAMLIST -B -f $REF -R commonVar.chr${TCHR}.bed -I -Ou | bcftools call -Ou -m | bcftools stats -i '%QUAL>20' -s - |  grep -B 1 ^PSC | cut -f3,4,5,6,14 >${BAMLIST}.callCounts.chr${TCHR}.txt &



