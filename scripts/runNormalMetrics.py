import sys, yaml
import subprocess

def die(msg):
  sys.stderr.write("%s: " % sys.argv[0])
  sys.stderr.write(msg)
  sys.stderr.write("\n")
  exit(1)

def chr_lists(autosomes, sexs, prefix):
  a_start, a_end = map(int, autosomes.split(":"))

  sex_list = list(map(lambda x: prefix + x, sexs.split(",")))

  autosomes_list = list(map(lambda x: prefix + str(x), range(a_start, a_end + 1)))

  return autosomes_list, sex_list

if __name__ == "__main__":
  if len(sys.argv) < 2:
    die("missing positional argument: config")

  if len(sys.argv) < 3:
    die("missing positional argument: chromosome")

  cfile=sys.argv[1]
  with open(cfile) as f:
    config = yaml.load(f)

  try:
    chro = config['contigPrefix'] + sys.argv[2]
    autosomes, sexs = chr_lists(config["autosomes"], config["sexChr"], config['contigPrefix'])
  except:
    chro = sys.argv[2]
    autosomes, sexs = chr_lists(config["autosomes"], config["sexChr"],'')
  chrs = autosomes + sexs

  sexlist = config["sexList"].split(",")
  nsample = len(sexlist)

  mlist = list(map(int, config["M"].split(",")))
  flist = list(map(int, config["F"].split(",")))

  try:
    index = sexs.index(chro)
    ploidy = [None] * nsample
    for i, g in enumerate(sexlist):
      if g == "F":
        ploidy[i] = str(flist[index])
      else:
        ploidy[i] = str(mlist[index])
    ploidyStr=''.join(ploidy)
  except ValueError:
    ploidyStr=('2' * nsample)

  command=config["gvmPath"] + "/gvm -c " + cfile + " -C " + chro + " -P -E -N --ploidystr=" + ploidyStr
  print("Running:\n" + command)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  stdout, stderr = process.communicate()
  print(stdout)
  if not process.returncode == 0:
    die("gvm failed:\n" + stderr)

  command1="sort -n -k2 " + config["outfile"] + chro + ".txt"
  command2="bgzip -c >" + config["outfile"] + chro + ".txt.gz"
  print("Running:\n" + command1 + " | " + command2)
  p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, shell=True)
  p2 = subprocess.Popen(command2, stdin=p1.stdout, shell=True)
  #stdout = p1.stdout() + "\n" + p2.stdout()
  #print(stdout)
  p1.wait()
  print(p1.returncode)
  if not p1.returncode == 0:
    print(p1.stderr)
    die("sort failed")
  p2.wait()
  if not p2.returncode == 0:
    print(p2.stderr)
    die("bgzip failed")

  command="tabix -b 2 -e 2 " + config["outfile"] + chro +  ".txt.gz"
  print("Running:\n" + command)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  stdout, stderr = process.communicate()
  print(stdout)
  if not process.returncode == 0:
    print(stderr)
    die("tabix failed")
