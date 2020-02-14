import sys
import subprocess
import pandas as pd

input_file = sys.argv[1]
output_dir = sys.argv[2]

df = pd.read_csv(input_file)

drugs =  ['R0-3306',
          'BMS-265246',
          'LEE011_Ribociclib',
          'Abemaciclib_LY2835219',
          'BSJ-03-124',
          'BSJ-03-123', 
          'SY-1365',
          'YKL-5-124',
          'THZ531',
          'THZ1',
          'senexinb',
          'ZZ133B',
          'BSJ-01-175',
          'MFH-2-90',
          'E17',
          'FMF041072',
          'Flavopiridol',
          'LY2606368',
          'LY3023414',
          'AZD2014',
          'AZD5363',
          'AZD1775',
          'CDK1_2iIII',
          'Torin2']

for drug in drugs:
    cmd = ("sbatch -o logs/%J.log -p long -t 15-00:00 --mem=10G").split(' ')
    cmd +=  ["--job-name=%s" % (drug)]
    cmd += ["--wrap=python run_fit_rep.py %s %s %s" % (drug, input_file, output_dir)]
    subprocess.call(cmd)
