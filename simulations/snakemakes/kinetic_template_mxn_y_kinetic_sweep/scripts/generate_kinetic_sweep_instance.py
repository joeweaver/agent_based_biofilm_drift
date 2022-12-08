import os
import sys
import shutil
import fileinput
import re

print('generating run')


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        if item != 'cell_rel_volumes.csv':
            s = os.path.join(src, item)
            d = os.path.join(dst, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, symlinks, ignore=shutil.ignore_patterns('cell_rel_volumes.csv'))
            else:
                shutil.copy2(s, d)


source_dir = sys.argv[1]
template_dir = os.path.dirname(sys.argv[1])
base_grid_dir = os.path.basename(template_dir)
seed_dir = os.path.basename(sys.argv[1])
ks = sys.argv[2]
mu = sys.argv[3]
bioyield = sys.argv[4]

target_dir = os.path.join('results', base_grid_dir, f'ks_{ks}-mu_{mu}-yield_{bioyield}', seed_dir)

if os.path.isdir(target_dir):
    copytree(source_dir, target_dir)
else:
    shutil.copytree(source_dir, target_dir,ignore=shutil.ignore_patterns('cell_rel_volumes.csv'))

ifname = os.path.join('.', target_dir, 'Inputscript.templated')
ofname = os.path.join('.', target_dir, 'Inputscript.lmp')
with open(ifname, 'r') as infile:
    with open(ofname, 'w') as outfile:
        for line in infile:
            if re.search('het_ks', line):
                print(line)
                res = re.sub('het_ks', f'{ks}', line)
                print(res)
                res = re.sub('het_mu', f'{mu}', res)
                print(res)
                res = re.sub('het_yield', f'{bioyield}', res)
                print(res)
                outfile.write(res)
            else:
                outfile.write(line)
