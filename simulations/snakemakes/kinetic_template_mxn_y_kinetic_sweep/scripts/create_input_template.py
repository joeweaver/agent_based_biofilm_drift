import sys
import os
import fileinput
import re

run_dir = sys.argv[1]
loser_file = os.path.join(run_dir, 'biggest_loser.txt')
original_input_file = os.path.join(run_dir, 'Inputscript.orig')
templated_input_file = os.path.join(run_dir, 'Inputscript.templated')

with open(loser_file, 'r') as f:
    biggest_loser = f.readline().strip()


m_str = f'fix monod_het_{biggest_loser} het_{biggest_loser} nufeb/monod/het sub \(.+\) o2'
with open(original_input_file, 'r') as oif:
    with open(templated_input_file, 'w') as tif:
        for line in oif:
            if re.match(f'fix monod_het_{biggest_loser} ', line):
                res = re.sub('sub .+ o2', 'sub het_ks o2', line)
                res = re.sub('growth .+ yield .+ decay', 
                            'growth het_mu yield het_yield decay', res)
                tif.write(res)
            else:
                tif.write(line)
