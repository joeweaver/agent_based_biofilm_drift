import sys
import os
import re
import numpy as np
import pandas as pd

indir = sys.argv[1]
grid_N = int(sys.argv[2])
deltaT = float(sys.argv[3])
basedir = os.path.basename(indir)
#print(f'cwd {os.getcwd()}')
#print(f'indir: {indir}')
#print(f'basedir {basedir}')
#results/3x3_10_default_mu_ks_yield_conc/ks_1.75e-05-mu_1.40e-04-yield_7.00e-01/rand1701/
rematch = re.search('conc\/ks_(.*)-mu_(.*)-yield_(.*)\/rand(.*)\/', indir)
#rematch = re.search('rand(.*)-ks_(.*)-mu_(.*)-yield_(.*)', indir)
#print(rematch)
#if(rematch):
#    print(f'ks: {rematch.group(1)} mu: {rematch.group(2)} yield: {rematch.group(3)} seed: {rematch.group(4)}') 
#else:
#    print('did not find parameters from pathname')
#    exit(1)

seed = rematch.group(4)
ks = rematch.group(1)
mu = rematch.group(2)
bioyield = rematch.group(3)

def site_fate(site_name, cell_volumes,n_bugs,timestep):
    #print(site_name)
    #timestep = 1000  # seconds
    char_diam = 1e-6
    char_vol = (char_diam/2)*(char_diam/2)*(char_diam/2)*4/3*np.pi
    rvm = cell_volumes[f'{site_name}v'].rolling(25).median()
    rrv = cell_volumes[f'{site_name}rv'].rolling(25).median()
    step_tail = cell_volumes["step"][-500:]
    het5_tail = rvm[-10:]
    y1 = het5_tail.iloc[0]
    y2 = het5_tail.iloc[9]
    x1 = step_tail.iloc[0]
    x2 = step_tail.iloc[9]
    dvol = (y2 - y1)  # cubic meters
    dtime = (x2 - x1)*timestep  # seconds
    het5_slope = dvol/dtime
    het5_slope_cph = het5_slope/char_vol*60*60
    cph_thresh = 1
    expected = 1/n_bugs
    het5_rv = rrv[-1:].iloc[0]
    fate = "Unidentified"
    if(het5_rv > 0.9*expected):
        fate = "Thriving"
    elif(het5_slope_cph > cph_thresh):
        fate = "Surviving"
    elif(het5_slope_cph <= cph_thresh):
        if(het5_rv > 0.9*expected):
            fate = "Surviving"
        else:
            fate = "Languishing"
    return((fate,expected,het5_rv,het5_slope_cph))


results_file = 'cell_rel_volumes.csv'
results = os.path.join(indir, results_file)
rv = pd.read_csv(results)

biggest_loser = 0
with open(os.path.join(indir, 'biggest_loser.txt')) as loser_file:
    for line in loser_file:
        biggest_loser = int(line.strip())

# TODO generate this using MxN        
#sites = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
sites = range(1,grid_N*grid_N+1)

header = 'seed,ks,mu,yield,colony,category,biggest_loser,expected,rv,cph\n'

out_file = 'colony_outcomes.csv'
output = os.path.join(indir, out_file)

with open(output, 'w') as f:
    f.write(header)
    for site in sites:
        site_prefix = f'het{site}_'
        # print(f'Site {site_prefix} is {site_fate(site_prefix,rv)}')
        category, expected, rvout, cph = site_fate(site_prefix, rv, grid_N*grid_N, deltaT)
        f.write(f'{seed},{ks},{mu},{bioyield},{site},{category},{biggest_loser==site},{expected:.2E},{rvout:.4E},{cph:.4E}\n')
