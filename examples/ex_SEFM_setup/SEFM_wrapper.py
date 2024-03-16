import numpy as np
import pandas as pd
import os
from functions import NSE, KGE, PBIAS, pkerr, call_func
from contextlib import contextmanager
import threading
import _thread
import subprocess
import shutil

dir = os.getcwd()

# Read specification file
vars = dict()
with open("obj_in.txt") as f:
    for line in f:
        line = line.partition('#')[0]
        eq_index = line.find('=')
        var_name = line[:eq_index].strip()
        value = line[eq_index + 1:].strip()
        vars[var_name] = value

# reservoir name
res = vars['res']

# event names
evs = vars['events'].replace(" ","").split(",")

# SEFM project folder name
folder = vars['proj_folder']

# objective functions
obj_fns = vars['obj_fns'].replace(" ","").split(",")

## Run SEFM, if model crashes re-run simulation until it succeeds
# delete existing simulation results
if os.path.isdir(f'{dir}\{folder}\Simulations'):
    shutil.rmtree(f'{dir}\{folder}\Simulations')

# if no simulations exist, run SEFM
print('Attempting SEFM simulation...')
while not os.path.isdir(f'{dir}\{folder}\Simulations'):
    class TimeoutException(Exception):
        def __init__(self, msg=''):
            self.msg = msg

# if SEFM crashes, return to top of loop, if successful exit
    @contextmanager
    def time_limit(seconds):
        timer = threading.Timer(seconds, lambda: _thread.interrupt_main())
        timer.start()
        try:
            subprocess.call("execute.bat")
        except KeyboardInterrupt:
            #raise TimeoutException("Timed out for operation {}".format(msg))
            print('SEFM failed...attempting simulation again...')
        finally:
            # if the action ends in specified time, timer is canceled
            print('SEFM simulation successful')
            timer.cancel()

    import time
    # ends after 5 seconds
    time_limit(5)


# map input objectives to functions
fn_mapper = {'NSE': NSE, 'KGE': KGE, 'pkerr': pkerr, 'PBIAS': PBIAS}

obj = np.zeros((len(evs),len(obj_fns)))
for i,e in enumerate(evs):
    # read observed data
    obs = pd.read_csv(f'observations/{e}_streamflow.csv')
    #TODO: add code to confirm time series coincide

    # find most recent simulation
    all_sims = []
    for d in os.listdir(f'{dir}\{folder}\Simulations\_EventCal\Storm{i+1}'):
        if os.path.isdir(f'{dir}\{folder}\Simulations\_EventCal\Storm{i+1}\{d}'):
            all_sims.append(f'{dir}\{folder}\Simulations\_EventCal\Storm{i+1}\{d}')
    sim_dir = max(all_sims, key=os.path.getmtime)

    #read model output
    sim = pd.read_csv(f'{sim_dir}\\{res}_Res.csv',skiprows=7,names=['Date','Flow','Outflow','Bypass','WSEL'])

    #calculate objective functions
    x = sim['Flow']
    y = obs['Flow']

    for j,fn in enumerate(obj_fns):
        val = call_func(x,y,fn,fn_mapper)

        obj[i,j] = val

output = pd.DataFrame(obj,columns=obj_fns,index=evs)
output.to_csv('objectives.csv')

if os.path.isdir(f'{dir}\{folder}\Simulations'):
    shutil.rmtree(f'{dir}\{folder}\Simulations')