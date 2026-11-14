"""Reproducible three-method comparison. See README for scope and limitations.

python benchmarks/compare.py reference
python benchmarks/compare.py calibrate
python benchmarks/compare.py validate
python benchmarks/compare.py report
"""
import argparse
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

from synthetic import Target,Surface,run_pde,run_particle,reprice,true_vol

HERE=Path(__file__).resolve().parent
WORK=HERE/'work'
RESULTS=HERE/'results'


def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def reference():
    for label,nx,nr,dt in [('reference',201,101,.002),('reference_check',161,81,.002)]:
        if (WORK/f'{label}.npz').exists():continue
        r=run_pde('density',None,nx,nr,dt,WORK/f'{label}.npz',reference=True,width=1.5)
        dump(WORK/f'{label}.json',r)


def calibrate():
    target=Target(WORK/'reference.npz');allresults=[]
    for method in ['density','moments']:
        for name,nx,nr,dt in [('coarse',61,31,.004),('medium',101,51,.004),('fine',141,71,.002)]:
            label=method+'_'+name
            if not (WORK/f'{label}.json').exists():
                r=run_pde(method,target,nx,nr,dt,WORK/f'{label}.npz',width=1.5)
                r['label']=label;dump(WORK/f'{label}.json',r)
            allresults.append(json.loads((WORK/f'{label}.json').read_text()))
    for particles in [20000,80000,320000]:
        for seed in [17,29,43]:
            label=f'particle_{particles}_{seed}'
            if not (WORK/f'{label}.json').exists():
                r=run_particle(target,particles,.002,seed,WORK/f'{label}.npz')
                r['label']=label;dump(WORK/f'{label}.json',r)
            allresults.append(json.loads((WORK/f'{label}.json').read_text()))
    dump(WORK/'runs.json',allresults)


def validate():
    runs=json.loads((WORK/'runs.json').read_text())
    surfaces={r['label']:Surface(WORK/(r['label']+'.npz')) for r in runs}
    # Three independent particle calibrations at each N, all repriced with the
    # same independent path sample. Calibration-seed variability remains visible.
    rows=reprice(surfaces,particles=120000,dt=.001,seed=98181)
    pd.DataFrame(rows).to_csv(RESULTS/'repricing.csv',index=False)


def report():
    runs=json.loads((WORK/'runs.json').read_text());target=Target(WORK/'reference.npz')
    price=pd.read_csv(RESULTS/'repricing.csv');out=[]
    for r in runs:
        surf=Surface(WORK/(r['label']+'.npz'))
        sq=[];ws=[];mx=0.
        for t in [.25,.5,1.,1.5,2.]:
            x=np.linspace((.025-.015)*t-2*.16*np.sqrt(t),(.025-.015)*t+2*.16*np.sqrt(t),101)
            w=target.density(t,x);w=w/w.sum();d=surf(t,x)-true_vol(t,x)[0]
            sq.extend((w*d*d).tolist());ws.extend(w.tolist());mx=max(mx,float(np.max(abs(d))))
        pr=price[price.label==r['label']]
        out.append(dict(label=r['label'],method=r['method'],runtime=r['runtime'],
                        sigma_rmse_vol_points=100*np.sqrt(np.sum(sq)/np.sum(ws)),
                        sigma_max_vol_points=100*mx,
                        repricing_mean_abs_vol_points=pr.error_vol_points.abs().mean(),
                        repricing_max_abs_vol_points=pr.error_vol_points.abs().max(),
                        repricing_max_se_vol_points=pr.se_vol_points.max(),
                        floor_nodes=r['floor_nodes'],cap_nodes=r['cap_nodes']))
    tab=pd.DataFrame(out);tab.to_csv(RESULTS/'comparison.csv',index=False)
    # Target sensitivity, on a common central strip (same interpolation in time).
    check=Target(WORK/'reference_check.npz');diff=[]
    for t in [.25,.5,1.,1.5,2.]:
        x=np.linspace(.01*t-2*.16*np.sqrt(t),.01*t+2*.16*np.sqrt(t),101)
        diff.extend((target.variance(t,x)-check.variance(t,x)).tolist())
    diagnostics=dict(machine=platform.platform(),python=platform.python_version(),numpy=np.__version__,
                     target_variance_max_change=float(np.max(abs(np.array(diff)))),
                     target_call_max_change=float(np.max(abs(target.calls-check.calls))),
                     reference=json.loads((WORK/'reference.json').read_text()),
                     reference_check=json.loads((WORK/'reference_check.json').read_text()))
    dump(RESULTS/'diagnostics.json',diagnostics)
    dump(RESULTS/'runs.json',runs)
    print(tab.to_string(index=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    colors={'density':'#2457a7','moments':'#ba3c2e','particle':'#21836c'}
    for method,g in tab.groupby('method'):
        axes[0].scatter(g.runtime,g.sigma_rmse_vol_points,label=method,color=colors[method],s=40)
    axes[0].set(xscale='log',yscale='log',xlabel='Calibration wall time (seconds)',ylabel='Local volatility RMSE (vol points)',title='Recovery of the hidden nonflat coefficient')
    axes[0].legend();axes[0].grid(alpha=.2)
    x=np.linspace(-.45,.45,201);axes[1].plot(x,true_vol(2,x)[0]*100,'k',lw=2,label='generating volatility')
    for label in ['density_fine','moments_fine','particle_320000_17']:
        axes[1].plot(x,Surface(WORK/f'{label}.npz')(2,x)*100,label=label,color=colors[label.split('_')[0]],ls='--')
    axes[1].set(xlabel='log spot / initial spot',ylabel='Local volatility (%)',title='Recovered coefficient at two years')
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    fig.savefig(RESULTS/'comparison.png',dpi=180)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['reference','calibrate','validate','report'])
    args=parser.parse_args();WORK.mkdir(exist_ok=True);RESULTS.mkdir(exist_ok=True)
    globals()[args.stage]()
