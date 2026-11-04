"""Run local restricted inputs. Generated artifacts must not be uploaded by default."""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.cluster import KMeans, kmeans_plusplus
from threadpoolctl import threadpool_limits
from .environment import environment_manifest
from .core import prepare, Evaluator, metrics, describe, violation, partition_key
from .optimizers import initial_population, mopso, nsga2, fronts, objective_front


def serial(value):
    if isinstance(value,np.ndarray): return serial(value.tolist())
    if isinstance(value,np.generic): return serial(value.item())
    if isinstance(value,float) and not np.isfinite(value): return None
    if isinstance(value,dict): return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [serial(v) for v in value]
    return value


def save(path,value):
    Path(path).write_text(json.dumps(serial(value),indent=2,allow_nan=False),encoding='utf-8')


def choose(pool,bounds):
    feasible=[p for p in pool if p['v']==0]
    if not feasible: return None
    nondominated=[feasible[i] for i in objective_front(feasible)]
    lo,hi=bounds; width=hi-lo; active=width>0
    def key(p):
        dist=float(np.sum(((p['f'][active]-lo[active])/width[active])**2))
        return (dist,-p['info']['coverage'],-p['info']['minimum_cluster_size'],tuple(p['decoded']))
    winner=min(nondominated,key=key).copy(); winner['distance_squared']=key(winner)[0]
    return winner


def report_metrics(z,labels,method,include_dbcv):
    try:
        score=metrics(z,labels,complementary=True)
    except ValueError as error:
        score={'SI':None,'DBI':None,'CHI':None,'metric_status':str(error)}
    if include_dbcv and method!='kmeans':
        from hdbscan.validity import validity_index
        assigned=labels>=0
        _,relabeled=np.unique(labels[assigned],return_inverse=True)
        if len(np.unique(relabeled))<2:
            score.update(DBCV=None,DBCV_status='undefined: fewer than two assigned clusters')
        else:
            try:
                value=float(validity_index(z[assigned].astype(np.float64),relabeled,metric='euclidean'))
                score.update(DBCV=value if np.isfinite(value) else None,
                    DBCV_status='valid' if np.isfinite(value) else 'nonfinite numerical result')
            except (ValueError,ZeroDivisionError,FloatingPointError) as error:
                score.update(DBCV=None,DBCV_status=str(error))
    return score


def run(input_path,config_path,output):
    config=json.loads(Path(config_path).read_text()); out=Path(output)
    out.mkdir(parents=True,exist_ok=False)  # protect previous results
    frame=pd.read_csv(input_path,dtype={'sample_id':str},float_precision='round_trip')
    if config.get('expected_samples') is not None and len(frame)!=config['expected_samples']:
        raise ValueError('Sample count differs from configuration; change only for an explicitly labelled experiment.')
    z,preprocessing=prepare(frame,config['units'])
    save(out/'preprocessing.json',preprocessing)
    save(out/'run_manifest.json',{'configuration':config,'input_sha256':hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
        'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,
        'pandas':pd.__version__,'sklearn':sklearn.__version__,'thread_limit':1,'installed_environment':environment_manifest(),
        'status':'REIMPLEMENTATION; manuscript output agreement not established'})
    summary=[]; audit=[]
    for method in config['methods']:
        full_scores={}
        pools={'baseline':[],'mopso':[],'nsga2':[]}; runs=[]
        ks=config['k_values'] if method=='kmeans' else [None]
        for k in ks:
            base_ev=Evaluator(z,method,k,config['silhouette_subset_size'])
            for seed in (config['seeds'] if method=='kmeans' else [config['seeds'][0]]):
                if method=='kmeans':
                    # Explicit 30-start baseline retains its selected *initial* centroids
                    # for the manuscript's fixed-initialization subsampling assessment.
                    baseline_rng=np.random.RandomState(seed)
                    best_fit=None; best_initial=None
                    for restart in range(config['baseline_n_init']):
                        start=kmeans_plusplus(z,n_clusters=k,random_state=baseline_rng)[0]
                        candidate=KMeans(n_clusters=k,init=start,n_init=1,random_state=seed,
                            max_iter=300,tol=1e-4,algorithm='lloyd').fit(z)
                        if best_fit is None or candidate.inertia_<best_fit.inertia_:
                            best_fit=candidate;best_initial=start.copy()
                    labels=best_fit.labels_; decoded=tuple(best_initial.ravel())
                else:
                    decoded=(base_ev.baseline_radius,8); labels=base_ev.fit(decoded)
                info=describe(labels); v=violation(info,method!='kmeans')
                try: score=metrics(z,labels); f=np.array([-score['SI'],score['DBI']])
                except ValueError: f=np.array([np.inf,np.inf]); v=np.inf
                rec={'labels':labels,'decoded':decoded,'f':f,'v':v,'info':info,'key':partition_key(labels),'K_search':k,'seed':seed}
                if v==0 or method!='kmeans': pools['baseline'].append(rec)
                audit.append({'method':method,'configuration':'baseline','K_search':k,'seed':seed,'info':info,'violation':v})
            for seed in config['seeds']:
                # Same initial candidates for both optimizers within each K / seed.
                initial=initial_population(base_ev,seed,config['population'])
                for name,optimizer in [('mopso',mopso),('nsga2',nsga2)]:
                    ev=Evaluator(z,method,k,config['silhouette_subset_size'])
                    optimizer(ev,initial.copy(),seed,config['generations'])
                    expected=config['population']*(1+config['generations'])
                    assert ev.requests==expected
                    candidates=[]
                    for rec in ev.partitions.values():
                        if rec['key'] not in full_scores:
                            full_scores[rec['key']]=metrics(z,rec['labels'])
                        score=full_scores[rec['key']]
                        candidates.append({**rec,'f':np.array([-score['SI'],score['DBI']]),'K_search':k,'seed':seed})
                    pools[name].extend(candidates)
                    runs.append((name,k,seed,candidates))
                    prefix=f'{method}_{name}_k{k}_seed{seed}'
                    save(out/f'{prefix}_requests.json',ev.logs)
                    save(out/f'{prefix}_audit.json',{'requests':ev.requests,'decoded_fits_or_extractions':ev.fits,
                        'fit_seconds':ev.fit_seconds,'unique_feasible_partitions':len(candidates),
                        'optics_ordering_fits':len(ev.orderings),'bounds_lower':ev.lo,'bounds_upper':ev.hi,
                        'silhouette_subset_indices':ev.subset,'initial_population':initial})
        # Common normalization; baselines included only for K-means as in Methods.
        all_candidates=pools['mopso']+pools['nsga2']+(pools['baseline'] if method=='kmeans' else [])
        unique={p['key']:p for p in all_candidates}
        if unique:
            values=np.array([p['f'] for p in unique.values()]); bounds=(values.min(axis=0),values.max(axis=0))
            save(out/f'{method}_normalization.json',{'min':bounds[0],'max':bounds[1],'unique_partition_count':len(unique)})
        else:
            bounds=None
            audit.append({'method':method,'status':'no feasible optimized partitions; no constraints relaxed'})
        for name,k,seed,candidates in runs:
            winner=choose(candidates,bounds) if bounds is not None else None
            audit.append({'method':method,'configuration':name,'K_search':k,'seed':seed,
                'selected':None if winner is None else {key:winner[key] for key in ['decoded','info','f','key']}})
        for name,pool in pools.items():
            # Fixed density baselines remain reportable references even when an
            # optimizer's feasibility constraints are not met.
            if name=='baseline' and method!='kmeans' and pool:
                winner={**pool[0],'distance_squared':None}
            else:
                winner=choose(pool,bounds) if bounds is not None else None
            if winner is None:
                summary.append({'algorithm':method,'configuration':name,'status':'infeasible'}); continue
            labels=winner['labels']; score=report_metrics(z,labels,method,config.get('dbcv',False))
            summary.append({'algorithm':method,'configuration':name,**winner['info'],**score,
                'distance_squared':winner['distance_squared'],
                'status':'feasible' if winner['v']==0 else 'baseline_constraint_violation',
                'constraint_violation':winner['v']})
            save(out/f'{method}_{name}_selected.json',{key:value for key,value in winner.items() if key!='labels'})
            pd.DataFrame({'sample_id':frame.sample_id,'cluster':labels}).to_csv(out/f'{method}_{name}_labels.csv',index=False)
            save(out/f'{method}_{name}_pareto.json',[
                {key:pool[i][key] for key in ['decoded','f','info','key','K_search','seed']}
                for i in objective_front(pool)])
    pd.DataFrame(summary).to_csv(out/'validation_summary.csv',index=False)
    save(out/'run_selections.json',audit)


def main():
    parser=argparse.ArgumentParser(description='Local reimplementation; requires verified positive processed concentrations.')
    parser.add_argument('--input',required=True); parser.add_argument('--config',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    with threadpool_limits(limits=1): run(args.input,args.config,args.output)

if __name__=='__main__': main()
