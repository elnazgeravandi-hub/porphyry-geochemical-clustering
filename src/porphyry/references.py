"""Interpretation-only enrichment, PCA selection and continuous C-N regression."""
from itertools import combinations
import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import rankdata, spearmanr
from sklearn.decomposition import PCA


def enrichment(concentrations):
    n=len(concentrations)
    return np.column_stack([(rankdata(concentrations[e],method='average')-.5)/n for e in ['Cu','Au','Mo']]).mean(axis=1)


def pca_reference(z,enrichment_score):
    model=PCA(svd_solver='full').fit(z); scores=model.transform(z)
    candidates=[]
    for j,var in enumerate(model.explained_variance_):
        if var<=np.finfo(float).eps*max(model.explained_variance_[0],1): continue
        correlation=float(spearmanr(scores[:,j],enrichment_score).statistic)
        if np.isfinite(correlation): candidates.append((abs(correlation),model.explained_variance_ratio_[j],-j,correlation))
    if not candidates: raise ValueError('No informative PCA component.')
    _,_,negative_index,correlation=max(candidates); index=-negative_index
    sign=1 if correlation>=0 else -1
    return {'component':index+1,'scores':sign*scores[:,index], 'coefficients':sign*model.components_[index],
            'explained_variance_ratio':model.explained_variance_ratio_[index], 'spearman':abs(correlation)}


def cn_reference(concentrations,min_segment=20,seed=42):
    values=np.asarray(concentrations,float)
    if not np.isfinite(values).all() or (values<=0).any(): raise ValueError('C-N requires positive finite concentrations.')
    thresholds,counts=np.unique(values,return_counts=True)
    x=np.log10(thresholds); y=np.log10(np.cumsum(counts[::-1])[::-1]); n=len(x)
    if n<min_segment: raise ValueError('Not enough distinct thresholds.')
    models=[]
    for segments in [1,2,3]:
        breaks_count=segments-1
        if n<segments*min_segment: break
        def fit(breakpoints):
            b=np.sort(np.asarray(breakpoints)); sizes=np.diff(np.r_[0,np.searchsorted(x,b,side='right'),n])
            if (sizes<min_segment).any(): return np.inf,None
            design=np.column_stack([np.ones(n),x]+[np.maximum(x-point,0) for point in b])
            beta=np.linalg.lstsq(design,y,rcond=None)[0]; residual=y-design@beta
            return float(residual@residual),beta
        if breaks_count:
            result=differential_evolution(lambda b:fit(b)[0],[(x[min_segment-1],x[-min_segment-1])]*breaks_count,
                seed=seed,popsize=20,maxiter=500,tol=1e-9,polish=False)
            breakpoints=np.sort(result.x); rss,beta=fit(breakpoints)
        else:
            breakpoints=np.array([]); rss,beta=fit(breakpoints)
        if not np.isfinite(rss): continue
        parameters=2+2*breaks_count
        bic=n*np.log(max(rss,np.finfo(float).eps*max(float(np.sum((y-y.mean())**2)),1.))/n)+parameters*np.log(n)
        models.append({'segments':segments,'breakpoints_log10':breakpoints,'thresholds':10**breakpoints,
                       'coefficients':beta,'rss':rss,'bic':bic})
    selected=models[0]
    for model in models[1:]:
        if model['segments']==selected['segments']+1 and selected['bic']-model['bic']>=6: selected=model
        else: break
    return {'selected':selected,'models':models,'distinct_thresholds':thresholds,'cumulative_counts':10**y,
            'note':'Numerical breakpoint search; global least-squares optimum not certified.'}
