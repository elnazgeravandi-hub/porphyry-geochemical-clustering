"""Partition agreement and whole-group cumulative target evaluation."""
import numpy as np
from sklearn.metrics import adjusted_rand_score
from .core import describe


def agreement(a,b):
    a=np.asarray(a); b=np.asarray(b)
    if a.shape!=b.shape: raise ValueError('Align sample IDs before comparing labels.')
    common=(a>=0)&(b>=0)
    def ari(mask):
        if len(np.unique(a[mask]))<2 or len(np.unique(b[mask]))<2: return None
        return float(adjusted_rand_score(a[mask],b[mask]))
    return {'ARI_all':ari(np.ones(len(a),dtype=bool)),'ARI_common_assigned':ari(common),
            'common_assigned_count':int(common.sum()),'common_assigned_fraction':float(common.mean())}


def sampling_stability(z,full_labels,refit):
    """refit(X_subset) must preserve full-data preprocessing and selected hyperparameters."""
    rows=[]
    for seed in range(42,62):
        indices=np.sort(np.random.default_rng(seed).choice(len(z),int(.8*len(z)),replace=False))
        labels=refit(z[indices])
        rows.append({'seed':seed,**agreement(np.asarray(full_labels)[indices],labels),**describe(labels)})
    return rows


def capture(groups,domain,occurrences,budgets=(.05,.10,.20)):
    """Groups: iterable (rank, geometries), highest ranks first; equal ranks merged.

    Inputs must already share one appropriate projected CRS. Domain and eligibility
    are caller-verified. Geometries are clipped and unioned; points counted once.
    """
    from shapely import union_all, GeometryCollection
    if domain.is_empty or not domain.is_valid or domain.area<=0: raise ValueError('Invalid evaluation domain.')
    points={(float(point.x),float(point.y)):point for point in occurrences if domain.covers(point)}
    if not points: raise ValueError('No eligible unique occurrences inside domain.')
    combined={}
    for rank,geometries in groups:
        if not np.isfinite(rank): raise ValueError('Invalid target rank.')
        combined.setdefault(rank,[]).extend(geometries)
    cumulative=GeometryCollection(); steps=[{'area_fraction':0.,'capture_percent':0.}]
    for rank in sorted(combined,reverse=True):
        geoms=combined[rank]
        if any(not g.is_valid for g in geoms): raise ValueError('Resolve invalid target geometries explicitly.')
        cumulative=union_all([cumulative,*geoms]).intersection(domain)
        steps.append({'rank':float(rank),'area_fraction':cumulative.area/domain.area,
            'capture_percent':100*sum(cumulative.covers(p) for p in points.values())/len(points)})
    return {'eligible_occurrences':len(points),'steps':steps,'budgets':[
        {'budget':budget,**max((s for s in steps if s['area_fraction']<=budget),key=lambda s:s['area_fraction'])}
        for budget in budgets]}
