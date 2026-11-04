"""Local profiles, PCA/C-N and optional catchment-union evaluation.

Figures are rendered from computed outputs, never transcribed Table 7 values.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .environment import environment_manifest
from .core import ELEMENTS,prepare
from .references import enrichment,pca_reference,cn_reference
from .evaluation import capture
from .cli import save
from .preprocess import UNITS


def read_indexed(path):
    data=pd.read_csv(path,dtype={'sample_id':str},float_precision='round_trip')
    if 'sample_id' not in data or data.sample_id.isna().any() or data.sample_id.str.strip().eq('').any() or data.sample_id.duplicated().any():
        raise ValueError('Unique nonempty sample_id column required.')
    return data.set_index('sample_id')


def load_labels(directory,index):
    models={}
    for algorithm in ['kmeans','dbscan','optics']:
        for method in ['baseline','mopso','nsga2']:
            name=f'{algorithm}_{method}'; path=Path(directory)/f'{name}_labels.csv'
            if not path.is_file(): raise FileNotFoundError(f'Required model labels missing: {path.name}')
            data=read_indexed(path)
            if set(data.index)!=set(index): raise ValueError(f'Label IDs do not match concentrations: {name}')
            values=pd.to_numeric(data.loc[index,'cluster'],errors='raise').to_numpy(float)
            if not np.isfinite(values).all() or not (values==np.floor(values)).all() or (values< -1).any():
                raise ValueError(f'Invalid cluster labels: {name}')
            models[name]=values.astype(int)
    return models


def elemental_profiles(concentrations,models):
    common=concentrations[ELEMENTS].median(); rows=[]
    for name,labels in models.items():
        for cluster in sorted(np.unique(labels)):
            part=concentrations.loc[labels==cluster,ELEMENTS]
            for element in ELEMENTS:
                values=part[element]; q25,median,q75=values.quantile([.25,.5,.75])
                rows.append({'model':name,'cluster':int(cluster),'group_type':'Noise' if cluster==-1 else 'Cluster',
                    'n':len(part),'element':element,'unit':UNITS[element],'q25':q25,'median':median,'q75':q75,
                    'common_median':common[element],'relative_q25':q25/common[element],
                    'relative_median':median/common[element],'relative_q75':q75/common[element]})
    return pd.DataFrame(rows)


def calculate_references(concentrations,z):
    score=enrichment(concentrations); pca=pca_reference(z,score)
    cn=cn_reference(concentrations.Cu.to_numpy(),20)
    # Values equal to a breakpoint enter the upper concentration class.
    classes=np.searchsorted(cn['selected']['thresholds'],concentrations.Cu,side='right')+1
    scores=pd.DataFrame({'enrichment_score':score,'selected_pc_score':pca['scores'],
        'cn_class':classes,'Cu':concentrations.Cu},index=concentrations.index)
    return scores,pca,cn


def plot_profiles(profiles,out):
    import matplotlib.pyplot as plt
    algorithms=['kmeans','dbscan','optics']; methods=['baseline','mopso','nsga2']
    fig,axes=plt.subplots(3,3,figsize=(14,11),sharex=True,sharey=True,layout='constrained')
    for r,algorithm in enumerate(algorithms):
        for c,method in enumerate(methods):
            ax=axes[r,c]; table=profiles[profiles.model==f'{algorithm}_{method}']
            for j,(cluster,group) in enumerate(table.groupby('cluster',sort=True)):
                g=group.set_index('element').loc[ELEMENTS]
                color='gray' if cluster==-1 else f'C{j%10}'
                ax.errorbar(np.arange(6)+(j-(table.cluster.nunique()-1)/2)*.025,g.relative_median,
                    yerr=np.vstack([g.relative_median-g.relative_q25,g.relative_q75-g.relative_median]),
                    marker='s' if cluster==-1 else 'o',linestyle='none',capsize=2,color=color,
                    label='Noise' if cluster==-1 else f'Cluster {cluster}')
            ax.axhline(1,color='black',ls=':',lw=.8); ax.set_yscale('log')
            ax.set_xticks(range(6),ELEMENTS); ax.set_title(f'{algorithm} / {method}')
            if c==0: ax.set_ylabel('Concentration / common median')
            ax.legend(fontsize=7)
    fig.savefig(out/'elemental_profiles.png',dpi=300); fig.savefig(out/'elemental_profiles.pdf'); plt.close(fig)


def plot_references(scores,pca,cn,out):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    axes[0].bar(ELEMENTS,pca['coefficients']);axes[0].axhline(0,color='black',lw=.6)
    axes[0].set_title(f'Selected PC{pca["component"]} coefficients');axes[0].set_ylabel('Oriented coefficient')
    thresholds=cn['distinct_thresholds']; counts=cn['cumulative_counts']; model=cn['selected']; x=np.log10(thresholds)
    design=np.column_stack([np.ones(len(x)),x]+[np.maximum(x-b,0) for b in model['breakpoints_log10']])
    axes[1].loglog(thresholds,counts,'.',ms=2,label='Observed cumulative counts')
    axes[1].loglog(thresholds,10**(design@model['coefficients']),label='Selected continuous fit')
    for threshold in model['thresholds']:axes[1].axvline(threshold,color='gray',ls='--')
    axes[1].set_xlabel('Cu (ppm)');axes[1].set_ylabel('N(>= concentration)');axes[1].legend(fontsize=7)
    fig.savefig(out/'pca_cn_reference.png',dpi=300);fig.savefig(out/'pca_cn_reference.pdf');plt.close(fig)


def read_spatial(config_path,index):
    import geopandas as gpd
    from pyproj import CRS
    from shapely import union_all
    config=json.loads(Path(config_path).read_text())
    crs=CRS.from_user_input(config['area_crs'])
    if not crs.is_projected or not crs.axis_info or not all(a.unit_name.lower() in ['metre','meter'] for a in crs.axis_info[:2]):
        raise ValueError('Area CRS must be a documented projected CRS in metres.')
    def load(path):
        layer=gpd.read_file(path)
        if layer.crs is None: raise ValueError('Missing spatial CRS.')
        if layer.empty or layer.geometry.isna().any() or layer.geometry.is_empty.any() or not layer.geometry.is_valid.all():
            raise ValueError('Resolve empty, missing or invalid geometries explicitly.')
        return layer.to_crs(crs)
    basins=load(config['catchments_path']); boundary=load(config['study_boundary_path']); points=load(config['occurrences_path'])
    if not basins.geom_type.isin(['Polygon','MultiPolygon']).all() or not boundary.geom_type.isin(['Polygon','MultiPolygon']).all():
        raise ValueError('Catchments and boundary must be polygons.')
    if not points.geom_type.eq('Point').all():raise ValueError('Occurrences must be point locations.')
    if not config.get('occurrence_inventory_preverified',False) or not str(config.get('eligibility_description','')).strip():
        raise ValueError('Supply a preverified eligible occurrence inventory and document eligibility criteria.')
    id_field=config['catchment_id_field']
    if basins[id_field].isna().any():raise ValueError('Missing catchment sample links.')
    basins[id_field]=basins[id_field].astype(str)
    if not set(index).issubset(set(basins[id_field])):raise ValueError('Some samples lack a catchment link.')
    basins=basins[basins[id_field].isin(index)].copy()
    # Multiple features per sample are dissolved, not duplicated in area/counts.
    basins=basins.dissolve(by=id_field).loc[index]
    roi=union_all(boundary.geometry.to_numpy())
    basins.geometry=basins.geometry.intersection(roi)
    if basins.geometry.is_empty.any() or not basins.geometry.is_valid.all():raise ValueError('Invalid or empty clipped linked catchment.')
    domain=union_all(basins.geometry.to_numpy())
    return basins,domain,list(points.geometry),config


def spatial_assessment(scores,models,basins,domain,points):
    rows=[]; all_results={}
    for name in list(models)+['pca','cn']:
        if name in models:
            labels=models[name]; ranks=pd.Series(scores.enrichment_score.to_numpy(),index=labels).groupby(level=0).median()
            sample_ranks=np.array([ranks[label] if label>=0 else np.nan for label in labels])
        elif name=='pca':sample_ranks=scores.selected_pc_score.to_numpy()
        else:
            ranks=scores.groupby('cn_class').Cu.median()
            sample_ranks=scores.cn_class.map(ranks).to_numpy()
        groups=[]
        for rank in sorted(set(sample_ranks[np.isfinite(sample_ranks)]),reverse=True):
            groups.append((float(rank),list(basins.geometry[sample_ranks==rank])))
        result=capture(groups,domain,points);all_results[name]=result
        for item in result['budgets']:
            rows.append({'model':name,'nominal_budget_percent':100*item['budget'],
                'actual_area_percent':100*item['area_fraction'],'capture_percent':item['capture_percent'],
                'eligible_occurrences':result['eligible_occurrences']})
    return pd.DataFrame(rows),all_results


def plot_spatial(scores,basins,domain,points,captures,out):
    import matplotlib.pyplot as plt
    import geopandas as gpd
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    layer=basins.copy(); layer['PCA']=scores.selected_pc_score.to_numpy();layer['CN']=scores.cn_class.to_numpy()
    for ax,column in zip(axes,['PCA','CN']):
        layer.plot(column=column,ax=ax,legend=True,cmap='viridis',edgecolor='none',alpha=.8)
        gpd.GeoSeries([domain],crs=layer.crs).boundary.plot(ax=ax,color='black',linewidth=.5)
        visible=[p for p in points if domain.covers(p)]
        if visible:gpd.GeoSeries(visible,crs=layer.crs).plot(ax=ax,color='red',markersize=7)
        ax.set_title('Oriented PCA score' if column=='PCA' else 'Cu C-N classes');ax.set_axis_off()
    fig.savefig(out/'reference_maps.png',dpi=300);fig.savefig(out/'reference_maps.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True,layout='constrained')
    for ax,algorithm in zip(axes,['kmeans','dbscan','optics']):
        for method in ['baseline','mopso','nsga2']:
            table=captures[captures.model==f'{algorithm}_{method}'].sort_values('nominal_budget_percent')
            ax.plot(table.nominal_budget_percent,table.capture_percent,marker='o',label=method)
        ax.set_title(algorithm);ax.set_xlabel('Nominal target-area budget (%)');ax.set_xticks([5,10,20]);ax.legend(fontsize=8)
    axes[0].set_ylabel('Eligible occurrence capture (%)');axes[0].set_ylim(0,100)
    fig.savefig(out/'occurrence_capture.png',dpi=300);fig.savefig(out/'occurrence_capture.pdf');plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description='Local computed profiles/reference plots and optional spatial capture.')
    parser.add_argument('--input',required=True);parser.add_argument('--labels-dir',required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--spatial-config')
    args=parser.parse_args();concentrations=read_indexed(args.input)[ELEMENTS]
    z,scaler=prepare(concentrations.reset_index(),UNITS)
    models=load_labels(args.labels_dir,concentrations.index)
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    save(out/'environment.json',environment_manifest())
    profiles=elemental_profiles(concentrations,models);profiles.to_csv(out/'elemental_profiles.csv',index=False)
    scores,pca,cn=calculate_references(concentrations,z);scores.to_csv(out/'reference_scores.csv')
    save(out/'reference_metadata.json',{'pca':{k:v for k,v in pca.items() if k!='scores'},'cn':cn,'scaler':scaler,
        'note':'Computed reimplementation outputs; agreement with manuscript results not established.'})
    sensitivity=[]
    for minimum in [10,20,30]:
        try:
            fit=cn_reference(concentrations.Cu.to_numpy(),minimum)['selected']
            sensitivity.append({'minimum_distinct_thresholds':minimum,**fit,'status':'fitted'})
        except ValueError as exc:sensitivity.append({'minimum_distinct_thresholds':minimum,'status':str(exc)})
    save(out/'cn_threshold_sensitivity.json',sensitivity)
    supporting={}
    for element in ['Au','Mo']:
        try:supporting[element]={'unit':UNITS[element],'fit':cn_reference(concentrations[element].to_numpy(),20)}
        except ValueError as exc:supporting[element]={'unit':UNITS[element],'status':str(exc)}
    save(out/'supporting_cn_references.json',supporting)

    plot_profiles(profiles,out);plot_references(scores,pca,cn,out)
    if args.spatial_config:
        basins,domain,points,config=read_spatial(args.spatial_config,concentrations.index)
        captures,detail=spatial_assessment(scores,models,basins,domain,points)
        captures.to_csv(out/'occurrence_capture.csv',index=False)
        save(out/'spatial_audit.json',{'configuration':config,'domain_area_m2':domain.area,'evaluations':detail,
            'ranking':'Whole cluster or C-N class groups; exact tied ranks merged. PCA ranked by oriented scores.',
            'noise':'Included in common domain; not independently selected as clustering targets.'})
        plot_spatial(scores,basins,domain,points,captures,out)

if __name__=='__main__':main()
