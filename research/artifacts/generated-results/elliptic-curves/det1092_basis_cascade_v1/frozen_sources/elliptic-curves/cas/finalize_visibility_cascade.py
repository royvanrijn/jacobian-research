#!/usr/bin/env python3
"""Audit finite-atlas minima and export the full stage table and scientific plot."""
import csv,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
OUT=ART/'visibility_cascade_summary_v1.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))
from research_runtime.store import checkpoint

def build():
    if OUT.exists():raise FileExistsError('preserve cascade summary')
    report=read(ART/'adaptive_visibility_cascade_v1.json')
    selector=read(LOCAL/'adaptive-visibility-cascade-v1/selection-replay.json')
    if selector['status']!='PASS_GENERIC_ONLY_SELECTION_REPLAY' or selector['charts']!=sum(r['charts'] for r in report['rows']):raise ArithmeticError('selector replay incomplete')
    inputs={rel(p):sha(p) for p in (Path(__file__),ART/'adaptive_visibility_cascade_v1.json',LOCAL/'adaptive-visibility-cascade-v1/selection-replay.json')}
    oldpaths={17:LOCAL/'curve302-focused-point-exposure-v2/curve302-generic17/maps.json',**{rank:LOCAL/f'curve302-recovered-subgroup-followup-v1/wave-{wave:02d}/maps.json' for rank,wave in ((19,1),(22,2),(24,3))}}
    old={rank:read(path)['centres'] for rank,path in oldpaths.items()};inputs.update({rel(path):sha(path) for path in oldpaths.values()})
    tables=[]
    for dirname,filename in [('visibility-cascade-matrix-v1','curve302_visibility_cascade_matrix_v1.json'),('visibility-cascade-support-matrix-v1','curve302_visibility_cascade_support_matrix_v1.json')]:
        data=read(ART/filename);inputs[rel(ART/filename)]=sha(ART/filename);kept={}
        for cell in data['cells']:
            name=cell['direction'];rank=cell['rank'];path=LOCAL/dirname/f'{name}-M{rank}.json';full=read(path);inputs[rel(path)]=sha(path)
            records=kept.get(name,[])+full['trials']
            unique={json.dumps([r['centre_word']+[0]*(rank-len(r['centre_word'])),r['translation_word']+[0]*(rank-len(r['translation_word']))]):r for r in records}
            kept[name]=list(unique.values());finite=[r for r in kept[name] if r['height'] is not None]
            winner=min(finite,key=lambda r:int(r['height'])) if finite else None
            cheap=sum(int(r['height'])<=125000 for r in finite)
            if winner!=cell['winner'] or cheap!=cell['cheap_representatives'] or len(kept[name])!=cell['retained_witness_count']:raise ArithmeticError('finite-atlas minimum/count does not replay')
            if cell['contained_by_exact_endpoint']!=any(r['status']=='EXACT_SUBGROUP_CONTAINMENT_ENDPOINT' for r in kept[name]):raise ArithmeticError('endpoint containment summary differs')
            w=winner;word=w['centre_word']+[0]*(rank-len(w['centre_word']));mask=sum((x%2)<<i for i,x in enumerate(word))
            tables.append({'atlas':'earliest_support' if 'support' in dirname else 'historical_introduction','direction':name,'subgroup_rank':rank,
                'height':w['height'],'coordinate':'/'.join(w['coordinate']),'cheap_representatives':cheap,'contained':cell['contained_by_exact_endpoint'],
                'base_orbit':json.dumps(w['orbit'],sort_keys=True),'parity':mask,'extension_bits':mask>>17,'quartic_bits':w['quartic_max_coefficient_bits'],
                'old_policy_parity_selected':any(c['parity']==mask for c in old[rank]) if rank in old else None,
                'old_policy_exact_centre_selected':any(c['representative']==word for c in old[rank]) if rank in old else None})
    tsv=ART/'visibility_cascade_table_v1.tsv'
    if tsv.exists():raise FileExistsError('preserve full matrix table')
    with tsv.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(tables[0]),delimiter='\t');writer.writeheader();writer.writerows(tables)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    support=[c for c in tables if c['atlas']=='earliest_support'];names=sorted({c['direction'] for c in support});ranks=sorted({c['subgroup_rank'] for c in support})
    image=np.full((len(names),len(ranks)),np.nan)
    for c in support:
        if not c['contained']:image[names.index(c['direction']),ranks.index(c['subgroup_rank'])]=math.log10(int(c['height']))
    fig,ax=plt.subplots(figsize=(12,6.8));cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#d9d9d9')
    plot=ax.imshow(image,aspect='auto',cmap=cmap,vmin=0,vmax=180)
    ax.set_xticks(range(len(ranks)),[str(r) for r in ranks]);ax.set_yticks(range(len(names)),names,fontsize=9)
    ax.set_xlabel('Historical certified subgroup rank (not the autonomous run)')
    ax.set_title('Curve 302: finite-atlas visibility at earliest centre availability')
    fig.colorbar(plot,ax=ax,label='log10 of minimum vetted projective coordinate height')
    for c in support:
        if not c['contained'] and int(c['height'])<=125000:
            ax.text(ranks.index(c['subgroup_rank']),names.index(c['direction']),c['height'],ha='center',va='center',color='white',fontsize=7)
    fig.text(.02,.015,'Grey: an exact subgroup-containment endpoint was found. Finite vetted atlas; no global coordinate minimum claimed.',fontsize=9)
    fig.tight_layout(rect=(0,.035,1,1));svg=ART/'visibility_cascade_v1.svg';png=ART/'visibility_cascade_v1.png'
    fig.savefig(svg);fig.savefig(png,dpi=160);plt.close(fig)
    inputs.update({rel(p):sha(p) for p in (tsv,svg,png)})
    payload={'status':'PASS_CASCADE_MINIMA_COUNTS_AND_SELECTION','inputs':inputs,'autonomous_path':report['rows'][0]['rank_path'],
      'control_paths':{r['id']:r['rank_path'] for r in report['rows'][1:]},'total_charts':sum(r['charts'] for r in report['rows']),
      'completed_charts':sum(r['completed_charts'] for r in report['rows']),'matrix_cells':len(tables),'selection_replay':selector,
      'boundary':'Historical subgroup and autonomous subgroup of the same rank are distinct. Blank old-policy fields mean that the old policy was not run at that stage. Chart timing is in the autonomous report; matrix cost measures are coordinate height and quartic coefficient bits.'}
    checkpoint(OUT,payload);print('FINAL CASCADE AUDIT PASS',payload['autonomous_path'],len(tables),'matrix cells')

if __name__=='__main__':build()
