"""Complete individual-model diagnostics after sealed, possibly incomplete arms.

An unsupported joint formula remains UNKNOWN. This does not change a selector,
repair an experimental arm, or make point-search calls. Reuse existing maps;
charge missing map construction and all exact diagnostic replays separately.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time


def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2)+'\n')


def run(folder):
    start=time.process_time(); sealed=read(folder/'all-arms-sealed.json')
    for relative,digest in sealed['seal_hashes'].items():
        if sha(folder/relative/'seal.json')!=digest:raise ArithmeticError('arm seal changed')
    runtime=folder/'runtime';sys.path.insert(0,str(runtime/'elliptic-curves/cas'))
    from height_model_benchmark import bound_input, selection_from_portfolio
    from pointed_height_bounds import build as bounds
    from pointed_minimal_neighbours import build as neighbours
    from verify_pointed_height_bounds import verify
    from verify_pointed_minimal_neighbours import verify as verify_neighbours
    from verify_pointed_height_portfolio import run as portfolio
    from memory_rank_certificate import checked_rank
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    from lean_preconditioned_map_receipts import obtain
    from search_observability import point_visibility
    from half_lattice_pointed_sieve import linear_combination
    from sage.all import EllipticCurve,QQ
    # The existing portable portfolio checker binds packets relative to its
    # research root. Keep successor diagnostics inside that frozen root;
    # preserve the earlier diagnostic and its path-binding failure unchanged.
    output=runtime/'height-model-accessibility-v2';output.mkdir(exist_ok=False)
    results=[]
    for armrow in sealed['rows']:
        if not armrow['success']:continue
        arm=folder/armrow['folder'];seed=read(arm/'input.json')['payload']
        gain=read(arm/'gain.json');protocol=read(arm/'protocol.json')
        point=tuple(map(F,gain['points'][-1]));ci=gain['centre']
        dest=output/seed['id'];dest.mkdir()
        model=tuple(map(F,seed['curve']));basis=tuple(tuple(map(F,p)) for p in seed['points'])
        proof=checked_rank(model,basis,seed['primes'],seed['torsion_prime'])
        state=certified_state(model,basis,proof);centre=seed['centres'][ci]
        data_by_key={}
        for policy in ('factor_free','preconditioned_full'):
            # Retained, already checked same-input maps are allowed after sealing.
            candidates=[arm/f'map-{ci:04d}-{policy}/result.json',
                folder/'diagnostics'/seed['id']/armrow['arm']/f'map-{ci:04d}-{policy}/result.json',
                folder/'individual-accessibility-v1'/seed['id']/f'map-{ci:04d}-{policy}/result.json']
            source=next((p for p in candidates if p.exists()),None)
            if source:
                mapping=read(source)['mapping']
                import pari_pointed_backend as backend
                search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},
                    coordinate_policy=mapping['coordinate_policy'])
                backend.validate_map(search,mapping)
                save(dest/f'{policy}-map-provenance.json',{'source':str(source.relative_to(folder)),'sha256':sha(source)})
            else:
                mapping,_,limited=obtain(dest,ci,policy,model,basis,centre,state,protocol)
                if limited:raise ArithmeticError('diagnostic map remains incomplete: '+policy)
            data=bound_input(seed,mapping,policy,protocol)
            save(dest/f'{policy}-input.json',data);data_by_key[policy]=data
        models,construction=neighbours(data_by_key['preconditioned_full'])
        save(dest/'neighbours/construction.json',construction)
        if len(models)!=2:raise ArithmeticError('individual neighbour vocabulary incomplete')
        for i,data in enumerate(models):save(dest/f'neighbours/neighbour-{i}-input.json',data)
        save(dest/'neighbours/minimality-verified.json',verify_neighbours(dest,runtime))
        all_data=[data_by_key['factor_free'],data_by_key['preconditioned_full']]+models
        names=['factor_free','preconditioned_full','neighbours/neighbour-0','neighbours/neighbour-1']
        packets=[]
        for name,data in zip(names,all_data):
            packet=bounds(data);save(dest/(name+'-bounds.json'),packet)
            save(dest/(name+'-verified.json'),verify(packet,runtime));packets.append(packet)
        joint=None;joint_error=None
        if len({m['neighbour_witness']['prime'] for m in models})!=2:
            joint_error='UNKNOWN: same-prime neighbours lie outside the frozen joint formula'
        else:
            try:
                result=portfolio(dest,runtime)
                save(dest/'portfolio.json',result)
                joint=selection_from_portfolio(result,packets[1:])
            except ArithmeticError as error:
                joint_error='UNKNOWN: '+str(error)
        anchor=tuple(map(F,packets[0]['anchor']))
        R=linear_combination(model,(point,anchor),(2,-1))
        E=EllipticCurve(QQ,list(map(str,model)))
        check=2*E(list(map(str,point)))-E(list(map(str,anchor)))
        if not ((R is None and check.is_zero()) or
                (R is not None and tuple(map(str,R))==tuple(str(c) for c in check.xy()))):
            raise ArithmeticError('Sage and rational 2P-Q transport disagree')
        Hx=1 if R is None else max(abs(R[0].numerator),R[0].denominator)
        rows=[]
        for data,packet in zip(all_data,packets):
            search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},
                coordinate_policy=data['mapping']['coordinate_policy'])
            location=point_visibility(search.chart_record(),point)
            n,d=map(int,location['coordinate']);H=max(abs(n),abs(d))
            ev=lambda key:sum(F(c)*n**i*d**(4-i) for i,c in enumerate(packet[key]))
            N,D=ev('numerator'),ev('denominator')
            if N.denominator!=1 or D.denominator!=1:raise ArithmeticError('nonintegral point evaluation')
            g=math.gcd(N.numerator,D.numerator);S=max(abs(N),abs(D))/H**4
            C=int(packet['finite_gcd_divisor']);L=F(packet['real']['lower']);bound=F(C)/L
            if C%g or max(abs(N),abs(D))/g!=Hx:raise ArithmeticError('height cancellation transport differs')
            ratio=F(H**4,Hx)/bound
            if not 0<ratio<=1 or ratio!=L*g/(C*S):raise ArithmeticError('height bound or residual decomposition failed')
            rows.append({'model':data['policy'],'coordinate':[str(n),str(d)],
                'parameter_height':str(H),'D':str(bound),'B2':math.log(C)/4-(math.log(L.numerator)-math.log(L.denominator))/4,
                'actual_distortion':math.log(H)-math.log(Hx)/4,
                'actual_over_sufficient_fourth_power':str(ratio),
                'real_slack_factor':str(S/L),'finite_slack_factor':str(F(C,g)),
                'actual_cancellation_gcd':str(g),'actual_real_norm':str(S),
                'inside_H125000':H<=125000,'guaranteed_at_H125000':F(Hx)<=F(125000**4)/bound})
        best_bound=min(range(4),key=lambda i:(F(rows[i]['D']),i))
        best_actual=min(range(4),key=lambda i:(int(rows[i]['parameter_height']),i))
        report={'case':seed['id'],'source_arm':armrow['folder'],'centre_index':ci,
            'winning_point':gain['points'][-1], 'anchor':list(map(str,anchor)),
            'R_2P_minus_Q':None if R is None else list(map(str,R)), 'x_height':str(Hx),'models':rows,
            'best_bound_model':rows[best_bound]['model'],'smallest_actual_height_model':rows[best_actual]['model'],
            'joint_selection':joint,'joint_unknown_reason':joint_error,
            'neighbour_primes':[m['neighbour_witness']['prime'] for m in models],
            'boundary':'Sealed literal winning point at its fixed winning anchor. Four individual bounds and exact '
                'coordinates replayed. No changed selector, additional point search, or time-to-gain measurement.'}
        save(dest/'report.json',report);results.append(report)
        save(output/'summary.json',{'status':'PARTIAL_DIAGNOSTICS','rows':results,'cpu_seconds':time.process_time()-start})
        print('INDIVIDUAL_ACCESSIBILITY',seed['id'],[r['parameter_height'] for r in rows],joint_error,flush=True)
    save(output/'summary.json',{'status':'PASS_INDIVIDUAL_ACCESSIBILITY_REPLAY','rows':results,
        'cpu_seconds':time.process_time()-start,'source_sha256':sha(Path(__file__)),
        'original_failed_diagnostic_preserved':True,'point_search_calls':0})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True)
    a=p.parse_args();run(a.folder.resolve())
