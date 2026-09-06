#!/usr/bin/env python3
"""Height-proposed, exactly verified coordinates for two retrospective quartets."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import subprocess
import traceback
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'PAIRED_QUARTET_RELATIONS_PROTOCOL.json'
COHORT=r.OUT/'rank_jump_completed_cohort_inputs_v1.json'
LIFTS=r.OUT/'rank_jump_solubility_first_v1.json'
INPUT=r.OUT/'rank_jump_paired_quartet_relations_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_paired_quartet_relations_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-paired-quartet-relations-v1'


def capture():
    cohort=r.read(COHORT);lifts=r.read(LIFTS);cases=[]
    for name in ('08234-003','08234-009'):
        c=next(row for row in cohort['rows'] if row['source_id']==name and row['phase']=='initial')
        l=next(row for row in lifts['rows'] if row['source_id']==name)
        model,basis=r.short(c['model'],c['points'])
        assert len(l['lifted_points'])==4 and l['generic_rank']==17
        r.short(model,[x['point'] for x in l['lifted_points']])
        cases.append({'id':name,'model':list(map(str,model)),'basis':[list(map(str,P)) for P in basis],
                      'rank_certificate':c['rank_certificate'],'lifts':l['lifted_points'],'old_quotient_lower':l['soluble_quotient_rank_lower_bound']})
    r.write_new(INPUT,{'schema':'rank-jump.paired-quartet-relations-inputs.v1','cases':cases,
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,COHORT,LIFTS)}})


def worker(index):
    from sage.all import QQ,RealField,EllipticCurve,matrix,vector,pari,lcm
    c=r.read(INPUT)['cases'][index];E=EllipticCurve(list(map(QQ,c['model'])))
    points=[E(list(map(QQ,P))) for P in c['basis']]+[E(list(map(QQ,P['point']))) for P in c['lifts']]
    n=len(c['basis']);pari.set_real_precision(100)
    e=pari.ellinit(E.a_invariants());H=pari.ellheightmatrix(e,[[P[0],P[1]] for P in points])
    RF=RealField(300);G=matrix(RF,[[RF(str(H[i,j])) for j in range(n+4)] for i in range(n+4)])
    coordinates=[];relations=[];proposal_residuals=[]
    for j in range(4):
        solution=G[:n,:n].solve_right(G[:n,n+j].column(0))
        coeff=[QQ(str(F(str(x)).limit_denominator(64))) for x in solution]
        den=int(lcm([x.denominator() for x in coeff]));ints=[int(den*x) for x in coeff]
        residual=den*points[n+j]-sum((a*P for a,P in zip(ints,points[:n])),E(0))
        assert residual.is_zero(), 'Numerical coordinates failed exact group addition'
        proposal_residuals.append(str(max(abs(RF(a)-b) for a,b in zip(coeff,solution))))
        coordinates.append(list(map(str,coeff)));relations.append({'lift_multiplier':den,'basis_coefficients':ints,'exact_residual':'O'})
    Q=matrix(QQ,[[QQ(x) for x in row[17:]] for row in coordinates]);kernel=Q.left_kernel().basis_matrix()
    integral_kernel=[];generic_relations=[]
    for v in kernel.rows():
        den=lcm([x.denominator() for x in v]);vv=[int(den*x) for x in v]
        generic=[sum(vv[j]*QQ(coordinates[j][i]) for j in range(4)) for i in range(17)]
        integral_kernel.append(vv);generic_relations.append(list(map(str,generic)))
    r.write_new(WORK/(c['id']+'_result.json'),{'id':c['id'],'status':'PASS','basis_rank':n,'generic_rank':17,
                 'labels':[x['label'] for x in c['lifts']],'coordinates_in_witness_basis':coordinates,
                 'exact_relations':relations,'numerical_proposal_max_errors':proposal_residuals,
                 'exact_quotient_rank':int(Q.rank()),'quotient_coordinates':[list(map(str,row)) for row in Q.rows()],
                 'kernel_integer_vectors':integral_kernel,'kernel_generic_coordinates':generic_relations})
    print(c['id'],'quotient rank',Q.rank(),'kernel',integral_kernel,flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i,c in enumerate(r.read(INPUT)['cases']):
        log=WORK/(c['id']+'.log');execution=WORK/(c['id']+'_execution.json')
        if not log.exists():
            with log.open('x') as out:
                try:
                    p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i)],stdout=out,stderr=out,timeout=60)
                    status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
                except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
            r.write_new(execution,status)
        path=WORK/(c['id']+'_result.json')
        rows.append({'id':c['id'],'execution':r.read(execution),'log':log.read_text(),
                     'result':r.read(path) if path.exists() else {'status':'UNKNOWN','quotient_rank_interval':[3,4]}})
    r.write_new(OUTPUT,{'schema':'rank-jump.paired-quartet-relations.v1','rows':rows,
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
                       'boundary':'Retrospective coordinates in an already independent witness subgroup. Heights propose; exact group addition certifies. Full curve ranks remain UNKNOWN.'})
    for row in rows:print(row['id'],row['execution'],row['result'].get('exact_quotient_rank'),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker']);p.add_argument('--case',type=int);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        try:worker(a.case)
        except Exception:traceback.print_exc();raise
