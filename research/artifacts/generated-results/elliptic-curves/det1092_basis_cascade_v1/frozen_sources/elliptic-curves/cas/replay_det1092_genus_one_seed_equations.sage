#!/usr/bin/env sage-python
"""Equation-only execution on the unchanged nine addresses; no point oracle.

Exact quartic-square evaluation, rational maps, and complete finite group
rank proof. Read no saved unlock coordinates, later points, or V3 artifacts.
The calibrated member is not relabelled as prospectively selected.25s cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves';DIR=ART/'det1092_norm8_seed_cover_v2'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(v,p):
    v=QQ(v);assert v.denominator()%p
    return int(v.numerator()*v.denominator().inverse_mod(p)%p)
def group_quotient(p,a,b,c):
    roots={}
    for y in range(p):roots.setdefault(y*y%p,[]).append(y)
    pts=[None]+[(x,y) for x in range(p) for y in roots.get((x**3+a*x*x+b*x+c)%p,[])]
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;xx,yy=Q
        if x==xx and (y+yy)%p==0:return None
        if x==xx:m=(3*x*x+2*a*x+b)*pow(2*y,-1,p)%p
        else:m=(yy-y)*pow(xx-x,-1,p)%p
        rx=(m*m-a-x-xx)%p;return (rx,(m*(x-rx)-y)%p)
    doubled={add(P,P) for P in pts};labels={P:0 for P in doubled};reps=[None];dim=0
    for P in pts:
        if P in labels:continue
        new=[]
        for j,Q in enumerate(reps):
            R=add(P,Q);new.append(R)
            for D in doubled:
                V=add(R,D);assert V not in labels;labels[V]=j|(1<<dim)
        reps+=new;dim+=1
    assert set(labels)==set(pts) and len(pts)==len(doubled)*2**dim
    return pts,doubled,labels,dim
def main():
    roster_path=ART/'det1092_rr_generic_point_controls_v2/protocol.json';roster=read(roster_path)
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    cases=[];inputs=[roster_path,Path(__file__)]
    for index in range(2):
        payload_path=DIR/f'equation-only-cover-{index:02d}.json';data=read(payload_path);inputs.append(payload_path)
        # Hash verification reads source bytes for binding, never their
        # point/outcome fields as an execution input.
        for path,digest in data['inputs'].items():assert sha(ROOT/path)==digest
        assert not {'historical_point','literal_point','control_outcomes','reconstructed_first_literal_point'}.intersection(data)
        f=R(data['quartic_coefficients']);ai=list(map(dec,data['parent_a_invariants']))
        generic=[[dec(v) for v in P] for P in data['generic_point_functions']]
        x0,x1,y0,y1=[dec(data['maps'][k]) for k in ['x0','x1','y0','y1']]
        exclusions=[f,R(EllipticCurve(K,ai).discriminant())]
        exclusions += [v.denominator() for v in [*ai,x0,x1,y0,y1,*[a for P in generic for a in P]]]
        for case in roster['cases']:
            tau=QQ(case['parameter']);assert all(v(tau) for v in exclusions)
            value=f(tau)
            row={'cover_index':index,'label':case['label'],'parameter':str(tau)}
            if value<0 or not value.numerator().is_square() or not value.denominator().is_square():
                row['status']='NONSPLIT_THIS_COVER_ONLY';cases.append(row);continue
            square=value.sqrt();E=EllipticCurve(QQ,[v(tau) for v in ai])
            P=E([x0(tau)+x1(tau)*square,y0(tau)+y1(tau)*square])
            points=[E([v(tau) for v in Q]) for Q in generic]+[P]
            rows=[];certificates=[]
            for p in data['finite_certificate_primes']:
                if any(v.denominator()%p==0 for v in E.a_invariants()) or E.discriminant().valuation(p)!=0:
                    certificates.append({'p':p,'status':'SKIP_EQUATION_NOT_GOOD_INTEGRAL'});continue
                a,b,c=[residue(v,p) for v in [E.b2(),8*E.b4(),16*E.b6()]]
                pts,doubled,labels,dim=group_quotient(p,a,b,c)
                reductions=[]
                for Q in points:
                    if Q.is_zero() or Q[0].valuation(p)<0:reductions.append(None)
                    else:reductions.append((residue(4*Q[0],p),residue(8*Q[1]+4*E.a1()*Q[0]+4*E.a3(),p)))
                assert all(Q in labels for Q in reductions)
                codes=[labels[Q] for Q in reductions]
                block=[[(code>>j)&1 for code in codes] for j in range(dim)];rows+=block
                certificates.append({'p':p,'status':'EXACT_FINITE_GROUP','points':pts,
                  'doubled':sorted(doubled,key=lambda P:(P is not None,P)),
                  'point_reductions':reductions,'rows':block})
            M=matrix(GF(2),rows);assert M.ncols()==18 and M.rank()==18 and M[:,:17].rank()==17
            p=data['no_two_torsion_prime'];assert E.discriminant().valuation(p)==0
            a,b,c=[residue(v,p) for v in [E.b2(),8*E.b4(),16*E.b6()]]
            assert all((x**3+a*x*x+b*x+c)%p for x in range(p))
            row.update(status='CERTIFIED_RANK18_FROM_EQUATIONS',rank_lower_bound=18,
                       generic_rank=17,no_two_torsion_prime=p,literal_point=list(map(str,P[:2])),
                       quartic_point=['0' if tau==0 else str(tau),str(square)],
                       finite_group_certificates=certificates)
            cases.append(row)
    gains=[r for r in cases if r['status']=='CERTIFIED_RANK18_FROM_EQUATIONS']
    assert len(gains)==1 and gains[0]['cover_index']==1 and gains[0]['parameter']=='0'
    result={'status':'PASS_EQUATION_ONLY_GENUS_ONE_FIRST_SEED_REPLAY',
      'classification':'verified equation-only execution of a retrospective calibration',
      'cases':cases,'execution_exceptional_point_inputs':0,'later_point_inputs':0,'point_searches':0,
      'boundary':'The known first point is not an execution input. It was used to choose cover01, so this is not a prospective discovery or a seed-existence exclusion on nonsplit controls.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in inputs}}
    retain(DIR/'equation-only-replay.json',result)
    print(result['status'],'rank18 at302; seventeen cover/address nonsplits',flush=True)
if __name__=='__main__':
    signal.alarm(25);main()
