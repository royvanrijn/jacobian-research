#!/usr/bin/env sage-python
"""Independent replay: no classifier/producer imports, no rational point search.

Reconstruct finite group quotients using Sage's elliptic group implementation;
check exact doubling/cycle equations. Nonhalving terminals receive an elementary
root-free modular quartic certificate at one prime <=257, a fixed proof cap.
Default checks the immutable replay; --write creates it. Run under timeout25s.
"""
import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, GF, EllipticCurve, PolynomialRing, matrix, vector, prime_range, gcd, lcm

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves/det1092_split_descent_v1'
OUT=ART/'independent-replay.json'


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(ok,message):
    if not ok:raise ArithmeticError(message)
def point(E,p):return E(0) if p is None else E(list(map(QQ,p)))
def word(E,basis,v):return sum((ZZ(n)*P for n,P in zip(v,basis)),E(0))
def key(P):return tuple(int(v) for v in P)


def replay():
    protocol=read(ART/'protocol.json')
    for name,h in {**protocol['sources'],**protocol['inputs']}.items():check(sha(ROOT/name)==h,'bound input changed')
    inputs={str((ART/'protocol.json').relative_to(ROOT)):sha(ART/'protocol.json')}
    rows=[]
    for family in ['302',*protocol['parameters']]:
        folder=ART/family;frame_path=folder/'generic-frame.json';frame=read(frame_path)
        inputs[str(frame_path.relative_to(ROOT))]=sha(frame_path)
        E=EllipticCurve(QQ,list(map(QQ,frame['curve'])))
        basis=[point(E,p) for p in frame['basis']]
        groups=[];base_rows=[];torsion=False
        for record in frame['records']:
            p=record['prime'];check(ZZ(p).is_prime(proof=True) and p>2,'prime')
            check(all(a.denominator()%p for a in E.a_invariants()) and E.discriminant()%p,'good reduction')
            ep=EllipticCurve(GF(p),E.a_invariants());points=sorted(ep.points(),key=key)
            doubles={key(2*P):2*P for P in points}
            labels={k:0 for k in doubles};reps=[ep(0)];dim=0
            for P in points:
                if key(P) in labels:continue
                new=[]
                for i,Q in enumerate(reps):
                    V=P+Q;new.append(V)
                    for D in doubles.values():
                        k=key(V+D);check(k not in labels,'coset overlap');labels[k]=i|(1<<dim)
                reps+=new;dim+=1
            check(len(labels)==len(points)==len(doubles)*2**dim,'quotient completeness')
            check(len(points)==record['order'] and dim==record['dimension'],'finite group size')
            if p==frame['no_two_torsion_prime']:
                check(len(points)%2==1,'torsion prime');torsion=True
            groups.append((p,ep,labels,dim))
        def code(P):
            out=[]
            for p,ep,labels,dim in groups:
                den=lcm([v.denominator() for v in P]);v=[ZZ(a*den) for a in P];g=gcd(v)
                red=ep([int(a/g % p) for a in v]);c=labels[key(red)]
                out.extend((c>>j)&1 for j in range(dim))
            return vector(GF(2),out)
        M=matrix(GF(2),[code(P) for P in basis]).transpose()
        check(torsion and M.rank()==17,'independent generic footprint')
        for path in sorted(folder.glob('*.json')):
            if path==frame_path:continue
            result=read(path);inputs[str(path.relative_to(ROOT))]=sha(path)
            check(result['protocol_sha256']==sha(ART/'protocol.json') and result['frame_sha256']==sha(frame_path),'case binding')
            omitted=result['omitted_basis_index'];B=[P for i,P in enumerate(basis) if i!=omitted]
            A=matrix(GF(2),[code(P) for P in B]).transpose();r=len(B)
            check(A.rank()==r,'retained basis rank')
            original=point(E,result['original']);current=original
            for step,h in enumerate(result['history']):
                check(step==h['step'] and current==point(E,h['point']),'history order')
                check(len(h['parity_word'])==r and set(h['parity_word'])<={0,1},'parity word')
                target=current-word(E,B,h['parity_word'])
                check(target==point(E,h['target']) and not any(code(target)),'unique footprint-compatible parity')
                if 'next' in h:
                    half=point(E,h['next']);check(2*half==target,'doubling identity');current=half
                else:
                    check(step==len(result['history'])-1 and result['reason']=='GLOBAL_NONHALVING_ESCAPE','incomplete nonterminal half')
            record=dict(family=family,label=result['label'],status=result['status'],reason=result.get('reason'),steps=result['steps'])
            if result['status']=='INHERITED_RATIONAL_SPAN':
                n=ZZ(result['relation_multiplier']);check(n>0 and n%2==1,'positive odd relation')
                check(n*original==word(E,B,result['relation_word']),'exact dependence identity')
                record.update(relation_multiplier=str(n),relation_word=result['relation_word'])
            elif result['status']=='NEW_INDEPENDENT_DIRECTION':
                check(current==point(E,result['terminal']),'terminal point')
                check(original==2**result['steps']*current+word(E,B,result['chain_word']),'original-to-terminal identity')
                if result['reason']=='FINITE_FOOTPRINT_ESCAPE':
                    check(A.augment(matrix(GF(2),len(code(current)),1,list(code(current)))).rank()==r+1,'finite independent column')
                else:
                    check(result['reason']=='GLOBAL_NONHALVING_ESCAPE','new-direction reason')
                    target=point(E,result['history'][-1]['target']);check(not target.is_zero(),'nonzero target')
                    R=PolynomialRing(QQ,'x');x=R.gen();a=target[0];a4,a6=E.a4(),E.a6()
                    poly=x**4-4*a*x**3-2*a4*x*x-(8*a6+4*a4*a)*x+a4*a4-4*a6*a
                    check(poly==R(result['history'][-1]['halving']['polynomial']),'duplication polynomial')
                    den=lcm([c.denominator() for c in poly]);integer=[ZZ(c*den) for c in poly];g=gcd(integer)
                    integer=[int(c/g) for c in integer]
                    obstruction=None
                    for p0 in prime_range(3,258):
                        p=int(p0)
                        if integer[-1]%p==0:continue
                        values=[sum(c*pow(z,i,p) for i,c in enumerate(integer))%p for z in range(p)]
                        if all(values):
                            obstruction=dict(prime=p,primitive_coefficients=list(map(str,integer)),
                                residues=[v%p for v in integer],values=values,no_infinity_root=True)
                            break
                    check(obstruction is not None,'UNKNOWN_INDEPENDENT_PROOF_PRIME_CAP')
                    record['rational_nonhalving_obstruction']=obstruction
            else:
                check(result['status']=='UNKNOWN_STEP_CAP','unknown status');record['unresolved']=True
            rows.append(record)
    check(len(rows)==38,'fixed five-fibre/38-decision roster')
    return dict(schema='independent-halving-cycle-replay-v1',status='PASS_INDEPENDENT_HALVING_CYCLE_REPLAY',
        cases=rows,case_count=len(rows),proof_prime_cap=257,inputs=inputs,checker_sha256=sha(Path(__file__)),
        candidate_generation_boundary='Candidate incidence is checked by the generic-equation producer and existing cover proofs; this independent replay checks rational curve membership, subgroup independence/dependence, and all doubling identities. No later point or search input.',
        argument='Generic mod2 injection proves2-saturation in the rational span. Exact odd relations certify inherited cases. A finite escape or unique-parity rational-nonhalving obstruction certifies independence; earlier doubling identities transport it to the original candidate.',
        point_searches=0,rank_upper_bounds=0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=replay();text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if a.write and not OUT.exists():
        with OUT.open('x') as f:f.write(text)
    else:check(OUT.read_text()==text,'immutable replay mismatch')
    print(result['status'],result['case_count'],'decisions',flush=True)
