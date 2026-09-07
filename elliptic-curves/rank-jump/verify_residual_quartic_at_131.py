#!/usr/bin/env python3
"""Independent exact replay using SymPy and integer Legendre symbols."""
import argparse
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_residual_quartic_incidence_inputs_v1.json'
SOURCE=r.OUT/'rank_jump_residual_quartic_at_131_v1.json'
OUTPUT=r.OUT/'rank_jump_residual_quartic_at_131_verification_v1.json'


def compute():
    import sympy as s
    data=r.read(INPUT);result=r.read(SOURCE);p=131
    for path,sha in result['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    t=s.Symbol('t')
    def poly(cs,modulus=None):
        return s.Poly(sum(s.Rational(c)*t**i for i,c in enumerate(cs)),t,modulus=modulus)
    gram=s.eye(17)*4;counts=0
    for modulus,key in [(None,'intersections_Q'),(p,'intersections_F131')]:
        A=poly(data['A'],modulus);B=poly(data['B'],modulus)
        D=-4*A**3-27*B**2
        assert D.degree()==24 and s.gcd(D,D.diff()).degree()==0 and s.gcd(A,D).degree()==0
        pts=[(poly(z['x'],modulus),poly(z['y'],modulus)) for z in data['sections']]
        for x,y in pts:assert y*y==x**3+A*x+B
        pairs=[];h=s.eye(17)*4
        for i,(x,y) in enumerate(pts):
            for j in range(i):
                dx=x-pts[j][0];dy=y-pts[j][1]
                finite=int(s.gcd(dx,dy).degree());infinity=int(min(4-dx.degree(),6-dy.degree()))
                pairs.append([i,j,finite,infinity]);h[i,j]=h[j,i]=2-finite-infinity;counts+=1
        assert pairs==result[key] and h.tolist()==result['gram']
        if modulus is None:gram=h
        else:assert h==gram
    assert gram.det()==948
    def ev(cs,n):
        value=0
        for c in reversed(cs):
            q=s.Rational(c);value=(value*n+int(q.p)*pow(int(q.q),-1,p))%p
        return value
    rows=[];skipped=[];irreducible=[]
    for n in range(p):
        a=ev(data['A'],n);b=ev(data['B'],n)
        if (-4*a**3-27*b*b)%p==0 or any(ev(z['y'],n)==0 for z in data['sections']):
            skipped.append(n);continue
        roots=[z for z in range(p) if (z**3+a*z+b)%p==0]
        if not roots:irreducible.append(n)
        for theta in roots:
            values=[(ev(z['x'],n)-theta)%p for z in data['sections']]
            assert all(values)
            mask=sum((pow(v,65,p)==130)<<i for i,v in enumerate(values))
            rows.append({'t':n,'theta':theta,'values':values,'mask':mask})
    assert rows==result['kummer_characters'] and skipped==result['skipped_base_values']
    assert irreducible==result['irreducible_cubic_base_values'] and irreducible
    # Separate bit elimination, not the producer's rank helper.
    pivots={};witness=[]
    for i,row in enumerate(rows):
        mask=row['mask']
        while mask:
            k=mask.bit_length()-1
            if k in pivots:mask^=pivots[k]
            else:pivots[k]=mask;witness.append(i);break
    assert len(pivots)==17 and witness==result['independent_character_row_indices']
    moments=r.read(r.OUT/'rank_jump_native_twist_frobenius_v1.json')
    traces=[]
    for n in (1,2):
        ledger=[z for z in moments['fibre_trace_ledger'] if z[0]==n]
        assert sum(z[3] for z in ledger)==p**n+1
        traces.append(-sum(z[3]*z[4] for z in ledger))
    assert traces==result['Frobenius_traces']==[1884,319520]
    X=s.Symbol('X');residual=(X+131)*(X*X+212*X+17161)
    coeff=s.Poly(residual,X).all_coeffs()
    assert -coeff[1]==traces[0]-17*p and coeff[1]**2-2*coeff[2]==traces[1]-17*p*p
    assert residual.subs(X,p)!=0
    assert s.Rational(residual.subs(X,p),p**3)==s.Rational(948,p)
    # A four-letter permutation has fixed-module dimension equal to its orbit count.
    from itertools import permutations
    cycle_types=set()
    for perm in permutations(range(4)):
        seen=set();cycles=[]
        for i in range(4):
            if i in seen:continue
            size=0;j=i
            while j not in seen:seen.add(j);size+=1;j=perm[j]
            cycles.append(size)
        if 16+len(cycles)==17:cycle_types.add(tuple(sorted(cycles)))
    assert cycle_types=={(4,)}
    return {'schema':'rank-jump.residual-quartic-at-131-verification.v1','status':'PASS',
        'height_intersections_recomputed':counts,'characters_recomputed':len(rows),'character_rank':17,
        'determinant':948,'original_fibre_ledger_orbits_reassembled':len(moments['fibre_trace_ledger']),
        'residual_polynomial':str(s.expand(residual)),'Artin_Tate_factor':'948/131',
        'quartic_cycle_types_if_two_torsion_dimension_17':[[4]],
        'bindings':{str(path.relative_to(r.ROOT)):r.digest(path.read_bytes()) for path in (Path(__file__),INPUT,SOURCE)},
        'boundary':'Arithmetic replay only. Selmer/Jacobian identification and Artin-Tate implications are proved in the companion note; no full point recount or exceptional coordinates.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);args=parser.parse_args()
    result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['characters_recomputed'],'characters; determinant 948; residual Frobenius polynomial')
