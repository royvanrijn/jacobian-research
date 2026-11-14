"""Reproduce the finite polarization attempts, without claiming minimality."""
from __future__ import annotations
from itertools import product
from fractions import Fraction
import contextlib
import io
import json
import sympy as s
from verify import X,ROOT,cubic12,homogeneous18
from polarize import terms,build,search


def main():
    masks=[]
    for mask in range(512):
        ts=terms(True,mask=mask)
        masks.append([mask,len(build(ts,[0]*len(ts))[0])])
    assert min(n for _,n in masks)==90
    _,k,_=cubic12(2);h,_=homogeneous18(k)
    basis=[]
    for bits in product([0,1],repeat=4):
        sub={X[6]:3*X[1]-X[6] if bits[0] else X[6],
             X[9]:X[9]-X[1] if bits[1] else X[9],
             X[8]:(X[8]+3*X[7])/2 if bits[2] else X[8]}
        sub[X[10]]=3*sub[X[9]]-X[10] if bits[3] else X[10]
        ht=[s.expand(f.subs(sub,simultaneous=True)) for f in h]
        vs=[tuple(int(i==j) for j in range(18)) for i in range(18)]
        ts=[]
        for row,f in enumerate(ht):
            for mon,co in s.Poly(f,*X).terms():
                if co:ts.append((row,Fraction(co),[v for v,e in zip(vs,mon) for _ in range(e)]))
        ts.append((17,Fraction(1),[vs[17]]*3))
        with contextlib.redirect_stdout(io.StringIO()):
            n,_,_,_=search(ts,restarts=5)
        basis.append({'bits':bits,'size':n})
    result={'grouping_scope':'all 512 subsets of nine published regroupings, fixed seven-cube identity',
            'grouping_results':masks,'grouping_minimum':90,
            'basis_scope':'16 specified invertible triangular rational changes; five seeded greedy restarts each',
            'basis_results':basis,'basis_minimum':min(r['size'] for r in basis),
            'boundary':'No global lower bound, no exhaustive search over decompositions or coordinate changes.'}
    (ROOT/'bounded_search.json').write_text(json.dumps(result,indent=2)+'\n')
    print('BOUNDED SEARCH: grouping minimum 90; basis minimum',result['basis_minimum'])


if __name__=='__main__':main()
