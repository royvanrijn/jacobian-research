#!/usr/bin/env python3
from pathlib import Path
import exact_core as ec

ROOT=Path(__file__).resolve().parent
B,E,C,F,np=ec.high_layers()
assert np==3
# J1: solve G_1..G_12 in 2 A G' + (B F'-B'F)-2 C'E = 0
rows3=[[ec.K(0) for _ in range(12)] for _ in range(20)]
for j in range(1,13):
    for i,Ai in enumerate(ec.Avec):
        deg=i+j-1
        if 0 <= deg < 20:
            rows3[deg][j-1] += Ai*(2*j)
known3=ec.tadd(ec.tmul(B,ec.tder(F)),ec.tscale(ec.tmul(ec.tder(B),F),ec.K(-1)))
known3=ec.tadd(known3,ec.tscale(ec.tmul(ec.tder(C),E),ec.K(-2)))
keep3=[]; rhs3=[]
for deg,row in enumerate(rows3):
    kn=known3[deg] if deg<len(known3) else ec.PP()
    if any(row) or kn:
        keep3.append(row); rhs3.append(-kn)
sol3,free3,comp3,np2=ec.linear_solve(keep3,rhs3,np)
assert (len(free3),np2)==(0,3)
G=[ec.PP() for _ in range(13)]
for j in range(1,13): G[j]=sol3[j-1]
# J0 = B G' - C' F
J0=ec.tadd(ec.tmul(B,ec.tder(G)),ec.tscale(ec.tmul(ec.tder(C),F),ec.K(-1)))
final=[]; seen=set()
for q in comp3+[q for q in J0 if q]:
    key=q.sing(3)
    if key not in seen:
        seen.add(key); final.append(q)
print(f'comp3={len(comp3)} J0_nonzero={sum(bool(q) for q in J0)} unique={len(final)}')
assert len(final)==25
inds=(0,1,7,9)
# write exact number-field Singular computation
H=ec.K(ec.MOD).sing() # zero after reduction, do not use
# MOD is monic fmpq_poly. Format directly.
def mod_sing(poly):
    terms=[]
    for i in range(len(poly)):
        c=poly[i]
        if not c: continue
        num=int(c.p); den=int(c.q)
        cs=str(num) if den==1 else f'({num}/{den})'
        mon='' if i==0 else ('u' if i==1 else f'u^{i}')
        if not mon: term=cs
        elif c==1: term=mon
        elif c==-1: term='-'+mon
        else: term=f'{cs}*{mon}'
        terms.append(term)
    return ('+'.join(terms).replace('+-','-'))
minpoly=mod_sing(ec.MOD)
compact=[final[i].sing(3) for i in inds]
full=[q.sing(3) for q in final]
for name,eqs in [('case2_compact4_exact.sing',compact),('case2_full25_exact.sing',full)]:
    p=ROOT/name
    p.write_text(
        'LIB "resources.lib";\nResources::setcores(1);\nLIB "nfmodstd.lib";\n'
        +'ring R=(0,u),(r,s,h),dp;\n'
        +'minpoly='+minpoly+';\noption(redSB);\n'
        +'ideal I=\n'+',\n'.join(eqs)+';\n'
        +'print("INPUT_SIZE="+string(size(I)));\n'
        +'ideal J=nfmodStd(I);\n'
        +'print("SIZE="+string(size(J)));\nJ;\n'
        +'if(size(J)==1 && J[1]==1){print("UNIT");}else{print("NONUNIT");}\nquit;\n'
    )
    print(p, p.stat().st_size)
# Save machine-readable residuals
(ROOT/'case2_residuals_exact.txt').write_text('\n'.join(f'R{i}={q.sing(3)}' for i,q in enumerate(final)))
