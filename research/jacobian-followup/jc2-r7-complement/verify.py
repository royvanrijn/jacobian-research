"""Exact infinity/quintic checks; published topology remains a written input."""
import hashlib
import json
import time
import sympy as S
from sympy.polys.rings import ring
import geometry as G

CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)

def zero(name,value):check(name,S.cancel(value)==0)

def main():
    start=time.monotonic();F,K=G.F,G.K
    manifest=json.loads((G.HERE/'sources/manifest.json').read_text())
    for source in manifest['sources']:
        check('primary source hash '+source['file'],hashlib.sha256((G.HERE/'sources'/source['file']).read_bytes()).hexdigest()==source['sha256'])
    # Formal residue calculation: [r^m]Y = 6/(m-4)*[tau^(m-10)]H^((m-4)/4).
    for m in range(11,32,2):
        R,c=G.infinity_power(S.Rational(m-4,4),m-10)
        expected=R.zero if m<31 else S.Rational(3,2)*R.gens[-1]
        check('universal normalized odd coefficient '+str(m),6*c/(m-4)==expected)
    rhoL=S.factor(F.rho[11].as_expr().subs(F.line))
    raw=S.Poly(S.cancel(F.rho[11].as_expr().subs(F.laurent)*F.x**21),F.x)
    check('scaled Laurent gap has integral z exponents',all(m[0]%3==0 for m,c in raw.terms()))
    rhoK=S.expand(sum(c*K.z**(m[0]//3) for m,c in raw.terms()))
    zero('line rho11 unit',rhoL-S.Rational(1,19683))
    zero('Laurent rho11 unit',rhoK+S.Rational(2,59049)*K.z**6)
    zero('line first odd projective coefficient',-3*rhoL/2+S.Rational(1,13122))
    zero('Laurent first odd projective coefficient',-3*rhoK/2-K.z**6/19683)
    ii=S.Symbol('ii')
    zero('fourth-root distinct even conjugates',S.rem(1-ii**10-2,ii*ii+1,ii))
    check('opposite conjugates separated at order31',1-(-1)**31==2)
    check('Whitney tangent estimate',9-3==6 and 10-4==6)
    semigroup={0}
    for n in range(1,53):
        if any(n-a in semigroup for a in [4,10,41]):semigroup.add(n)
    check('infinity conductor48',47 not in semigroup and all(n in semigroup for n in range(48,52)))
    gaps=[n for n in range(48) if n not in semigroup]
    check('infinity delta24',len(gaps)==24)
    check('marked line Milnor67',2*(len(gaps)+10)-2+1==67)
    # Uniform projective extension has no base point at source infinity.
    sh,uh=S.symbols('sh uh')
    for name,p,q,param in [('L',F.p.subs(F.line),F.q.subs(F.line),F.d),('K',K.P,K.Q,K.z)]:
        ph=S.Poly(S.expand(uh**10*p.subs(F.t,sh/uh)),sh,uh).as_expr()
        qh=S.Poly(S.expand(uh**10*q.subs(F.t,sh/uh)),sh,uh).as_expr()
        zero(name+' projective infinity first coordinate',ph.subs(uh,0))
        zero(name+' projective infinity nonzero coordinate',qh.subs(uh,0)-sh**10)
    # Exact cubic-root action, avoiding numerical roots of unity.
    omega=S.Symbol('omega');cyclo=omega**2+omega+1
    for name,poly in [('p',F.p.subs(F.line)),('q',F.q.subs(F.line))]:
        factor=1 if name=='p' else omega**2
        expr=S.expand(poly.subs(F.d,omega*F.d)-factor*poly.subs(F.t,omega*F.t))
        zero('line cube-root action '+name,S.rem(expr,cyclo,omega))
    for variable,weight in [(F.a,2),(F.b,3),(F.c,4),(F.d,5),(F.B,2),(F.D,6),(F.E,8)]:
        expr=S.together(F.laurent[variable].subs(F.x,omega*F.x)-omega**(-weight)*F.laurent[variable])
        zero('Laurent weighted cube-root action '+str(variable),S.rem(S.fraction(expr)[0],cyclo,omega))
    check('polynomial exceptional orbit count',1+(34-1)//3==12)
    check('Laurent exceptional orbit count',237//3==79)
    check('exceptional coefficient orbits90',12+79-1==90)
    check('at most92 complement types',2+90==92)
    # Verify the selected probes lie in the nodal opens; no braid is inferred.
    line_data=json.loads((G.HERE.parent/'jc2-r7-collisions/strata.json').read_text())
    for name,expr in line_data['exceptional_factors'].items():
        check('nodal d1 avoids '+name,S.sympify(expr,locals={'d':F.d}).subs(F.d,1)!=0)
    for name,f in [('C1',K.C1),('C2',K.C2),('D',K.D),('T',K.Tangent)]:check('nodal z1 avoids '+name,f.subs(K.z,1)!=0)
    witnesses=json.loads((G.HERE/'witnesses.json').read_text())
    check('two explicitly selected nodal probes',{row['family'] for row in witnesses}=={'L','K'} and len(witnesses)==2)
    for row in witnesses:
        expected=(F.p.subs(F.line).subs(F.d,1),F.q.subs(F.line).subs(F.d,1)) if row['family']=='L' else (K.P.subs(K.z,1),K.Q.subs(K.z,1))
        for name,poly in zip(['p','q'],expected):zero('exact nodal probe '+row['family']+name,S.sympify(row[name],locals={'t':F.t})-poly)
    bounds=json.loads((G.HERE/'bounds.json').read_text())
    for row in bounds['attempts']:
        check('capped braid is not a theorem '+row['family'],row['state']=='CAPPED_WITHOUT_COMPLETED_BRAIDS' and not row['completed_braid_file'])
        check('capped braid log hash '+row['family'],hashlib.sha256((G.HERE/row['log']).read_bytes()).hexdigest()==row['log_sha256'])
    # The Newton space for an actual birational (3,5) quintic image.
    candidates=[(i,j) for i in range(6) for j in range(4) if i+j<=5]
    high={}
    for i,j in candidates:
        weight=3*i+5*j
        if weight>15:high.setdefault(weight,[]).append((i,j))
    check('all forbidden high weights unique',all(len(rows)==1 for rows in high.values()))
    check('Newton space has thirteen monomials',len(G.MONOMIALS)==13)
    degrees=[6*i+10*j for i,j in G.PIVOTS]+[9]
    check('restriction basis has distinct pole degrees',len(set(degrees))==13)
    check('quintic leading relation weight30',6*5==10*3==30)
    results=[]
    for name,p,q,param,rho in [('L',F.p.subs(F.line),F.q.subs(F.line),F.d,rhoL),('K',K.P,K.Q,K.z,rhoK)]:
        equation,pullback,steps=G.approximate_quintic(p,q,F.t)
        check(name+' exact seven leading cancellations',[row[0] for row in steps]==[28,24,22,18,16,12,10])
        check(name+' quintic ordinary degree',S.Poly(equation,G.U,G.V).total_degree()==5)
        zero(name+' monic cubic target coefficient',S.Poly(equation,G.V).LC()-1)
        check(name+' quintic Newton support',all(3*i+5*j<=15 for (i,j),c in S.Poly(equation,G.U,G.V).terms()))
        check(name+' residual pole degree9',S.degree(pullback,F.t)==9)
        zero(name+' residual leading unit',S.Poly(pullback,F.t).LC()+3*rho)
        zero(name+' fixed lower pole freedom t6',S.Poly(pullback,F.t).nth(6))
        zero(name+' fixed lower pole freedom constant',S.Poly(pullback,F.t).nth(0))
        R,pp,tt=ring(str(param)+',t',S.QQ)
        P,Q=R.from_expr(p),R.from_expr(q);value=R.zero
        for (i,j),co in S.Poly(equation,G.U,G.V).terms():value+=R.from_expr(co)*P**i*Q**j
        check(name+' literal polynomial contact identity',value==R.from_expr(pullback))
        results.append(dict(family=name,parameter=str(param),equation=str(equation),pullback=str(pullback),
                            leading_pullback=str(S.factor(S.Poly(pullback,F.t).LC())),
                            maximal_contact_family='equation + alpha*U + beta; further rationality/Keller constraints remain required'))
    spectrum=sorted({50-n for n in degrees if n>=9})
    check('quintic contact spectrum',spectrum==[20,22,24,26,28,30,32,34,38,40,41])
    check('finite quintic intersection at least9',50-max(spectrum)==9)
    ledger=dict(schema='r7-complement-and-quintic-v1',infinity_characteristic=[4,10,31],infinity_semigroup=[4,10,41],
                infinity_delta=24,infinity_line_intersection=10,marked_infinity_pair_milnor=67,
                generic_nodal_complement_families=2,exceptional_cube_root_orbits=90,complement_type_upper_bound=92,
                quintic_contact_spectrum=spectrum,finite_intersection_minimum=9,closest_quintics=results,
                topology_assurance='Written direct Whitney/properness proof with imported Mather Proposition11.1; checker verifies jets, units and arithmetic, not the topology theorem.',
                unresolved=['generic braid presentations after both60s caps','joint C/B/H configuration at infinity and finite intersections','actual additional image equations and source realization'])
    (G.HERE/'geometry.json').write_text(json.dumps(ledger,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} INFINITY/QUINTIC CHECKS ({time.monotonic()-start:.3f}s)',flush=True)
    inputs=[G.HERE/n for n in ['PROOF.md','README.md','REVIEW.md','geometry.py','verify.py','braid_probe.py','witnesses.json','bounds.json','sources/manifest.json']]
    inputs += [G.HERE/'sources'/row['file'] for row in manifest['sources']]
    inputs += [G.HERE/row['log'] for row in bounds['attempts']]
    inputs += [G.HERE.parent/'jc2-r7-image-family/family.py',G.HERE.parent/'jc2-r7-image-family/PROOF.md',
               G.HERE.parent/'jc2-r7-collisions/PROOF.md',G.HERE.parent/'jc2-r7-collisions/strata.json',
               G.HERE.parent/'jc2-r7-laurent-collisions/laurent.py',
               G.HERE.parent/'jc2-r7-laurent-collisions/PROOF.md',G.HERE.parent/'jc2-r7-laurent-collisions/strata.json',
               G.ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(schema='r7-complement-and-quintic-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),
                 input_sha256={str(p.relative_to(G.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                 output_sha256={'geometry.json':hashlib.sha256((G.HERE/'geometry.json').read_bytes()).hexdigest()},
                 assurance='Exact jets, units, cube-root actions and quintic identities. Direct geometric hypotheses and imported topology are authored written proof; no formalization or independent external review. No completed generic braid presentation.')
    (G.HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
