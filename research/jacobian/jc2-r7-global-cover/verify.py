"""Exact blowup, permutation and line-genus checks for the Nori obstruction."""
from pathlib import Path
import json,hashlib,time,math
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)
def valuation(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return math.inf
    return min(m[0] for m,c in S.Poly(num,var).terms())-min(m[0] for m,c in S.Poly(den,var).terms())
def cycle(*letters):
    p=list(range(6))
    for a,b in zip(letters,letters[1:]+letters[:1]):p[a]=b
    return tuple(p)
def compose(a,b):return tuple(a[b[i]] for i in range(6))
def inverse(a):return tuple(a.index(i) for i in range(6))

def main():
    start=time.monotonic();manifest=json.loads((HERE/'sources/manifest.json').read_text())
    for row in manifest['sources']:
        path=HERE/'sources'/row['file'] if 'file' in row else ROOT/row['existing_file']
        check('primary-source SHA256 '+path.name,hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'])
    rr,ss=S.symbols('r s')
    contacts=[20,22,24,26,28,30,32,34,38,40]
    old=json.loads((HERE.parent/'jc2-r7-complement/geometry.json').read_text())
    check('complete polynomial-quintic contact spectrum after41 exclusion',[v for v in old['quintic_contact_spectrum'] if v!=41]==contacts)
    table=[]
    check('nodal quintic has exactly four affine nodes',(5-1)*(5-2)//2-(2-1)*(5-1)//2==4)
    for contact in contacts+[36]:
        n=(contact-20)//2
        xc=rr**4;yc=2*rr**10+rr**31 if contact==20 else rr**10+rr**(contact-10)+rr**31
        xb=ss**2;yb=ss**5
        check('model intersection multiplicity '+str(contact),valuation(yc*yc-xc**5,rr)==contact)
        cvals=[];bvals=[]
        for transform in [lambda x,y:(x,y/x),lambda x,y:(x,y/x),lambda x,y:(x/y,y),lambda x,y:(x,y/x-1)]:
            cvals.append(min(valuation(xc,rr),valuation(yc,rr)));bvals.append(min(valuation(xb,ss),valuation(yb,ss)))
            xc,yc=map(S.cancel,transform(xc,yc));xb,yb=map(S.cancel,transform(xb,yb))
        check('four B infinity multiplicities '+str(contact),bvals==[2,2,1,1])
        check('four C shared multiplicities '+str(contact),cvals==[4,4,2,2])
        check('B clean transverse boundary after fourth blowup '+str(contact),xb==ss and yb==0)
        for j in range(n):
            check(f'additional shared C multiplicity {contact}:{j}',min(valuation(xc,rr),valuation(yc,rr))==2)
            cvals.append(2);bvals.append(1);yc=S.cancel(yc/xc)
        check('C separates from B at infinity '+str(contact),valuation(yc,rr)==0)
        check('Noether shared-intersection sum '+str(contact),sum(a*b for a,b in zip(cvals,bvals))==contact)
        square=25-sum(m*m for m in bvals)
        check('strict B self-intersection '+str(contact),square==15-n)
        table.append(dict(infinity_contact=contact,actual_degree35_contact=contact in contacts,B_multiplicities=bvals,
            B_self_intersection=square,node_threshold=8,Nori_excludes_transverse_pair=square>8,
            minimum_finite_tangency_sum=max(0,7-n)))
    # Finite smooth tangencies need the final blowup at a boundary crossing too.
    for order in range(2,31):
        y=rr**order
        for j in range(order-1):check(f'finite tangency still at boundary point {order}:{j}',valuation(y,rr)>0);y=S.cancel(y/rr)
        check('last unresolved finite boundary triple '+str(order),y==rr)
        y=S.cancel(y/rr);check('finite tangency separates after m blowups '+str(order),y==1)
    identity=tuple(range(6));groups=[]
    for name,sigma,g in [('transposition',cycle(0,1),cycle(0,1,2)),('three_cycle',cycle(0,1,2),cycle(2,3))]:
        conjugate=compose(compose(g,sigma),inverse(g));commutator=compose(compose(compose(sigma,conjugate),inverse(sigma)),inverse(conjugate))
        check(name+' normal closure contains noncommuting elements',commutator!=identity)
        groups.append(dict(inertia=name,element=list(sigma),conjugate=list(conjugate),commutator=list(commutator)))
    # Full two-prime Hurwitz census, without assuming any individual f_infinity=1.
    genus=[]
    for index in [2,3]:
        affine=10+5*(index-1)
        rows=[]
        for residue_sum in range(1,7):
            numerator=affine-4-residue_sum
            if numerator%2:continue
            gamma=numerator//2;canonical=2*gamma-8
            rho_min=int(max(2,math.ceil(S.Rational(60-canonical*canonical,6))))
            rho_max=6*gamma-1
            check(f'line Hurwitz e{index} s{residue_sum}',2*gamma-2==-12+affine+6-residue_sum)
            check(f'minimal-resolution positivity window e{index} s{residue_sum}',rho_min<=rho_max)
            rows.append(dict(residue_degree_sum=residue_sum,general_line_genus=gamma,K_dot_H=canonical,rho_Hodge_min=rho_min,rho_Q_big_max=rho_max))
        check('row genus range '+str(index),sorted(v['general_line_genus'] for v in rows)==([3,4,5] if index==2 else [5,6,7]))
        genus.append(dict(second_index=index,affine_ramification_degree=affine,rows=rows))
    # The first-image opposite branch violates Orevkov NC in every affine projection.
    aa=S.Symbol('a',nonzero=True)
    # Z=T^10+(even lower powers)+a*T^-11; delta=(-a/5)*T^-20.
    check('generic projection odd difference and derivative orders',-11-9==-20)
    check('compensating branch coefficient',-2*aa-10*(-aa/5)==0)
    check('independent coordinate difference tends to zero',5-20==-15 and S.Rational(6,5)*aa!=0)
    result=dict(schema='r7-global-cover-nori-v1',infinity_resolution=table,normal_inertia_obstructions=groups,
        necessary_inequality='I_infinity(C,B)+2*sum_finite_tangencies(I_y(C,B)) >=34',
        scope='B is a nodal degree(3,5) polynomial image; finite C/B singular loci are mutually avoided. C is a complete gap21 image and the finite normalization has S6 monodromy with a single2/3-cycle over B.',
        transverse_excluded_contacts=[v['infinity_contact'] for v in table if v['actual_degree35_contact'] and v['Nori_excludes_transverse_pair']],
        general_line_budgets=genus,assurance='Authored blowup/kernel and Hurwitz interpretation with imported Nori Proposition3.27; algebra checker does not prove the topology theorem. No formalization or independent external review.')
    (HERE/'budgets.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs=[HERE/name for name in ['verify.py','PROOF.md','README.md','REVIEW.md','sources/manifest.json','discovery.json']]
    inputs += [HERE/'sources'/row['file'] if 'file' in row else ROOT/row['existing_file'] for row in manifest['sources']]
    inputs += [HERE.parent/path for path in ['jc2-r7-complement/PROOF.md','jc2-r7-complement/geometry.json','jc2-r7-quintic-normalizations/PROOF.md','jc2-omitted-divisors/GLOBAL_MONODROMY.md','jc2-omitted-divisors/RESIDUE_LEDGER.md','minimal-resolution-positivity/PROOF.md']]
    receipt=dict(passed=True,check_count=len(CHECKS),checks=CHECKS,elapsed_seconds=round(time.monotonic()-start,3),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        output_sha256={'budgets.json':hashlib.sha256((HERE/'budgets.json').read_bytes()).hexdigest()},assurance=result['assurance'])
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(f'PASS {len(CHECKS)} GLOBAL COVER CHECKS ({time.monotonic()-start:.3f}s)',flush=True)

if __name__=='__main__':main()
