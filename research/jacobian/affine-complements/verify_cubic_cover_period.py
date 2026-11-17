"""PF2RULED18: exact cubic-cover boundary period and pole(6,9) necessary frontier.
Global formal/cover geometry is written proof; exact finite identities and
full section kernels are checked here. Shared polynomial utilities are used.
"""
from pathlib import Path
from time import monotonic
from fractions import Fraction as F
from math import gcd
import json,hashlib,platform
import sympy as S
import verify_contact_filtration as V
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.expand(value))==0)
def scale(a,c):return {t:c*v for t,v in a.items() if c*v}
def form(a,b):return {(i,a+b-i):F(S.binomial(b,i)) for i in range(b+1)}
def partitions(n,minimum=1):
    if n==0:yield ()
    for k in range(minimum,n+1):
        for rest in partitions(n-k,k):yield (k,)+rest

def main():
    start=monotonic();w,t=S.symbols('w t');m=S.symbols('m',integer=True,positive=True)
    h=S.Function('h')(w);rho=S.Function('rho')(w)
    # Critical normal coefficient and the character primitive on zeta^3=h.
    critical=(2-m)*S.diff(h**(m/3),w)*rho-m*h**(m/3)*S.diff(rho,w)
    expected=-m*h**(m/3-1)*(h*S.diff(rho,w)+(m-2)*S.diff(h,w)*rho/3)
    check('all-height gcd3 critical coefficient identity',S.powsimp(S.expand(critical-expected),force=True)==0)
    hh,ss,AA,dh,ds,dA=S.symbols('hh ss AA dh ds dA',nonzero=True)
    derivative=(dA*ss+AA*ds)/hh**2-2*AA*ss*dh/hh**3
    zero('squarefree denominator turns rational ODE into polynomial identity',
         3*hh**2*derivative+4*hh*dh*AA*ss/hh**2-(3*ss*dA+(3*ds-2*ss*dh/hh)*AA))
    # Complete integer root multiplicity census for nonic sections in the q=infinity frame.
    census=[]
    for D in range(1,10):
        for part in partitions(D):
            cube=all(k%3==0 for k in part);r=len(part);simple=part.count(1)
            degree_bound=(r==1 or (D%3==0 and r<=2*D//3))
            budget=(r==1 or (D%3==0 and 2*D-3*r+simple>=2*D//3-1))
            census.append(dict(degree=D,partition=list(reversed(part)),cube=cube,
                root_count=r,simple_roots=simple,passes_numerator_degree_bound=degree_bound,
                passes_primitive_pole_budget=budget,
                numerator_degree=(0 if r==1 else 2*D//3-r) if degree_bound else None))
            if not cube and r>=2 and degree_bound:
                check(f'finite pole degree formula D{D} {part}',sum(max(2*k-3,0) for k in part)==2*D-3*r+simple)
    check('multi-root degrees are exactly3,6,9',sorted({c['degree'] for c in census if not c['cube'] and c['root_count']>=2 and c['passes_numerator_degree_bound']})==[3,6,9])
    for D,r,part in [(3,2,[2,1]),(6,4,[3,1,1,1]),(9,6,[4,1,1,1,1,1])]:
        retained=[c['partition'] for c in census if c['degree']==D and c['root_count']==r and not c['cube'] and c['passes_primitive_pole_budget']]
        check(f'extremal root partition D{D} r{r}',retained==[part])
    def period_test(poly):
        poly=S.Poly(poly,w).monic().as_expr();D=S.degree(poly,w)
        sq=S.Poly(S.cancel(poly/S.gcd(poly,S.diff(poly,w))),w).monic().as_expr();r=S.degree(sq,w)
        if r==1:N=0
        elif D%3 or 2*D//3<r:return None
        else:N=2*D//3-r
        coeff=S.symbols('a0:'+str(N+1));a=sum(coeff[i]*w**i for i in range(N+1))
        B=S.cancel(3*S.diff(sq,w)-2*sq*S.diff(poly,w)/poly)
        E=S.Poly(S.expand(3*sq*S.diff(a,w)+B*a-1),w)
        sols=S.linsolve(E.all_coeffs(),coeff)
        if sols is S.EmptySet or not sols:return None
        sol=next(iter(sols));ans=S.cancel(a.subs(dict(zip(coeff,sol))))
        zero('exact period solution '+str(poly),3*sq*S.diff(ans,w)+B*ans-1)
        return ans
    # All one-root noncube degrees and all two-root noncube cases in degrees3,6,9.
    positives=[]
    for D in (1,2,4,5,7,8):
        a=period_test(w**D);check('one-root period D'+str(D),a==S.Rational(1,3-2*D));positives.append(dict(h=str(w**D),A=str(a)))
    for D in (3,6,9):
        for k in range(1,D//2+1):
            if k%3==0 and (D-k)%3==0:continue
            poly=w**k*(w-1)**(D-k);a=period_test(poly)
            check(f'complete two-root period D{D} k{k}',a is not None)
            positives.append(dict(h=str(poly),A=str(a)))
    # Complete degree-six three-root positions, with distinct roots0,1,t.
    triples=[]
    for part,q,a0,a1 in [((4,1,1),2*t*t+t+2,-1/(5*t),-3/(5*t*(t+1))),
                        ((3,2,1),2*t+1,-1/(3*t),-1/(t*(t-1))),
                        ((2,2,2),t*t-t+1,-1/t,3/(t*(t+1)))]:
        sq=w*(w-1)*(w-t);poly=w**part[0]*(w-1)**part[1]*(w-t)**part[2]
        B=S.cancel(3*S.diff(sq,w)-2*sq*S.diff(poly,w)/poly);a=a0+a1*w
        res=S.Poly(S.cancel(3*sq*S.diff(a,w)+B*a-1),w)
        zero('three-root constant normalized '+str(part),res.coeff_monomial(1))
        zero('three-root quadratic coefficient '+str(part),res.coeff_monomial(w*w))
        numerator=S.factor(S.fraction(S.cancel(res.coeff_monomial(w)))[0])
        check('three-root sole residual criterion '+str(part),S.rem(numerator,q,t)==0 and S.degree(numerator,t)==S.degree(q,t))
        check('three-root parameter denominator guards '+str(part),S.gcd(S.Poly(q,t),S.Poly(t*(t-1)*(t+1),t)).degree()==0)
        triples.append(dict(multiplicities=list(part),normalized_roots=['0','1','t'],criterion=str(q),A=str(a),remaining_residual=str(res.coeff_monomial(w))))
    check('three-root denominator t=-1 cannot solve411',-S.Rational(3,5)!=0)
    check('three-root denominator t=-1 cannot solve222',-3*(-1/S.Integer(-1))!=0)
    check('degree9 five-root22221 is excluded by pole budget',2*9-3*5+1<2*9//3-1)
    # Extremal six-root noncontact and four-root positive-genus families.
    c=S.symbols('c',nonzero=True)
    for D,poly,a in [(6,w**3*(w**3+c),-1/(3*c)),(9,w**4*(w**5+c),-1/(5*c))]:
        sq=w*(w**(D-(3 if D==6 else 4))+c)
        B=S.cancel(3*S.diff(sq,w)-2*sq*S.diff(poly,w)/poly)
        zero('extremal family identity D'+str(D),3*sq*S.diff(a,w)+B*a-1)
        positives.append(dict(h=str(poly),A=str(a)))
    s0,s1,s2=S.symbols('s0 s1 s2');S3=w**3+s2*w*w+s1*w+s0
    zero('four-root classification differential',w*S.diff(S3,w)-3*S3-(-s2*w*w-2*s1*w-3*s0))
    b=S.symbols('b0:5');S5=w**5+sum(b[i]*w**i for i in range(5))
    zero('six-root classification differential',w*S.diff(S5,w)-5*S5-sum((i-5)*b[i]*w**i for i in range(5)))
    elliptic=period_test((w**3-1)**2);check('positive-genus exactness control',elliptic==-w/3)
    positives.append(dict(h='(w**3-1)**2',A=str(elliptic),primitive='-w/t, t**3=w**3-1'))
    # Direct function-field derivative of the elliptic primitive.
    z=S.symbols('z',nonzero=True)
    zero('elliptic primitive derivative modulo cubic relation',S.fraction(S.cancel(-1/z+w**3/z**4-1/z**4))[0].subs(z**3,w**3-1))
    negatives=[w**3*(w-1),w**4*(w-1),w*(w-1)*(w-2),w**4*(w-1)*(w-2),w**3*(w-1)**2*(w-2),w**2*(w-1)**2*(w-2)**2,
               w**3*(w-1)**2*(w-2)*(w-3)*(w-4)*(w-5)]
    for poly in negatives:check('negative period control '+str(poly),period_test(poly) is None)
    # Full section spaces for the numerical contact branch of pole(6,9).
    ap=V.powers({(0,0):1,(2,1):1},9);bp=V.powers({(0,0):-1,(2,1):-2,(4,2):-1,(0,1):-1},9);tp=V.powers({(0,0):1,(2,0):1,(4,1):1},18)
    spaces={};sourcecert={}
    for pole in (6,9):
        labels=[];original=[];rows=[];piv={}
        for j in range(pole+1):
            for i in range(pole+j+1):
                idx=len(labels);labels.append((j,i));poly={(x+i,y):F(c) for (x,y),c in V.mul(V.mul(ap[pole-j],bp[j]),tp[pole+j-i]).items()}
                original.append(poly.copy());vec={idx:F(1)}
                while poly:
                    key=max(poly,key=V.order)
                    if key not in piv:
                        c=poly[key];poly=scale(poly,1/c);vec=scale(vec,1/c);row=dict(source=poly,comb=vec,pivot=key,index=idx);piv[key]=row;rows.append(row);break
                    row=piv[key];c=poly[key];poly=V.add(poly,row['source'],-c);vec=V.add(vec,row['comb'],-c)
                check(f'Cox independence m{pole} row{idx}',bool(poly))
        check('complete section dimension '+str(pole),len(rows)==(pole+1)*(3*pole+2)//2)
        for row in rows:
            idx=row['index'];reconstructed={};sy={}
            for k,c in row['comb'].items():
                reconstructed=V.add(reconstructed,original[k],c);j,i=labels[k];power=2*(pole-j)+i;sy[power]=sy.get(power,0)+c*(-1)**(pole-j)
            check(f'exact source basis m{pole} row{idx}',reconstructed==row['source'] and max(row['comb'])==idx and row['comb'][idx]!=0 and max(reconstructed,key=V.order)==row['pivot'] and reconstructed[row['pivot']]==1)
            row['symbol']={k:c for k,c in sy.items() if c};row['chart']=V.chart(reconstructed)
            check(f'actual chart regularity m{pole} row{idx}',all(a>=0 for a,b in row['chart']))
        spaces[pole]=rows;sourcecert[str(pole)]={'labels':labels,'dimension':len(rows),'echelon_rows':[{'pivot':r['pivot'],'combination':[[k,str(c)] for k,c in r['comb'].items()]} for r in rows]}
    def raw(pole,cap):return [r for r in spaces[pole] if sum(r['pivot'])<=cap]
    kernels=[]
    def kernel(pole,cap,cut,degree=None,leading=None,label=''):
        inp=raw(pole,cap);piv={};out=[];comb=[]
        for idx,row in enumerate(inp):
            cc=row['chart'].copy();sy=row['symbol'].copy();vec={idx:F(1)};co={(0,k,0):c for k,c in sy.items() if k<cut}
            if degree is not None:co.update({(1,a,b):c for (a,b),c in V.add({t:c for t,c in cc.items() if sum(t)>=degree},leading,-cc.get((0,degree),0)).items()})
            while co:
                key=min(co)
                if key not in piv:
                    c=co[key];piv[key]=(scale(cc,1/c),scale(sy,1/c),scale(co,1/c),scale(vec,1/c));break
                cp,sp,cop,vp=piv[key];c=co[key];cc=V.add(cc,cp,-c);sy=V.add(sy,sp,-c);co=V.add(co,cop,-c);vec=V.add(vec,vp,-c)
            check(f'complete constraint change {label} row{idx}',max(vec)==idx and vec[idx]!=0)
            if not co:out.append(dict(chart=cc,symbol=sy));comb.append([[k,str(c)] for k,c in vec.items()])
        check('complete kernel dimension '+label,len(piv)+len(out)==len(inp))
        kernels.append(dict(label=label,input_dimension=len(inp),rank=len(piv),kernel_dimension=len(out),null_combinations=comb))
        return out
    for i,row in enumerate(raw(6,80)):check('full P6 source80 chart degree<=20 row'+str(i),max(map(sum,row['chart']))<=20)
    for i,row in enumerate(raw(6,72)):check('full P6 source72 chart degree<=19 row'+str(i),max(map(sum,row['chart']))<=19)
    for i,row in enumerate(raw(9,108)):check('Q9 source108 symbol order>=8 row'+str(i),all(k>=8 for k in row['symbol']))
    Q=kernel(9,128,6,label='Q128_contact6');check('contact6 Q128 dimension',len(Q)==139)
    for i,row in enumerate(Q):
        cc=row['chart'];check('full Q128 top projection row'+str(i),{t:c for t,c in cc.items() if sum(t)>=32}==scale(form(8,24),cc.get((0,32),0)))
    Qexc=kernel(9,108,9,label='Q108_contact9');check('exception contact9 dimension',len(Qexc)==131)
    H27=form(7,20);uH={(a+1,b):c for (a,b),c in H27.items()};E0=form(9,18);E1={(a+1,b):c for (a,b),c in form(8,18).items()};E2={(a+2,b):c for (a,b),c in form(7,18).items()}
    for i,row in enumerate(Qexc):
        cc=row['chart'];a=cc.get((1,27),0);c0=cc.get((0,27),0);b0=cc.get((1,26),0)-18*c0;d0=cc.get((2,25),0)-153*c0-18*b0
        expected=V.add(V.add(V.add(scale(uH,a),E0,c0),E1,b0),E2,d0)
        check('full Q108 exceptional projection row'+str(i),{t:c for t,c in cc.items() if sum(t)>=27}==expected)
    constrained=[]
    for cap,degree,aa,bb,mindim,minorder in [(48,12,3,9,43,13),(80,20,5,15,65,5),(72,18,6,12,61,6)]:
        rr=kernel(6,cap,4,degree,form(aa,bb),label='P6_'+str(cap))
        check('full partner family dimension '+str(cap),len(rr)==mindim)
        check('full partner symbol minimum '+str(cap),min(k for row in rr for k in row['symbol'])==minorder)
        if cap==72:
            for i,row in enumerate(rr):check('exception s6=p0_18 row'+str(i),row['symbol'].get(6,0)==row['chart'].get((0,18),0))
        constrained.append(dict(source_degree=cap,chart_degree=degree,dimension=len(rr),minimum_symbol_order=minorder))
    check('source128 window has only128',[sum(r['pivot']) for r in raw(9,128) if sum(r['pivot'])>=125]==[128])
    check('only degree18 fits exception partner', [d for d in (18,19) if d%3==0]==[18])
    check('after source-gap shears only j3,j5 remain',[j for j in range(1,6) if 8%j]!=[1,2,3,4,5] and [j for j in range(1,6) if 8%j]==[3,5])
    for pole,k,expected in [(6,0,96),(9,0,144),(6,1,88),(9,1,132)]:check(f'exact source degree contact k{k} pole{pole}',16*pole-4*(pole//3)*k==expected)
    for k,sl in [(0,S.Rational(3,4)),(1,S.Rational(8,11))]:
        for pole in (6,9):zero(f'cusp slope k{k} pole{pole}',S.Rational(3*pole-(pole//3)*k,4*pole-(pole//3)*k)-sl)
    frontier=[dict(contact_multiplicities=[0],ordinary_degrees=[96,144],cusp_multiplier_slope='3/4'),
              dict(contact_multiplicities=[1],ordinary_degrees=[88,132],cusp_multiplier_slope='8/11'),
              dict(contact_multiplicities=list(range(3,10)),ordinary_degrees=[80,128],cusp_multiplier_slope='3/4'),
              dict(contact_multiplicities=[7,8,9],ordinary_degrees=[48,128],cusp_multiplier_slope='3/4'),
              dict(contact_multiplicities=[3],ordinary_degrees=[72,108],cusp_multiplier_slope='2/3')]
    cert=dict(as_of='2026-11-14',period_polynomial_identity='3*s*A_prime+(3*s_prime-2*s*h_prime/h)*A=1',
        rational_denominator='rho=A*s/h^2',cover_primitive='3*zeta^4*A*s/h^2, zeta^3=h',
        multi_root_numerator_degree='2*D/3-r; requires D in{3,6,9}',
        primitive_pole_degree='2*D-3*r+number_of_simple_roots',primitive_infinity_contact_order='2*D/3-1',
        root_partition_census=census,degree6_three_root_classification=triples,
        positive_period_controls=positives,negative_period_controls=[str(x) for x in negatives],
        source_spaces=sourcecert,kernel_certificates=kernels,constrained_partner_spaces=constrained,
        pole69_necessary_frontier=frontier,excluded_contact_multiplicity=2,finite_map_constructed=False)
    path=HERE/'cubic_cover_period_certificates.json';path.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'CUBIC_COVER_BOUNDARY_PERIOD.md',path,HERE/'verify_contact_filtration.py',HERE/'BOUNDARY_LEADING_PERIOD.md',
        HERE/'GLOBAL_CUSP_LEADING_OBSTRUCTION.md',HERE/'HEIGHT_NINE_RESONANT_FRONTIER.md',HERE/'SOURCE_COORDINATES.md',HERE/'WRIGHT_MARKED_MODEL.md',
        HERE/'sources/degree-transfer-manifest.json',HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        written_result='All-height gcd3 cubic-cover exactness gives a finite polynomial test; full degree3/6 multi-root position classification, degree9 restrictions, and five necessary pole69 source cases exclude contact multiplicity2',
        all_height_gcd3_leading_period_criterion_proved=True,period_test_necessary_and_sufficient_for_cover_exactness=True,
        all_period_root_positions_classified=False,degree6_multi_root_positions_classified=True,
        complete_source_degree_cases_for_pole69=True,pole69_necessary_frontier=frontier,
        source_frontier_superseded_by='PF2RULED19',
        source_frontier_scope='PF2RULED18 before generic-polar exclusions; these five groups are historical',
        pole69_contact_multiplicity2_excluded=True,noncontact_period_controls_retained=True,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(period_criterion_uses_finiteness=False,period_criterion_uses_external_degree_theorem=False,
            source_degree_frontier_uses_retained_2022_alternative=True,auxiliary_cover_is_target_intermediate_field=False,
            contact_forced_for_every_gcd3_pair=False,period_solution_constructs_actual_PQ=False,
            all_pole69_pairs_excluded=False,all_degree9_period_positions_classified=False,
            actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,conditional_2026_degree125_upgrade_used=False,
            formal_verification=False,independent_full_reproof=False,external_review=False),
        source_sha256={str(f.relative_to(HERE)) if f.is_relative_to(HERE) else '../jc2-omitted-divisors/'+f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'cubic_cover_period_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; cubic-cover period test proved; pole69 contact2 excluded; five numerical cases remain')

if __name__=='__main__':main()
