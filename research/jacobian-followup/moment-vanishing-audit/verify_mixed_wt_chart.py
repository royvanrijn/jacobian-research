#!/usr/bin/env python3
"""Exact mixed-chart reduction, complete spine certificate and prefix trap.

The full parent remains UNKNOWN. The large generic-spine cofactor identity
is multiplied independently with Fraction arithmetic, without Singular.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from collections import defaultdict
import hashlib,json,gzip,re,time
import sympy as sp
from verify_moment_transfers import add,mul,scale,gaussian
from compile_mixed_wt_parent import moments as parent_moment,VARIABLES,compositions


def trunc(poly,n):return {mon:c for mon,c in poly.items() if mon[0]<=n}


def power(poly,k,n):
    out={(0,0):F(1)}
    for _ in range(k):out=trunc(mul(out,poly),n)
    return out


def phase_series(a,b,c,d,G,n):
    Y={(0,1):F(1)};Z={(1,0):a,(1,1):c}
    # Z=t*(a+cY).  The two corrections are -(b/c)*Z and -(d/c)*Z^2.
    Ft=add(Y,scale(Z,-b/c),scale(mul(Z,Z),-d/c))
    exponent=scale(add(mul(Y,Y),scale(mul(Ft,Ft),-1)),F(1,2))
    for (i,j),coeff in G.items():
        term=mul(power(Z,i,n),power(Ft,j,n))
        exponent=add(exponent,{(t+1,y):coeff*v for (t,y),v in term.items() if t+1<=n})
    exponent=trunc(exponent,n);assert all(t>=1 for t,y in exponent)
    out=[F(1)]+[F(0)]*n;pk={(0,0):F(1)}
    for k in range(1,n+1):
        pk=trunc(mul(pk,exponent),n)
        for (t,y),coef in pk.items():
            if y%2==0:out[t]+=coef*F(factorial(y),2**(y//2)*factorial(y//2)*factorial(k))
    return out


def direct(a,b,c,d,G,n):
    P=add({(0,1,0):a,(1,1,0):b,(0,1,1):c,(2,1,0):d},
          {(i,0,j):v for (i,j),v in G.items()})
    P={m:v for m,v in P.items() if v};pk={(0,0,0):F(1)};out=[]
    for m in range(1,n+1):pk=mul(pk,P);out.append(gaussian(pk))
    return out


def numeric_checks(parent):
    inputs=[(F(2),F(0),F(1),F(1),{(0,0):F(-1),(0,1):F(2),(0,2):F(-3),(0,3):F(1),(1,0):F(4),(1,1):F(-2),(1,2):F(3),(2,0):F(1),(2,1):F(-1),(3,0):F(2)}),
            (F(0),F(0),F(1),F(1),{(3,0):F(2),(2,1):F(-1),(1,2):F(3),(0,3):F(1)}),
            (F(2),F(3),F(1),F(1),{(3,0):F(1),(2,1):F(-2),(0,3):F(1),(1,0):F(2),(0,0):F(-1)}),
            (F(1),F(-2),F(2),F(0),{(2,1):F(1),(0,3):F(2),(0,1):F(-1)})]
    rows=[]
    names_to_powers={'r0':(0,0),'u0':(0,1),'v0':(0,2),'s':(0,3),'r1':(1,0),'u1':(1,1),'v1':(1,2),'r2':(2,0),'u2':(2,1),'r3':(3,0)}
    for a,b,c,d,G in inputs:
        values=direct(a,b,c,d,G,8);series=phase_series(a,b,c,d,G,8)
        assert all(values[m-1]==factorial(m)*series[m] for m in range(1,9))
        if b==0 and c==d==1:
            params=[a]+[G.get(names_to_powers[name],F(0)) for name in VARIABLES[1:]]
            for m in range(1,9):
                compiled=parent_moment(m);stored={tuple(mon):coef for mon,coef in parent['moments'][m-1]['terms']};assert compiled==stored
                value=sum((F(coef)*sp_product([v**e for v,e in zip(params,mon)]) for mon,coef in compiled.items()),F(0))
                assert value==values[m-1]
        rows.append({'a_b_c_d':list(map(str,(a,b,c,d))),'pure_moments_through_8':list(map(str,values))})
    return rows


def sp_product(values):
    out=F(1)
    for v in values:out*=v
    return out


def symbolic_spine():
    a,A,B,C,S=sp.symbols('a A B C S')
    terms=[(a,1,0,0),(1,1,2,0),(1,1,0,1),(A,0,3,0),(B,0,2,1),(C,0,1,2),(S,0,0,3)]
    out={}
    for m in range(1,9):
        poly=0
        for counts in compositions(m,7):
            w=sum(n*x[1] for n,x in zip(counts,terms));z=sum(n*x[2] for n,x in zip(counts,terms));t=sum(n*x[3] for n,x in zip(counts,terms))
            if w!=z or t%2:continue
            num=factorial(m)*factorial(w)*factorial(t)
            den=2**(t//2)*factorial(t//2)*sp_product([factorial(n) for n in counts]);assert num%den==0
            poly+=(num//den)*sp.prod(x[0]**n for n,x in zip(counts,terms))
        out[m]=sp.expand(poly)
    return (a,A,B,C,S),out


def fraction_certificate(cert):
    w,a,R=sp.symbols('w a R');result=defaultdict(F);products=0
    assert cert['variables']==['w','a','R'] and len(cert['generators'])==len(cert['cofactors'])==5
    assert sp.expand(sp.sympify(cert['generators'][4],locals={'w':w,'a':a,'R':R})-(w*a*(a*a+3)*R-1))==0
    for expr,rows in zip(cert['generators'],cert['cofactors']):
        assert re.fullmatch(r'[0-9aRw+\-*/^()\s]+',expr)
        poly=sp.Poly(sp.sympify(expr.replace('^','**'),locals={'w':w,'a':a,'R':R}),w,a,R,domain=sp.QQ)
        gen=[(mon,F(str(c))) for mon,c in poly.terms()]
        for mon,coefficient in rows:
            cf=F(coefficient)
            for gm,gc in gen:
                result[tuple(i+j for i,j in zip(mon,gm))]+=cf*gc;products+=1
    result={mon:c for mon,c in result.items() if c}
    assert result=={(0,0,0):F(1)},result
    return {'cofactor_monomials':sum(map(len,cert['cofactors'])),'rational_products':products,'identity':'1'}


def proof_checks(cert,syms,M):
    a,A,B,C,S=syms;R=sp.symbols('R');q4=5*C*C+6*C+4;q8=33*C**4+60*C**3+60*C*C+40*C+16
    assert sp.expand(M[2].subs(a,0)-15*S*S)==0
    assert sp.expand(M[4].subs({a:0,S:0})-36*q4)==0
    assert sp.expand(M[5].subs({a:0,S:0})-360*B*(5*C+4))==0
    assert sp.expand(M[6].subs({a:0,S:0})-5400*(2*A*(C+1)+B*B))==0
    assert sp.expand(M[8].subs({a:0,S:0,A:0,B:0})-529200*q8)==0
    assert sp.expand(5*q4-(5*C+2)*(5*C+4))==12
    assert sp.expand(q4-(5*C+1)*(C+1))==3
    assert sp.expand((6567*C**3+7452*C*C+3912*C+3472)*q4+(-995*C-514)*q8)==5664
    Cs=-2-15*S*S/(2*a);Bs=sp.factor(-S*(15*Cs+6)/(2*a))
    As=sp.factor(-M[4].subs({A:0,B:Bs,C:Cs})/(24*a*(a*a+3)))
    assert sp.cancel(M[2].subs(C,Cs))==0
    assert sp.cancel(M[3].subs({C:Cs,B:Bs}))==0
    assert sp.cancel(M[4].subs({A:As,B:Bs,C:Cs}))==0
    for m,stored in zip(range(5,9),cert['generators'][:4]):
        expr=sp.factor(M[m].subs({A:As,B:Bs,C:Cs}));num,den=sp.fraction(expr)
        factors=sp.factor_list(den)[1]
        assert {f for f,e in factors}<={a,a*a+3}
        if m%2:num=sp.cancel(num/S)
        poly=sp.Poly(num,S);assert all(mon[0]%2==0 for mon,c in poly.terms())
        f=sum(c*R**(mon[0]//2) for mon,c in poly.terms());f=sp.Poly(f,a,R).primitive()[1].as_expr()
        assert sp.expand(f-sp.sympify(stored,locals={'a':a,'R':R}))==0
    x=sp.symbols('x');p6=7*x*x+72*x-45;p8=2*x**3-135*x*x-160*x-45
    assert sp.expand((71198*x*x-4754655*x-9221545)*p6+(-249193*x-2742363)*p8)==538375860
    witness={S:0,B:0,C:-2,A:-6/a}
    assert sp.cancel(M[6].subs(witness)+1440*(7*a**4+72*a*a-45)/a)==0
    assert sp.cancel(M[8].subs(witness)+241920*(2*a**6-135*a**4-160*a*a-45)/(a*a))==0
    # Exceptional a^2+3=0: all calculations in an exact rational quotient.
    rem=lambda e:sp.rem(sp.Poly(sp.expand(e),a),sp.Poly(a*a+3,a)).as_expr()
    assert sp.expand(rem(M[6].subs({S:0,B:0,C:-2}))-(-4320*(-5*A+32*a)))==0
    assert rem(M[8].subs({S:0,B:0,C:-2,A:sp.Rational(32,5)*a}))==-sp.Rational(601316352,5)
    num,den=sp.fraction(sp.factor(M[4].subs({A:0,B:Bs,C:Cs})))
    assert sp.factor(rem(num)/(S*S*(115*S*S+128*a)))==81
    assert sp.gcd(den,a*a+3)==1
    GB=sp.groebner([S*S+sp.Rational(128,115)*a,a*a+3],S,a,A,order='lex')
    Ce=sp.Rational(146,23);Be=sp.Rational(388,23)*a*S
    assert sp.expand(GB.reduce(sp.expand(M[5].subs({C:Ce,B:Be})))[1]+1440*S*(529*A+155545*a)/529)==0
    assert GB.reduce(sp.expand(M[6].subs({C:Ce,B:Be,A:-sp.Rational(155545,529)*a})))[1]==-sp.Rational(3120931595808,60835)*a
    fa=7*a**4+72*a*a-45
    for m in range(1,8):
        num,den=sp.fraction(sp.cancel(M[m].subs(witness)));assert sp.rem(num,fa,a)==0 and sp.gcd(den,fa)==1
    num,den=sp.fraction(sp.cancel(M[8].subs(witness)));assert sp.gcd(num,fa)==1 and sp.gcd(den,fa)==1
    return {'branches':['a=0','a*(a^2+3)!=0,S=0','a*(a^2+3)*S!=0','a^2+3=0,S=0','a^2+3=0,S!=0'],'sharp_prefix':8,'witness':'a root of 7a^4+72a^2-45; A=-6/a,B=S=0,C=-2','witness_zero_moments':7,'eighth_nonzero_at_every_root':True}


def parent_control_checks(parent,p):
    syms=sp.symbols(' '.join(VARIABLES)); env=dict(zip(VARIABLES,syms));a=env['a']
    M={m:sp.Add(*(int(c)*sp.prod(v**e for v,e in zip(syms,mon)) for mon,c in parent['moments'][m-1]['terms'])) for m in range(1,6)}
    substitutions={};pivots=[1,2*a,6*(a*a+1),24*a*(a*a+3)]
    for m,name,pivot in zip(range(1,5),['r0','r1','r2','r3'],pivots):
        v=env[name];eq=sp.cancel(M[m].subs(substitutions,simultaneous=True))
        assert sp.cancel(sp.diff(eq,v)-pivot)==0
        value=sp.cancel(-eq.subs(v,0)/pivot);assert sp.cancel(eq.subs(v,value))==0
        substitutions[v]=value
    J5=sp.cancel(M[5].subs(substitutions,simultaneous=True))
    saved=(p/'mixed_wt_parent_12_reduced.ideal').read_text().split(',')[0]
    retained=sp.sympify(saved.replace('^','**'),locals=env);assert sp.cancel(J5-retained)==0
    u1,v1,u2=[env[n] for n in ['u1','v1','u2']];K=sp.cancel(sp.diff(J5,u2))
    assert sp.diff(K,u2)==0
    ku=120*((a*a+1)**3-4)/(a*(a*a+1));kv=720*(a**4+a*a+2)/(a*a+1)
    ru=a*(a*a-6)*(a*a+1)/5280;rv=(-a**4+4*a*a+13)*(a*a+1)/31680
    assert sp.cancel(sp.diff(K,u1)-ku)==0 and sp.cancel(sp.diff(K,v1)-kv)==0
    K0=K.subs({u1:0,v1:0});assert sp.cancel(K-ku*u1-kv*v1-K0)==0
    assert sp.cancel(ru*ku+rv*kv)==1
    k,l=sp.symbols('k l');inverse={u1:ru*(k-K0)-kv*l,v1:rv*(k-K0)+ku*l}
    assert sp.cancel(K.subs(inverse,simultaneous=True)-k)==0
    assert sp.cancel((-rv*u1+ru*v1).subs(inverse,simultaneous=True)-l)==0
    profile=json.loads((p/'mixed_wt_moment5_profile.json').read_text())
    assert sp.cancel(sp.sympify(profile['u2_pivot'],locals=env)-K)==0
    return {'triangular_pivots':list(map(str,pivots)),'fifth_matches_retained_input':True,
            'fifth_linear_in_u2':True,'control_determinant':'1',
            'boundary':'Split K!=0 and K=0; a*(a^2+1)*(a^2+3)!=0 only. No parent exclusion.'}


def exceptional_parent_checks(parent,p):
    vs=sp.symbols(' '.join(VARIABLES));e=dict(zip(VARIABLES,vs));a=e['a']
    M={m:sp.Add(*(int(c)*sp.prod(v**k for v,k in zip(vs,mon)) for mon,c in parent['moments'][m-1]['terms'])) for m in range(1,7)}
    profiles=json.loads((p/'mixed_wt_exceptional_quadratic_profiles.json').read_text())
    def reducer(k):
        mod=a*a+k
        def red(expr):
            n,d=sp.fraction(sp.cancel(expr));dr=sp.rem(d,mod,a);assert dr!=0
            return sp.rem(sp.rem(n,mod,a)*sp.invert(dr,mod,a),mod,a).expand()
        return red
    for k in [1,3]:
        red=reducer(k);sol={}
        for row in profiles[str(k)]['triangular']:
            v=e[row['variable']];eq=red(M[row['moment']].subs(sol,simultaneous=True))
            assert red(sp.diff(eq,v)-sp.sympify(row['pivot'],locals=e))==0
            value=sp.sympify(row['solution'],locals=e);assert red(eq.subs(v,value))==0;sol[v]=value
        J5=red(M[5].subs(sol,simultaneous=True));assert red(J5-sp.sympify(profiles[str(k)]['J5'],locals=e))==0
        if k==1:
            K=240*(6*a*e['s']*e['v0']-a*e['u0']*e['v0']+3*a*e['v1']+2*a+6*e['s']**2-e['u1']-2*e['v0']**2)
            assert red(sp.diff(J5,e['r2'])-K)==0
            assert sp.diff(K,e['u1'])==-240
        else:
            J4=red(M[4].subs(sol,simultaneous=True));assert red(J4-sp.sympify(profiles['3']['J4'],locals=e))==0
            H=e['u0']+2*e['s']-a*e['v0']
            assert red(sp.diff(J4,e['u2'])+144*H)==0 and sp.diff(J4,e['r3'])==0
            assert red(sp.diff(J5,e['r3'])+720*H)==0
    data=json.loads((p/'mixed_wt_alpha0_first4_profile.json').read_text());sol={}
    for row in data['first_three']:
        v=e[row['variable']];eq=sp.expand(M[row['moment']].subs(a,0).subs(sol,simultaneous=True))
        assert sp.diff(eq,v)==int(row['pivot'])
        value=sp.sympify(row['solution'],locals=e);assert sp.expand(eq.subs(v,value))==0;sol[v]=value
    J={m:sp.expand(M[m].subs(a,0).subs(sol,simultaneous=True)) for m in [4,5,6]}
    for m in [4,5]:assert sp.expand(J[m]-sp.sympify(data['J'+str(m)],locals=e))==0
    H=e['u0']+5*e['s'];assert sp.diff(J[4],e['u2'])==72*H
    assert sp.diff(J[4],e['r3'])==0 and sp.diff(J[5],e['r3'])==360*H
    boundary={e['u0']:-5*e['s']};D=e['r1']+10*e['s']*e['v0']+5*e['v1']+4
    assert sp.expand(sp.diff(J[5],e['u2']).subs(boundary)-360*D)==0
    assert sp.expand(sp.diff(J[6],e['r3']).subs(boundary)-2160*(D+1))==0
    return {'a=0':'r0,u1,r2 pivots 1,2,6; H=u0+5s controls u2 and r3; boundary D,D+1 retained',
            'a^2+1=0':'r0,r1,u2,r3 pivots 1,2a,12a,48a; fifth r2 control has u1 slope -240',
            'a^2+3=0':'r0,r1,r2 pivots 1,2a,-12; H=u0+2s-a*v0 controls fourth u2 and fifth r3',
            'status':'Exact reductions only; no exceptional full-parent face excluded.'}


def main():
    start=time.monotonic();p=Path(__file__).parent
    parent=json.loads((p/'mixed_wt_parent_moments_12.json').read_text())
    assert parent['compiler_sha256']==hashlib.sha256((p/'compile_mixed_wt_parent.py').read_bytes()).hexdigest()
    cert=json.loads(gzip.decompress((p/'mixed_wt_spine_prefix8_certificate.json.gz').read_bytes()))
    assert cert['log_sha256']==hashlib.sha256((p/'mixed_wt_spine_prefix8_lift.log').read_bytes()).hexdigest()
    control=parent_control_checks(parent,p);exceptional=exceptional_parent_checks(parent,p);rows=numeric_checks(parent);syms,M=symbolic_spine();proof=proof_checks(cert,syms,M);identity=fraction_certificate(cert)
    run=json.loads((p/'mixed_wt_parent_12_exact_run.json').read_text());assert run['process_state']=='terminal_timeout'
    assert run['input_sha256']==hashlib.sha256((p/'mixed_wt_parent_12.sing').read_bytes()).hexdigest()
    assert run['log_sha256']==hashlib.sha256((p/'mixed_wt_parent_12_exact.log').read_bytes()).hexdigest()
    assert (p/'mixed_wt_parent_12_reduced.ideal').is_file()
    rec={'scope':'Complete homogeneous-cubic spine and exact full-parent reductions; no full-parent closure claimed.',
         'exceptional_parent_checks':exceptional,'parent_control_checks':control,'numeric_Wick_phase_checks':rows,'spine_proof':proof,'independent_cofactor_identity':identity,
         'parent_status':'UNKNOWN: exact generic Q(a) standard-basis stage timed out; reduced input retained.',
         'parent_moment_term_counts':[len(x['terms']) for x in parent['moments']],
         'sympy_version':sp.__version__,'elapsed_seconds':time.monotonic()-start,
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_hashes':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['verify_moment_transfers.py','compile_mixed_wt_parent.py','mixed_wt_parent_moments_12.json','mixed_wt_spine_prefix8_certificate.json.gz','mixed_wt_parent_12_reduced.ideal','mixed_wt_moment5_profile.json','mixed_wt_parent_12_exact_run.json','mixed_wt_alpha0_first4_profile.json','mixed_wt_exceptional_quadratic_profiles.json']}}
    (p/'mixed_wt_chart_verification.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
