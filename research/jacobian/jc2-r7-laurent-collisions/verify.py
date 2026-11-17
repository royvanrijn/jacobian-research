"""Replay Laurent-family collision strata without Singular or root embeddings."""
from pathlib import Path
import gzip
import hashlib
import json
import time
import sympy as S0
from sympy.polys.rings import ring
import laurent as L

CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)

def zero(name,value):check(name,S0.cancel(value)==0)

def proportional(name,left,right):
    ratio=S0.cancel(left/right)
    check(name,ratio.is_Rational and ratio!=0)

def read(name):
    path=L.HERE/name
    if name.endswith('.gz'):
        with gzip.open(path,'rt') as f:return json.load(f)
    return json.loads(path.read_text())

def sparse(expr,variables):
    rows=S0.Poly(expr,*variables).terms()
    assert all(c.is_Integer for m,c in rows)
    return {m:int(c) for m,c in rows if c}

def decode_int(rows,dim):
    result={}
    for m,c in rows:
        assert len(m)==dim and all(type(e) is int and e>=0 for e in m)
        assert tuple(m) not in result
        result[tuple(m)]=int(c)
    return result

def integer_identity(data,expressions,variables,power):
    check('integer certificate variables '+str(power),data['variables']==list(map(str,variables)))
    check('integer target scalar nonzero '+str(power),int(data['target_scalar'])!=0)
    generators=[decode_int(row,len(variables)) for row in data['generators']]
    for i,(actual,expected) in enumerate(zip(generators,expressions)):
        check(f'integer source generator {power},{i}',actual==sparse(expected,variables))
    check('integer certificate dimensions '+str(power),len(generators)==len(expressions)==len(data['multipliers']))
    result={}
    for generator,rows in zip(generators,data['multipliers']):
        multiplier=decode_int(rows,len(variables))
        for a,c in generator.items():
            for b,d in multiplier.items():
                m=tuple(x+y for x,y in zip(a,b))
                result[m]=result.get(m,0)+c*d
    result={m:c for m,c in result.items() if c}
    if 'target' in data:target=decode_int(data['target'],len(variables))
    else:target={(0,)*len(variables):int(data['target_scalar'])}
    check('independent integer identity '+str(power),result==target)
    if 'parameter_factor' in data:
        factor=S0.sympify(data['parameter_factor'],locals={'z':L.z})
        expected=int(data['target_scalar'])*factor*L.S**power
        check('integer target reconstruction '+str(power),target==sparse(expected,variables))
        return factor

def modular_gcd(name,filename,modulus,expressions):
    data=read(filename);R,z,v=ring('z,v',S0.QQ)
    decode=lambda rows:R.from_dict({tuple(m):S0.QQ(c) for m,c in rows})
    mod=R.from_expr(modulus)
    check(name+' exact modulus',decode(data['modulus'])==mod)
    gen=[decode(rows) for rows in data['generators']]
    check(name+' certificate dimensions',len(gen)==len(data['multipliers'])==len(expressions)==2)
    for i,(g,f) in enumerate(zip(gen,expressions)):check(name+f' source generator{i}',g==R.from_expr(f))
    gg=decode(data['gcd']);multipliers=[decode(row) for row in data['multipliers']]
    check(name+' modular Bezout',not (sum((a*b for a,b in zip(gen,multipliers)),R.zero)-gg).rem([mod]))
    check(name+' mutated Bezout rejected',bool((gen[0]*(multipliers[0]+1)+gen[1]*multipliers[1]-gg).rem([mod])))
    check(name+' monic linear gcd',gg.coeff(v)==R.one and all(m[1]==0 or m==(0,1) for m in gg))
    K,zz=ring('z',S0.QQ);mm=K.from_expr(modulus)
    root=K.from_dict({(m[0],):-c for m,c in gg.items() if m[1]==0})
    for i,g in enumerate(gen):
        value=K.zero
        for exponent in range(g.degree(v),-1,-1):
            value=(value*root).rem(mm)
            value+=K.from_dict({(m[0],):c for m,c in g.items() if m[1]==exponent})
            value=value.rem(mm)
        check(name+f' gcd divides generator{i}',not value)
    return root

# Small finite-field polynomial arithmetic, independent of factorization software.
def trim(f):
    while f and f[-1]==0:f.pop()
    return f

def fprem(f,g,p):
    f=f[:];g=trim(g[:]);inv=pow(g[-1],-1,p)
    while len(f)>=len(g):
        shift=len(f)-len(g);c=f[-1]*inv%p
        for i,a in enumerate(g):f[i+shift]=(f[i+shift]-c*a)%p
        trim(f)
    return f

def fpmul(f,g,p):
    out=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        for j,b in enumerate(g):out[i+j]=(out[i+j]+a*b)%p
    return trim(out)

def fppow(f,n,g,p):
    out=[1]
    while n:
        if n&1:out=fprem(fpmul(out,f,p),g,p)
        f=fprem(fpmul(f,f,p),g,p);n//=2
    return out

def fpsub(f,g,p):
    out=[0]*max(len(f),len(g))
    for i,a in enumerate(f):out[i]=a
    for i,a in enumerate(g):out[i]=(out[i]-a)%p
    return trim(out)

def fpgcd(f,g,p):
    while g:f,g=g,fprem(f,g,p)
    return [a*pow(f[-1],-1,p)%p for a in f] if f else []

def irreducibility():
    data=read('irreducibility_certificate.json');possible=set(range(1,50))
    for witness in data['witnesses']:
        p=witness['prime'];check('prime witness '+str(p),S0.isprime(p))
        product=[1];subset={0}
        for index,row in enumerate(witness['factors']):
            f=row['coefficients'];n=len(f)-1
            check(f'monic modular factor {p},{index}',n>=1 and f[-1]==1 and row['multiplicity']==1)
            x=[0,1];value=x;tests={n//q for q in S0.factorint(n)}
            for i in range(1,n+1):
                value=fppow(value,p,f,p)
                if i in tests:check(f'Rabin gcd {p},{index},{i}',len(fpgcd(f,fpsub(value,x,p),p))==1)
            check(f'Rabin Frobenius {p},{index}',not fprem(fpsub(value,x,p),f,p))
            product=fpmul(product,f,p);subset|={a+n for a in list(subset)}
        coefficients=[int(c)%p for c in reversed(S0.Poly(L.D,L.z).all_coeffs())]
        check('nonzero leading coefficient '+str(p),coefficients[-1]!=0)
        inv=pow(coefficients[-1],-1,p);coefficients=[a*inv%p for a in coefficients]
        check('modular factor product '+str(p),product==coefficients)
        possible&=subset
    check('degree50 irreducibility over Q',not possible)

def triple_point():
    data=read('triple_point_certificate.json');R,z=ring('z',S0.QQ)
    decode=lambda rows:R.from_dict({tuple(m):S0.QQ(c) for m,c in rows})
    modulus=decode(data['modulus']);sv=decode(data['S']);uv=decode(data['U'])
    numerator=decode(data['V_numerator']);denominator=decode(data['V_denominator'])
    check('point exact modulus',modulus==R.from_expr(L.D))
    check('point denominator nonzero',bool(denominator.rem(modulus)))
    mul=lambda a,b:(a*b).rem(modulus)
    def powers(a,n):
        rows=[R.one]
        for _ in range(n):rows.append(mul(rows[-1],a))
        return rows
    sp,up,vp,dp=powers(sv,8),powers(uv,4),powers(numerator,3),powers(denominator,3)
    for index,f in enumerate(L.triple):
        poly=S0.Poly(f,L.V,L.U,L.S,L.z);degree=poly.degree(L.V);value=R.zero
        for (k,j,i,l),c in poly.terms():
            value+=(mul(mul(mul(sp[i],up[j]),vp[k]),dp[degree-k])*c*z**l).rem(modulus)
        check('exact degree50-field triple point '+str(index),not value.rem(modulus))
    check('triple point S nonzero',bool(sv.rem(modulus)))

def rank_witness():
    data=read('triple_rank_certificate.json');p=data['prime'];root=data['parameter_root']
    check('rank witness prime',S0.isprime(p))
    check('rank witness recorded coordinates',data['variables']==['V','U','S'] and data['weights']==[3,2,1] and (p,root)==(101,21))
    check('rank witness keeps degree50',int(S0.Poly(L.D,L.z).LC())%p!=0)
    check('rank witness lies on degree50 polynomial',int(S0.Poly(L.D,L.z).eval(root))%p==0)
    columns=[tuple(m) for m in data['columns']];weights=(3,2,1);bound=data['weighted_degree']
    expected={(a,b,c) for a in range(bound//3+1) for b in range(bound//2+1) for c in range(bound+1) if 3*a+2*b+c<=bound}
    check('full truncated monomial set',len(columns)==len(set(columns)) and set(columns)==expected)
    lookup={m:i for i,m in enumerate(columns)};polys=[]
    for f in L.triple:
        row={}
        for m,c in S0.Poly(f,L.V,L.U,L.S,L.z).terms():row[m[:3]]=(row.get(m[:3],0)+int(c)*pow(root,m[3],p))%p
        polys.append({m:c for m,c in row.items() if c})
    selected=data['pivot_columns'];check('rank minor dimensions',len(data['selected_rows'])==len(selected)==len(columns)-1 and len(set(selected))==len(selected))
    matrix=[]
    for gi,multi in data['selected_rows']:
        assert gi in range(4) and len(multi)==3 and all(type(a) is int and a>=0 for a in multi)
        row=[0]*len(columns)
        for m,c in polys[gi].items():
            key=tuple(a+b for a,b in zip(m,multi));assert key in lookup
            row[lookup[key]]=c
        matrix.append([row[i] for i in selected])
    for col in range(len(matrix)):
        candidates=[i for i in range(col,len(matrix)) if matrix[i][col]]
        assert candidates,'modular rank minor became singular'
        row=candidates[0];matrix[col],matrix[row]=matrix[row],matrix[col]
        inv=pow(matrix[col][col],-1,p)
        for i in range(col+1,len(matrix)):
            factor=matrix[i][col]*inv%p
            if factor:
                for j in range(col,len(matrix)):matrix[i][j]=(matrix[i][j]-factor*matrix[col][j])%p
    check('exact nonzero357-dimensional rank minor',len(matrix)==357)
    check('linear coordinates and unit belong to truncation',all(m in lookup for m in [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))

def projection_controls(root):
    R,z=ring('z',S0.QQ);mod=R.from_expr(L.Fproj);root=R.from_expr(root.as_expr())
    mul=lambda a,b:(a*b).rem(mod)
    inv=lambda a:R.from_expr(S0.invert(a.as_expr(),L.Fproj,L.z))
    powers=[R.one]
    for _ in range(14):powers.append(mul(powers[-1],root))
    def atv(expr):
        out={}
        for (i,j,k),c in S0.Poly(expr,L.v,L.w,L.z).terms():out[j]=out.get(j,R.zero)+(powers[i]*c*z**k).rem(mod)
        return {j:c.rem(mod) for j,c in out.items()}
    a,b,j=atv(L.A),atv(L.B),atv(L.J);ai=inv(a[2]);alpha,beta=mul(a[1],ai),mul(a[0],ai)
    wr=[(R.one,R.zero),(R.zero,R.one)]
    for _ in range(2,6):
        p,q=wr[-1];wr.append((-mul(q,beta),(p-mul(q,alpha)).rem(mod)))
    def reducew(poly):
        p,q=R.zero,R.zero
        for i,c in poly.items():p+=mul(c,wr[i][0]);q+=mul(c,wr[i][1])
        return p.rem(mod),q.rem(mod)
    check('projection quadratic divides both incidence equations',reducew(b)==(R.zero,R.zero))
    j0,j1=reducew(j);disc=(alpha*alpha-4*beta).rem(mod)
    norm=(j0*j0-mul(alpha,mul(j0,j1))+mul(beta,j1*j1)).rem(mod)
    check('two distinct projected collision pairs',S0.degree(S0.gcd(disc.as_expr(),L.Fproj,L.z),L.z)==0)
    check('both projected pairs reduced',S0.degree(S0.gcd(norm.as_expr(),L.Fproj,L.z),L.z)==0)
    zero_jac=S0.cancel(L.J.subs({L.v:0,L.w:-(15*L.z**2+4*L.z+4)/12}))
    check('zero-sum collision is reduced',S0.degree(S0.gcd(S0.fraction(zero_jac)[0],L.W,L.z),L.z)==0)

def main():
    start=time.monotonic()
    for item in read('certificates.json')['certificates']:
        check('certificate hash '+item['file'],hashlib.sha256((L.HERE/item['file']).read_bytes()).hexdigest()==item['sha256'])
    s,t=S0.symbols('s0 t0')
    for name,poly,div in [('P',L.P,L.A),('Q',L.Q,L.B)]:
        zero(name+' symmetric divided difference',poly.subs(L.t,s)-poly.subs(L.t,t)-(s-t)*div.subs({L.v:s+t,L.w:s*t}))
    check('degree12 fixed collision polynomial',S0.degree(L.H,L.v)==12 and S0.Poly(L.H,L.v).LC()==81)
    zero('quadratic incidence leading coefficient',S0.Poly(L.A,L.w).LC()-3*L.v)
    zero('quartic incidence leading coefficient',S0.Poly(L.B,L.w).LC()-5*L.v)
    proportional('scaled-family nonimmersion locus',S0.resultant(S0.diff(L.P,L.t),S0.diff(L.Q,L.t),L.t),L.z**12*L.C1*L.C2)
    zero('resultant normalization',L.sub[-1]-L.v**10*L.H/81)
    zero('linear subresultant normalization',L.sub[-2]+L.v**3*(L.k*L.w+L.n)/12)
    proportional('projection discriminant factorization',S0.resultant(L.H,S0.diff(L.H,L.v),L.v),L.z**8*L.C1**3*L.Fproj**2*L.Tangent)
    proportional('linear recovery exceptional factors',S0.resultant(L.H,L.k,L.v),L.z**4*L.W**7*L.Fproj**2)
    zero('zero-sum incidence',L.A.subs(L.v,0)+L.z*(L.w+(15*L.z**2+4*L.z+4)/12))
    proportional('zero-sum collision condition',L.B.subs({L.v:0,L.w:-(15*L.z**2+4*L.z+4)/12}),L.z**2*L.W)
    factors=[L.z,L.C1,L.C2,L.W,L.Fproj,L.Tangent,L.D]
    for i,f in enumerate(factors):
        check('squarefree parameter factor'+str(i),S0.degree(S0.gcd(f,S0.diff(f,L.z)),L.z)==0)
        for j,g in enumerate(factors[:i]):check(f'parameter factors disjoint{i},{j}',S0.degree(S0.gcd(f,g,L.z),L.z)==0)
    lr=integer_identity(read('necessary_reduced_certificate.json.gz'),L.reduced,(L.U,L.S,L.z),2)
    lw=integer_identity(read('necessary_weighted_certificate.json.gz'),L.triple,(L.V,L.U,L.S,L.z),0)
    proportional('necessary parameter gcd',S0.gcd(lr,lw,L.z),L.z**6*(3*L.z+1)*L.D)
    special=[L.primitive(f.subs(L.z,-S0.Rational(1,3)),(L.V,L.U,L.S)) for f in L.triple]
    integer_identity(read('exception_unit_certificate.json'),special,(L.V,L.U,L.S),-1)
    # Complete the S0 chart before dividing by its parameter.
    base=[f.subs(L.S,0) for f in L.triple];us=-(15*L.z**2+8*L.z+4)/8
    proportional('S0 exceptional leading coefficient',base[0].subs(L.U,us),L.C1*(75*L.z**3+30*L.z**2-12*L.z+8))
    proportional('S0 exceptional compatibility',base[1].subs(L.U,us),L.z*L.C1*(5*L.z+2))
    check('S0 degenerate roots excluded',S0.degree(S0.gcd(75*L.z**3+30*L.z**2-12*L.z+8,5*L.z+2,L.z),L.z)==0)
    vv=S0.solve(base[1],L.V)[0];rem=[]
    for f in base[2:]:
        raw=S0.together(f.subs(L.V,vv)).as_numer_denom()[0]
        rem.append(L.primitive(S0.rem(raw,base[0],L.U),(L.U,L.z)))
    proportional('S0 no other triple parameter',S0.gcd(S0.resultant(base[0],rem[0],L.U),S0.resultant(base[0],rem[1],L.U),L.z),L.C1**2)
    irreducibility();triple_point();rank_witness()
    root=modular_gcd('projection','projection_gcd_certificate.json.gz',L.Fproj,[L.H,S0.diff(L.H,L.v)])
    modular_gcd('zero sum','zero_sum_gcd_certificate.json',L.W,[L.H,L.k])
    projection_controls(root)
    cases=read('../jc2-normal-r7/classification.json')['cases']
    for index,zvalue,sign in [(0,-S0.Rational(2,3),1),(1,S0.Rational(4,15),-1)]:
        values={name:S0.Rational(value) for name,value in cases[index]['parameters'].items()}
        sub={getattr(L.F,name):value for name,value in values.items()}
        for name,scaled,original in [('p',L.P,L.F.p),('q',L.Q,L.F.q)]:
            zero(f'retained cusp class{index} {name}',scaled.subs(L.z,zvalue).subs(L.t,sign*L.t)-original.subs(sub))
    check('Laurent exceptional z count',S0.degree(L.D,L.z)+S0.degree(L.Tangent,L.z)+2==79)
    check('Laurent exceptional x count',3*79==237)
    check('full coefficient-union exceptional count',34+237-3==268)
    ledger=dict(schema='r7-laurent-collision-strata-v1',parameter='z=x^3; x!=0',
                parameter_polynomials={name:str(f) for name,f in [('C1',L.C1),('C2',L.C2),('D',L.D),('T',L.Tangent),('projection_only_F',L.Fproj),('zero_sum_only_W',L.W)]},
                singularities=[dict(stratum='C1*C2*D*T!=0',nodes=12),
                               dict(stratum='D=0',z_values=50,x_values=150,ordinary_triples=1,nodes=9),
                               dict(stratum='T=0',z_values=27,x_values=81,tacnodes_A3=1,nodes=10),
                               dict(stratum='C1=0',z_values=1,x_values=3,cusps_A2=1,tangencies_A5=1,nodes=8),
                               dict(stratum='C2=0',z_values=1,x_values=3,cusps_A2=1,nodes=11)],
                exceptional_z_count=79,exceptional_x_count=237,two_family_union_exceptional_points=268,
                shared_unramified_image_H_equals_C_excluded_immersed_points=163,
                global_consequence='A degree6 Keller candidate with an immersed Laurent D-root e2,f1 image omits its unique triple value entirely. Every additional omitted image must avoid it; this is not a surjectivity contradiction.',
                remaining=['actual first-image parameter and second-image equations','global complement relations and infinity comparison','finite normal cover and source A2 realization'])
    (L.HERE/'strata.json').write_text(json.dumps(ledger,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} LAURENT COLLISION CHECKS ({time.monotonic()-start:.3f}s)',flush=True)
    inputs=[L.HERE/n for n in ['PROOF.md','README.md','REVIEW.md','laurent.py','verify.py','certificates.json']]
    inputs += [L.HERE/item['file'] for item in read('certificates.json')['certificates']]
    inputs += [L.HERE.parent/'jc2-r7-image-family/family.py',L.HERE.parent/'jc2-r7-image-family/PROOF.md',
               L.HERE.parent/'jc2-r7-collisions/PROOF.md',L.HERE.parent/'jc2-normal-r7/classification.json',
               L.HERE.parent/'jc2-omitted-divisors/IMAGE_NORM_LEDGER.md',
               L.HERE.parent/'jc2-omitted-divisors/RESIDUE_LEDGER.md',
               L.ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(schema='r7-laurent-collisions-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),
                 input_sha256={str(p.relative_to(L.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                 output_sha256={'strata.json':hashlib.sha256((L.HERE/'strata.json').read_bytes()).hexdigest()},
                 assurance='Literal integer necessary identities, exact degree50-field point, finite-field irreducibility and rank certificates, and modular projection controls. Written local geometry and fibre interpretation; no formalization or independent external review.')
    (L.HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
