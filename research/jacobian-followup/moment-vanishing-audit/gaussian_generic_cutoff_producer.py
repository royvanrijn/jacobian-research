"""Actual Gaussian certificates over Frac(QQ[c]/prime) and family cutoffs.

Generic operators are produced over the exact fraction field by Macaulay2.
Finite Weyl/IBP/indicial/actual-jet identities replay in Python modulo the
provided prime. Primality and the family decomposition retain their CAS
boundary. Original moment bounds extend to every point of a safe component.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

import sympy as sp
import gaussian_family_producer as F
import gaussian_source_hom_producer as G

CORE=Path(__file__).with_name('gaussian_generic_field_core.m2')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


class Fraction:
    """An exact quotient-domain fraction; equality uses cross multiplication."""
    def __init__(self,field,numerator,denominator=1):
        self.field=field
        n,d=field.normal(numerator),field.normal(denominator)
        if d==0:raise ZeroDivisionError('zero denominator modulo coefficient prime')
        if n==0:n,d=sp.S.Zero,sp.S.One
        else:
            common=sp.gcd(n,d)
            if common!=1:n,d=sp.cancel(n/common),sp.cancel(d/common)
            leading=sp.Poly(d,*field.parameters,domain=sp.QQ).LC()
            n,d=sp.expand(n/leading),sp.expand(d/leading)
        self.n,self.d=n,d
    def _other(self,x):return self.field.convert(x)
    def __bool__(self):return self.n!=0
    def __neg__(self):return Fraction(self.field,-self.n,self.d)
    def __add__(self,x):
        x=self._other(x);return Fraction(self.field,self.n*x.d+x.n*self.d,self.d*x.d)
    __radd__=__add__
    def __sub__(self,x):return self+-self._other(x)
    def __rsub__(self,x):return self._other(x)+-self
    def __mul__(self,x):
        x=self._other(x);return Fraction(self.field,self.n*x.n,self.d*x.d)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self._other(x);return Fraction(self.field,self.n*x.d,self.d*x.n)
    def __rtruediv__(self,x):return self._other(x)/self
    def __eq__(self,x):
        x=self._other(x);return self.field.normal(self.n*x.d-x.n*self.d)==0
    def __pow__(self,k):
        if type(k) is not int:raise TypeError('integer fraction power required')
        return Fraction(self.field,self.n**k,self.d**k) if k>=0 else Fraction(self.field,self.d**(-k),self.n**(-k))


class GenericField:
    """The input prime is a field-presentation contract, checked by the CAS."""
    def __init__(self,parameter_count,prime):
        self.parameters=sp.symbols('c:'+str(parameter_count))
        self.prime=prime
        polynomials=[self.polynomial(g) for g in prime]
        self.gb=sp.groebner(polynomials,*self.parameters,domain=sp.QQ)
        if any(p.as_expr()==1 for p in self.gb.polys):raise ValueError('generic prime must be proper')
        self.K=self
        self.zero=Fraction(self,0);self.one=Fraction(self,1)
    def polynomial(self,terms):
        return sp.expand(sum(sp.Rational(c)*sp.prod(x**k for x,k in zip(self.parameters,e)) for e,c in terms))
    @lru_cache(maxsize=16384)
    def normal(self,value):
        return sp.expand(self.gb.reduce(sp.expand(value))[1])
    def convert(self,value):
        if isinstance(value,Fraction):
            if value.field is not self:raise ValueError('different coefficient fields')
            return value
        value=sp.Rational(value)
        return Fraction(self,value)
    def parse(self,value):
        if isinstance(value,dict):
            if set(value)!={'numerator','denominator'}:raise ValueError('invalid generic coefficient')
            return Fraction(self,self.polynomial(value['numerator']),self.polynomial(value['denominator']))
        return self.convert(value)
    def encode_polynomial(self,value):
        return [[list(e),str(c)] for e,c in sp.Poly(value,*self.parameters,domain=sp.QQ).terms() if c]
    def encode(self,value):
        value=self.convert(value)
        return dict(numerator=self.encode_polynomial(value.n),denominator=self.encode_polynomial(value.d))
    def integer_roots(self,coefficients):
        coefficients={int(k):self.convert(v) for k,v in coefficients.items() if v}
        if not coefficients:raise ValueError('zero generic indicial polynomial')
        denominator=sp.prod(v.d for v in coefficients.values())
        polynomials={}
        z=sp.Symbol('z')
        items=list(coefficients.items())
        for i,(degree,c) in enumerate(items):
            cleared=self.normal(c.n*sp.prod(v.d for j,(_,v) in enumerate(items) if j!=i))
            for mon,scalar in sp.Poly(cleared,*self.parameters,domain=sp.QQ).terms():
                polynomials[mon]=polynomials.get(mon,sp.S.Zero)+scalar*z**degree
        content=sp.Poly(0,z,domain=sp.QQ)
        for poly in polynomials.values():content=sp.gcd(content,sp.Poly(poly,z,domain=sp.QQ))
        if content.is_zero:raise ValueError('zero rational generic content')
        content=content.monic()
        roots=sorted(int(v) for v in content.ground_roots() if v.is_Integer)
        return roots,dict(generic_coefficients=[[k,self.encode(v)] for k,v in sorted(coefficients.items())],
                          clearing_denominator=self.encode_polynomial(self.normal(denominator)),
                          rational_content_terms=[[int(e[0]),str(c)] for e,c in content.terms()])


def indicial(operator,field,which):
    weight=(max if which=='infinity' else min)(a-b for (a,b),c in operator if field.parse(c))
    coefficients={};z=sp.Symbol('z')
    for (a,b),c in operator:
        if a-b!=weight:continue
        for (k,),v in sp.Poly(sp.ff(z,b),z).terms():
            coefficients[k]=coefficients.get(k,field.zero)+field.parse(c)*field.convert(v)
    roots,certificate=field.integer_roots(coefficients)
    nonnegative=[v for v in roots if v>=0]
    certificate.update(weight=weight,nonnegative_integer_roots=nonnegative,bound=max(nonnegative,default=-1))
    return certificate


def normalize(data):
    family={k:v for k,v in data.items() if k!='generic_prime'}
    n,r,phase,base,precision=F.normalize(family)
    prime=data['generic_prime']
    for g in prime:F.terms(g,r)
    field=GenericField(r,prime)
    for g in base:
        poly=sum(c*sp.prod(x**k for x,k in zip(field.parameters,e)) for e,c in g.items())
        if field.normal(poly)!=0:raise ValueError('generic prime does not contain parameter ideal')
    generic={}
    for es,c in phase.items():
        coefficient=Fraction(field,c*sp.prod(x**k for x,k in zip(field.parameters,es[n:])))
        G.sparse_add(generic,{tuple(es[:n]):coefficient},field)
    return family,n,r,phase,precision,field,generic


def make_script(data):
    family,n,r,phase,precision,field,generic=normalize(data)
    cs=['c'+str(i) for i in range(r)]
    xs=['x'+str(i) for i in range(n)]
    coords=xs+['t'];derivatives=['d'+x for x in coords]
    prime=','.join(F.expression(F.terms(g,r),cs) for g in data['generic_prime']) or '0_B'
    gens=[]
    for i,x in enumerate(xs):
        linear='+'.join('('+str(precision[i,j])+')*'+y for j,y in enumerate(xs) if precision[i,j])
        gens.append('d'+x+'+('+linear+')-t*diff('+x+',P)')
    gens.append('dt-P')
    return f'''load "{CORE.as_posix()}";
B=QQ[{','.join(cs)}];
PP=ideal matrix(B,{{{{{prime}}}}});
if not isPrime PP then error "generic coefficient ideal must be PRIME";
gCoefficientBase=B;gCoefficientQuotient=B/PP;
K=frac gCoefficientQuotient;
W=K[{','.join(coords+derivatives)},WeylAlgebra=>{{{','.join(x+'=>'+d for x,d in zip(coords,derivatives))}}}];
P=sub({F.expression(phase,xs+cs)},W);
I=ideal({','.join(gens)});
R=gGenericSource(I,{n});
R=new HashTable from (pairs R | {{"generic_prime_checked"=>true}});
stdio << "RESULT_JSON|" << toJSON R << endl;
print("GENERIC_COMPLETE");exit 0;
'''


def connections(n,phase,precision,field):
    m=n+1;P={e+(0,):c for e,c in phase.items()};hs=[]
    for i in range(n):
        h={}
        for j in range(n):
            e=[0]*m;e[j]=1
            G.sparse_add(h,{tuple(e):-field.convert(precision[i,j])},field)
        for e,c in G.derivative(P,i,field).items():
            ee=list(e);ee[n]+=1;G.sparse_add(h,{tuple(ee):c},field)
        hs.append(h)
    return hs+[P]


def original_generators(n,phase,precision,field):
    m=n+1;zero=(0,)*(2*m)
    hs=connections(n,phase,precision,field)
    gens=[]
    for i in range(m):
        e=list(zero);e[m+i]=1;g={tuple(e):field.one}
        G.sparse_add(g,{a+(0,)*m:-c for a,c in hs[i].items()},field)
        gens.append(g)
    return gens,hs


def verify_source(source,n,phase,precision,field):
    m=n+1;zero=(0,)*(2*m)
    parse=lambda op:{tuple(e):field.parse(c) for e,c in op}
    encode=lambda op:[[list(e),field.encode(c)] for e,c in sorted(op.items()) if c]
    roots,root_certificate=field.integer_roots({k:field.parse(c) for k,c in source['integration_b_terms']})
    if sorted(source['integration_integer_roots'])!=roots or source['integration_cutoff']!=max(roots):
        raise ValueError('wrong generic integration threshold')
    B=source['integration_cutoff'];basis=[tuple(e) for e in source['amplitude_exponents']]
    if B<0 or len(set(basis))!=len(basis) or len(basis)!=math.comb(n+B,n) or any(len(e)!=n or min(e)<0 or sum(e)>B for e in basis):
        raise ValueError('incomplete original generic amplitude coverage')
    rank=len(basis);cols=source['relation_columns']
    if source['rank']!=rank:raise ValueError('wrong generic module rank')
    if type(cols) is not int or cols<0:raise ValueError('invalid generic relation count')
    if len(source['relations'])!=rank or any(len(row)!=cols for row in source['relations']):
        raise ValueError('incomplete generic relation matrix')
    if len(source['coordinate_annihilators'])!=rank or len(source['coordinate_annihilator_left_witnesses'])!=rank:
        raise ValueError('missing generic coordinate certificate')
    if any(len(row)!=cols for row in source['coordinate_annihilator_left_witnesses']):
        raise ValueError('incomplete generic coordinate witness')
    gens,hs=original_generators(n,phase,precision,field)
    if len(gens)!=len(source['original_recurrence_left_witnesses']):raise ValueError('incomplete original recurrence witness')
    total={}
    for w,g in zip(source['original_recurrence_left_witnesses'],gens):
        G.sparse_add(total,G.sparse_product(parse(w),g,field),field)
    recurrence=parse(source['original_recurrence_operator'])
    if total!=recurrence:raise ValueError('invalid ORIGINAL generic recurrence identity')
    theta={zero:field.convert(-n)}
    for i in range(n):
        e=list(zero);e[i]=e[m+i]=1;theta[tuple(e)]=-field.one
    bterms={k:field.parse(c) for k,c in source['integration_b_terms']};power={zero:field.one};beuler={}
    for k in range(max(bterms)+1):
        G.sparse_add(beuler,power,field,bterms.get(k,field.zero))
        power=G.sparse_product(power,theta,field)
    remainder=dict(recurrence);G.sparse_add(remainder,beuler,field,-field.one)
    C=sp.Rational(0)
    for e in remainder:
        drop=sum(e[m:m+n])-sum(e[:n])
        if drop<1:raise ValueError('generic recurrence fails to lower original degree')
        C=max(C,sp.Rational(e[n]-e[-1],drop))
    relation_primitives=[];one={(0,)*m:field.one}
    cols=source['relation_columns']
    if len(source['unprojected_relation_operators'])!=cols:raise ValueError('missing generic source columns')
    for j in range(cols):
        op={}
        for i,e in enumerate(basis):
            time=G.lift_time_operator(source['relations'][i][j],n,field)
            G.sparse_add(op,G.sparse_product(time,{e+(0,)+(0,)*m:field.one},field),field)
        difference=parse(source['unprojected_relation_operators'][j]);G.sparse_add(difference,op,field,-field.one)
        qs,rest=G.divide_spatial_derivatives_left(difference,n,field)
        if rest:raise ValueError('wrong generic original amplitude projection')
        primitives=[{e:-c for e,c in G.exponential_action(q,one,hs,field).items()} for q in qs]
        rhs={}
        for k,p in enumerate(primitives):G.sparse_add(rhs,G.connection_apply(p,k,hs,field),field)
        if rhs!=G.exponential_action(op,one,hs,field):raise ValueError('invalid ORIGINAL generic relation IBP')
        relation_primitives.append(primitives)
    scalar_primitives=[];degrees=[]
    for i,ann in enumerate(source['coordinate_annihilators']):
        if not ann:raise ValueError('zero generic scalar operator')
        for k in range(len(basis)):
            left={}
            for j,w in enumerate(source['coordinate_annihilator_left_witnesses'][i]):
                G.sparse_add(left,G.weyl_product(w,source['relations'][k][j],field),field)
            if left!=(parse(ann) if k==i else {}):raise ValueError('wrong generic LEFT coordinate identity')
        primitives=[{} for _ in range(n)]
        for j,w in enumerate(source['coordinate_annihilator_left_witnesses'][i]):
            for k in range(n):G.sparse_add(primitives[k],G.exponential_action(G.lift_time_operator(w,n,field),relation_primitives[j][k],hs,field),field)
        rhs={}
        for k,p in enumerate(primitives):G.sparse_add(rhs,G.connection_apply(p,k,hs,field),field)
        marker={basis[i]+(0,):field.one}
        if rhs!=G.exponential_action(G.lift_time_operator(ann,n,field),marker,hs,field):
            raise ValueError('invalid ORIGINAL generic scalar IBP')
        scalar_primitives.append(primitives)
        degrees.append(dict(infinity=indicial(ann,field,'infinity'),origin=indicial(ann,field,'origin')))
    denominators={}
    def retain_denominators(value):
        if isinstance(value,dict):
            if set(value)=={'numerator','denominator'}:
                denominator=field.parse(value).d
                if denominator!=1:
                    encoded=field.encode_polynomial(denominator)
                    denominators[json.dumps(encoded)]=encoded
            else:
                for v in value.values():retain_denominators(v)
        elif isinstance(value,list):
            for v in value:retain_denominators(v)
    retain_denominators(source)
    return dict(integration_content_certificate=root_certificate,original_recurrence_identity_verified=True,
                original_degree_lowering_remainder=encode(remainder),recurrence_slope=str(C),
                original_relation_IBP_primitives=[[encode(p) for p in qs] for qs in relation_primitives],
                original_coordinate_IBP_primitives=[[encode(p) for p in qs] for qs in scalar_primitives],
                scalar_indicial_certificates=degrees,
                nonconstant_coefficient_denominators=[denominators[k] for k in sorted(denominators)],
                specialization_contract='Operators live in the generic field and specialize where their retained denominators are nonzero; ACTUAL moment cutoffs extend to EVERY point of the component by regular coefficient identities.')


def wick_engine(covariance):
    @lru_cache(maxsize=None)
    def wick(alpha):
        if not any(alpha):return sp.S.One
        if sum(alpha)%2:return sp.S.Zero
        i=next(i for i,k in enumerate(alpha) if k)
        beta=list(alpha);beta[i]-=1;out=sp.S.Zero
        for j,k in enumerate(beta):
            if k:
                gamma=list(beta);gamma[j]-=1;out+=k*covariance[i,j]*wick(tuple(gamma))
        return out
    return wick


def actual_periods(source,verified,n,phase,precision,field):
    wick=wick_engine(precision.inv());periods=[];certificates=[]
    for alpha,ann,cert in zip(source['amplitude_exponents'],source['coordinate_annihilators'],verified['scalar_indicial_certificates']):
        D,R=cert['infinity']['bound'],cert['origin']['bound']
        power={(0,)*n:field.one};coefficients=[]
        for k in range(max(D,R)+1):
            value=sum((c*field.convert(wick(tuple(a+b for a,b in zip(alpha,e)))) for e,c in power.items()),field.zero)
            coefficients.append(value/field.convert(math.factorial(k)))
            power=G.commutative_product(power,phase,field)
        polynomial={k:coefficients[k] for k in range(D+1) if coefficients[k]}
        residual={}
        for k,c in polynomial.items():
            G.sparse_add(residual,{e:c*v for e,v in G.operator_action(ann,k,field).items()},field)
        matches=not residual and all(not coefficients[k] for k in range(D+1,R+1))
        periods.append([[[k],field.encode(c)] for k,c in polynomial.items()] if matches else None)
        certificates.append(dict(original_Wick_coefficients=[field.encode(c) for c in coefficients],
            coefficient_bound=max(D,R),candidate_degree_bound=D,origin_jet_bound=R,
            candidate_operator_residual=[[e,field.encode(c)] for e,c in sorted(residual.items())],
            original_coefficients_beyond_degree_bound_zero=all(not coefficients[k] for k in range(D+1,R+1)),
            actual_period_polynomial=matches))
    safe=all(p is not None for p in periods)
    cutoff=None
    if safe:
        C=sp.Rational(verified['recurrence_slope']);E=sp.Rational(0);pureC=C
        for alpha,p in zip(source['amplitude_exponents'],periods):
            degree=sum(alpha);time_degree=max((e[0] for e,c in p),default=-1)
            E=max(E,time_degree-C*degree)
            if degree:pureC=max(pureC,sp.Rational(time_degree,degree))
        unit=next(p for alpha,p in zip(source['amplitude_exponents'],periods) if not any(alpha))
        pure=(len(unit)==1 and unit[0][0]==[0] and field.parse(unit[0][1])==field.one)
        cutoff=dict(slope=str(C),intercept=str(E),conditional_pure_linear_slope=str(C+E),
                    entire_component_pure=pure,pure_linear_slope=str(pureC) if pure else None,
                    convention='E(Q P^m)=0 for every integer m>slope*deg(Q)+intercept',
                    extends_to_ALL_specializations_of_component=True)
    return dict(safe=safe,actual_generator_periods=periods,actual_jet_certificates=certificates,cutoff=cutoff,
                original_bad_marker=next((list(a) for a,p in zip(source['amplitude_exponents'],periods) if p is None),None))


def verify_result(result):
    if result['status']!='DECIDED':raise ValueError('generic result is not decided')
    family,n,r,phase,precision,field,generic=normalize(result['input'])
    verified=verify_source(result['source'],n,generic,precision,field)
    selected=actual_periods(result['source'],verified,n,generic,precision,field)
    if verified!=result['verified_original_certificate'] or selected!=result['actual_period_selection']:
        raise ValueError('generic finite certificate mismatch')
    if result['source']['generic_prime_checked'] is not True:raise ValueError('missing generic prime contract')
    return True


def run(data,workdir,timeout=30):
    if type(timeout) is not int or timeout<0:raise ValueError('nonnegative timeout required')
    family,n,r,phase,precision,field,generic=normalize(data)
    out=Path(workdir).resolve();out.mkdir(parents=True,exist_ok=True)
    script=out/'generic.m2';script.write_text(make_script(data))
    try:
        process=subprocess.run(['M2','--no-readline','--silent','--stop',str(script)],capture_output=True,text=True,timeout=timeout or None)
    except subprocess.TimeoutExpired as exc:
        (out/'generic.stdout').write_bytes((exc.stdout or b'')+(exc.stderr or b''))
        return dict(status='UNKNOWN',reason='generic source CAS timeout',input=data)
    output=process.stdout+process.stderr;(out/'generic.stdout').write_text(output)
    lines=[l.split('|',1)[1] for l in output.splitlines() if l.startswith('RESULT_JSON|')]
    if process.returncode or len(lines)!=1 or 'GENERIC_COMPLETE' not in output:
        raise RuntimeError('generic source failed closed; see retained stdout')
    source=json.loads(lines[0]);verified=verify_source(source,n,generic,precision,field)
    selected=actual_periods(source,verified,n,generic,precision,field)
    result=dict(status='DECIDED',input=data,source=source,verified_original_certificate=verified,actual_period_selection=selected,
                core_sha256=sha(CORE),cas_files_sha256={p.name:sha(p) for p in [script,out/'generic.stdout']},
                assurance='finite original identities/IBP/indicial/actual-jet replay modulo supplied prime; primality is checked by Macaulay2,not independently re-proved here',
                actual_image_kernel_computed=False,source_relation_completeness_needed_for_decision=False)
    verify_result(json.loads(json.dumps(result)))
    return result


def verify_family_torsion(data,family):
    """Check the finite ordinary annihilation needed for numerical ring bounds."""
    n,r,phase,base,precision=F.normalize(data)
    xs=sp.symbols('c:'+str(r))
    poly=lambda ts:sum(sp.Rational(c)*sp.prod(x**k for x,k in zip(xs,e)) for e,c in ts)
    original=[sum(c*sp.prod(x**k for x,k in zip(xs,e)) for e,c in g.items()) for g in base]
    gb=sp.groebner(original,*xs,domain=sp.QQ)
    nf=lambda f:gb.reduce(sp.expand(f))[1]
    J=[poly(g) for g in family['safe_radical_ideal']]
    N=[poly(g) for g in family['polynomial_Hom_normalization_ideal_lift']]
    ngb=sp.groebner(N,*xs,domain=sp.QQ)
    if any(ngb.reduce(f)[1]!=0 for f in original):raise ValueError('normalization lift lost ordinary parameter relations')
    L=family['normalization_torsion_exponent']
    if type(L) is not int or L<1:raise ValueError('invalid normalization torsion exponent')
    power=[sp.S.One]
    for k in range(L):
        product=[nf(f*g) for f in power for g in J]
        pg=sp.groebner(original+product,*xs,domain=sp.QQ)
        power=[nf(g.as_expr()) for g in pg.polys if nf(g.as_expr())!=0]
        if not power:break
    if any(nf(a*b)!=0 for a in N for b in power):raise ValueError('claimed normalization torsion power is false')
    unit=ngb.reduce(sp.S.One)[1]==0
    if unit!=family['coefficient_ring_phase_safe']:raise ValueError('wrong normalization-one decision')
    return dict(ordinary_power_annihilation_verified=True,normalization_one_verified=True,
                input_parameter_ideal_retained_ordinary=True,
                exact_normalization_ideal_from_joint_family_solver=True)


def family_bounds(data,workdir,timeout=30,family_result=None):
    out=Path(workdir).resolve();out.mkdir(parents=True,exist_ok=True)
    family=F.run(data,out/'family',timeout) if family_result is None else family_result
    if family['input']!=data:raise ValueError('family result belongs to another input')
    if family['status']!='DECIDED':return dict(status='UNKNOWN',reason='family decomposition not decided',family_result=family)
    if family['core_sha256']!=sha(F.CORE):raise ValueError('family result uses another core revision')
    torsion_verified=verify_family_torsion(data,family)
    components=[]
    for i,prime in enumerate(family['safe_closed_components']):
        generic=run(dict(data,generic_prime=prime),out/('component_'+str(i)),timeout)
        if generic['status']!='DECIDED':return dict(status='UNKNOWN',reason='generic component certificate not decided',family_result=family,completed_components=components,pending_component=generic)
        if not generic['actual_period_selection']['safe']:raise ValueError('family safe component failed actual generic certificate')
        components.append(generic)
    C=max((sp.Rational(c['actual_period_selection']['cutoff']['slope']) for c in components),default=sp.S.Zero)
    E=max((sp.Rational(c['actual_period_selection']['cutoff']['intercept']) for c in components),default=sp.S.Zero)
    n,r,phase,base,precision=F.normalize(data)
    D=max((sum(e[:n]) for e in phase),default=0)
    L=family['normalization_torsion_exponent']
    thick_E=E+(C*D+1)*(r*L-1) if L>1 else E
    return dict(status='DECIDED',input=data,family_result=family,component_certificates=components,
        verified_ordinary_torsion=torsion_verified,
        all_safe_field_fibers_cutoff=dict(slope=str(C),intercept=str(E),conditional_pure_linear_slope=str(C+E)),
        presented_ring_normalization_cutoff=dict(slope=str(C),intercept=str(thick_E),conditional_pure_linear_slope=str(C+thick_E),
            parameter_count=r,ambient_spatial_phase_degree=D,torsion_exponent=L,
            applies_to='a*E(Q P^m) for EVERY normalization scalar a; a=1 only when the presented-ring phase is safe'),
        ordinary_cofinal_forward_index_formula='For thickness L>=2: ceil(max(C,E+(C*D+1)*(r*L-1))); for L=1 use ceil(max(C,E)).',
        whole_family_coverage_contract='safe component decomposition from the exact joint family solver; finite original generic certificates prove their bounds on every specialization',
        numerical_cutoffs_sufficient_not_least=True,ordinary_full_cutoff_ideal_producer_implemented=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--family',action='store_true');parser.add_argument('--family-result',type=Path)
    parser.add_argument('--workdir',type=Path);parser.add_argument('--timeout',type=int,default=30)
    args=parser.parse_args();data=json.loads(args.input.read_text());work=args.workdir or args.output.with_suffix('.cas')
    if args.family:
        previous=json.loads(args.family_result.read_text()) if args.family_result else None
        result=family_bounds(data,work,args.timeout,previous)
    else:result=run(data,work,args.timeout)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
