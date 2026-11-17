"""Ordinary Gaussian coefficient ideals with primary-thickening stopping tests.

An exact scalar rational-residue IBP search supplies the
general stopping oracle; finite nilpotence and isotropic support are exact
shortcuts. A timeout returns UNKNOWN, never equality of consecutive prefixes.
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
import time

import sympy as sp
import gaussian_family_producer as F
import gaussian_generic_cutoff_producer as G

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
dump=lambda d:json.dumps(d,indent=2)+'\n'


class BudgetExpired(Exception):pass


class Budget:
    def __init__(self,seconds):
        if type(seconds) is not int or seconds<0:raise ValueError('nonnegative integer timeout required')
        self.end=time.monotonic()+seconds if seconds else None
    def check(self):
        if self.end is not None and time.monotonic()>=self.end:raise BudgetExpired
    def remaining(self):
        self.check();return max(.01,self.end-time.monotonic()) if self.end is not None else None


def polynomial(terms,variables):
    return sp.expand(sum(sp.Rational(c)*sp.prod(x**k for x,k in zip(variables,e)) for e,c in terms))


def encode_polynomial(expr,variables):
    return [[list(e),str(c)] for e,c in sp.Poly(expr,*variables,domain=sp.QQ).terms() if c]


def primary_components(ideal,r,out,budget):
    """Trusted primary decomposition, with exact intersection and prime checks."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cs=['c'+str(i) for i in range(r)]
    expression=','.join(F.expression(F.terms(f,r),cs) for f in ideal) or '0_R'
    script=f'''needsPackage "PrimaryDecomposition";needsPackage "JSON";
R=QQ[{','.join(cs)}];I=ideal matrix(R,{{{{{expression}}}}});
qs=primaryDecomposition I;
if intersect qs != I then error "primary intersection mismatch";
encode=f->apply(listForm f,u->{{u#0,toString u#1}});
rows=apply(qs,Q->(
 P:=radical Q;
 if not isPrime P then error "nonprime radical";
 U:=support ((independentSets(P,Limit=>1))_0);
 if #U != dim(R/P) then error "wrong independent-set cardinality";
 new HashTable from {{"primary_ideal"=>apply(first entries gens Q,encode),
 "radical_ideal"=>apply(first entries gens P,encode),
 "independent_indices"=>select(toList(0..numgens R-1),i->member(R_i,U)),
 "dimension"=>dim(R/P),"prime_radical_checked"=>true}}
 ));
stdio << "PRIMARY_JSON|" << toJSON rows << endl;
print("PRIMARY_COMPLETE");exit 0;
'''
    (out/'primary.m2').write_text(script)
    try:proc=subprocess.run(['M2','--no-readline','--silent','--stop',str((out/'primary.m2').resolve())],capture_output=True,text=True,timeout=budget.remaining())
    except subprocess.TimeoutExpired as exc:
        (out/'primary.stdout').write_bytes((exc.stdout or b'')+(exc.stderr or b''));raise BudgetExpired
    output=proc.stdout+proc.stderr;(out/'primary.stdout').write_text(output)
    lines=[s.split('|',1)[1] for s in output.splitlines() if s.startswith('PRIMARY_JSON|')]
    if proc.returncode or len(lines)!=1 or 'PRIMARY_COMPLETE' not in output:raise RuntimeError('primary decomposition failed closed')
    return json.loads(lines[0])


def exponents(n,d):
    if n==0:
        yield ();return
    if d<0:return
    for k in range(d+1):
        for tail in exponents(n-1,d-k):yield (k,)+tail


class FiniteAlgebra:
    """QQ(c_independent)[c_dependent]/Q, extended by independent ray parameters."""
    def __init__(self,r,component,n):
        self.cs=sp.symbols('c:'+str(r));self.component=component
        indices=component['independent_indices']
        if len(indices)!=len(set(indices)) or any(type(i) is not int or i<0 or i>=r for i in indices):raise ValueError('invalid independent variables')
        if len(indices)!=component['dimension'] or component['prime_radical_checked'] is not True:raise ValueError('missing primary field contract')
        self.independent=tuple(self.cs[i] for i in indices)
        self.dependent=tuple(c for i,c in enumerate(self.cs) if i not in indices)
        self.rays=sp.symbols('a:'+str(n+1));self.parameters=self.independent+self.rays
        self.K=sp.QQ.frac_field(*self.parameters);self.zero=self.K.zero;self.one=self.K.one
        eqs=[polynomial(f,self.cs) for f in component['primary_ideal']]
        if self.dependent:
            self.gb=sp.groebner(eqs,*self.dependent,domain=self.K)
            if not self.gb.is_zero_dimensional:raise ValueError('localized primary algebra is not finite')
            leading=[tuple(g.LM(order=self.gb.order).exponents) for g in self.gb.polys]
            if (0,)*len(self.dependent) in leading:raise ValueError('localized primary algebra vanished')
            limits=[]
            for j in range(len(self.dependent)):
                limits.append(min(e[j] for e in leading if e[j] and all(k==0 for i,k in enumerate(e) if i!=j)))
            self.basis=[tuple(e) for e in itertools.product(*(range(k) for k in limits))
                        if not any(all(a>=b for a,b in zip(e,lm)) for lm in leading)]
        else:
            self.gb=None;self.basis=[()]
            if any(eqs):raise ValueError('independent field quotient is nonzero')
        self.rank=len(self.basis);self.positions={e:i for i,e in enumerate(self.basis)}
        self.unit=self.vector(1)
        self.table=[[self.vector(sp.prod(x**(a+b) for x,a,b in zip(self.dependent,e,f))) for f in self.basis] for e in self.basis]
        # Commutative multiplication tables and the original ordinary relations.
        for f in eqs:
            if any(self.vector(f)):raise ValueError('lost primary relation')
        for i in range(self.rank):
            if self.mul(self.unit,self.basis_vector(i))!=self.basis_vector(i):raise ValueError('wrong algebra unit')
    def basis_vector(self,i):return tuple(self.one if j==i else self.zero for j in range(self.rank))
    def vector(self,expr):
        if not self.dependent:return (self.K.convert(expr),)
        remainder=self.gb.reduce(sp.expand(expr))[1]
        poly=sp.Poly(remainder,*self.dependent,domain=self.K);out=[self.zero]*self.rank
        for e,c in poly.rep.to_dict().items():out[self.positions[e]]=c
        return tuple(out)
    def mul(self,x,y):
        out=[self.zero]*self.rank
        for i,a in enumerate(x):
            if not a:continue
            for j,b in enumerate(y):
                if not b:continue
                for k,c in enumerate(self.table[i][j]):
                    if c:out[k]+=a*b*c
        return tuple(out)
    def add(self,x,y):return tuple(a+b for a,b in zip(x,y))
    def scale(self,x,c):return tuple(a*c for a in x)
    def encode_coefficient(self,c):
        expr=self.K.to_sympy(c);num,den=sp.fraction(expr)
        return dict(numerator=encode_polynomial(num,self.parameters),denominator=encode_polynomial(den,self.parameters))
    def parse_coefficient(self,c):
        return self.K.convert(polynomial(c['numerator'],self.parameters)/polynomial(c['denominator'],self.parameters))
    def encode_vector(self,v):return [self.encode_coefficient(c) for c in v]
    def manifest(self):
        return dict(component=self.component,scalar_parameters=[str(x) for x in self.parameters],
                    dependent_variables=[str(x) for x in self.dependent],standard_monomials=[list(e) for e in self.basis],
                    rank=self.rank,unit=self.encode_vector(self.unit),multiplication_table=[[self.encode_vector(v) for v in row] for row in self.table],
                    injectivity_contract='Q is primary and nonzero polynomials in its maximal independent variables avoid radical(Q); localization is injective. Primary status is trusted CAS.')


def phase_in_algebra(data,A):
    n,r,phase,base,precision=F.normalize(data);out={}
    for e,c in phase.items():
        a=tuple(e[:n]);v=A.vector(c*sp.prod(x**k for x,k in zip(A.cs,e[n:])))
        out[a]=A.add(out.get(a,(A.zero,)*A.rank),v)
    return {e:v for e,v in out.items() if any(v)}


def multiply_polynomials(P,Q,A):
    out={}
    for a,x in P.items():
        for b,y in Q.items():
            e=tuple(i+j for i,j in zip(a,b));v=A.mul(x,y)
            out[e]=A.add(out.get(e,(A.zero,)*A.rank),v)
    return {e:v for e,v in out.items() if any(v)}


class Moments:
    def __init__(self,data,A):
        self.data=data;self.A=A;self.n=data['dimension'];self.phase=phase_in_algebra(data,A)
        self.powers=[{(0,)*self.n:A.unit}]
        self.precision=F.normalize(data)[4];self.wick=G.wick_engine(self.precision.inv())
    def power(self,m):
        while len(self.powers)<=m:self.powers.append(multiply_polynomials(self.powers[-1],self.phase,self.A))
        return self.powers[m]
    def moment(self,alpha,m):
        out=(self.A.zero,)*self.A.rank
        for e,v in self.power(m).items():out=self.A.add(out,self.A.scale(v,self.A.K.convert(self.wick(tuple(a+b for a,b in zip(alpha,e))))))
        return out


def forbidden_indices(n,N,E,total):
    for m in range(1,total+1):
        for alpha in exponents(n,total-m):
            if m>N*sum(alpha)+E:yield tuple(alpha),m


def first_nonzero(moments,indices):
    for alpha,m in indices:
        value=moments.moment(alpha,m)
        if any(value):return dict(marker=list(alpha),power=m,value=moments.A.encode_vector(value))
    return None


def indicial(operator,A):
    w=min(a-b for (a,b),c in operator.items() if c);z=sp.Symbol('z');poly={}
    for (a,b),c in operator.items():
        if a-b!=w:continue
        for (j,),v in sp.Poly(sp.ff(z,b),z).terms():poly[j]=poly.get(j,A.zero)+c*A.K.convert(v)
    poly={j:c for j,c in poly.items() if c}
    if not poly:raise ValueError('zero scalar indicial polynomial')
    common=sp.S.One
    fractions=[]
    for j,c in poly.items():
        num,den=sp.fraction(A.K.to_sympy(c));fractions.append((j,num,den));common*=den
    polys={}
    for i,(j,num,den) in enumerate(fractions):
        q=sp.Poly(sp.expand(num*sp.prod(d for k,(_,_,d) in enumerate(fractions) if k!=i)),*A.parameters,domain=sp.QQ)
        for e,c in q.terms():polys[e]=polys.get(e,0)+c*z**j
    content=sp.Poly(0,z,domain=sp.QQ)
    for q in polys.values():content=sp.gcd(content,sp.Poly(q,z,domain=sp.QQ))
    if content.is_zero:raise ValueError('zero rational indicial content')
    content=content.monic();roots=sorted(int(x) for x in content.ground_roots() if x.is_Integer and x>=0)
    return dict(weight=w,nonnegative_integer_roots=roots,bound=max(roots,default=-1),
                coefficients=[[j,A.encode_coefficient(c)] for j,c in sorted(poly.items())],
                rational_content=[[int(e[0]),str(c)] for e,c in content.terms()])


class LinearWitness:
    """Sparse exact column reduction, retaining a dependence in original columns."""
    def __init__(self,K,budget):self.K=K;self.pivots={};self.budget=budget
    def add(self,column,tag):
        v=dict(column);w={tag:self.K.one}
        while v:
            self.budget.check();row=min(v);c=v[row]
            if row not in self.pivots:
                self.pivots[row]=({e:a/c for e,a in v.items()},{e:a/c for e,a in w.items()});return None
            q,certificate=self.pivots[row]
            for e,a in q.items():
                x=v.get(e,self.K.zero)-c*a
                if x:v[e]=x
                elif e in v:del v[e]
            for e,a in certificate.items():
                x=w.get(e,self.K.zero)-c*a
                if x:w[e]=x
                elif e in w:del w[e]
        return w


def make_ibp_problem(data,A,N=0,E=0,pure=False):
    n=data['dimension'];Xs=sp.symbols('x:'+str(n));u=sp.Symbol('u');vs=sp.symbols('v:'+str(n+1))
    variables=(u,)+Xs if pure else (u,)+Xs+vs
    K=A.K;poly=lambda x:sp.Poly(x,*variables,domain=K)
    zero=poly(0);one=poly(1)
    phase=phase_in_algebra(data,A)
    P=[zero for _ in range(A.rank)]
    for e,v in phase.items():
        mon=sp.prod(x**k for x,k in zip(Xs,e))
        for j,c in enumerate(v):P[j]+=poly(K.to_sympy(c)*mon)
    precision=F.normalize(data)[4]
    directions=list(Xs) if pure else list(Xs)+list(vs)
    connections=[]
    for var in directions:
        columns=[]
        for i in range(A.rank):
            if var in Xs:
                j=Xs.index(var);multipliers=[q.diff(var) for q in P]
                factor=u if pure else vs[0]
                scalar=-sum(precision[j,k]*Xs[k] for k in range(n))+(0 if pure else vs[j+1])
            elif var==vs[0]:multipliers=P;factor=1;scalar=0
            else:multipliers=[zero]*A.rank;factor=0;scalar=Xs[vs.index(var)-1]
            result=[zero for _ in range(A.rank)]
            for j,q in enumerate(multipliers):
                for k,c in enumerate(A.table[j][i]):
                    if c:result[k]+=poly(K.to_sympy(c)*factor)*q
            result[i]+=poly(scalar);columns.append(result)
        connections.append(columns)
    if pure:
        f=one;hg=0;gnum=one
    else:
        a=A.rays;t=a[0]*u/vs[0]
        rational=t**(E+1)/(1-t)
        for j in range(n):rational/=1-a[j+1]*u/vs[j+1]*t**N
        rational/=sp.prod(vs)
        f=poly(sp.prod(vs)*(vs[0]-a[0]*u)*sp.prod(vs[j+1]*vs[0]**N-a[j+1]*a[0]**N*u**(N+1) for j in range(n)))
        hg=E+2;gnum=poly(sp.cancel(rational*f.as_expr()**hg))
    return dict(variables=variables,u=u,Xs=Xs,directions=directions,connections=connections,P=P,f=f,hg=hg,gnum=gnum,pure=pure)


def ibp_column(problem,A,box,tag):
    K=A.K;zs=problem['variables'];f=problem['f'];hg=problem['hg'];pure=problem['pure']
    poly=lambda e:sp.Poly(e,*zs,domain=K);zero=poly(0)
    degree=box if pure else box*(int(f.total_degree())+1)
    if type(box) is not int or box<1:raise ValueError('invalid scalar IBP box')
    if tag[0]=='P':
        direction,i=tag[1:3];e=tag[3:]
        if len(e)!=len(zs) or any(type(k) is not int or k<0 for k in e) or sum(e)>degree:
            raise ValueError('primitive outside declared finite box')
        if type(direction) is not int or type(i) is not int or not 0<=direction<len(problem['directions']) or not 0<=i<A.rank:
            raise ValueError('invalid primitive coordinate')
        mon=poly(sp.prod(x**k for x,k in zip(zs,e)));var=problem['directions'][direction]
        post=poly(1) if pure else f**(hg-1)
        diagonal=mon.diff(var) if pure else f*mon.diff(var)-box*f.diff(var)*mon
        vector=[post*((mon if pure else f*mon)*q) for q in problem['connections'][direction][i]]
        vector[i]+=post*diagonal
        return vector
    if tag[0]!='L' or len(tag)!=3:raise ValueError('unknown scalar identity column')
    k,j=tag[1:]
    if any(type(v) is not int or v<0 or v>box for v in [k,j]):raise ValueError('operator outside declared finite box')
    if pure:
        if 'phase_powers' not in problem:
            phase={}
            for i,q in enumerate(problem['P']):
                for e,c in q.rep.to_dict().items():
                    v=list(phase.get(e,(A.zero,)*A.rank));v[i]=c;phase[e]=tuple(v)
            problem['phase_powers']=[{(0,)*len(zs):A.unit}];problem['phase_sparse']=phase
        powers=problem['phase_powers']
        while len(powers)<=j:powers.append(multiply_polynomials(powers[-1],problem['phase_sparse'],A))
        vector=[zero for _ in range(A.rank)]
        for e,v in powers[j].items():
            mon=poly(sp.prod(x**d for x,d in zip(zs,e))*problem['u']**k)
            for i,c in enumerate(v):
                if c:vector[i]+=poly(K.to_sympy(c))*mon
        return vector
    if 'rational_derivative_numerators' not in problem:problem['rational_derivative_numerators']=[problem['gnum']]
    numerators=problem['rational_derivative_numerators']
    while len(numerators)<=j:
        ell=len(numerators)-1;q=numerators[-1]
        numerators.append(f*q.diff(problem['u'])-(hg+ell)*f.diff(problem['u'])*q)
    value=numerators[j]*f**(box-j)*poly(problem['u']**k)
    return [value*poly(K.to_sympy(c)) for c in A.unit]


def ibp_columns(problem,A,box):
    degree=box if problem['pure'] else box*(int(problem['f'].total_degree())+1)
    for e in exponents(len(problem['variables']),degree):
        for direction in range(len(problem['directions'])):
            for i in range(A.rank):
                tag=('P',direction,i)+e
                yield tag,ibp_column(problem,A,box,tag)
    for j in range(box+1):
        for k in range(box+1):
            tag=('L',k,j)
            yield tag,ibp_column(problem,A,box,tag)


def ibp_search(data,A,budget,N=0,E=0,pure=False):
    problem=make_ibp_problem(data,A,N,E,pure);box=0
    while True:
        budget.check();box+=1;linear=LinearWitness(A.K,budget)
        for tag,vector in ibp_columns(problem,A,box):
            budget.check();column={(i,)+tuple(e):c for i,q in enumerate(vector) for e,c in q.rep.to_dict().items() if c}
            witness=linear.add(column,tag)
            if witness is not None and any(t[0]=='L' for t in witness):
                first=next(c for t,c in witness.items() if t[0]=='L');witness={t:c/first for t,c in witness.items()}
                operator={(t[1],t[2]):c for t,c in witness.items() if t[0]=='L'}
                certificate=dict(kind='polynomial_Gaussian_IBP' if pure else 'rational_cone_residue_IBP',box=box,
                    directions=[str(x) for x in problem['directions']],variables=[str(x) for x in problem['variables']],
                    primitive_denominator_power=0 if pure else box,
                    denominator_polynomial=[[list(e),A.encode_coefficient(c)] for e,c in problem['f'].rep.to_dict().items()],
                    identity_terms=[[list(t),A.encode_coefficient(c)] for t,c in witness.items()],
                    operator=[[list(e),A.encode_coefficient(c)] for e,c in sorted(operator.items())])
                verify_ibp(data,A,certificate,N,E,pure)
                return operator,certificate


def verify_ibp(data,A,certificate,N=0,E=0,pure=False):
    problem=make_ibp_problem(data,A,N,E,pure);wanted={tuple(t):A.parse_coefficient(c) for t,c in certificate['identity_terms']}
    total={}
    # Retained identities are replayed against the original rational connection,
    # using only the columns that occur, not a rank/completeness assertion.
    if len(wanted)!=len(certificate['identity_terms']):raise ValueError('duplicate retained identity column')
    if certificate['variables']!=[str(x) for x in problem['variables']] or certificate['directions']!=[str(x) for x in problem['directions']]:
        raise ValueError('wrong original connection coordinates')
    expected_denominator=[[list(e),A.encode_coefficient(c)] for e,c in problem['f'].rep.to_dict().items()]
    if certificate['denominator_polynomial']!=expected_denominator:raise ValueError('wrong rational localization polynomial')
    if certificate['primitive_denominator_power']!=(0 if pure else certificate['box']):raise ValueError('wrong primitive localization power')
    for tag,scalar in list(wanted.items()):
        vector=ibp_column(problem,A,certificate['box'],tag);del wanted[tag]
        for i,q in enumerate(vector):
            for e,c in q.rep.to_dict().items():
                key=(i,)+tuple(e);v=total.get(key,A.zero)+scalar*c
                if v:total[key]=v
                elif key in total:del total[key]
    if wanted or total:raise ValueError('invalid original ordinary-ideal IBP identity')
    op={tuple(e):A.parse_coefficient(c) for e,c in certificate['operator']}
    from_terms={tuple(t[1:]):A.parse_coefficient(c) for t,c in certificate['identity_terms'] if t[0]=='L'}
    if not op or op!=from_terms:raise ValueError('wrong retained scalar operator')
    return True


def subtract_unit_operator(operator,A):
    p0={a:c for (a,b),c in operator.items() if b==0}
    if not p0:return operator
    left={}
    for a,c in p0.items():
        left[(a,1)]=left.get((a,1),A.zero)+c
        if a:left[(a-1,0)]=left.get((a-1,0),A.zero)-a*c
    out={}
    for (a,b),c in left.items():
        for (i,j),d in operator.items():
            for ell in range(min(b,i)+1):
                e=(a+i-ell,b+j-ell);v=c*d*math.comb(b,ell)*(math.factorial(i)//math.factorial(i-ell))
                out[e]=out.get(e,A.zero)+v
    return {e:c for e,c in out.items() if c}


def component_test(data,component,N,E,pure,budget,force_ibp=False):
    A=FiniteAlgebra(data['parameter_count'],component,data['dimension']);moments=Moments(data,A);n=data['dimension']
    result=dict(algebra=A.manifest())
    if not force_ibp:
        indices=[((0,)*n,m) for m in range(1,5)] if pure else forbidden_indices(n,N,E,4)
        bad=first_nonzero(moments,indices)
        if bad:return dict(result,zero=False,method='original_nonzero_coefficient',witness=bad)
        if not moments.power(A.rank):
            indices=[((0,)*n,m) for m in range(1,A.rank)] if pure else (
                (alpha,m) for m in range(1,A.rank) for alpha in exponents(n,max(-1,(m-E-1)//N)) if m>N*sum(alpha)+E)
            bad=first_nonzero(moments,indices)
            return dict(result,zero=bad is None,method='finite_nilpotent_phase',nilpotence_exponent=A.rank,witness=bad)
        active=[i for i in range(n) if any(e[i] for e in moments.phase)]
        covariance=moments.precision.inv()
        dmin=min((sum(e) for e in moments.phase),default=0)
        if dmin>0 and all(covariance[i,j]==0 for i in active for j in active) and (pure or N*dmin>=1):
            return dict(result,zero=True,method='totally_isotropic_support',active_variables=active,minimum_phase_degree=dmin)
    operator,certificate=ibp_search(data,A,budget,N,E,pure)
    zero_operator=subtract_unit_operator(operator,A) if pure else operator
    origin=indicial(zero_operator,A);bound=origin['bound']
    indices=[((0,)*n,m) for m in range(1,bound+1)] if pure else forbidden_indices(n,N,E,bound)
    checked=[];bad=None
    for alpha,m in indices:
        value=moments.moment(alpha,m);checked.append(dict(marker=list(alpha),power=m,value=A.encode_vector(value)))
        if any(value) and bad is None:bad=checked[-1]
    return dict(result,zero=bad is None,method='scalar_IBP_origin_jet',IBP_certificate=certificate,
        zero_series_operator=[[list(e),A.encode_coefficient(c)] for e,c in sorted(zero_operator.items())],
        origin_indicial=origin,checked_original_moments=checked,witness=bad)


def original_moment(data,alpha,m):
    n,r,phase,base,precision=F.normalize(data);xs=sp.symbols('x:'+str(n));cs=sp.symbols('c:'+str(r))
    P=sum(c*sp.prod(x**k for x,k in zip(xs+cs,e)) for e,c in phase.items())
    wick=G.wick_engine(precision.inv());power=sp.Poly(sp.expand(P**m),*xs)
    return sp.expand(sum(c*wick(tuple(a+b for a,b in zip(alpha,e))) for e,c in power.terms()))


def normalized_input(data,N,E,pure):
    n,r,phase,base,precision=F.normalize(data)
    if type(N) is not int or N<0 or type(E) is not int or E<0:raise ValueError('nonnegative integer cutoff slope and intercept required')
    cs=sp.symbols('c:'+str(r));original=[polynomial(g,cs) for g in data.get('parameter_ideal',[])]
    return n,r,cs,original


def constant_phase(data):
    n,r,phase,base,precision=F.normalize(data)
    if any(any(e[:n]) for e in phase):return None
    cs=sp.symbols('c:'+str(r))
    return sp.expand(sum(c*sp.prod(x**k for x,k in zip(cs,e[n:])) for e,c in phase.items()))


def zero_slope_certificate(data,E):
    n,r,cs,original=normalized_input(data,0,E,False)
    xs=sp.symbols('x:'+str(n));phase=F.normalize(data)[2];precision=F.normalize(data)[4]
    P=sum(c*sp.prod(x**k for x,k in zip(xs+cs,e)) for e,c in phase.items())
    power=sp.Poly(sp.expand(P**(E+1)),*xs)
    degree=max((sum(e) for e,c in power.terms() if c),default=0)
    basis=list(exponents(n,degree));wick=G.wick_engine(precision.inv())
    gram=sp.Matrix([[wick(tuple(i+j for i,j in zip(a,b))) for b in basis] for a in basis])
    inverse=gram.inv()
    coefficients=sp.Matrix([power.coeff_monomial(e) for e in basis]);moments=gram*coefficients
    if inverse*moments!=coefficients:raise ValueError('invalid Gaussian pairing inversion')
    return dict(phase_power=E+1,spatial_degree=degree,marker_exponents=[list(e) for e in basis],
                Wick_Gram=[[str(c) for c in row] for row in gram.tolist()],
                coefficient_from_moment_matrix=[[str(c) for c in row] for row in inverse.tolist()],
                original_moments=[dict(marker=list(e),power=E+1,polynomial=encode_polynomial(sp.expand(m),cs)) for e,m in zip(basis,moments)],
                phase_power_coefficients=[encode_polynomial(c,cs) for c in coefficients])


def zero_test(data,ideal,out,N=0,E=0,pure=False,timeout=30,budget=None):
    n,r,cs,original=normalized_input(data,N,E,pure);budget=budget or Budget(timeout)
    polys=[polynomial(g,cs) for g in ideal];gb=sp.groebner(polys,*cs,domain=sp.QQ)
    if any(gb.reduce(f)[1]!=0 for f in original):raise ValueError('proposed ideal omits original ordinary constraints')
    try:
        if gb.reduce(sp.S.One)[1]==0:return dict(status='DECIDED',zero=True,method='unit_ideal',components=[])
        constant=constant_phase(data)
        if constant is not None:
            m=1 if pure else E+1;value=sp.expand(constant**m);zero=gb.reduce(value)[1]==0
            return dict(status='DECIDED',zero=zero,method='constant_phase',components=[],
                witness=None if zero else dict(marker=[0]*n,power=m,polynomial=encode_polynomial(value,cs)))
        # Nondegeneracy of the Gaussian polynomial pairing makes slope zero
        # equivalent to P^(E+1)=0 coefficientwise, over EVERY QQ-algebra.
        if not pure and N==0:
            xs=sp.symbols('x:'+str(n));phase=F.normalize(data)[2]
            P=sum(c*sp.prod(x**k for x,k in zip(xs+cs,e)) for e,c in phase.items())
            power=sp.Poly(sp.expand(P**(E+1)),*xs);coefficients=power.coeffs()
            zero=all(gb.reduce(sp.expand(c))[1]==0 for c in coefficients);witness=None
            degree=max((sum(e) for e,c in power.terms() if c),default=0)
            if not zero:
                for alpha in exponents(n,degree):
                    value=original_moment(data,alpha,E+1)
                    if gb.reduce(value)[1]!=0:
                        witness=dict(marker=list(alpha),power=E+1,polynomial=encode_polynomial(value,cs));break
                if witness is None:raise ValueError('Gaussian pairing failed to detect a nonzero coefficient')
            return dict(status='DECIDED',zero=zero,
                method='nondegenerate_Gaussian_pairing',phase_power=E+1,
                coefficient_remainders=[encode_polynomial(gb.reduce(sp.expand(c))[1],cs) for c in coefficients],
                witness=witness,components=[])
        components=primary_components(ideal,r,Path(out)/'primary',budget);tested=[]
        for component in components:
            budget.check();result=component_test(data,component,N,E,pure,budget);tested.append(result)
            if not result['zero']:
                witness=result['witness'];alpha=witness['marker'];m=witness['power']
                original_polynomial=original_moment(data,alpha,m)
                if gb.reduce(original_polynomial)[1]==0:raise ValueError('primary witness belongs to proposed ideal')
                return dict(status='DECIDED',zero=False,components=tested,primary_components=components,
                    witness=dict(marker=alpha,power=m,polynomial=encode_polynomial(original_polynomial,cs)))
        return dict(status='DECIDED',zero=True,components=tested,primary_components=components)
    except BudgetExpired:return dict(status='UNKNOWN',reason='ordinary coefficient oracle timeout')


def run(data,out,N=0,E=0,pure=False,timeout=30):
    n,r,cs,original=normalized_input(data,N,E,pure);out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    budget=Budget(timeout);ideal=list(original);trace=[];witnesses=[]
    constant=constant_phase(data)
    if constant is not None:
        m=1 if pure else E+1;value=sp.expand(constant**m)
        gb=sp.groebner(ideal+[value],*cs,domain=sp.QQ)
        return dict(status='DECIDED',input=data,cutoff_slope=N,cutoff_intercept=E,pure_moments_only=pure,
                    ideal_lift=[encode_polynomial(g.as_expr(),cs) for g in gb.polys],method='constant_phase',
                    original_generating_moment=dict(marker=[0]*n,power=m,polynomial=encode_polynomial(value,cs)),
                    retains_ordinary_parameter_ideal=True)
    if not pure and N==0:
        xs=sp.symbols('x:'+str(n));phase=F.normalize(data)[2]
        P=sum(c*sp.prod(x**k for x,k in zip(xs+cs,e)) for e,c in phase.items())
        ideal+=sp.Poly(sp.expand(P**(E+1)),*xs).coeffs()
        gb=sp.groebner(ideal,*cs,domain=sp.QQ)
        return dict(status='DECIDED',input=data,cutoff_slope=N,cutoff_intercept=E,pure_moments_only=pure,
            ideal_lift=[encode_polynomial(g.as_expr(),cs) for g in gb.polys],method='nondegenerate_Gaussian_pairing',
            phase_power=E+1,pairing_certificate=zero_slope_certificate(data,E),retains_ordinary_parameter_ideal=True)
    try:
        while True:
            budget.check();gb=sp.groebner(ideal,*cs,domain=sp.QQ)
            encoded=[encode_polynomial(g.as_expr(),cs) for g in gb.polys]
            oracle=zero_test(data,encoded,out/('step_'+str(len(trace))),N,E,pure,budget=budget)
            trace.append(dict(proposed_ideal=encoded,oracle=oracle))
            if oracle['status']!='DECIDED':raise BudgetExpired
            if oracle['zero']:
                return dict(status='DECIDED',input=data,cutoff_slope=N,cutoff_intercept=E,pure_moments_only=pure,
                    ideal_lift=encoded,original_generating_moments=witnesses,stopping_trace=trace,
                    retains_ordinary_parameter_ideal=True,actual_image_kernel_computed=False,
                    assurance='Primary decomposition/injectivity contract uses trusted CAS; original moment witnesses and scalar IBP identities are finite exact certificates. No universal independent assurance.')
            witness=oracle['witness'];f=polynomial(witness['polynomial'],cs)
            if gb.reduce(f)[1]==0:raise ValueError('ordinary ideal ascent is not strict')
            ideal.append(f);witnesses.append(witness)
    except BudgetExpired:
        return dict(status='UNKNOWN',reason='ordinary coefficient ideal timeout',input=data,cutoff_slope=N,cutoff_intercept=E,
                    pure_moments_only=pure,completed_strict_ascent=trace,original_generating_moments=witnesses)


def verify_component(data,result,N,E,pure):
    A=FiniteAlgebra(data['parameter_count'],result['algebra']['component'],data['dimension'])
    if A.manifest()!=result['algebra']:raise ValueError('wrong finite primary algebra')
    moments=Moments(data,A);n=data['dimension'];method=result['method']
    if method=='original_nonzero_coefficient':
        w=result['witness'];value=moments.moment(tuple(w['marker']),w['power'])
        if not any(value) or A.encode_vector(value)!=w['value'] or result['zero']:raise ValueError('false original primary witness')
    elif method=='finite_nilpotent_phase':
        if result['nilpotence_exponent']!=A.rank or moments.power(A.rank):raise ValueError('false polynomial nilpotence')
        indices=[((0,)*n,m) for m in range(1,A.rank)] if pure else (
            (alpha,m) for m in range(1,A.rank) for alpha in exponents(n,max(-1,(m-E-1)//N)) if m>N*sum(alpha)+E)
        bad=first_nonzero(moments,indices)
        if result['witness']!=bad or result['zero']!=(bad is None):raise ValueError('wrong finite nilpotent tail')
    elif method=='totally_isotropic_support':
        active=[i for i in range(n) if any(e[i] for e in moments.phase)]
        dmin=min((sum(e) for e in moments.phase),default=0);covariance=moments.precision.inv()
        if not dmin>0 or not all(covariance[i,j]==0 for i in active for j in active) or not (pure or N*dmin>=1):raise ValueError('invalid isotropic support bound')
        if not result['zero'] or active!=result['active_variables'] or dmin!=result['minimum_phase_degree']:raise ValueError('wrong isotropic manifest')
    elif method=='scalar_IBP_origin_jet':
        cert=result['IBP_certificate'];verify_ibp(data,A,cert,N,E,pure)
        operator={tuple(e):A.parse_coefficient(c) for e,c in cert['operator']}
        zero_operator=subtract_unit_operator(operator,A) if pure else operator
        if [[list(e),A.encode_coefficient(c)] for e,c in sorted(zero_operator.items())]!=result['zero_series_operator']:raise ValueError('wrong scalar zero-series operator')
        origin=indicial(zero_operator,A)
        if origin!=result['origin_indicial']:raise ValueError('wrong scalar origin threshold')
        indices=[((0,)*n,m) for m in range(1,origin['bound']+1)] if pure else forbidden_indices(n,N,E,origin['bound'])
        checked=[dict(marker=list(a),power=m,value=A.encode_vector(moments.moment(a,m))) for a,m in indices]
        bad=next((w for w in checked if any(A.parse_coefficient(c) for c in w['value'])),None)
        if checked!=result['checked_original_moments'] or bad!=result['witness'] or result['zero']!=(bad is None):raise ValueError('wrong original scalar jets')
    else:raise ValueError('unknown primary stopping certificate')
    if not result['zero']:
        w=result['witness']
        if pure:
            if any(w['marker']) or w['power']<1:raise ValueError('not a pure witness')
        elif w['power']<=N*sum(w['marker'])+E:raise ValueError('witness is not forbidden')
    return True


def verify_result(result):
    if result['status']!='DECIDED':raise ValueError('ordinary ideal not decided')
    data=result['input'];N=result['cutoff_slope'];E=result['cutoff_intercept'];pure=result['pure_moments_only']
    n,r,cs,original=normalized_input(data,N,E,pure)
    expected=list(original)
    def basis(fs):return [encode_polynomial(f.as_expr(),cs) for f in sp.groebner(fs,*cs,domain=sp.QQ).polys]
    if result.get('method')=='constant_phase':
        constant=constant_phase(data)
        if constant is None:raise ValueError('phase is not constant')
        m=1 if pure else E+1;value=sp.expand(constant**m)
        if result['original_generating_moment']!=dict(marker=[0]*n,power=m,polynomial=encode_polynomial(value,cs)):raise ValueError('wrong constant-phase generator')
        expected.append(value)
    elif result.get('method')=='nondegenerate_Gaussian_pairing':
        if pure or N!=0:raise ValueError('invalid slope-zero method')
        cert=zero_slope_certificate(data,E)
        if cert!=result['pairing_certificate']:raise ValueError('wrong Gaussian pairing certificate')
        expected += [polynomial(f,cs) for f in cert['phase_power_coefficients']]
    else:
        witnesses=result['original_generating_moments'];traces=result['stopping_trace']
        if len(traces)!=len(witnesses)+1:raise ValueError('incomplete ordinary ideal trace')
        for i,step in enumerate(traces):
            if basis(expected)!=step['proposed_ideal']:raise ValueError('wrong ordinary ascending ideal')
            oracle=step['oracle']
            if oracle['status']!='DECIDED' or oracle['zero']!=(i==len(witnesses)):raise ValueError('invalid ideal stopping decision')
            verify_zero_test(data,step['proposed_ideal'],oracle,N,E,pure)
            if oracle.get('method')=='unit_ideal':
                if basis(expected)!=[[[[0]*r,'1']]]:raise ValueError('false unit ideal')
            else:
                components=oracle['primary_components'];tested=oracle['components']
                if len(tested)>len(components) or (oracle['zero'] and len(tested)!=len(components)):raise ValueError('missing primary component')
                for j,component in enumerate(tested):
                    if component['algebra']['component']!=components[j]:raise ValueError('wrong primary component order')
                    verify_component(data,component,N,E,pure)
                if oracle['zero']!=all(c['zero'] for c in tested):raise ValueError('wrong primary intersection decision')
            if i<len(witnesses):
                w=witnesses[i]
                if w!=oracle['witness']:raise ValueError('lost original generating witness')
                value=original_moment(data,tuple(w['marker']),w['power'])
                if encode_polynomial(value,cs)!=w['polynomial']:raise ValueError('false original generating moment')
                gb=sp.groebner(expected,*cs,domain=sp.QQ)
                if gb.reduce(value)[1]==0:raise ValueError('non-strict ordinary ascent')
                expected.append(value)
    if basis(expected)!=result['ideal_lift']:raise ValueError('wrong full ordinary ideal')
    return True


def intersection(ideals,cs):
    z=sp.Symbol('intersection_aux');current=[sp.S.One]
    for following in ideals:
        gb=sp.groebner([z*f for f in current]+[(1-z)*g for g in following],z,*cs,domain=sp.QQ)
        current=[f.as_expr() for f in gb.polys if not f.as_expr().has(z)]
    return [g.as_expr() for g in sp.groebner(current,*cs,domain=sp.QQ).polys]


def verify_zero_test(data,ideal,result,N=0,E=0,pure=False):
    if result['status']!='DECIDED':raise ValueError('coefficient zero test not decided')
    n,r,cs,original=normalized_input(data,N,E,pure)
    polys=[polynomial(f,cs) for f in ideal];gb=sp.groebner(polys,*cs,domain=sp.QQ)
    if any(gb.reduce(f)[1]!=0 for f in original):raise ValueError('zero test lost ordinary base constraints')
    method=result.get('method')
    if method=='unit_ideal':
        if gb.reduce(sp.S.One)[1]!=0 or not result['zero']:raise ValueError('false unit-ideal zero test')
    elif method=='constant_phase':
        c=constant_phase(data)
        if c is None:raise ValueError('not a constant phase')
        m=1 if pure else E+1;value=sp.expand(c**m);zero=gb.reduce(value)[1]==0
        witness=None if zero else dict(marker=[0]*n,power=m,polynomial=encode_polynomial(value,cs))
        if result['zero']!=zero or result['witness']!=witness:raise ValueError('false constant-phase zero test')
    elif method=='nondegenerate_Gaussian_pairing':
        if N!=0 or pure:raise ValueError('invalid slope-zero pairing test')
        xs=sp.symbols('x:'+str(n));phase=F.normalize(data)[2]
        P=sum(c*sp.prod(x**k for x,k in zip(xs+cs,e)) for e,c in phase.items())
        coefficients=sp.Poly(sp.expand(P**(E+1)),*xs).coeffs()
        remainders=[encode_polynomial(gb.reduce(sp.expand(c))[1],cs) for c in coefficients]
        if remainders!=result['coefficient_remainders'] or result['zero']!=all(not f for f in remainders):raise ValueError('false Gaussian coefficient-pairing test')
    else:
        components=result['primary_components'];tested=result['components']
        if not components or len(tested)>len(components) or (result['zero'] and len(tested)!=len(components)):raise ValueError('missing primary stopping component')
        ideals=[[polynomial(f,cs) for f in c['primary_ideal']] for c in components]
        if sp.groebner(intersection(ideals,cs),*cs,domain=sp.QQ)!=gb:raise ValueError('wrong ordinary primary intersection')
        for j,component in enumerate(tested):
            if component['algebra']['component']!=components[j]:raise ValueError('wrong primary stopping component')
            verify_component(data,component,N,E,pure)
        if result['zero']!=all(c['zero'] for c in tested):raise ValueError('wrong ordinary primary zero decision')
    if not result['zero']:
        w=result['witness'];alpha=tuple(w['marker']);m=w['power']
        if m<1 or (pure and any(alpha)) or (not pure and m<=N*sum(alpha)+E):raise ValueError('invalid forbidden moment index')
        value=original_moment(data,alpha,m)
        if encode_polynomial(value,cs)!=w['polynomial'] or gb.reduce(value)[1]==0:raise ValueError('false ordinary nonmembership witness')
    return True


def radical_ideal(ideal,r,out,budget):
    cs=sp.symbols('c:'+str(r));z=sp.Symbol('intersection_aux')
    gb=sp.groebner([polynomial(f,cs) for f in ideal],*cs,domain=sp.QQ)
    if gb.reduce(sp.S.One)[1]==0:return [[[[0]*r,'1']]],[]
    components=primary_components(ideal,r,out,budget)
    current=[sp.S.One]
    for component in components:
        following=[polynomial(f,cs) for f in component['radical_ideal']]
        gb=sp.groebner([z*f for f in current]+[(1-z)*g for g in following],z,*cs,domain=sp.QQ)
        current=[f.as_expr() for f in gb.polys if not f.as_expr().has(z)]
    return [encode_polynomial(f.as_expr(),cs) for f in sp.groebner(current,*cs,domain=sp.QQ).polys],components


def verify_numerical_family(data,numerical):
    if numerical.get('input')!=data or numerical['status']!='DECIDED':raise ValueError('wrong numerical family input')
    family=numerical['family_result']
    if family['input']!=data or family['core_sha256']!=G.sha(F.CORE):raise ValueError('wrong joint family source')
    G.verify_family_torsion(data,family)
    components=numerical['component_certificates'];primes=family['safe_closed_components']
    if len(components)!=len(primes):raise ValueError('incomplete numerical safe-component coverage')
    for component,prime in zip(components,primes):
        if component['input']!=dict(data,generic_prime=prime):raise ValueError('wrong generic safe component')
        G.verify_result(component)
        if not component['actual_period_selection']['safe']:raise ValueError('unsafe numerical component')
    C=max((sp.Rational(c['actual_period_selection']['cutoff']['slope']) for c in components),default=sp.S.Zero)
    E=max((sp.Rational(c['actual_period_selection']['cutoff']['intercept']) for c in components),default=sp.S.Zero)
    n,r,phase,base,precision=F.normalize(data);D=max((sum(e[:n]) for e in phase),default=0)
    L=family['normalization_torsion_exponent'];ringE=E if L==1 else E+(C*D+1)*(r*L-1)
    for bounds,e in [(numerical['all_safe_field_fibers_cutoff'],E),(numerical['presented_ring_normalization_cutoff'],ringE)]:
        if bounds['slope']!=str(C) or bounds['intercept']!=str(e):raise ValueError('wrong numerical upper certificate')
    return True


def least_bound(data,out,pure_linear=False,field_safe=False,timeout=30,numerical_result=None):
    """Least nonnegative integer N for m>N*degQ (+N in affine mode).

    field_safe restricts to all safe field fibers, and additionally pure
    fibers in linear mode. Otherwise the identity must hold in the input ring.
    """
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True);budget=Budget(timeout)
    numerical=G.family_bounds(data,out/'numerical',timeout) if numerical_result is None else numerical_result
    if numerical.get('input')!=data:raise ValueError('numerical bound belongs to another family')
    if numerical['status']!='DECIDED':return dict(status='UNKNOWN',reason='numerical family bound unavailable',numerical_result=numerical)
    verify_numerical_family(data,numerical)
    target=numerical['family_result']['safe_radical_ideal'] if field_safe else data.get('parameter_ideal',[])
    try:
        if not field_safe and not numerical['family_result']['coefficient_ring_phase_safe']:
            return dict(status='NO_FINITE_BOUND',reason='presented-ring phase is not eventually safe',numerical_result=numerical)
        pure_result=None;pure_safe_decomposition=None
        if pure_linear:
            if field_safe:
                pure_result=run(data,out/'pure_ideal',pure=True,timeout=0 if budget.end is None else max(1,math.ceil(budget.remaining())))
                if pure_result['status']!='DECIDED':raise BudgetExpired
                target,pure_safe_decomposition=radical_ideal(target+pure_result['ideal_lift'],data['parameter_count'],out/'pure_safe_radical',budget)
            else:
                pure_result=zero_test(data,target,out/'ring_purity',pure=True,budget=budget)
                if pure_result['status']!='DECIDED':raise BudgetExpired
                if not pure_result['zero']:return dict(status='NO_LINEAR_BOUND',reason='positive pure moments do not all vanish',purity_result=pure_result,numerical_result=numerical)
        bounds=numerical['all_safe_field_fibers_cutoff'] if field_safe else numerical['presented_ring_normalization_cutoff']
        C=sp.Rational(bounds['slope']);E=sp.Rational(bounds['intercept'])
        upper=int(sp.ceiling(C+E if pure_linear else max(C,E)));tests=[]
        for N in range(upper+1):
            budget.check();test=zero_test(data,target,out/('cutoff_'+str(N)),N,0 if pure_linear else N,budget=budget)
            tests.append(dict(N=N,oracle=test))
            if test['status']!='DECIDED':raise BudgetExpired
            if test['zero']:
                return dict(status='DECIDED',input=data,least_integer_coefficient=N,upper_certificate=upper,
                    pure_linear=pure_linear,field_safe_only=field_safe,target_ideal=target,numerical_result=numerical,
                    pure_result=pure_result,pure_safe_decomposition=pure_safe_decomposition,ordered_cutoff_tests=tests,
                    convention='m>N*deg(Q)' if pure_linear else 'm>N*(deg(Q)+1)',
                    least_real_slope_claimed=False)
        raise ValueError('numerical upper bound failed exact prescribed-cutoff test')
    except BudgetExpired:return dict(status='UNKNOWN',reason='least-bound calculation timeout',numerical_result=numerical)


def verify_least_result(result):
    numerical=result['numerical_result'];data=numerical['input'];verify_numerical_family(data,numerical)
    if result['status']=='NO_FINITE_BOUND':
        if numerical['family_result']['coefficient_ring_phase_safe']:raise ValueError('false non-safety decision')
        return True
    if result['status']=='NO_LINEAR_BOUND':
        target=data.get('parameter_ideal',[])
        verify_zero_test(data,target,result['purity_result'],pure=True)
        if result['purity_result']['zero']:raise ValueError('false purity obstruction')
        return True
    if result['status']!='DECIDED':raise ValueError('least bound not decided')
    pure=result['pure_linear'];fields=result['field_safe_only'];target=result['target_ideal'];r=data['parameter_count'];cs=sp.symbols('c:'+str(r))
    expected=numerical['family_result']['safe_radical_ideal'] if fields else data.get('parameter_ideal',[])
    if pure:
        if fields:
            verify_result(result['pure_result']);combined=expected+result['pure_result']['ideal_lift']
            combinedgb=sp.groebner([polynomial(f,cs) for f in combined],*cs,domain=sp.QQ)
            parts=result['pure_safe_decomposition']
            if combinedgb.reduce(sp.S.One)[1]==0:expected=[[[[0]*r,'1']]]
            else:
                if not parts:raise ValueError('missing pure-safe primary decomposition')
                primary=intersection([[polynomial(f,cs) for f in c['primary_ideal']] for c in parts],cs)
                if sp.groebner(primary,*cs,domain=sp.QQ)!=combinedgb:raise ValueError('wrong pure-safe primary intersection')
                radical=intersection([[polynomial(f,cs) for f in c['radical_ideal']] for c in parts],cs)
                expected=[encode_polynomial(f,cs) for f in radical]
        else:
            verify_zero_test(data,expected,result['pure_result'],pure=True)
            if not result['pure_result']['zero']:raise ValueError('missing ring purity')
    def gb(fs):return sp.groebner([polynomial(f,cs) for f in fs],*cs,domain=sp.QQ)
    if gb(expected)!=gb(target):raise ValueError('wrong least-bound target')
    tests=result['ordered_cutoff_tests'];least=result['least_integer_coefficient']
    if type(least) is not int or least<0 or [t['N'] for t in tests]!=list(range(least+1)):raise ValueError('missing smaller integer tests')
    for t in tests:
        N=t['N'];verify_zero_test(data,target,t['oracle'],N,0 if pure else N)
        if t['oracle']['zero']!=(N==least):raise ValueError('false least integer coefficient')
    return True


def least_generating_prefix(result):
    """Least total-index prefix; for pure moments this is the usual moment prefix."""
    verify_result(result);data=result['input'];N=result['cutoff_slope'];E=result['cutoff_intercept'];pure=result['pure_moments_only']
    n,r,cs,original=normalized_input(data,N,E,pure)
    witnesses=result.get('original_generating_moments',[])
    if result.get('method')=='constant_phase':witnesses=[result['original_generating_moment']]
    if result.get('method')=='nondegenerate_Gaussian_pairing':witnesses=result['pairing_certificate']['original_moments']
    high=max((sum(w['marker'])+w['power'] for w in witnesses),default=0);low=0;trials=[]
    target=sp.groebner([polynomial(f,cs) for f in result['ideal_lift']],*cs,domain=sp.QQ)
    @lru_cache(maxsize=None)
    def contains(k):
        indices=[((0,)*n,m) for m in range(1,k+1)] if pure else forbidden_indices(n,N,E,k)
        gb=sp.groebner(original+[original_moment(data,a,m) for a,m in indices],*cs,domain=sp.QQ)
        answer=gb==target;trials.append(dict(total_index_prefix=k,full_ideal_generated=answer));return answer
    if not contains(high):raise ValueError('retained generators do not give a prefix upper bound')
    while low<high:
        middle=(low+high)//2
        if contains(middle):high=middle
        else:low=middle+1
    if not contains(low) or (low and contains(low-1)):raise ValueError('incorrect least generating prefix')
    return dict(least_total_index_prefix=low,pure_moment_prefix=pure,finite_membership_trials=trials,
                convention='all forbidden (alpha,m) with |alpha|+m<=prefix; pure mode uses alpha=0')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--slope',type=int,default=0);parser.add_argument('--intercept',type=int,default=0)
    parser.add_argument('--pure',action='store_true');parser.add_argument('--timeout',type=int,default=30);parser.add_argument('--workdir',type=Path)
    parser.add_argument('--least',action='store_true');parser.add_argument('--safe-fields',action='store_true')
    parser.add_argument('--prefix',action='store_true')
    parser.add_argument('--numerical-result',type=Path)
    args=parser.parse_args();data=json.loads(args.input.read_text())
    work=args.workdir or args.output.with_suffix('.cas')
    if args.least:
        numerical=json.loads(args.numerical_result.read_text()) if args.numerical_result else None
        result=least_bound(data,work,args.pure,args.safe_fields,args.timeout,numerical)
    else:result=run(data,work,args.slope,args.intercept,args.pure,args.timeout)
    if args.prefix and not args.least and result['status']=='DECIDED':result['least_generating_prefix']=least_generating_prefix(result)
    args.output.write_text(dump(result))


if __name__=='__main__':main()
