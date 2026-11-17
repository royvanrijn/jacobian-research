"""Exact Gaussian graph certificates, rational marker slopes and real regions.

QQ graph enumeration is general but can be impractical. Finite-primary jet
certificate construction over general presented rings is a written algorithm.
The region routine requires a valid order-law certificate and exact zero oracle.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations,count
import copy,json,math,time
import sympy as s
import algebraic_bernstein_generator as W
import bernstein_slope_producer as S
import gaussian_ordinary_cutoff_producer as O

def monomials(n,d):
    if n==0:
        yield ();return
    for k in range(d+1):
        for e in compositions(n,k):yield e
def compositions(n,k):
    if n==1:yield (k,);return
    for j in range(k+1):
        for rest in compositions(n-1,k-j):yield (j,)+rest

def source(data):
    if set(data)-{'dimension','phase_terms','precision_matrix'}:raise ValueError('unknown QQ source input field')
    n=data['dimension'];phase=W.rational_terms(data['phase_terms'],n)
    if type(n)!=int or n<1 or not phase or max(sum(e) for e,c in phase)<1:raise ValueError('positive-degree nonzero QQ phase required')
    H=s.Matrix(data.get('precision_matrix',s.eye(n).tolist()))
    if H.shape!=(n,n) or H!=H.T or not H.det() or any(not v.is_Rational for v in H):raise ValueError('fixed nonsingular rational precision required')
    A=W.Weyl(n+1);v=[A.var(i) for i in range(2*n+2)];u=v[:n];du=v[n+1:2*n+1];tau=v[n];dtau=v[-1]
    f=A.substitute(phase,du);gens=[]
    for i in range(n):
        diff=[]
        for e,c in phase:
            if e[i]:
                j=list(e);j[i]-=1;diff.append([j,str(Q(c)*e[i])])
        gens.append(A.add(A.scale(u[i],-1),*[A.scale(du[j],Q(str(H[i,j]))) for j in range(n)],A.mul(A.substitute(diff,du),dtau)))
    gens.append(A.sub(tau,f))
    ES=A.scale(A.mul(dtau,tau),-1);E=A.add(*[A.mul(a,b) for a,b in zip(u,du)])
    return A,gens,u+du,ES,E,f,H,phase

def sparse_relation(columns,deadline=None):
    basis={}
    for index,(tag,column) in enumerate(columns):
        if deadline and time.monotonic()>deadline:raise TimeoutError('exact graph enumeration cap')
        row=dict(column);comb={index:Q(1)}
        while row:
            lead=max(row)
            if lead not in basis:break
            b,c=basis[lead];scale=row[lead]
            for key,value in b.items():
                row[key]=row.get(key,Q(0))-scale*value
                if not row[key]:del row[key]
            for key,value in c.items():
                comb[key]=comb.get(key,Q(0))-scale*value
                if not comb[key]:del comb[key]
        if not row:
            if any(columns[j][0][0]=='p' and c for j,c in comb.items()):return comb
            continue
        scale=row[max(row)];basis[max(row)]=({k:v/scale for k,v in row.items()},{k:v/scale for k,v in comb.items()})
    return None

def search_identity(data,kind,timeout=30):
    if kind not in ['ordinary','Euler']:raise ValueError('unknown identity')
    A,gens,base,ES,E,f,H,phase=source(data);n=data['dimension'];deadline=time.monotonic()+timeout if timeout else None
    if timeout<0:raise ValueError('nonnegative search timeout required')
    for box in count(1):
        columns=[]
        for exponent in monomials(2*(n+1),box):
            if deadline and time.monotonic()>deadline:raise TimeoutError('exact graph enumeration cap')
            term={exponent:Q(1)}
            for i,g in enumerate(gens):columns.append((('w',i,exponent),A.scale(A.mul(term,g),-1)))
        for degree in range(box+1):
            for exponent in compositions(2*n,degree):
                if kind=='Euler' and sum(exponent[:n])-sum(exponent[n:])<1:continue
                term=A.substitute([[exponent,'1']],base)
                for i in range(box-degree+1):
                    value=A.mul(term,A.power(ES,i))
                    if kind=='ordinary':value=A.mul(value,f)
                    columns.append((('r',exponent,i),A.scale(value,-1)))
        for i in range(box+1):
            for j in range((box-i+1) if kind=='Euler' else 1):
                columns.append((('p',i,j),A.mul(A.power(ES,i),A.power(E,j))))
        relation=sparse_relation(columns,deadline)
        if relation is None:continue
        scalar={};operator={};witness=[{} for g in gens]
        for index,c in relation.items():
            tag=columns[index][0]
            if tag[0]=='p':scalar[(tag[1],tag[2])]=scalar.get((tag[1],tag[2]),Q(0))+c
            elif tag[0]=='r':
                exp=tag[1]+(tag[2],);operator[exp]=operator.get(exp,Q(0))+c
            else:witness[tag[1]]=A.add(witness[tag[1]],{tag[2]:c})
        if not scalar:raise ValueError('missing nonzero scalar identity')
        return dict(kind=kind,box=box,scalar_terms=[[list(e),str(c)] for e,c in sorted(scalar.items()) if c],
                    operator_terms=[[list(e),str(c)] for e,c in sorted(operator.items()) if c],left_witnesses=[A.serialize(w) for w in witness])
    raise AssertionError('unreachable exhaustive search limit')

def verify(data,certificate):
    A,gens,base,ES,E,f,H,phase=source(data);n=data['dimension']
    for kind in ['ordinary','Euler']:
        c=certificate[kind];p=A.substitute(c['scalar_terms'],[ES,E]);R=A.substitute(c['operator_terms'],base+[ES])
        if not p:raise ValueError('zero scalar polynomial')
        if kind=='ordinary' and any(e[1] for e,co in c['scalar_terms']):raise ValueError('ordinary scalar depends on Euler')
        if kind=='Euler' and any(sum(e[:n])-sum(e[n:2*n])<1 for e,co in c['operator_terms'] if Q(co)):raise ValueError('Euler remainder does not strictly raise order')
        difference=A.sub(p,A.mul(R,f) if kind=='ordinary' else R)
        if len(c['left_witnesses'])!=len(gens):raise ValueError('wrong witness count')
        total=A.add(*[A.mul(A.parse(w),g) for w,g in zip(c['left_witnesses'],gens)])
        if total!=difference:raise ValueError('original Gaussian Fourier graph identity failed')
    if certificate.get('graph_generators')!=[A.serialize(g) for g in gens]:raise ValueError('changed original graph presentation')
    if certificate.get('specialization_start')!=0:raise ValueError('graph identities specialize from zero')
    return True

def generate(data,timeout=30):
    A,gens,*_=source(data)
    certificate=dict(ordinary=search_identity(data,'ordinary',timeout),Euler=search_identity(data,'Euler',timeout),
                     graph_generators=[A.serialize(g) for g in gens],specialization_start=0,
                     torsion='Absorbed in the exact UNLOCALIZED graph scalar polynomials;no minimal localized b-function asserted.')
    verify(data,certificate);return certificate

def heat_order(data,m):
    if type(m)!=int or m<0:raise ValueError('nonnegative original power required')
    n,r,phase,ideal,H=O.F.normalize(data);xs=s.symbols('x:'+str(n));cs=s.symbols('c:'+str(r));C=H.inv()
    P=sum(v*s.prod(z**e for z,e in zip(xs+cs,exponent)) for exponent,v in phase.items())
    term=s.expand(P**m);value=term
    for k in range(1,(s.Poly(term,*xs).total_degree()//2)+1):
        term=s.expand(sum(C[i,j]*s.diff(term,xs[i],xs[j]) for i in range(n) for j in range(n))/(2*k))
        if term==0:break
        value+=term
    gb=s.groebner([O.polynomial(g,cs) for g in data.get('parameter_ideal',[])],*cs,domain=s.QQ)
    exponents=[sum(e) for e,c in s.Poly(s.expand(value),*xs).terms() if gb.reduce(s.expand(c))[1]!=0]
    return min(exponents) if exponents else None

def field_family(data):
    n=data['dimension'];phase=W.rational_terms(data['phase_terms'],n)
    return dict(dimension=n,parameter_count=1,phase_terms=[[e+[0],c] for e,c in phase],parameter_ideal=[[[[1],'1']]],
                precision_matrix=data.get('precision_matrix',[[str(v) for v in row] for row in s.eye(n).tolist()]))

def order_law(data,certificate):
    verify(data,certificate);n=data['dimension'];D=max(sum(e) for e,c in data['phase_terms']);D=max(1,D)
    X,Y=S.X,S.Y
    p=sum(s.Rational(c)*X**e[0]*Y**e[1] for e,c in certificate['Euler']['scalar_terms'])
    shifted=s.Poly(s.expand(p.subs(Y,Y-D*X)),X,Y,domain=s.QQ)
    B=s.Poly(sum(s.Rational(c)*X**e[0] for e,c in certificate['ordinary']['scalar_terms']),X,domain=s.QQ)
    J=max((sum(e[n:2*n]) for e,c in certificate['ordinary']['operator_terms']),default=0)
    start=max([0]+S.integer_roots(B))+1
    raw=S.make_plan(dict(polynomial_terms=S.terms(shifted),validity_start=start,step_bound=D+J))
    K=raw['query_operator_power'];d=heat_order(field_family(data),K)
    if d is None:raise ValueError('nonnull field phase unexpectedly nilpotent')
    matches=[q for q in raw['strips'] if abs(q['b']*(d+D*K)-q['a']*K)<q['normal_radius']]
    if len(matches)!=1:raise ValueError('actual order did not select exactly one certified rational strip')
    strip=matches[0];rho=s.Rational(strip['a'],strip['b'])-D
    if rho<0:raise ValueError('negative Gaussian marker slope')
    offsets=[];threshold=K
    for c in range(-strip['normal_radius']+1,strip['normal_radius']):
        line=s.Poly(s.expand(shifted.as_expr().subs(Y,(strip['a']*X+c)/strip['b'])),X,domain=s.QQ)
        if line.is_zero:offsets.append(str(s.Rational(c,strip['b'])))
        else:threshold=max(threshold,int(s.floor(S.root_bound(line)))+1)
    return dict(marker_slope=str(rho),safe=bool(rho>0),asymptotic_coefficient=str(1/rho) if rho else None,
        exact_line_threshold=threshold,late_offsets=offsets,shift_D=D,ordinary_upper_step=J,
        query_power=K,original_order_at_query=d,shifted_strip=strip,
        arithmetic_scope='Only root-free-collar arithmetic is reused;formal-operator safety labels are NOT imported.')

def clock_input(data,r):
    n,k,phase,ideal,H=O.F.normalize(data)
    lifted=[[[0,1]+[0]*(n+k),'1']]
    lifted.extend([[[r-1,0]+list(e),str(c)] for e,c in phase.items()])
    return dict(dimension=n+2,parameter_count=k,phase_terms=lifted,parameter_ideal=copy.deepcopy(data.get('parameter_ideal',[])),
                precision_matrix=[[str(v) for v in row] for row in s.diag(s.Matrix([[0,1],[1,0]]),H).tolist()])

def rational_zero(data,C,E,out,timeout=30):
    C,E=Q(C),Q(E)
    if min(C,E)<0:raise ValueError('nonnegative rational cutoff required')
    n,k,phase,ideal,H=O.F.normalize(data)
    cs=s.symbols('c:'+str(k));gb=s.groebner([O.polynomial(g,cs) for g in data.get('parameter_ideal',[])],*cs,domain=s.QQ)
    # Exact early failures,never a positive prefix inference. Nilpotence is
    # separately certified by the invertible heat map,then all powers are finite.
    orders=[]
    for m in range(1,4):
        d=heat_order(data,m);orders.append(d)
        if d is None:
            return dict(status='DECIDED',zero=True,method='exact finite nilpotent powers',nilpotence_index=m,
                        original_orders=orders,C=str(C),E=str(E))
        if m>C*d+E:
            for alpha in compositions(n,d):
                value=O.original_moment(data,alpha,m)
                if gb.reduce(value)[1]!=0:
                    return dict(status='DECIDED',zero=False,method='exact original moment',C=str(C),E=str(E),
                                witness=dict(marker=list(alpha),power=m,polynomial=O.encode_polynomial(value,cs)))
            raise ValueError('heat order did not supply a detecting original marker')
    active=[i for i in range(n) if any(e[i] for e in phase)]
    degree=min((sum(e[:n]) for e in phase),default=0)
    if degree and C*degree>=1 and all(H.inv()[i,j]==0 for i in active for j in active):
        return dict(status='DECIDED',zero=True,method='exact isotropic support inequality',minimum_phase_order=degree,C=str(C),E=str(E))
    r=1 if C.denominator==E.denominator==1 or C==0 else math.lcm(C.denominator,E.denominator)
    lifted=data if r==1 else clock_input(data,r)
    N=int(C*r);M=math.floor(E) if C==0 else int(E*r)
    result=O.zero_test(lifted,lifted.get('parameter_ideal',[]),Path(out),N,M,timeout=timeout)
    answer=dict(status=result['status'],C=str(C),E=str(E),clock_degree=r,underlying_input=lifted,underlying_oracle=result)
    if result['status']!='DECIDED':return answer
    answer['zero']=result['zero']
    if not result['zero']:
        witness=result['witness'];alpha=witness['marker'];m=witness['power']
        if r!=1:
            numerator=m-alpha[0]+alpha[1]
            if numerator<0 or numerator%r:raise ValueError('clock witness does not specialize to an original power')
            m=numerator//r;alpha=alpha[2:]
        if not m>C*sum(alpha)+E:raise ValueError('mapped witness outside original forbidden cone')
        value=O.original_moment(data,alpha,m)
        cs=s.symbols('c:'+str(k));gb=s.groebner([O.polynomial(g,cs) for g in data.get('parameter_ideal',[])],*cs,domain=s.QQ)
        if gb.reduce(value)[1]==0:raise ValueError('mapped original moment is zero')
        answer['witness']=dict(marker=list(alpha),power=m,polynomial=O.encode_polynomial(value,cs))
    return answer

def rational_ideal(data,C,E,out,timeout=30):
    """Full ordinary ideal through the sealed integer producer and clock law."""
    C,E=Q(C),Q(E)
    if min(C,E)<0:raise ValueError('nonnegative rational cutoff required')
    r=1 if C.denominator==E.denominator==1 or C==0 else math.lcm(C.denominator,E.denominator)
    lifted=data if r==1 else clock_input(data,r)
    N=int(C*r);M=math.floor(E) if C==0 else int(E*r)
    result=O.run(lifted,Path(out),N,M,timeout=timeout)
    return dict(status=result['status'],input=data,C=str(C),E=str(E),clock_degree=r,underlying_result=result,
                ideal_lift=result.get('ideal_lift'),retains_ordinary_nilpotents=True,
                equality_basis='exact monomial-clock ordinary ideal law,not radical equality')

def vertices(lower,points):
    lower=Q(lower);lines=[(Q(1),Q(0),lower),(Q(0),Q(1),Q(0))]
    lines.extend((Q(q),Q(1),Q(m)) for m,q in points)
    out=set()
    for (a,b,c),(d,e,f) in combinations(lines,2):
        det=a*e-b*d
        if not det:continue
        C,E=(c*e-b*f)/det,(a*f-c*d)/det
        if all(x*C+y*E>=z for x,y,z in lines):out.add((C,E))
    return sorted(out)

def refine_region(lower,order,oracle):
    """Exact separation with a certified eventual finite-parallel-line law.

    That law proves termination;arbitrary unverified callbacks do not.
    Each false oracle answer must contain an ORIGINAL nonzero moment.
    """
    points=[];trace=[]
    while True:
        vs=vertices(lower,points);passed=[];changed=False
        for C,E in vs:
            result=oracle(C,E);trace.append(dict(C=str(C),E=str(E),result=result))
            if result.get('status')!='DECIDED':return dict(status='UNKNOWN',trace=trace,retained_original_points=points)
            if result['zero']:passed.append((C,E));continue
            witness=result['witness'];m=witness['power'];q=order(m)
            if q is None or not m>C*q+E or q>sum(witness['marker']):raise ValueError('invalid exact order/nonmembership witness')
            if (m,q) in points:raise ValueError('cutting plane failed strict refinement')
            points.append((m,q));changed=True;break
        if not changed:return dict(status='DECIDED',asymptotic_coefficient=str(lower),vertices=[[str(C),str(E)] for C,E in vs],
            supporting_original_points=[dict(power=m,least_marker_degree=q) for m,q in points],trace=trace,
            recession_generators=[[1,0],[0,1]],strict_moment_convention='vanishing for m>C*deg(Q)+E',
            termination_scope='Requires certified eventual finite parallel lines;no periodicity assumed.')

def solve_field(data,out,timeout=30):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    n=data['dimension'];phase=W.rational_terms(data['phase_terms'],n)
    assurance=dict(independent_replay=False,formal_verification=False,external_review=False)
    if not phase:return dict(status='DECIDED',input=data,method='zero phase',region=dict(vertices=[['0','0']],recession_generators=[[1,0],[0,1]]),**assurance)
    if max(sum(e) for e,c in phase)==0:return dict(status='UNSAFE',input=data,method='nonzero constant phase',region_empty=True,**assurance)
    try:certificate=generate(data,timeout)
    except TimeoutError as exc:return dict(status='UNKNOWN',reason=str(exc),**assurance)
    law=order_law(data,certificate)
    if not law['safe']:return dict(status='UNSAFE',certificate=certificate,order_law=law,region_empty=True,**assurance)
    family=field_family(data);counter=[0]
    def oracle(C,E):
        counter[0]+=1;return rational_zero(family,C,E,out/('vertex_'+str(counter[0])),timeout)
    region=refine_region(Q(law['asymptotic_coefficient']),lambda m:heat_order(family,m),oracle)
    return dict(status=region['status'],input=data,certificate=certificate,order_law=law,region=region,
        general_QQ_graph_search_implemented=True,general_presented_ring_jet_certificate_backend_implemented=False,
        independent_replay=False,formal_verification=False,external_review=False)

def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--timeout',type=int,default=30)
    parser.add_argument('--rational-zero',nargs=2,metavar=('C','E'));parser.add_argument('--rational-ideal',nargs=2,metavar=('C','E'))
    args=parser.parse_args();data=json.loads(args.input.read_text());work=args.output.with_suffix('.work')
    if args.rational_zero:result=rational_zero(data,*args.rational_zero,work,args.timeout)
    elif args.rational_ideal:result=rational_ideal(data,*args.rational_ideal,work,args.timeout)
    else:result=solve_field(data,work,args.timeout)
    args.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
