"""Explicit source coordinates for every forced F1 boundary pair.

Produces universal polynomial identities over Q[a,b,b^-1], b!=0.
The geometric coverage and pole-filtration interpretation are written.
"""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,shutil
HERE=Path(__file__).resolve().parent
x,y,s,t,alpha,beta,a,b=S.symbols('x y s t alpha beta a b')
start=monotonic()

def degree(f):return S.Poly(S.expand(f),x,y).total_degree()
def top(f):
    p=S.Poly(S.expand(f),x,y);d=p.total_degree()
    return S.expand(sum(c*x**i*y**j for (i,j),c in p.terms() if i+j==d))
def certificate(expr,R,var):
    # Division over the coefficient fraction field is followed by the
    # explicit check that no parameter other than b was inverted.
    others=sorted((expr.free_symbols|R.free_symbols)-{var},key=str)
    field=S.QQ.frac_field(*others)
    q,rem=S.div(S.Poly(expr,var,domain=field),S.Poly(R,var,domain=field))
    assert rem.is_zero
    q=S.cancel(q.as_expr())
    den=S.denom(q)
    assert (den.free_symbols-{b})==set(),den
    assert S.cancel(expr-q*R)==0
    return q

def case(marked):
    if marked:
        B=alpha*s+beta*t*t
        E=alpha*s*(a*s+b*t)+beta*(s**3+t*t*(a*s+b*t))
        BB=b*b+x*y;AA=a+BB*x*x
        TT=(1-x*AA)/b
        AL=(BB*(2*AA-x*AA*AA)-y)/(b*b)
        YI=beta*(2*(a+beta*s*s)-s*(a+beta*s*s)**2)-b*b*alpha
        tv=(1-a*s-beta*s**3)/b
        R=b*b*alpha*s+beta*(1-a*s-beta*s**3)**2-b*b
        section=[alpha*t,beta*t*t,beta*s*t,beta*s*s]
        normal_B=(s,t*t)
        normal_E=(s*(a*s+b*t),s**3+t*t*(a*s+b*t))
        variable=alpha
    else:
        AL=b+a*b*x+x*x*y
        TT=(1-a*x+AL*x*x)/b
        BB=(a*a*b-(1-a*x)*y-AL*AL)/b
        B=alpha*t+beta*s*s
        E=alpha*(-s*s+t*(a*s+b*t))+beta*s*s*(a*s+b*t)
        YI=a*a*alpha-(1+a*s)*(alpha*alpha+b*beta)
        tv=(1-a*s+alpha*s*s)/b
        R=(1-a*s)*alpha+s*s*(alpha*alpha+b*beta)-b
        section=[beta*t*t,beta*s*t,alpha*t,alpha*s]
        normal_B=(t,s*s)
        normal_E=(-s*s+t*(a*s+b*t),s*s*(a*s+b*t))
        variable=beta
    sub={s:x,t:TT,alpha:AL,beta:BB}
    assert S.cancel(B.subs(sub,simultaneous=True)-1)==0
    assert S.cancel(E.subs(sub,simultaneous=True)-1)==0
    assert S.cancel(YI.subs(sub,simultaneous=True)-y)==0
    inverse={x:s,y:YI}
    reverse={}
    for name,value,original in [('alpha',AL,alpha),('beta',BB,beta),('t',TT,t)]:
        expr=S.cancel((value.subs(inverse,simultaneous=True)-original).subs(t,tv))
        reverse[name]=str(certificate(expr,R,variable))
    gens=[S.expand(v.subs(sub,simultaneous=True)) for v in section]
    assert [degree(v) for v in gens]==[16,12,8,4]
    assert S.expand(normal_B[0]*normal_E[1]-normal_B[1]*normal_E[0]-s**4)==0
    return dict(name='contact_at_marked' if marked else 'contact_away_marked',
        boundary_columns=dict(B=list(map(str,normal_B)),E=list(map(str,normal_E))),
        cox_equations=dict(B=str(B),E=str(E)),
        forward=dict(s=str(x),t=str(S.expand(TT)),alpha=str(S.expand(AL)),beta=str(S.expand(BB))),
        inverse_y=str(S.expand(YI)),eliminated_t=str(tv),single_relation=str(S.expand(R)),
        reverse_identity_multipliers=reverse,
        ordered_section_basis=list(map(str,section)),source_sections=list(map(str,gens)),
        section_degrees=[degree(v) for v in gens],section_top_forms=[str(top(v)) for v in gens],
        section_vanishing_orders_at_contact=[0,1,2,3])
rows=[case(True),case(False)]
payload=dict(as_of='2026-11-14',outcome='EXACT_COORDINATES_COMPUTED',
    parameter_domain='Q[a,b,b^-1], b!=0; every complex specialization with b!=0',cases=rows,
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    assurance=dict(normal_form_coverage='written proof',cox_slice_interpretation='written proof',
                   all_height_section_generation='written monomial-factorization proof',
                   imported_degree_theorem_independently_reproved=False),
    elapsed_seconds=round(monotonic()-start,6))
(HERE/'source_coordinate_certificates.json').write_text(json.dumps(payload,indent=2)+'\n')
# Preserve the exact primary source already retained in this repository.
root=HERE.parents[1]
source=root/'foundations/plane-jc/external/arxiv-2204.14178'
manifest=[]
for old,new in [('Increasing_Lower_Bound_20_04_2022.tex','2204.14178v1.tex'),
                ('2204.14178.pdf','2204.14178v1.pdf')]:
    destination=HERE/'sources'/new
    shutil.copyfile(source/old,destination)
    data=destination.read_bytes()
    manifest.append(dict(file=new,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),
        url='https://arxiv.org/'+('pdf/' if new.endswith('.pdf') else 'src/')+'2204.14178v1',
        copied_from_research_relative=str((source/old).relative_to(root))))
(HERE/'sources/degree-transfer-manifest.json').write_text(json.dumps(dict(
    as_of='2026-11-14',sources=manifest,version='arXiv:2204.14178v1, 29 April 2022',
    imported_result='Theorem2.1: every noninvertible plane Keller pair has maximum degree>=125 or degrees(72,108) up to exchange',
    external_theorem=True,independent_whole_proof_replay=False,
    conditional_2026_degree_125_exclusion_used=False),indent=2)+'\n')
for row in rows:print(row['name'],row['section_degrees'],row['section_top_forms'])
print('seconds',monotonic()-start)
