"""PF2RULED15: complete exact source/pole filtration and contact cusp-jet controls.
The finite-dimensional section-space coverage is PF2COORD1; this is not
an independent implementation or a formal proof assistant verification.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from time import monotonic
import hashlib,json,platform
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]

def check(name,value):
    assert value,name
    checks.append(name)

def add(a,b,k=1):
    c=a.copy()
    for t,v in b.items():
        c[t]=c.get(t,0)+k*v
        if not c[t]:del c[t]
    return c

def mul(a,b):
    c={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            c[i+k,j+l]=c.get((i+k,j+l),0)+x*y
    return {t:v for t,v in c.items() if v}

def powers(a,n):
    out=[{(0,0):1}]
    for _ in range(n):out.append(mul(out[-1],a))
    return out

def order(t):return sum(t),t[0]

def chart(poly):
    out={}
    # X=u^-1, Y=u^2*(u*v-1), exactly the actual marked inclusion.
    for (a,b),c in poly.items():
        for k in range(b+1):
            t=(2*b-a+k,k)
            out[t]=out.get(t,0)+c*comb(b,k)*(-1)**(b-k)
    return {t:v for t,v in out.items() if v}

def hull(points):
    ps=sorted(set(points))
    def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    lo=[]
    for p in ps:
        while len(lo)>1 and cross(lo[-2],lo[-1],p)<=0:lo.pop()
        lo.append(p)
    hi=[]
    for p in reversed(ps):
        while len(hi)>1 and cross(hi[-2],hi[-1],p)<=0:hi.pop()
        hi.append(p)
    return lo[:-1]+hi[:-1]

def rational_rank(columns):
    piv={}
    for a in columns:
        a={i:F(x) for i,x in enumerate(a) if x}
        while a:
            k=min(a)
            if k not in piv:
                c=a[k];piv[k]={i:x/c for i,x in a.items()};break
            a=add(a,piv[k],-a[k])
    return len(piv)

def coeff(row,u,v):return row['chart'].get((u,v),0)

def symbol(row,k):return row['symbol'].get(k,0)

def main():
    start=monotonic()
    alpha={(0,0):1,(2,1):1}
    beta={(0,0):-1,(2,1):-2,(4,2):-1,(0,1):-1}
    t={(0,0):1,(2,0):1,(4,1):1}
    ap=powers(alpha,8);bp=powers(beta,8);tp=powers(t,16)
    # Bcal=alpha*t+beta*X^2=1 on the actual source slice.
    check('actual Cox slice B=1',add(mul(alpha,t),{(a+2,b):c for (a,b),c in beta.items()})=={(0,0):1})
    spaces={};certspaces={}
    for m in (5,6,7,8):
        original=[];labels=[];piv={};rows=[]
        for j in range(m+1):
            for i in range(m+j+1):
                idx=len(original);labels.append((j,i))
                source={(x+i,y):F(v) for (x,y),v in mul(mul(ap[m-j],bp[j]),tp[m+j-i]).items()}
                original.append(source);poly=source.copy();vector={idx:F(1)}
                while poly:
                    pivot=max(poly,key=order)
                    if pivot not in piv:
                        c=poly[pivot]
                        poly={t:v/c for t,v in poly.items()};vector={t:v/c for t,v in vector.items()}
                        row=dict(pivot=pivot,poly=poly,combination=vector,index=idx)
                        piv[pivot]=row;rows.append(row);break
                    row=piv[pivot];c=poly[pivot]
                    poly=add(poly,row['poly'],-c);vector=add(vector,row['combination'],-c)
                check(f'Cox section independent m{m} index{idx}',bool(poly))
        check(f'complete dimension V{m}',len(rows)==(m+1)*(3*m+2)//2)
        for row in rows:
            idx=row['index'];vec=row['combination'];pivot=row['pivot'];poly=row['poly']
            reconstructed={}
            for k,c in vec.items():reconstructed=add(reconstructed,original[k],c)
            check(f'exact source reconstruction m{m} row{idx}',reconstructed==poly)
            check(f'invertible triangular basis change m{m} row{idx}',max(vec)==idx and vec[idx]!=0)
            check(f'source echelon m{m} row{idx}',poly[pivot]==1 and max(poly,key=order)==pivot)
            sym={}
            for k,c in vec.items():
                j,i=labels[k];power=2*(m-j)+i
                sym[power]=sym.get(power,0)+c*(-1)**(m-j)
            row['symbol']={k:v for k,v in sym.items() if v}
            row['chart']=chart(poly)
            check(f'whole marked chart polynomial m{m} row{idx}',all(a>=0 for a,b in row['chart']))
            row['degree']=sum(pivot)
        spaces[m]=rows
        certspaces[str(m)]={'dimension':len(rows),'cox_basis_ji':labels,
            'echelon_rows':[{'pivot':r['pivot'],'combination':[[k,str(v)] for k,v in r['combination'].items()],
                'boundary_symbol':{str(k):str(v) for k,v in r['symbol'].items()}} for r in rows]}
    def filtered(m,d):return [r for r in spaces[m] if r['degree']<=d]
    configurations=[(6,72,63,18,1),(7,72,77,18,2),(8,108,112,27,1)]
    summaries=[]
    for m,d,dimension,n,vface_degree in configurations:
        rows=filtered(m,d);top=[r for r in rows if r['degree']==d]
        check(f'exact filtered dimension V{m} source{d}',len(rows)==dimension)
        check(f'unique ordinary homogeneous form V{m} source{d}',len(top)==1 and
              {t:v for t,v in top[0]['poly'].items() if sum(t)==d}=={(3*n,n):1})
        points={t for row in rows for t in row['chart']}
        check(f'maximum v degree V{m} source{d}',max(b for a,b in points)==n)
        check(f'maximum u degree on v face V{m} source{d}',max(a for a,b in points if b==n)==vface_degree)
        for row in rows:
            # The v^n u^0 term equals the coefficient of the unique source top monomial.
            check(f'source top equals marked restriction lead V{m} row{row["index"]}',
                  coeff(row,0,n)==row['poly'].get((3*n,n),0))
        summaries.append(dict(pole=m,ordinary_bound=d,dimension=dimension,
            top_form=f'X^{3*n}*Y^{n}',chart_v_degree=n,chart_top_v_u_degree_bound=vface_degree,
            source_support_hull=hull({t for row in rows for t in row['poly']}),chart_support_hull=hull(points)))
    for row in filtered(7,72):
        check('pole7 source72 symbol vanishes to order8 row'+str(row['index']),all(k>=8 for k in row['symbol']))
    for row in filtered(6,72):
        check('pole6 first top jet identity row'+str(row['index']),coeff(row,1,18)==-symbol(row,5))
        check('pole6 next jet identity row'+str(row['index']),coeff(row,1,17)-13*coeff(row,0,18)==-symbol(row,6))
        check('pole6 third jet identity row'+str(row['index']),coeff(row,1,16)-12*coeff(row,0,17)==5*symbol(row,5)-symbol(row,7))
    for row in filtered(8,108):
        check('pole8 first top jet identity row'+str(row['index']),coeff(row,1,27)==-symbol(row,4))
        check('pole8 next jet identity row'+str(row['index']),coeff(row,1,26)-20*coeff(row,0,27)==-symbol(row,5))
        check('pole8 third jet identity row'+str(row['index']),coeff(row,1,25)-19*coeff(row,0,26)==7*symbol(row,4)-symbol(row,6))
        check('pole8 fourth jet identity row'+str(row['index']),coeff(row,1,24)-18*coeff(row,0,25)+7*coeff(row,0,27)==7*symbol(row,5)-symbol(row,7))
    # Complete lower-pole first-jet freedom is a deliberate negative control.
    jet_ranks=[]
    for m,d,n,expected_dim,expected_rank in [(5,72,18,49,35),(6,72,18,63,38),(7,108,27,91,52)]:
        rows=filtered(m,d);mon=[(a,b) for a in (0,1) for b in range(n+1)]
        columns=[[r['chart'].get(t,0) for t in mon] for r in rows]
        rank=rational_rank(columns)
        check(f'lower pole jet rank V{m} degree{d}',len(rows)==expected_dim and rank==expected_rank)
        jet_ranks.append(dict(pole=m,source_bound=d,jet_monomials=len(mon),dimension=len(rows),rank=rank))
    check('companion fixes multiplier slope',F(20,27)*27==20)
    check('pole6 higher contact slopes contradict companion',F(13,18)!=F(20,27))
    check('pole6 order6 coefficient relation',F(13)-18*F(20,27)==-F(1,3))
    check('next companion jet fixes cusp unit',F(19)-26*F(20,27)==-F(7,27))
    c,d,e=S.symbols('c d e',nonzero=True)
    delta=-7*d/(729*c)
    cusp=(2*e+27*c)/d-26*d/(27*c)
    check('cusp parameter determined by third companion coefficient',
          S.cancel(18*e-7*c-S.Rational(20,27)*25*e-delta*(26*d+27*cusp*c))==0)
    cert=dict(as_of='2026-11-14',exact_field='Q',spaces=certspaces,filtered_summaries=summaries,
        jet_ranks=jet_ranks,
        pole6_jet_identities=['p1_18=-s5','p1_17-13*p0_18=-s6','p1_16-12*p0_17=5*s5-s7'],
        pole8_jet_identities=['q1_27=-t4','q1_26-20*q0_27=-t5','q1_25-19*q0_26=7*t4-t6','q1_24-18*q0_25+7*q0_27=7*t5-t7'],
        remaining_contact_multiplicities={'6,8':[2],'7,8':[2,3]},
        cusp_multiplier={'lambda_slope':'20/27','delta':'-7*q0_26/(729*q0_27)','delta_nonzero':True,
                         'pole6_symbol_order6':'-p0_18/3','cusp_parameter':'(2*q0_25+27*q0_27)/q0_26-26*q0_26/(27*q0_27)'},finite_map_constructed=False)
    certpath=HERE/'contact_filtration_certificates.json';certpath.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'CONTACT_SOURCE_FILTRATION.md',certpath,HERE/'SOURCE_COORDINATES.md',
            HERE/'WRIGHT_MARKED_MODEL.md',HERE/'GLOBAL_CUSP_LEADING_OBSTRUCTION.md',
            HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        written_result='Complete exact contact source filtration excludes multiplicities3,4,6 for(6,8) and multiplicity1 for(7,8); all actual height-eight pairs have a linear cusp multiplier of slope20/27',
        frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],remaining_contact_multiplicities={'6,8':[2],'7,8':[2,3]},
        ordinary_degrees=[72,108],boundary_restriction_degrees=[18,27],
        common_cusp_multiplier_degree_upper_bound=1,common_cusp_multiplier_slope='20/27',
        pole7_first_jet_obstruction_from_unrestricted_lower_pole_lifts=False,
        assurance=dict(complete_finite_section_spaces=True,exact_rational_basis_change=True,
            finite_map_and_cusp_inputs='PF2OMIT9 written theorem',source_section_coverage='PF2COORD1 written theorem',
            source_top_terms_assigned_to_published_Newton_corner=False,conditional_2026_degree125_upgrade_used=False,
            complete_height_eight_excluded=False,all_pole_four_coordinates_excluded=False,
            actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False),
        source_sha256={str(f.relative_to(HERE)) if f.is_relative_to(HERE) else '../jc2-omitted-divisors/'+f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),elapsed_seconds=round(monotonic()-start,6))
    (HERE/'contact_filtration_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; contact(6,8) only k2; contact(7,8) only k2/k3; lambda slope20/27')

if __name__=='__main__':main()
