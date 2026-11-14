"""PF2RULED16: a complete chart-leading identity excludes all quartic height eight.
Reconstructs the retained PF2RULED15 rational basis using shared exact
polynomial utilities. This is not an independent whole-proof implementation.
"""
from pathlib import Path
from time import monotonic
from math import gcd
import hashlib,json,platform
import sympy as S
import verify_contact_filtration as V
HERE=Path(__file__).resolve().parent
checks=[]

def check(name,value):
    assert value,name
    checks.append(name)

def terms(poly):return [[a,b,str(c)] for (a,b),c in sorted(poly.items())]

def main():
    start=monotonic()
    prior=json.loads((HERE/'contact_filtration_certificates.json').read_text())
    receipt=json.loads((HERE/'contact_filtration_receipt.json').read_text())
    check('retained complete filtration input passed',receipt['outcome']=='PASS')
    for path,digest in receipt['source_sha256'].items():
        check('retained filtration input seal '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    alpha={(0,0):1,(2,1):1}
    beta={(0,0):-1,(2,1):-2,(4,2):-1,(0,1):-1}
    t={(0,0):1,(2,0):1,(4,1):1}
    ap=V.powers(alpha,8);bp=V.powers(beta,8);tp=V.powers(t,16)
    spaces={}
    for m in (7,8):
        data=prior['spaces'][str(m)]
        labels=[(j,i) for j in range(m+1) for i in range(m+j+1)]
        check(f'complete Cox labels V{m}',list(map(tuple,data['cox_basis_ji']))==labels)
        originals=[{(x+i,y):c for (x,y),c in V.mul(V.mul(ap[m-j],bp[j]),tp[m+j-i]).items()} for j,i in labels]
        check(f'complete certificate dimension V{m}',len(data['echelon_rows'])==len(labels)==(m+1)*(3*m+2)//2)
        pivots=set();rows=[]
        for idx,entry in enumerate(data['echelon_rows']):
            poly={};symbol={};vector={k:V.F(v) for k,v in entry['combination']}
            check(f'invertible triangular change V{m} row{idx}',max(vector)==idx and vector[idx]!=0)
            for k,c in vector.items():
                poly=V.add(poly,originals[k],c)
                j,i=labels[k];power=2*(m-j)+i
                symbol[power]=symbol.get(power,0)+c*(-1)**(m-j)
            symbol={k:v for k,v in symbol.items() if v}
            pivot=tuple(entry['pivot'])
            check(f'exact source echelon row V{m} row{idx}',max(poly,key=V.order)==pivot and poly[pivot]==1 and pivot not in pivots)
            pivots.add(pivot)
            check(f'exact symbol reconstruction V{m} row{idx}',symbol=={int(k):V.F(v) for k,v in entry['boundary_symbol'].items()})
            chart=V.chart(poly)
            check(f'actual chart regularity V{m} row{idx}',all(a>=0 for a,b in chart))
            rows.append(dict(pivot=pivot,source=poly,chart=chart,symbol=symbol,index=idx))
        spaces[m]=rows
    P=[r for r in spaces[7] if sum(r['pivot'])<=72]
    Q=[r for r in spaces[8] if sum(r['pivot'])<=108]
    check('complete lower-coordinate intersection dimension',len(P)==77)
    check('complete companion intersection dimension',len(Q)==112)
    for row in P:
        check('lower-coordinate chart total degree<=20 row'+str(row['index']),max(map(sum,row['chart']))<=20)
    check('bound20 is attained in the full P space',max(sum(t) for row in P for t in row['chart'])==20)
    # H=v^7*(u+v)^20. The whole >=27 part is a universal linear
    # combination of H and only the boundary coefficients s4,s5.
    H={(k,27-k):V.F(S.binomial(20,k)) for k in range(21)}
    uH={(a+1,b):c for (a,b),c in H.items()}
    uK={(k+1,26-k):V.F(S.binomial(19,k)) for k in range(20)}
    top_source_rows=[]
    for row in Q:
        idx=row['index'];chart=row['chart'];c=chart.get((0,27),0)
        s4=row['symbol'].get(4,0);s5=row['symbol'].get(5,0)
        expected=V.add(V.add({t:c*v for t,v in H.items() if c},uH,-s4),uK,-s5)
        actual={t:v for t,v in chart.items() if sum(t)>=27}
        check('universal companion high-total-degree identity row'+str(idx),actual==expected)
        check('source top coefficient equals v27 restriction coefficient row'+str(idx),c==row['source'].get((81,27),0))
        check('companion symbol has no coefficients below4 row'+str(idx),all(k>=4 for k in row['symbol']))
        if sum(row['pivot'])==108:top_source_rows.append(row)
    check('exact source108 forces nonzero H coefficient',len(top_source_rows)==1 and
          {t:c for t,c in top_source_rows[0]['source'].items() if sum(t)==108}=={(81,27):1})
    w,v=S.symbols('w v')
    check('primitive exponents',gcd(7,20)==1 and gcd(27,20)==1)
    # In w=u+v the leading bracket equals
    # v^6*w^19*(7*w*F_w-20*v*F_v).
    eigenvalues={}
    for d in range(1,21):
        eigenvalues[str(d)]=[]
        for i in range(d+1):
            mon=w**i*v**(d-i);eigen=27*i-20*d
            check(f'exact diagonal bracket operator d{d} i{i}',
                  S.expand(7*w*S.diff(mon,w)-20*v*S.diff(mon,v)-eigen*mon)==0)
            check(f'nonresonance d{d} i{i}',eigen!=0)
            eigenvalues[str(d)].append(eigen)
    check('leading bracket degree exceeds prescribed u degree for every d',all(d+27-2>1 for d in range(1,21)))
    # Explicit resonant and nonresonant controls prevent a false universal
    # centralizer claim: degree27 itself commutes, lower degrees do not.
    primitive=v**7*w**20
    check('degree27 resonance retained',S.expand(7*w*S.diff(primitive,w)-20*v*S.diff(primitive,v))==0)
    check('coordinate degree1 is nonresonant',S.expand(7*w*S.diff(w,w)-20*v*S.diff(w,v))==7*w)
    pairs=[[m,9] for m in range(4,9) if 9%m and not(m==4 and 9%2)]
    check('first unexcluded height9 reduced pairs',pairs==[[5,9],[6,9],[7,9],[8,9]])
    cert=dict(as_of='2026-11-14',exact_field='Q',
        input_basis='contact_filtration_certificates.json; reconstructed V7,V8',
        dimensions={'W_7_72':77,'W_8_108':112},lower_coordinate_chart_degree_upper_bound=20,
        universal_companion_high_total_identity='Q_[>=27]=c*v^7*(u+v)^20-s4*u*v^7*(u+v)^20-s5*u*v^7*(u+v)^19',
        coefficient_c='[X^81 Y^27]Q=[v^27]Q0; nonzero at exact ordinary degree108',
        high_polynomial_H=terms(H),symbol_coefficient4_polynomial=terms({t:-c for t,c in uH.items()}),
        symbol_coefficient5_polynomial=terms({t:-c for t,c in uK.items()}),
        contact_hypothesis='s4=s5=0; all previously retained contact Q symbols have order>=8',
        diagonal_operator='7*w*d/dw-20*v*d/dv',eigenvalue_formula='27*i-20*d',eigenvalues=eigenvalues,
        resonant_degree27_control=True,excluded_complete_pole_pairs=[[6,8],[7,8]],
        minimum_compatible_quartic_pole_height=9,first_open_reduced_pairs=pairs,finite_map_constructed=False)
    path=HERE/'height_eight_certificates.json';path.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'POLE_HEIGHT_EIGHT_EXCLUSION.md',path,HERE/'verify_contact_filtration.py',
            HERE/'contact_filtration_certificates.json',HERE/'contact_filtration_receipt.json',HERE/'CONTACT_SOURCE_FILTRATION.md',
            HERE/'SOURCE_COORDINATES.md',HERE/'WRIGHT_MARKED_MODEL.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        written_result='A primitive chart-leading form excludes all remaining contact(6,8),(7,8); complete compatible quartic pole height8 is excluded',
        excluded_complete_pole_pairs=[[6,8],[7,8]],all_three_PF2RULED15_contact_strata_excluded=True,
        complete_pole_height_eight_excluded=True,minimum_compatible_quartic_pole_height=9,
        first_open_reduced_pairs=pairs,minimum_target_coordinate_pole=4,
        leading_chart_companion='c*v^7*(u+v)^20',lower_coordinate_chart_degree_upper_bound=20,
        assurance=dict(core_chart_obstruction_uses_finiteness=False,core_chart_obstruction_uses_cusp_Euler_input=False,
            core_chart_obstruction_uses_external_degree_theorem=False,
            quartic_height_eight_corollary_uses_retained_2022_degree_alternative=True,
            uses_shared_polynomial_utilities_from_PF2RULED15=True,
            finite_section_basis_coverage='PF2COORD1 written theorem plus exact reconstructed rational basis',
            polynomial_target_shear_preserves_or_lowers_maximum_pole='PF2RULED2 written theorem',
            published_72_108_exception_excluded_in_all_pole_heights=False,
            source_top_terms_assigned_to_published_Newton_corner=False,conditional_2026_degree125_upgrade_used=False,
            height_nine_excluded=False,all_pole_four_coordinates_excluded=False,
            actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'height_eight_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; complete height8 excluded; minimum height9; first pairs{pairs}')

if __name__=='__main__':main()
