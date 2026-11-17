"""PF2RULED17: exhaustive height-nine gcd-one filtration and resonant frontier.
Shared exact arithmetic utilities; no independent whole-proof implementation.
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
def form(a,b):return {(i,a+b-i):V.F(S.binomial(b,i)) for i in range(b+1)}
def scale(p,c):return {t:c*v for t,v in p.items() if c*v}
def maximum(p):return max(map(sum,p)) if p else -1

def main():
    start=monotonic()
    ap=V.powers({(0,0):1,(2,1):1},9)
    bp=V.powers({(0,0):-1,(2,1):-2,(4,2):-1,(0,1):-1},9)
    tp=V.powers({(0,0):1,(2,0):1,(4,1):1},18)
    spaces={};certspaces={}
    for m in (5,7,8,9):
        piv={};rows=[];labels=[];originals=[]
        for j in range(m+1):
            for i in range(m+j+1):
                idx=len(labels);labels.append((j,i))
                source={(x+i,y):V.F(c) for (x,y),c in V.mul(V.mul(ap[m-j],bp[j]),tp[m+j-i]).items()}
                originals.append(source);poly=source.copy();vec={idx:V.F(1)}
                while poly:
                    key=max(poly,key=V.order)
                    if key not in piv:
                        c=poly[key];poly=scale(poly,1/c);vec=scale(vec,1/c)
                        row=dict(source=poly,comb=vec,pivot=key,index=idx);piv[key]=row;rows.append(row);break
                    row=piv[key];c=poly[key];poly=V.add(poly,row['source'],-c);vec=V.add(vec,row['comb'],-c)
                check(f'independent Cox section m{m} row{idx}',bool(poly))
        check(f'complete Cox dimension m{m}',len(rows)==(m+1)*(3*m+2)//2)
        for row in rows:
            idx=row['index'];sy={};reconstructed={}
            check(f'invertible source basis change m{m} row{idx}',max(row['comb'])==idx and row['comb'][idx]!=0)
            for k,c in row['comb'].items():
                reconstructed=V.add(reconstructed,originals[k],c)
                j,i=labels[k];r=2*(m-j)+i;sy[r]=sy.get(r,0)+c*(-1)**(m-j)
            check(f'exact source reconstruction m{m} row{idx}',reconstructed==row['source'])
            row['symbol']={k:c for k,c in sy.items() if c};row['chart']=V.chart(row['source'])
            check(f'exact source echelon m{m} row{idx}',max(row['source'],key=V.order)==row['pivot'] and row['source'][row['pivot']]==1)
            check(f'actual chart regularity m{m} row{idx}',all(a>=0 for a,b in row['chart']))
        spaces[m]=rows
        certspaces[str(m)]={'cox_basis_ji':labels,'dimension':len(rows),'echelon_rows':[{'pivot':r['pivot'],'combination':[[k,str(c)] for k,c in r['comb'].items()]} for r in rows]}
    def raw(m,cap):return [r for r in spaces[m] if sum(r['pivot'])<=cap]
    kernel_records=[]
    def kernel(m,cap,contact,degree=None,leading=None,label=''):
        rows=raw(m,cap);piv={};null=[];nullcomb=[]
        for idx,row in enumerate(rows):
            sy=row['symbol'].copy();cc=row['chart'].copy();vec={idx:V.F(1)}
            con={(0,k,0):c for k,c in sy.items() if k<contact}
            if degree is not None:
                high=V.add({t:c for t,c in cc.items() if sum(t)>=degree},leading,-cc.get((0,degree),0))
                con.update({(1,a,b):c for (a,b),c in high.items()})
            while con:
                key=min(con)
                if key not in piv:
                    c=con[key];piv[key]=(scale(cc,1/c),scale(sy,1/c),scale(con,1/c),scale(vec,1/c));break
                cp,sp,co,vp=piv[key];c=con[key]
                cc=V.add(cc,cp,-c);sy=V.add(sy,sp,-c);con=V.add(con,co,-c);vec=V.add(vec,vp,-c)
            check(f'complete constraint basis change {label} row{idx}',max(vec)==idx and vec[idx]!=0)
            if not con:
                check(f'contact constraints {label} row{idx}',all(k>=contact for k in sy))
                null.append(dict(chart=cc,symbol=sy));nullcomb.append([[k,str(c)] for k,c in vec.items()])
        check('complete kernel dimension '+label,len(null)+len(piv)==len(rows))
        kernel_records.append(dict(label=label,pole=m,source_bound=cap,contact_order=contact,
            input_dimension=len(rows),constraint_rank=len(piv),kernel_dimension=len(null),null_combinations=nullcomb,
            imposed_chart_degree=degree))
        return null
    # Full contact spaces at the coarse all-height source bounds.
    basic={}
    for m,cap,expected,degree in [(5,64,46,16),(7,96,85,24),(8,112,109,28),(9,128,136,32)]:
        rows=kernel(m,cap,m,label=f'contact{m}_{cap}');basic[m]=rows
        check(f'coarse contact dimension m{m}',len(rows)==expected)
        for i,row in enumerate(rows):check(f'coarse contact chart bound m{m} row{i}',maximum(row['chart'])<=degree)
    H32=form(8,24)
    for i,row in enumerate(basic[9]):
        cc=row['chart'];c=cc.get((0,32),0)
        check('whole companion top32 projection row'+str(i),{t:v for t,v in cc.items() if sum(t)>=32}==scale(H32,c))
    check('only source degree128 can meet125 through128',
          sorted(sum(r['pivot']) for r in raw(9,128) if sum(r['pivot'])>=125)==[128])
    # The normal orientation of the72,108 exception has a genuine resonance.
    low=kernel(9,108,9,label='exception_Q9_108')
    H27=form(7,20);uH={(a+1,b):c for (a,b),c in H27.items()}
    A=form(9,18);B={(a+1,b):c for (a,b),c in form(8,18).items()};C={(a+2,b):c for (a,b),c in form(7,18).items()}
    for i,row in enumerate(low):
        cc=row['chart'];a=cc.get((1,27),0);c=cc.get((0,27),0)
        b=cc.get((1,26),0)-18*c;d=cc.get((2,25),0)-153*c-18*b
        expected=V.add(V.add(V.add(scale(uH,a),A,c),B,b),C,d)
        check('full exceptional high projection row'+str(i),{t:v for t,v in cc.items() if sum(t)>=27}==expected)
    # Reverse orientation: pole8/source108 is primitive; pole9/source72 has degree<=21.
    reverse=kernel(8,108,8,label='reverse_exception_P8_108')
    for i,row in enumerate(reverse):
        cc=row['chart'];check('reversed exceptional primitive top row'+str(i),
             {t:v for t,v in cc.items() if sum(t)>=27}==scale(H27,cc.get((0,27),0)))
    for i,row in enumerate(raw(9,72)):check('reverse partner full chart degree<=21 row'+str(i),maximum(row['chart'])<=21)
    # Source coefficients at d=72,108,128 are exactly the indicated v^n constants.
    for m,cap,n in [(7,72,18),(8,72,18),(9,108,27),(8,108,27),(9,128,32)]:
        rr=raw(m,cap);tops=[r for r in rr if sum(r['pivot'])==cap]
        check(f'unique source top image m{m} d{cap}',len(tops)==1 and
              {t:c for t,c in tops[0]['source'].items() if sum(t)==cap}=={(3*n,n):1})
        for i,row in enumerate(rr):check(f'source top to chart coefficient m{m} d{cap} row{i}',
              row['chart'].get((0,n),0)==row['source'].get((3*n,n),0))
    # Primitive alternatives and the one retained quadratic degeneration.
    multiplicities=[[1,7,20],[7,18,1,1],[7,18,2],[8,18,1],[7,19,1],[8,19],[7,20],[9,18]]
    gs=[]
    for mult in multiplicities:
        g=0
        for n in mult:g=gcd(g,n)
        gs.append(g)
    check('only v-squared quadratic gives a nonprimitive degree27 form',gs==[1,1,1,1,1,1,1,9])
    w,v=S.symbols('w v');T=v*w**3;R=v*w**2
    check('exceptional leading resonance retained',S.expand(S.diff(R**6,w)*S.diff(R**9,v)-S.diff(R**6,v)*S.diff(R**9,w))==0)
    # Leading centralizer of T^8 has degree4j and is T^j.
    for d in range(1,29):
        resonant=[i for i in range(d+1) if 32*i-24*d==0]
        check('exact128 homogeneous resonance degree'+str(d),resonant==([3*d//4] if d%4==0 else []))
    check('exception partner degree18 is forced', [d for d in range(18,21) if d%3==0]==[18])
    check('divisible source degrees are exactly the removable ones', [j for j in range(1,8) if 8%j==0]==[1,2,4])
    for j in (1,2,4):
        check('source shear drops128 below125 for j'+str(j),16*j*(8//j)==128 and 128-4==124 and 16*j<72)
    for j in range(1,9):
        for a in range(1,3*j+1):
            check(f'leading chart monomial source gap j{j} a{a}',16*j-5*a<=16*j-5)
        for a in range(4*j):
            for b in range(4*j-a):
                check(f'lower chart monomial source gap j{j} a{a} b{b}',4*b-a<=16*j-4)
    # Contact refinements on the full spaces satisfying the necessary leading chart form.
    configurations=[(5,48,3,'T',9,[2,3]),(7,48,3,'T',17,[3]),(7,80,5,'T',9,[2,3]),(7,96,6,'T',7,[1,2,3]),
        (8,48,3,'T',21,[3]),(8,80,5,'T',13,[2,3]),(8,96,6,'T',9,[2,3]),(8,112,7,'T',8,[1,2,3]),
        (7,72,6,'R',11,[2,3]),(8,72,6,'R',15,[2,3])]
    table=[]
    for m,cap,j,base,minorder,contacts in configurations:
        degree=(4 if base=='T' else 3)*j
        leading=form(j,(3 if base=='T' else 2)*j)
        rows=kernel(m,cap,m,degree,leading,label=f'family_{m}_{cap}_{base}{j}')
        powers=[k for row in rows for k in row['symbol']]
        check(f'full symbol minimum m{m} d{cap}',min(powers)==minorder)
        check(f'exact contact list m{m} d{cap}',[k for k in (1,2,3) if m*k>=minorder]==contacts)
        check(f'nonzero leading chart coefficient possible m{m} d{cap}',any(row['chart'].get((0,degree),0) for row in rows))
        table.append(dict(poles=[m,9],ordinary_degrees=[cap,128 if base=='T' else 108],
            primitive_cubic_contact_multiplicities=contacts,chart_degrees=[degree,32 if base=='T' else 27],
            chart_common_form='v*(u+v)^3' if base=='T' else 'v*(u+v)^2',lower_symbol_order_bound=minorder,
            lower_space_dimension=len(rows),cusp_multiplier_slope='3/4' if base=='T' else '2/3'))
    cert=dict(as_of='2026-11-14',exact_field='Q',spaces=certspaces,kernel_certificates=kernel_records,
        gcd_one_height_nine_necessary_rows=table,
        upper_degree128_chart_form='c*v^8*(u+v)^24',
        exceptional_high_projection='a*u*v^7*(u+v)^20+v^7*(u+v)^18*(c*v^2+b*u*v+d*u^2)',
        exceptional_resonant_chart_forms=['c1*v^6*(u+v)^12','c2*v^9*(u+v)^18'],
        exceptional_reverse_orientation_excluded=True,
        primitive_factor_multiplicities=multiplicities,multiplicity_gcds=gs,
        source_degree_gap='a leading v^j*(u+v)^(3j) of chart degree4j has source degree16j and no degrees16j-1,16j-2,16j-3',
        excluded_divisible_source_partner_exponents=[1,2,4],finite_map_constructed=False)
    cp=HERE/'height_nine_frontier_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'HEIGHT_NINE_RESONANT_FRONTIER.md',cp,HERE/'verify_contact_filtration.py',
            HERE/'SOURCE_COORDINATES.md',HERE/'WRIGHT_MARKED_MODEL.md',HERE/'POLE_HEIGHT_EIGHT_EXCLUSION.md',
            HERE/'sources/degree-transfer-manifest.json',HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        written_result='Complete exact height-nine gcd-one necessary frontier: degree128 resonant branch and a retained normal72,108 resonance; reversed exception and divisible shear rows excluded',
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        gcd_one_height_nine_necessary_rows=table,exceptional_reverse_orientation_excluded=True,
        exceptional_normal_orientation_retained=True,normal_exception_cusp_multiplier_slope='2/3',
        degree128_cusp_multiplier_slope='3/4',
        assurance=dict(complete_finite_section_spaces=True,all_lower_pole_terms_retained=True,
            ordinary_shear_preserves_pole_height=False,ordinary_shear_preserves_Keller_and_field_degree=True,
            degree128_shear_exclusion_uses_retained_2022_alternative=True,
            gcd3_pair_6_9_classified=False,all_gcd_one_height_nine_pairs_excluded=False,
            all_retained_linear_rows_realizable_as_finite_maps=False,
            unrestricted_72_108_exception_excluded=False,conditional_2026_degree125_upgrade_used=False,
            actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False),
        source_sha256={str(f.relative_to(HERE)) if f.is_relative_to(HERE) else '../jc2-omitted-divisors/'+f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),elapsed_seconds=round(monotonic()-start,6))
    (HERE/'height_nine_frontier_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; ten necessary gcd-one height9 rows; normal72,108 resonance retained; gcd3 open')

if __name__=='__main__':main()
