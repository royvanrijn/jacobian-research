#!/usr/bin/env python3
"""Exact bounded controls for the all-index digit-weight proof; not a rank test."""
import hashlib, json, math
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
tests=[]
def check(name,ok,**data):
    if not ok: raise AssertionError((name,data))
    tests.append(name)
def s7(n):
    s=0
    while n: s,n=s+n%7,n//7
    return s
def vp(n):
    if not n: return 10**9
    n=abs(n); v=0
    while n%7==0:n//=7;v+=1
    return v
def trim(a):
    while a and a[-1]==0:a.pop()
    return a
def add(a,b):
    o=[0]*max(len(a),len(b))
    for i,x in enumerate(a):o[i]+=x
    for i,x in enumerate(b):o[i]+=x
    return trim([x%7 for x in o])
def sc(a,c): return trim([x*c%7 for x in a])
def deriv_v(a):return trim([i*x%7 for i,x in enumerate(a)])
def vtimes(a):return [0]+a if a else []
def raise_poly(a,c):
    o=[0]*(len(a)+1)
    for k,x in enumerate(a):
        o[k]+=(c-k)*x
        o[k+1]+=x
    return o
def actual_B(n):
    a=[1]
    for m in range(n): a=raise_poly(a,n+6*m+5)
    return a
orders=list(range(1,51))+[55,97,118,145,146,196,197,294,295,342,343,344,399,400]
examples=[]
actual_records={}
tight=0
cut_rejected=False
extra_division_rejected=False
for n in orders:
    ss=s7(n); ee=[(ss-a+4)//6 for a in range(6)]
    den=6**n*math.factorial(n); vden=vp(den)
    denunit=den//7**vden%7; invden=pow(denunit,-1,7)
    actual=actual_B(n)
    columns=[]
    for q in range(7):
        rows=[]
        for a in range(6):
            row=[]
            for j in range(a,len(actual),6):
                ell=(j-a)//6
                numer=actual[j]*math.factorial(a)
                value=vp(numer)+ell-vden
                check("actual_coefficient_bound",value>=ee[a],
                      N=n,Qpower=q,row=a,j=j,value=value,bound=ee[a])
                if value==ee[a]:
                    tight+=1
                    unit=(numer//7**vp(numer))%7
                    row.append(unit*invden%7)
                    if not extra_division_rejected:
                        check("reject_unjustified_extra_division",value<ee[a]+1)
                        extra_division_rejected=True
                else:row.append(0)
            rows.append(trim(row))
        columns.append(rows)
        actual=raise_poly(actual,n+5)
    actual_records[n]=columns
    cut=(ss-1)%6
    for q in range(6):
        for a in range(6):
            f=columns[q][a]
            pred=add(sc(f,n+5-a),deriv_v(f))
            prev=columns[q][a-1 if a else 5]
            if a!=cut:
                cross=sc(prev,a) if a else vtimes(prev)
                pred=add(pred,cross)
            check("actual_normalized_raising",pred==columns[q+1][a],
                  N=n,Qpower=q,row=a,cut=cut)
            if not cut_rejected and a==cut and prev:
                wrong=add(pred,sc(prev,a) if a else vtimes(prev))
                if wrong!=pred:
                    check("reject_restored_cut_link",wrong!=columns[q+1][a])
                    cut_rejected=True
    if ss>=8:
        examples.append({"N":n,"digit_sum":ss,"row_bounds":ee,
                         "seed_nonzero_rows":[a for a in range(6) if columns[0][a]]})
for ss in range(1,201):
    ee=[(ss-a+4)//6 for a in range(6)]
    links=[1+ee[5]-ee[0]]+[ee[a-1]-ee[a] for a in range(1,6)]
    check("cut_arithmetic",set(links)<={0,1} and sum(links)==1
          and links[(ss-1)%6]==1, digit_sum=ss)
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,0)+x*y
    return {i:x for i,x in out.items() if x}
f={0:1,1:1,7:2,8:7,20:49,49:3,342:343}
check("test_series_in_F1",all(6*vp(x)>=s7(i)-1 for i,x in f.items()))
power={0:1}
for _ in range(7):power=mul(power,f)
for i,x in power.items():
    check("seventh_power_F1",6*vp(x)>=s7(i)-1, degree=i)
bad=mul({1:1},{1:1})
check("reject_ordinary_square_F1",any(6*vp(x)<s7(i)-1 for i,x in bad.items()))
check("negative_controls_executed",cut_rejected and extra_division_rejected)
def vadd(a,b):return [add(x,y) for x,y in zip(a,b)]
def vscale(a,c):return [sc(x,c) for x in a]
def link(a,cut):
    out=[[] for _ in range(6)]
    for row in range(1,6):
        if row!=cut:out[row]=sc(a[row-1],row)
    if cut!=0:out[0]=vtimes(a[5])
    return out
def diagonal(a,kappa):
    return [add(sc(x,kappa-row),deriv_v(x)) for row,x in enumerate(a)]
def raising(a,kappa,cut):return vadd(diagonal(a,kappa),link(a,cut))
def gauge(a,cut,sign=1):
    out=[x[:]for x in a]; term=a
    for j in range(1,6):
        term=link(term,cut)
        factor=pow(math.factorial(j),-1,7)*pow(sign,j,7)%7
        out=vadd(out,vscale(term,factor))
    return out
for kappa in range(7):
    for cut in range(6):
        for row in range(6):
            for degree in [0,1,6,7,8,13]:
                vec=[[]for _ in range(6)]
                vec[row]=[0]*degree+[1]
                nn=vec
                for _ in range(6):nn=link(nn,cut)
                check("finite_link_nilpotence",not any(nn))
                check("finite_gauge_inverse",gauge(gauge(vec,cut,-1),cut)==vec)
                left=raising(vec,kappa,cut)
                right=gauge(diagonal(gauge(vec,cut,-1),kappa),cut)
                check("finite_operator_gauge",left==right)
                power=vec
                for _ in range(7):power=raising(power,kappa,cut)
                check("artin_schreier_operator",power==left)
for n,columns in actual_records.items():
    cut=(s7(n)-1)%6; kappa=(n+5)%7
    seed=gauge(columns[0],cut,-1)
    eig=[]
    for lam in range(7):
        section=[]
        for row,poly in enumerate(seed):
            section.append(trim([x if (degree+kappa-row-lam)%7==0 else 0
                                 for degree,x in enumerate(poly)]))
        eig.append(section)
    for j in range(5):
        middle=[[]for _ in range(6)]
        for lam in range(7):
            middle=vadd(middle,vscale(eig[lam],pow(lam,j,7)))
        check("actual_eigensection_factorization",
              gauge(middle,cut)==columns[j],N=n,Qpower=j)
def evalpoly(poly,x):return sum(a*pow(x,i,7)for i,a in enumerate(poly))%7
literal=[[6,1,6,1,4],[4,5,5,3,2],[4,5,3,1,1],
         [0,5,1,4,4],[0,0,6,5,3],[0,0,0,2,0]]
observed=[[evalpoly(actual_records[8][q][a],6)for q in range(5)]for a in range(6)]
check("actual_rank_countercontrol_matrix",observed==literal)
from itertools import combinations
def det7(mat):
    mat=[r[:]for r in mat]; out=1
    for j in range(len(mat)):
        row=next((i for i in range(j,len(mat))if mat[i][j]%7),None)
        if row is None:return 0
        if row!=j:mat[j],mat[row]=mat[row],mat[j];out=-out
        unit=mat[j][j]%7; out=out*unit%7; inv=pow(unit,-1,7)
        for i in range(j+1,len(mat)):
            factor=mat[i][j]*inv%7
            for k in range(j,len(mat)):mat[i][k]=(mat[i][k]-factor*mat[j][k])%7
    return out%7
for rows in combinations(range(6),5):
    check("countercontrol_all_maximal_minors_zero",
          det7([literal[i]for i in rows])==0,rows=rows)
check("countercontrol_rank_exactly_four",det7([r[:4]for r in literal[:4]])==1)
H8at1=7
for k in range(1,9):
    H8at1+=math.comb(8,k)*math.factorial(8+6*k-1)//math.factorial(7)//7**(k-1)
check("countercontrol_local_eisenstein",H8at1%2401==266 and vp(H8at1)==1)
H8leading=math.factorial(55)//math.factorial(7)//7**7
check("countercontrol_leading_unit",vp(H8leading)==0)

pins={}
for p in [Path(__file__),HERE/"GAP_SIX_DIGIT_NORMALIZATION_AND_FINITE_GAUGE.tex",
          BASE/"research/factorial/image-rigidity/GAP_SIX_TWO_POWER_PLUS_ONE_SFC.tex",
          BASE/"research/factorial/image-rigidity/GAP_SIX_BOUNDARY_COUPLED_FROBENIUS.tex",
          BASE/"research/factorial/image-rigidity/ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
          BASE/"research/factorial/image-rigidity/GAP_SIX_SMALL_OFFSET_FAMILIES.tex",
          BASE/"research/factorial/image-rigidity/PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex"]:
    pins[str(p.relative_to(BASE))]=hashlib.sha256(p.read_bytes()).hexdigest()
receipt={"status":"PASS","controls":len(tests),"orders":orders,
         "Q_powers":list(range(7)),"tight_actual_coefficients":tight,
         "examples":examples,
         "scope":"Bounded exact actual-coefficient and normalized-raising controls; the unbounded theorem is proved in the accompanying TeX.",
         "independent_whole_proof_replay":False,"formal_verification":False,
         "external_review":False,"source_sha256":pins}
(HERE/"gap_six_digit_normalization_and_finite_gauge_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({k:receipt[k] for k in ["status","controls","tight_actual_coefficients"]}))
