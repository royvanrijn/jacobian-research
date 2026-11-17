#!/usr/bin/env python3
"""Exact bounded controls; the companion note supplies all-valuation proofs."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import platform

HERE = Path(__file__).resolve().parent
A, B = -10758382, 6600759744
S0 = F(1816347,819200)
KAPPA = F(42649,16384)
checks = Counter()
largest_root_modulus_exponent = 0


def check(kind, condition):
    if not condition:
        raise AssertionError(kind)
    checks[kind] += 1


def val(q):
    q = F(q)
    if not q:
        return float("inf")
    a,b = abs(q.numerator),q.denominator
    n = 0
    while a%3 == 0:
        a//=3
        n+=1
    while b%3 == 0:
        b//=3
        n-=1
    return n


def p3(n):
    return F(3**n) if n>=0 else F(1,3**(-n))


def modq(q, modulus):
    q=F(q)
    return q.numerator*pow(q.denominator,-1,modulus)%modulus


def elliptic_f(X):
    return X**3+A*X+B


def on_curve(P):
    return P is None or P[1]**2 == elliptic_f(P[0])


def add(P,Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x,y=P
    u,v=Q
    if x==u and y==-v:
        return None
    slope=(3*x*x+A)/(2*y) if P==Q else (v-y)/(u-x)
    xx=slope*slope-x-u
    return xx,-y+slope*(x-xx)


def mul(P,n):
    if n<0:
        return mul((P[0],-P[1]),-n)
    out=None
    while n:
        if n%2:
            out=add(out,P)
        P=add(P,P)
        n//=2
    return out


def cubic(X):
    v=F(3,6400)*(X+3984)
    return (-F(12,5),S0-v,F(10,9)*v*v-F(10,9)*S0*v+F(125,1024))


def discriminant(a,b,c):
    return a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c


def independent_root_test(coefficients):
    """Translate/scale independently, then perform the complete FMA3a test."""
    global largest_root_modulus_exponent
    a,b,c=coefficients
    choices=[]
    for shift in [F(0),F(1),F(2),F(4,5)]:
        aa=a+3*shift
        bb=b+2*a*shift+3*shift*shift
        cc=c+b*shift+a*shift*shift+shift**3
        bounds=[F(val(q),degree) for q,degree in [(aa,1),(bb,2),(cc,3)] if q]
        exponent=min(bounds).__floor__()
        choices.append((exponent,aa,bb,cc))
    exponent,aa,bb,cc=max(choices,key=lambda row:row[0])
    aa/=p3(exponent)
    bb/=p3(2*exponent)
    cc/=p3(3*exponent)
    if min(val(aa),val(bb),val(cc))<0:
        raise AssertionError("nonintegral independent cubic model")
    dd=val(discriminant(aa,bb,cc))
    if not isinstance(dd,int) or dd<0:
        raise AssertionError("bad cubic discriminant")
    N=dd//2+1
    largest_root_modulus_exponent=max(largest_root_modulus_exponent,N)
    roots=[0]
    modulus=1
    for level in range(N):
        next_modulus=3*modulus
        ca,cb,cc0=[modq(q,next_modulus) for q in [aa,bb,cc]]
        roots=[r+j*modulus for r in roots for j in range(3)
               if ((r+j*modulus)**3+ca*(r+j*modulus)**2
                   +cb*(r+j*modulus)+cc0)%next_modulus==0]
        if not roots:
            return False
        modulus=next_modulus
    return True


def local_point_data(X):
    f=elliptic_f(X)
    if not f:
        return None
    vv=val(f)
    if vv%2 or modq(f/p3(vv),3)!=1:
        return None
    return vv//2


def predicted_test(X):
    m=local_point_data(X)
    if m is None:
        raise AssertionError("input not a non-torsion local elliptic point")
    if val(X)<0:
        n=-val(X)//2
        unit_y_square=modq(elliptic_f(X)/p3(2*m),9)
        unit_y=next(r for r in [1,2,4,5,7,8] if r*r%9==unit_y_square)
        unit_t=modq(-X*p3(2*n),9)*pow(unit_y,-1,9)%9
        return n%3==0 and unit_t in [1,8]
    reduction=modq(X,3)
    if reduction==0:
        return False
    if reduction==2:
        return m>=2
    if m>=2:
        return m>=3
    x=X+3552
    z=F(3,819200)*(128*x-25865)
    k=val(z)
    if k<4:
        raise AssertionError("branch-disc valuation bound failed")
    return k%3==0 and modq(KAPPA*z/p3(k),9) in [1,8]


def poly_add(*ps):
    return [sum((p[i] if i<len(p) else 0) for p in ps)
            for i in range(max(map(len,ps)))]


def poly_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


def quadratic_mul(u,v,d,modulus=9):
    a,b=u
    c,e=v
    return ((a*c+d*b*e)%modulus,(a*e+b*c)%modulus)


def quadratic_pow(u,n,d):
    out=(1,0)
    while n:
        if n%2:
            out=quadratic_mul(out,u,d)
        u=quadratic_mul(u,u,d)
        n//=2
    return out


def quotient_count(N,small):
    modulus=3**(N-1)
    domain=range(0,modulus,3 if small else 1)
    total=allowed=0
    for a in range(2):
        for b in range(2):
            for k in domain:
                total+=1
                if (a,b)==(1,0):
                    continue                         # X=0
                if (a,b)==(0,1):                     # X=2
                    allowed+=int(k%3==0)
                elif (a,b)==(1,1):                  # X=1
                    allowed+=int(k%9==0) if small else 1
                elif k==0:                         # unresolved formal tail
                    allowed+=1
                else:
                    depth=val(k)+1
                    unit=(k//3**(depth-1))%9
                    allowed+=int(depth%3==0 and unit in [1,8])
    return total,allowed


def hensel_quadratic_root(start,modulus):
    root=start
    current=3
    while current<modulus:
        root=next(root+j*current for j in range(3)
                  if ((root+j*current)**2-10656*(root+j*current)+27091730)%(3*current)==0)
        current*=3
    return root


def finite_add(P,Q,p):
    """Translated model y^2=x^3-10656*x^2+27091730*x."""
    if P is None:
        return Q
    if Q is None:
        return P
    x,y=P
    u,v=Q
    if x==u and (y+v)%p==0:
        return None
    if P==Q:
        slope=(3*x*x-2*10656*x+27091730)*pow(2*y,-1,p)%p
    else:
        slope=(v-y)*pow(u-x,-1,p)%p
    xx=(slope*slope+10656-x-u)%p
    return xx,(-y+slope*(x-xx))%p


def finite_mul(P,n,p):
    out=None
    while n:
        if n%2:
            out=finite_add(out,P,p)
        P=finite_add(P,P,p)
        n//=2
    return out


def saturation_certificate():
    P1=(F(373248,49),F(59459616,343))
    P2=(F(3442688,529),F(304670016,12167))
    expected={41:{"count":48,"P1":(39,37),"P2":(0,0),
                  "A":(37,37),"B":None},
              53:{"count":60,"P1":(21,16),"P2":(33,51),
                  "A":(13,45),"B":(13,8)}}
    data={}
    for p in [41,53]:
        check("saturation prime good reduction",
              (27091730**2*((-10656)**2-4*27091730))%p!=0)
        points=[(x,y) for x in range(p) for y in range(p)
                if (y*y-x**3+10656*x*x-27091730*x)%p==0]
        row=expected[p]
        check("complete saturation finite-field point count",1+len(points)==row["count"])
        pp1=tuple(modq(q,p) for q in P1)
        pp2=tuple(modq(q,p) for q in P2)
        check("exact saturation reduction coordinates",pp1==row["P1"] and pp2==row["P2"])
        scale=row["count"]//3
        aa,bb=finite_mul(pp1,scale,p),finite_mul(pp2,scale,p)
        check("exact saturation three-primary images",aa==row["A"] and bb==row["B"])
        check("saturation nonzero image has exact order three",aa is not None and finite_mul(aa,3,p) is None)
        if p==53:
            check("second saturation images are negatives",finite_add(aa,bb,p) is None)
        for P in points:
            check("complete finite-field group exponent control",finite_mul(P,row["count"],p) is None)
        data[str(p)]={"cardinality":row["count"],"P1":pp1,"P2":pp2,
                      "scale":scale,"image_P1":aa,"image_P2":bb}
    for a in range(3):
        for b in range(3):
            vanishes=[]
            for p in [41,53]:
                aa=expected[p]["A"]
                bb=expected[p]["B"]
                vanishes.append(finite_add(finite_mul(aa,a,p),finite_mul(bb,b,p),p) is None)
            check("complete independent mod-three saturation forms",
                  all(vanishes)==(a==0 and b==0))
    check("saturation forms determinant is a unit",(-1)%3==2)
    return data


def replay():
    saturation=saturation_certificate()
    z=[-F(3*25865,819200),F(3*128,819200)]
    ss=[F(5200),F(-1)]
    qh=[F(27091730),F(-10656),F(1)]
    check("exact Cardano norm identity",
          poly_add(poly_mul(ss,ss),[-c for c in qh],[-2*c for c in [-F(25865),F(128)]])==[0,0,0])
    lhs=poly_add([KAPPA],[-F(10,9)*c for c in z])
    check("exact Cardano sum identity",lhs==[F(5200,1920),-F(1,1920)])
    check("Cardano product normalization",F(2*819200,3*3840**2)==F(1,27))
    check("quadratic-root residue for X1",(67**2-10656*67+27091730)%3**5==0)
    check("X1 unit congruences",modq(F(128*67-25865,9*819200),9)==7
          and F(5200-67,3)%9==1 and modq(F(7,1280)*F(5200-67,3),9)==8)
    check("X2 unit congruence",modq(F((128*5-25865)*(5200-5),1048576000),9)==1)
    check("kappa and formal-leading unit",modq(KAPPA,9)==4 and 160**3%9==1)
    units=[u for u in range(9) if u%3]
    check("complete Q3 cube-unit residues",{u**3%9 for u in units}=={1,8})
    for d in [2,5,8]:
        field_units=[(a,b) for a in range(9) for b in range(9) if a%3 or b%3]
        cubes={quadratic_pow(u,3,d) for u in field_units}
        criterion={u for u in field_units if quadratic_pow(u,8,d)==(1,0)}
        check("complete unramified quadratic cube units d%d"%d,cubes==criterion and len(cubes)==8)
    for a in [0,9,18]:
        for b in [0,9,18]:
            for c in range(1,27):
                if c%3:
                    check("normalized cubic root modulus nine",
                          val(discriminant(a,b,c))==3 and
                          independent_root_test((F(a),F(b),F(c)))==(c%9 in [1,8]))

    profile=Counter()
    seen=set()
    def control(X,kind):
        X=F(X)
        if X in seen or local_point_data(X) is None:
            return
        co=cubic(X)
        if not discriminant(*co):
            return
        seen.add(X)
        predicted=predicted_test(X)
        actual=independent_root_test(co)
        check("independent complete local cubic test "+kind,predicted==actual)
        disc="O" if val(X)<0 else "X%d"%modq(X,3)
        profile[(disc,"tested")]+=1
        profile[(disc,"soluble")]+=int(actual)

    # Declared rational X domain, including non-rational Q3 Y coordinates.
    for n in range(4):
        for a in range(-60,61):
            if a:
                control(F(a,3**(2*n)),"rational-X")
    # Both finite simple-root discs at all these controlled depths.
    for residue in [1,2]:
        for m in range(1,6):
            root=hensel_quadratic_root(residue,3**(2*m+2))
            for u in range(1,9):
                if u%3:
                    control(F(root+3**(2*m)*u-3552),"finite-root-disc")
    # Every indicated branch depth, including soluble m=1 examples.
    for k in range(4,11):
        for u in range(1,9):
            if u%3:
                control(F(25865,128)+6400*p3(k-1)*u-3552,"branch-disc")
    for n in range(1,7):
        for a in [1,4,7,10,13,16,19,22,25]:
            control(F(a,3**(2*n)),"formal-disc")

    P1=(F(199200,49),F(59459616,343))
    P2=(F(1563680,529),F(304670016,12167))
    T0=(F(-3552),F(0))
    check("exact rational witnesses",all(on_curve(P) for P in [P1,P2,T0]))
    translated_P1=add(P1,T0)
    check("P1 generator belongs to the proven local subgroup",
          val(P1[0]+3552)==6 and
          translated_P1[0]+3552==F(27091730)/(P1[0]+3552) and
          val(-translated_P1[0]/translated_P1[1])==3)
    check("P2 generator local two-torsion displacement depth",val(P2[1])==2)
    for m in range(-2,3):
        for n in range(-2,3):
            if m>0 or (m==0 and n>0):
                for k in [0,1]:
                    P=add(mul(P1,m),mul(P2,n))
                    if k:
                        P=add(P,T0)
                    if P is not None:
                        check("bounded rational group law",on_curve(P))
                        control(P[0],"rational-MW")

    tables=[]
    for J in [1,2,3]:
        for small in [True,False]:
            total,retained=quotient_count(3*J,small)
            expected=(F(55,156)+F(45,26*27**J) if small
                      else F(53,156)+F(15,26*27**J))
            check("complete refined finite quotient",F(retained,total)==expected)
            tables.append({"J":J,"smaller_image":small,"total":total,
                           "retained":retained,"fraction":str(expected)})
    check("first refined upper coefficient",F(55,156)+F(45,26*27)==F(5,12))
    check("smaller-image local measure",(F(1,13)+F(1,3)+1)/4==F(55,156))
    check("full-image branch-disc local measure",F(1,9)+F(2,3)*F(1,39)==F(5,39))
    check("full-image total local measure",(F(1,39)+F(5,39)+F(1,3))/4==F(19,156))
    for J in range(1,20):
        small=F(1,3)+sum((F(1,4)*F(2,9)*F(1,3**(3*k-2))
                         for k in range(1,J)),F(0))+F(1,4*3**(3*J-2))
        check("geometric-sum refined finite bound",small==F(55,156)+F(45,26*27**J))
    check("three-adic class corruption rejected",not independent_root_test((F(0),F(0),F(2))))
    check("three-adic unit corruption rejected",quadratic_pow((1,3),8,5)!=(1,0))
    return {"passed":True,"exact_controls":sum(checks.values()),
            "checks_by_kind":dict(checks),"complete_finite_quotients":tables,
            "local_control_profile":{d:{k:profile[(d,k)] for k in ["tested","soluble"]}
                                     for d in ["O","X0","X1","X2"]},
            "largest_independent_cubic_root_modulus_exponent":largest_root_modulus_exponent,
            "bounded_domains":{"rational_X":"X=a/3^(2n), 0<=n<=3, -60<=a<=60, a!=0; retained only locally square elliptic RHS",
               "finite_root_discs":"residue1or2, 1<=m<=5, u=1..8 not divisible3",
               "branch_discs":"4<=k<=10, u=1..8 not divisible3",
               "formal_discs":"1<=n<=6, a in{1,4,7,10,13,16,19,22,25}",
               "rational_MW":"|m|,|n|<=2, sign-canonical, rational T0 added or not; no basis assumption"},
            "finite_upper_fraction":"5/12","limsup_upper_fraction":"55/156",
            "local_only_leading_fraction":"55/156",
            "ambient_full_group_local_measure":"19/156",
            "finite_three_saturation_certificate":saturation,
            "actual_MW_closure_case_determined":True,
            "actual_MW_closure_index":3,
            "saturation_at_other_primes_claimed":False,
            "all_valuation_quantifiers_supplied_by_written_proof":True,
            "bounded_tests_prove_unbounded_scope":False,
            "Hasse_infinitude_claimed":False,"uniform_J_error_claimed":False,
            "formal_verification":False,"external_review":False,"independent_replay":False,
            "python":platform.python_version(),"arithmetic":"standard-library exact integers and fractions"}


if __name__=="__main__":
    report=replay()
    files=["ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md",Path(__file__).name]
    report["input_sha256"]={f:sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    destination=HERE/"elliptic_exception_complete_three_adic_classification_checks.json"
    destination.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"passed":True,"exact_controls":report["exact_controls"],
                      "local_control_profile":report["local_control_profile"],
                      "input_sha256":report["input_sha256"]},indent=2))
