#!/usr/bin/env python3
"""Complete F29/mod29² certificate for the reciprocal degree26 split (4,8).

The companion note supplies the uniform integrality and whole-closure proof.
This checker independently solves the four linear leading equations for Q,
reconstructs every retained integer equation, runs a complete C++ field census,
and verifies a literal left-null obstruction at its sole strong residue point.
"""
from __future__ import annotations
import argparse,hashlib,json,platform,shutil,subprocess,tempfile,time
from pathlib import Path
import sympy as s
from verify_mixed_reciprocal_coupled_ode import from_critical,chebyshev,w,r

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CPP=HERE/'verify_degree_twenty_six_mixed_4_8_finite_field.cpp'
SYSTEM=HERE/'degree26_mixed_4_8_system.json'
CERT=HERE/'degree_twenty_six_mixed_4_8_global_certificate.json'
VARIABLES=s.symbols('v3 v2 v1 v0 c d r')
PRIME=29
RAW_POINTS=[[10,13,16,4,20,7,15],[18,3,2,17,6,9,12],
            [20,27,2,16,15,20,9],[22,0,8,3,4,4,16]]
SURVIVOR=[20,27,2,16,15,20,9]
LEFT_NULL=[18,8,27,21,0,0,23,17,11]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def primitive(expr):
    numerator,denominator=s.fraction(s.cancel(expr))
    assert not denominator.free_symbols,'Variable denominator in integer equations'
    _,integer_poly=s.Poly(numerator,*VARIABLES,domain=s.QQ).clear_denoms(convert=True)
    return integer_poly.primitive()[1]

def literal_system():
    v3,v2,v1,v0,c,d,rr=VARIABLES
    assert rr==r
    P=w*w+c*w+d;V=w**4+v3*w**3+v2*w*w+v1*w+v0
    C,D=chebyshev(6)
    target=s.Poly(s.expand(25*(C[6]-r*C[5])*P),w)
    q=s.symbols('q0:4');unknown_Q=w**4+sum(q[j]*w**j for j in range(4))
    def operator(Q):
        return s.Poly(s.expand(4*(w-2)*V*s.diff(Q,w)
             +(V+2*(w-2)*s.diff(V,w))*Q),w)
    error=operator(unknown_Q)-target
    matrix,rhs=s.linear_eq_to_matrix([error.nth(j)for j in range(4,8)],q)
    assert s.expand(matrix.det())==9*13*17*21
    solved=[s.expand(x)for x in matrix.inv()*rhs]
    assert all(not s.fraction(s.cancel(x))[1].free_symbols for x in solved)
    Q=w**4+sum(solved[j]*w**j for j in range(4))
    error=operator(Q)-target
    assert error.degree()<=3
    _,F,X,Y=from_critical(6,s.expand(P*Q))
    assert s.expand(X+2*(w+2)*s.diff(X,w)-25*(D[6]-r*D[5])*P*Q)==0
    assert s.expand(Y+2*(w-2)*s.diff(Y,w)-25*(C[6]-r*C[5])*P*Q)==0
    domain=s.QQ.poly_ring(*VARIABLES)
    Xpoly=s.Poly(X,w,domain=domain)
    weak=Xpoly.rem(s.Poly(P,w,domain=domain))
    strong=Xpoly.rem(s.Poly(P*P,w,domain=domain))
    G=s.cancel(sum(F[j]*r**j for j in range(14))/(2*r))
    raw=[primitive(error.nth(j))for j in range(4)]+[primitive(weak.nth(j))for j in range(2)]+[primitive(G)]
    extra=[primitive(strong.nth(j))for j in range(4)]
    assert all(all(ex[4]+ex[5]<=1 and ex[6]<=1 for ex,_ in poly.terms())for poly in raw[:4])
    assert [poly.total_degree()for poly in raw]==[5,5,5,5,11,11,12]
    assert [len(poly.terms())for poly in raw]==[73,90,88,60,652,712,240]
    record=dict(schema_version=1,variables=list(map(str,VARIABLES)),P=str(P),V=str(V),
        Q_coefficients=[str(s.expand(Q).coeff(w,j))for j in range(5)],
        constant_top_down_divisors=[21,17,13,9],
        raw_equation_order=['Y_compatibility_w'+str(j)for j in range(4)]+['X_mod_P_w0','X_mod_P_w1','F(r)/(2r)'],
        raw_equations=[str(poly.as_expr())for poly in raw],
        strong_extra_order=['X_mod_P_squared_w'+str(j)for j in range(4)],
        strong_extra_equations=[str(poly.as_expr())for poly in extra],
        strong_order=['Y_compatibility_w'+str(j)for j in range(4)]+['X_mod_P_squared_w'+str(j)for j in range(4)]+['F(r)/(2r)'],
        raw_term_counts=[len(poly.terms())for poly in raw],
        strong_equations_are_necessary_without_residue_nondegeneracy=True)
    return record,raw,extra

def value(poly,point):return int(poly.eval(dict(zip(VARIABLES,point))))

def rank_mod(matrix):
    matrix=[[int(a)%PRIME for a in row]for row in matrix]
    row=0
    for column in range(len(matrix[0])):
        pivot=next((i for i in range(row,len(matrix))if matrix[i][column]),None)
        if pivot is None:continue
        matrix[row],matrix[pivot]=matrix[pivot],matrix[row]
        scale=pow(matrix[row][column],-1,PRIME)
        matrix[row]=[a*scale%PRIME for a in matrix[row]]
        for i in range(len(matrix)):
            if i==row:continue
            multiplier=matrix[i][column]
            matrix[i]=[(a-multiplier*b)%PRIME for a,b in zip(matrix[i],matrix[row])]
        row+=1
        if row==len(matrix):break
    return row

def obstruction(raw,extra):
    equations=raw[:4]+extra+raw[-1:]
    evaluations=[value(poly,SURVIVOR)for poly in equations]
    assert all(a%PRIME==0 for a in evaluations)
    jacobian=[[value(poly.diff(variable),SURVIVOR)%PRIME for variable in VARIABLES]for poly in equations]
    rhs=[(-a//PRIME)%PRIME for a in evaluations]
    assert rhs==[1,14,4,15,11,23,20,3,22]
    assert all(sum(LEFT_NULL[i]*jacobian[i][j]for i in range(9))%PRIME==0 for j in range(7))
    assert sum(a*b for a,b in zip(LEFT_NULL,rhs))%PRIME==1
    assert rank_mod(jacobian)==6 and rank_mod([row+[b]for row,b in zip(jacobian,rhs)])==7
    checks=[]
    for point in RAW_POINTS:
        assert all(value(poly,point)%PRIME==0 for poly in raw)
        residues=[value(poly,point)%PRIME for poly in extra]
        assert (all(a==0 for a in residues))==(point==SURVIVOR)
        checks.append(dict(point=point,strong_extra_residues=residues))
    return dict(schema_version=1,prime=PRIME,raw_points=RAW_POINTS,survivor=SURVIVOR,
        strong_filter=checks,Jacobian=jacobian,rhs=rhs,left_null=LEFT_NULL,
        Jacobian_rank=6,augmented_rank=7,left_null_rhs=1,
        conclusion='The unique strong F29 point has no lift modulo29².')

def census(raw,extra):
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('The complete field census requires g++')
    payload=[str(PRIME)]
    for poly in raw+extra:
        payload.append(str(len(poly.terms())))
        for ex,co in poly.terms():
            assert abs(int(co))<2**63 and max(ex)<32
            payload.append(' '.join(map(str,[int(co),*ex])))
    scratch=ROOT/'.scratch';scratch.mkdir(exist_ok=True)
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='p2q4-census-',dir=scratch)as temporary:
        binary=Path(temporary)/'census'
        subprocess.run([compiler,'-O3','-std=c++17',str(CPP),'-o',str(binary)],check=True,capture_output=True,text=True)
        run=subprocess.run([str(binary)],input='\n'.join(payload)+'\n',capture_output=True,text=True,check=True)
    result=json.loads(run.stdout)
    assert result['outer_count']==PRIME**4==707281
    assert result['first_four_candidates']==25764
    assert result['consistent_rank_counts']==[0,0,0,36,706443]
    assert result['raw_count']==4 and sorted(result['raw_points'])==sorted(RAW_POINTS)
    assert result['strong_count']==1 and result['strong_points']==[SURVIVOR]
    result['compiler']=subprocess.run([compiler,'--version'],text=True,capture_output=True,check=True).stdout.splitlines()[0]
    result['elapsed_seconds']=round(time.monotonic()-started,3)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-system',type=Path)
    parser.add_argument('--write-certificate',type=Path)
    parser.add_argument('--write-receipt',type=Path)
    args=parser.parse_args()
    record,raw,extra=literal_system()
    if args.write_system:args.write_system.write_text(json.dumps(record,indent=2)+'\n')
    else:assert json.loads(SYSTEM.read_text())==record,'Retained integer system differs from reconstruction'
    certificate=obstruction(raw,extra)
    if args.write_certificate:args.write_certificate.write_text(json.dumps(certificate,indent=2)+'\n')
    else:assert json.loads(CERT.read_text())==certificate,'Retained no-lift certificate differs'
    report=dict(schema_version=1,result='PASS',python=platform.python_version(),sympy=s.__version__,
        checker_sha256=sha(__file__),cpp_sha256=sha(CPP),input_sha256=sha(args.write_system or SYSTEM),
        certificate_sha256=sha(args.write_certificate or CERT),
        coefficient_reconstruction='Exact independent linear solve for Q and Euler reconstruction',
        complete_census=census(raw,extra),no_lift_certificate=certificate,
        proof_inputs={name:sha(HERE/name)for name in ['verify_mixed_reciprocal_coupled_ode.py','MIXED_RECIPROCAL_COUPLED_ODE.md','RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md']},
        scope='Complete F29 strong-system census and literal mod29² obstruction. Uniform raw-parameter integrality and monic quotient division in the written proof exclude every Q29 point of the strong closure, including denominators and degenerate reduction. This excludes both reciprocal (4,8) orientations overQ; balanced and nonreciprocal cases remain open.',
        independent_replay=False,formal_verification=False,external_review=False)
    bad=LEFT_NULL.copy();bad[0]=(bad[0]+1)%PRIME
    assert any(sum(bad[i]*certificate['Jacobian'][i][j]for i in range(9))%PRIME for j in range(7))
    report['corrupted_left_null_rejected']=True
    if args.write_receipt:args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
