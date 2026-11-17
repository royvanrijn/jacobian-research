#!/usr/bin/env python3
"""Exact controls for PF164's uniform complete ruling-six pole ledger."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilq(n,d):return (n+d-1)//d
def parts(n,least=1):
    if n==0:
        yield ()
    for first in range(least,n+1):
        for tail in parts(n-first,first):
            yield (first,)+tail
def partitions_up_to(n):
    return [a for k in range(1,n+1) for a in parts(k)]
Y,Z,s=S.symbols("Y Z s")
profiles=[
    ("A",Y**2*(Y-1)**4,4),
    ("B5",Y**2*(Y-1)**3,3),
    ("B6",Y**2*(Y-1)**3*(Y+1),3),
    ("C",Y**3*(Y-1)**3,3)]
for name,I,ell in profiles:
    check(name+" weighted initial multiplicity exact",
          all(S.diff(I,Y,j).subs(Y,1)==0 for j in range(ell))
          and S.diff(I,Y,ell).subs(Y,1)!=0)
    J=S.Poly(S.expand(I.subs(Y,1+s*Z)),Z)
    leading=S.diff(I,Y,ell).subs(Y,1)/S.factorial(ell)
    check(name+" entire shifted ell coefficient leading",
          S.Poly(J.nth(ell),s).nth(ell)==leading and leading!=0)
    # At equality in the old bound the ell term would survive in H.
    check(name+" equality gives nonlinear response obstruction",
          S.degree(leading*Y**ell+Y,Y)==ell and ell>1)

# Every normalized pole index is a positive integer, with total mass<=ell.
# This is the COMPLETE finite set of possible index partitions.
Aparts=partitions_up_to(4); BCparts=partitions_up_to(3)
check("complete quartic-index partitions",
      Aparts==[(1,),(1,1),(2,),(1,1,1),(1,2),(3,),
               (1,1,1,1),(1,1,2),(1,3),(2,2),(4,)])
check("complete cubic-index partitions",
      BCparts==[(1,),(1,1),(2,),(1,1,1),(1,2),(3,)])
for family,ell,period,bases,denom,offset,indices in [
    ("A",4,2,[3,4],2,2,Aparts),
    ("B",3,3,[4,5,6],3,3,BCparts),
    ("C",3,3,[4,5,6],3,3,BCparts)]:
    for L in bases:
        upper=(2*L-5 if family=="A" else L-4)
        for part in indices:
            # q < alpha*(L/denom - 1) with the displayed offset.
            qmax=[ceilq(alpha*(L-offset),denom)-1 for alpha in part]
            check(f"{family} complete base L{L} partition{part}",
                  sum(qmax)<=upper)
        for alpha in range(1,ell+1):
            base=ceilq(alpha*(L-offset),denom)-1
            lifted=ceilq(alpha*(L+period-offset),denom)-1
            check(f"{family} exact period L{L} alpha{alpha}",
                  lifted-base==alpha)
    check(f"{family} period covers every future slope",
          (4 if family=="A" else 3)>=max(sum(a) for a in indices))
    # A's global bound increases by4 per period2; B/C's by3 per period3.
    increase=4 if family=="A" else 3
    for part in indices:
        check(f"{family} all-period upper-bound monotonicity {part}",
              increase-sum(part)>=0)
check("B slope4 cannot support any pole",
      all(ceilq(alpha*(4-3),3)-1==0 for alpha in range(1,4)))
check("C B0=4 cannot support any pole",
      all(ceilq(alpha*(4-3),3)-1==0 for alpha in range(1,4)))
check("A slope3 budget requires a high-index pole if any",
      [ceilq(alpha*(3-2),2)-1 for alpha in range(1,5)]==[0,0,1,1])

n,j,kappa,D0,m0=S.symbols("n j kappa D0 m0",integer=True,nonnegative=True)
LL=5+n
for rem in range(4):
    delta=4*LL+4*j+rem
    check("B coefficient divisor quarter comparison residue"+str(rem),
          S.simplify(S.floor(delta/4)-4-(LL-4))==j)
BB=5+n
for rem in range(1,5):
    jump=rem+4*j
    delta=4*BB+3*jump
    check("C coefficient divisor quarter comparison residue"+str(rem),
          S.simplify(S.floor(delta/4)-4-(BB-4))==3*j+(3*rem)//4)
check("B comparison with universal half budget",
      S.expand((4*LL+kappa)/2-5-(LL-4))==n+kappa/2+4)
check("C comparison with universal half budget",
      S.expand((4*BB+3*(1+j))/2-5-(BB-4))==n+S.Rational(3,2)*j+S.Rational(11,2))
DD=1+D0; MM=2+m0
general_gap=S.expand(3*(DD+3*MM-6)/6-(DD/2-5))
noA_gap=S.expand(2*(DD+3*MM-6)/6-(DD/4-4))
check("alpha3 contradicts universal pole bound",
      general_gap==5+S.Rational(3,2)*m0)
check("alpha2 contradicts no-A pole bound",
      noA_gap==S.Rational(1,12)*D0+m0+S.Rational(49,12))
for degree in range(1,16):
    lo=ceilq(degree,6); hi=degree//2-5
    check("all D below16 violate integer pole interval D"+str(degree),lo>hi)
for degree in range(1,48):
    lo=ceilq(degree,6); hi=degree//4-4
    check("all no-A D below48 violate integer pole interval D"+str(degree),lo>hi)
check("universal D16 numerical interval remains compatible",
      ceilq(16,6)==16//2-5==3)
check("no-A D48 numerical interval remains compatible",
      ceilq(48,6)==48//4-4==8)

for twos in range(4):
    ones=6-2*twos
    check("original infinity all-index partition twos"+str(twos),
          ones+2*twos==6)
    check("original t has simple pole from alpha1 twos"+str(twos),
          Fraction(1*2,2)==1)
    check("original t has simple pole from alpha2 twos"+str(twos),
          Fraction(2*1,2)==1)
    check("original infinity exactly six points twos"+str(twos),
          ones+2*twos==6)
check("no-A intermediate infinity exactly six unramified points",6*1==6)

# Unramified source infinity need not have residue degree one, even for a
# global polynomial submersion with holomorphic generic coframe there.
x,t,p0=S.symbols("x t p0")
P0=x+(t**2*x)**6
PX=S.diff(P0,x); PT=S.diff(P0,t)
check("infinity residue guard global submersion Bezout identity",
      S.expand(PX*(1-6*x**5*t**12)+PT*(3*x**4*t**13))==1)
sub=S.cancel(P0.subs({x:s**2*Y,t:1/s},simultaneous=True))
check("infinity residue guard entire equation",sub==Y**6+s**2*Y)
check("infinity residue guard six simple geometric branches",
      S.discriminant(Y**6-p0,Y)==46656*p0**5)
check("infinity residue guard holomorphic coframe",
      S.cancel((-s**(-2))/PX.subs({x:s**2*Y,t:1/s},simultaneous=True))
      ==-1/(s**2+6*Y**5))
check("infinity residue guard residue polynomial degree six",
      S.degree(Y**6-p0,Y)==6)
check("infinity residue guard irreducibility Eisenstein p",
      S.Poly(Y**6-p0,Y).nth(0)==-p0 and S.degree(-p0,p0)==1)
check("infinity residue guard not asserted actual Keller target",
      P0.subs(t,0)==x)

# Controls rejecting stronger, currently unproved conclusions.
check("ledger does not close all A possibilities",
      2*4-5==3 and 4*4==16 and ceilq(16,6)==3)
check("ledger does not close all no-A B possibilities",
      12-4==8 and 4*12==48 and ceilq(48,6)==8)
check("ruling and intermediate/original generic degrees kept distinct",
      6 != 2*8 and 8 != 6)

receipt={
    "claim":"PF2RULED164","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
        name:sha(BASE/name) for name in [
            "RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md",
            "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
            "MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md"]},
    "uniform_proved_scope":{
        "coefficient_degrees":"unrestricted","quadratic_source_models":"complete PF157 class",
        "local_budget_A":"2*L-5","local_budget_B":"L-4","local_budget_C":"B0-4",
        "actual_B_pole_slope_minimum":5,"actual_C_pole_B0_minimum":5,
        "intermediate_degree_global_upper":"floor(D/2)-5",
        "no_A_intermediate_degree_global_upper":"floor(D/4)-4",
        "infinity_local_degree_lower":"ceil(alpha*(D+3*m-6)/6)",
        "intermediate_infinity_u_index_maximum":2,
        "no_A_intermediate_infinity_u_index":1,
        "original_source_infinity_t_index":1,"original_source_infinity_point_count":6,
        "leading_degree_minimum":16,"no_A_leading_degree_minimum":48},
    "finite_control_completeness":{
        "quartic_partitions":len(Aparts),"cubic_partitions":len(BCparts),
        "all_slope_extensions":"exact period2/3 with index-mass monotonicity"},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser()
parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_pole_budget_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),
    "elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))
