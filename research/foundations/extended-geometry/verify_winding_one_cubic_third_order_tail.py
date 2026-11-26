#!/usr/bin/env python3
"""Replay the exact third-order tail exclusion with Python and Singular."""

from fractions import Fraction
from math import comb
import shutil
import subprocess


VARIABLES = "a3,a4,a5,b3,b4,b5,b6,b7,t,x"
GUARD_FACTORS = ("7*d^2+1", "3*d^2+1", "143*d^4+78*d^2+3", "13*d^4+10*d^2+1")
EXCEPTIONAL_FACTORS = ("d", "d-1", "d+1", *GUARD_FACTORS)


def interval_moment(n, ell):
    terms = []
    for j in range(n + 1):
        if (j + ell) % 2:
            continue
        weight = sum(
            Fraction((-1) ** (n - h) * comb(n, j) * comb(n, h), 2 * h + j + ell + 1)
            for h in range(n + 1)
        )
        if weight:
            term = f"({weight})"
            if j != n:
                term += f"*d^{n - j}"
            terms.append(term)
    return "+".join(terms) or "0"


def run(singular, name, source):
    result = subprocess.run(
        [singular, "-q"],
        input=source + '\nprint("THIRD-ORDER TAIL CHECK PASSED");\nquit;\n',
        text=True,
        capture_output=True,
        check=True,
    )
    diagnostics = result.stdout + result.stderr
    if (
        "?" in diagnostics
        or "FAIL" in diagnostics
        or "THIRD-ORDER TAIL CHECK PASSED" not in result.stdout.splitlines()
    ):
        raise RuntimeError(f"{name} failed:\n{result.stdout}\n{result.stderr}")
    for line in result.stdout.splitlines():
        if line.startswith("certificate "):
            print(line)
    print(f"passed {name}")


LIFT_CODE = """
proc checkedLift(ideal H, ideal target)
{
    matrix liftMatrix=lift(H,target);
    matrix discrepancy=matrix(H)*liftMatrix-matrix(target);
    for(int i=1;i<=nrows(discrepancy);i++){
        for(int j=1;j<=ncols(discrepancy);j++){
            if(discrepancy[i,j]!=0){ERROR("FAIL lift identity");}
        }
    }
    return(liftMatrix);
}
"""


def moment_code(minpoly=None):
    source = f"ring R=(0,d),({VARIABLES}),dp;\n"
    if minpoly is not None:
        source += f"minpoly={minpoly};\n"
    source += "matrix interval[10][4];\n"
    for n in range(10):
        for ell in range(4):
            source += f"interval[{n + 1},{ell + 1}]={interval_moment(n, ell)};\n"
    return source + LIFT_CODE + """
proc choose(int n, int k)
{
    number q=1;
    for(int j=1;j<=k;j++){q=q*(n-j+1)/j;}
    return(q);
}

proc integrate(poly P, int n)
{
    matrix C=coef(P,t);
    poly value=0;
    intvec ex;
    for(int i=1;i<=ncols(C);i++){
        if(C[1,i]!=0){
            ex=leadexp(C[1,i]);
            value=value+C[2,i]*interval[n+1,ex[9]+1];
        }
    }
    return(value);
}

proc sphereMoment(poly tail, int m)
{
    poly value=0;
    matrix C;
    int j;
    for(int l=1;3*l<=m;l++){
        C=coeffs(tail^l,x);
        for(j=3*l;j<=m && j<nrows(C);j++){
            value=value+choose(m,l)*choose(m-1,j-1)*integrate(C[j+1,1],m-l);
        }
    }
    return(value);
}
"""


GENERIC_CODE = """
if(interval[3,1]!=8/105*(7*d^2+1)
   || interval[4,1]!=-16/105*d*(3*d^2+1)
   || interval[5,1]!=128/45045*(143*d^4+78*d^2+3)
   || interval[6,1]!=-256/9009*d*(13*d^4+10*d^2+1)){
    ERROR("FAIL interval guard factors");
}

poly A3=8*a4/9-22*a5/27;
poly C3=-A3*interval[3,2]/interval[3,1];
poly C4=-3*C3-(3*A3+a4)*interval[4,2]/interval[4,1];
poly C5=-6*C3-4*C4-(6*A3+4*a4+a5)*interval[5,2]/interval[5,1];
poly C7=-(2/9*(C3-d*A3/3)-14/81*(C4-d*a4/3)
              +35/243*(C5-d*a5/3)-91/729*b6)/(728/6561);
poly tail=x^3*(A3*t+C3)+x^4*(a4*t+C4)+x^5*(a5*t+C5)+x^6*b6+x^7*C7;
poly M6=sphereMoment(tail,6);
poly C6=-subst(M6,b6,0)/(6*interval[6,1]);
if(subst(M6,b6,C6)!=0){ERROR("FAIL sixth moment substitution");}
tail=subst(tail,b6,C6);

poly M7=sphereMoment(tail,7);
poly M8=sphereMoment(tail,8);
poly M9=sphereMoment(tail,9);
poly M10=sphereMoment(tail,10);

proc denominators(matrix C)
{
    list factors;
    poly q;
    int j;
    for(int i=1;i<=nrows(C);i++){
        for(j=1;j<=ncols(C);j++){
            q=C[i,j];
            while(q!=0){
                factors[size(factors)+1]=denominator(leadcoef(q));
                q=q-lead(q);
            }
        }
    }
    return(factors);
}

ideal target=a4,a5;
list liftDenominators;
liftDenominators[1]=denominators(checkedLift(ideal(M7,M8,M9),target));
liftDenominators[2]=denominators(checkedLift(ideal(M7,M8,M10),target));
liftDenominators[3]=denominators(checkedLift(ideal(M7,M9,M10),target));
list transportDenominators=denominators(matrix(ideal(tail,M7,M8,M9,M10)));

ring Q=0,(d),dp;
list allDenominators=imap(R,liftDenominators);
list transport=imap(R,transportDenominators);
proc squarefreeLcm(list factors)
{
    poly value=1;
    poly den;
    for(int i=1;i<=size(factors);i++){
        den=factors[i];
        value=value*den/gcd(value,den);
    }
    return(value/gcd(value,diff(value,d)));
}
poly guard=d*(7*d^2+1)*(3*d^2+1)*(143*d^4+78*d^2+3)*(13*d^4+10*d^2+1);
poly transportLcm=squarefreeLcm(transport);
if(deg(gcd(transportLcm,guard))!=deg(transportLcm)){
    ERROR("FAIL transport includes a denominator outside the guard");
}
poly D1=squarefreeLcm(allDenominators[1]);
poly D2=squarefreeLcm(allDenominators[2]);
poly D3=squarefreeLcm(allDenominators[3]);
poly D=gcd(gcd(D1,D2),D3);
poly expected=d*(d-1)*(d+1)*(7*d^2+1)/7;
if(D/leadcoef(D)!=expected){ERROR("FAIL lift denominator gcd");}
print("certificate squarefree lift denominator degrees "
    +string(deg(D1))+" "+string(deg(D2))+" "+string(deg(D3)));
print("certificate monic common denominator divisor "+string(D/leadcoef(D)));
"""


EXCEPTIONAL_CODE = """
poly tail=x^3*(a3*t+b3)+x^4*(a4*t+b4)+x^5*(a5*t+b5)+x^6*b6+x^7*b7;
poly N2=a3-8*a4/9+22*a5/27;
poly N3=2/9*(b3-d*a3/3)-14/81*(b4-d*a4/3)+35/243*(b5-d*a5/3)
         -91/729*b6+728/6561*b7;
ideal H=N2,N3;
for(int m=3;m<=9;m++){H=H+ideal(sphereMoment(tail,m));}
ideal target=a3^2,a4^2,a5^2,b3^2,b4^2,b5^2,b6^2,b7^2;
matrix identity=checkedLift(H,target);
"""


def main():
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Install Singular 4 and place it on PATH.")

    irreducibility = "ring Q=0,(d),dp;\nideal factors;\n"
    for polynomial in GUARD_FACTORS:
        irreducibility += (
            f"factors=factorize({polynomial},1);\n"
            f'if(size(factors)!=1 || deg(factors[1])!=deg({polynomial}))'
            '{ERROR("FAIL guard factor irreducibility");}\n'
        )
    run(singular, "guard factor irreducibility", irreducibility)
    run(singular, "generic lift certificates", moment_code() + GENERIC_CODE)
    for polynomial in EXCEPTIONAL_FACTORS:
        run(singular, f"exceptional lift modulo {polynomial}", moment_code(polynomial) + EXCEPTIONAL_CODE)


if __name__ == "__main__":
    main()
