#!/usr/bin/env python3
"""Verify weighted specialization and tangent equations with Singular."""

from fractions import Fraction
import shutil
import subprocess

from verify_winding_one_cubic_third_order_tail import interval_moment


PRIME = 101
SUCCESS = "WEIGHTED SPECIALIZATION CHECK PASSED"


def binomial(alpha, k):
    value = Fraction(1)
    for j in range(k):
        value *= (alpha - j) / (j + 1)
    return value


def run(singular, label, source):
    result = subprocess.run(
        [singular, "-q"],
        input=source + f'\nprint("{SUCCESS}");\nquit;\n',
        text=True,
        capture_output=True,
        check=True,
    )
    diagnostics = result.stdout + result.stderr
    if "?" in diagnostics or "FAIL" in diagnostics or SUCCESS not in result.stdout.splitlines():
        raise RuntimeError(f"{label} failed:\n{diagnostics}")
    print(f"passed {label}")


def show_target_inclusion(target, generators="H"):
    return f"""
ideal gb=std({generators});
if(size(reduce(ideal({",".join(target)}),gb))!=0){{ERROR("FAIL ideal membership");}}
"""


def sector(d, mmax, *, infinity=False, rational=False):
    a_max, b_max = (4, 6) if infinity else (5, 7)
    a_names = [f"a{j}" for j in range(1, a_max + 1)]
    b_names = [f"b{j}" for j in range(2, b_max + 1)]
    variables = a_names + b_names
    ring_names = ",".join([*variables, "t", "x"])
    angular_weights = ",".join(["0"] * (len(variables) + 1) + ["1"])
    coefficient_weights = ",".join(["1"] * len(variables) + ["0", "0"])
    char = 0 if rational else PRIME
    source = f"ring R={char},({ring_names}),dp;\nnumber d={d};\n"
    source += f"matrix interval[{mmax + 1}][{mmax + 1}];\n"
    source += "".join(
        f"interval[{n + 1},{ell + 1}]={interval_moment(n, ell)};\n"
        for n in range(mmax + 1)
        for ell in range(mmax + 1 - n)
    )
    source += f"""
proc choose(int n, int k)
{{
    number value=1;
    for(int j=1;j<=k;j++){{value=value*(n-j+1)/j;}}
    return(value);
}}

proc integrate(poly p, int n)
{{
    matrix C=coef(p,t);
    poly value=0;
    intvec exponent;
    for(int i=1;i<=ncols(C);i++){{
        if(C[1,i]!=0){{
            exponent=leadexp(C[1,i]);
            value=value+C[2,i]*interval[n+1,exponent[{len(variables) + 1}]+1];
        }}
    }}
    return(value);
}}

intvec coefficientWeights={coefficient_weights};
proc moment(poly tail, int m, int jetDegree)
{{
    intvec angularWeights={angular_weights};
    tail=jet(tail,m,angularWeights);
    poly power=1;
    poly value=0;
    matrix C;
    int j;
    for(int l=1;l<=m;l++){{
        power=jet(power*tail,m,angularWeights);
        if(jetDegree>0){{power=jet(power,jetDegree,coefficientWeights);}}
        C=coeffs(power,x);
"""
    if infinity:
        source += "        if(m<nrows(C)){value=value+choose(m,l)*integrate(C[m+1,1],m-l);}\n"
    else:
        source += f"""
        for(j=l;j<=m && j<nrows(C);j++){{
            value=value+choose(m,l)*choose(m-1,j-1)*integrate(C[j+1,1],m-l);
        }}
"""
    source += """
    }
    return(value);
}
"""
    terms = []
    for j in range(1, b_max + 1):
        row = []
        if j <= a_max:
            row.append(f"a{j}*t")
        if j >= 2:
            row.append(f"b{j}")
        terms.append(f"x^{j}*({'+'.join(row)})")
    source += f"poly tail={'+'.join(terms)};\n"
    if infinity:
        source += """
if(110/243==0 || 728/6561==0){ERROR("FAIL inverse equations at infinity");}
"""
    else:
        first = [
            f"({binomial(Fraction(-2, 3), j - 1)})*a{j}" for j in range(1, a_max + 1)
        ]
        second = []
        for j in range(1, b_max + 1):
            b_term = f"b{j}" if j >= 2 else "0"
            a_term = f"-d*a{j}/3" if j <= a_max else ""
            second.append(f"({binomial(Fraction(-1, 3), j - 1)})*({b_term}{a_term})")
        source += f"poly N2={'+'.join(first)};\npoly N3={'+'.join(second)};\n"
    return source, variables


def reduced_fiber(d, mmax, infinity):
    source, variables = sector(d, mmax, infinity=infinity)
    source += "ideal H;\n" if infinity else "ideal H=N2,N3;\n"
    source += f"for(int m=1;m<={mmax};m++){{H=H+ideal(moment(tail,m,0));}}\n"
    return source + show_target_inclusion(
        [f"{variable}^20" if infinity else variable for variable in variables]
    )


def d_one_residual(mmax, *, rational=False):
    source, _ = sector(1, mmax, rational=rational)
    jet_degree = 3 if rational else 0
    source += """
poly A5=-(a1-2*a2/3+5*a3/9-40*a4/81)*243/110;
if(subst(N2,a5,A5)!=0){ERROR("FAIL second-deficit elimination");}
tail=subst(tail,a5,A5);
poly equation;
number pivot;
"""
    for j in range(2, 8):
        unit_test = (
            "bigint(numerator(pivot))%101==0 || bigint(denominator(pivot))%101==0"
            if rational else "pivot==0"
        )
        source += f"""
pivot=leadcoef({j}*interval[{j},1]);
if({unit_test}){{ERROR("FAIL moment pivot");}}
equation=moment(tail,{j},{jet_degree});
poly C{j}=-subst(equation,b{j},0)/pivot;
if(equation-pivot*b{j}+pivot*C{j}!=0 || jet(C{j},0)!=0){{
    ERROR("FAIL triangular elimination");
}}
tail=subst(tail,b{j},C{j});
N3=subst(N3,b{j},C{j});
"""
        if rational:
            source += "tail=jet(tail,3,coefficientWeights);\n"
    source += f"""
N3=subst(N3,a5,A5);
ideal H=N3;
for(int m=8;m<={mmax};m++){{H=H+ideal(moment(tail,m,{jet_degree}));}}
"""
    return source


TANGENT_CODE = """
ring LocalCoefficients=0,(a1,a2,a3,a4),dp;
ideal rows=imap(R,H);
ideal orderTwo=a1^2,a1*a2,a1*a3,a1*a4,a2^2,a2*a3,a2*a4,a3^2,a3*a4,a4^2;
ideal initial;
int n=size(rows);

proc coeffAt(poly p, poly rowMon)
{
    poly q=p/rowMon;
    return(leadcoef(subst(subst(subst(subst(q,a1,0),a2,0),a3,0),a4,0)));
}

int i,j,pivot;
number c;
poly r;
for(i=1;i<=n;i++){
    if(jet(rows[i],1)!=0){ERROR("FAIL exact residual linear term");}
    r=rows[i];
    while(r!=0){
        if(bigint(denominator(leadcoef(r)))%101==0){
            ERROR("FAIL input coefficient denominator");
        }
        r=r-lead(r);
    }
}
for(j=1;j<=size(orderTwo);j++){
    pivot=0;
    for(i=1;i<=n;i++){
        c=coeffAt(rows[i],orderTwo[j]);
        if(bigint(numerator(c))%101!=0){pivot=i;break;}
    }
    if(pivot!=0){
        r=rows[pivot]/coeffAt(rows[pivot],orderTwo[j]);
        initial=initial+ideal(jet(r,2));
        rows[pivot]=0;
        for(i=1;i<=n;i++){rows[i]=rows[i]-coeffAt(rows[i],orderTwo[j])*r;}
    }
}
if(size(initial)!=4){ERROR("FAIL quadratic pivot count");}
for(i=1;i<=n;i++){
    if(jet(rows[i],2)!=0){ERROR("FAIL unresolved exact quadratic");}
    initial=initial+ideal(rows[i]);
}
ring F=101,(a1,a2,a3,a4),dp;
ideal H=imap(LocalCoefficients,initial);
"""


def main():
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Install Singular 4 and place it on PATH.")
    for d in (0, 2):
        for infinity in (False, True):
            chart = "infinity" if infinity else "affine"
            run(
                singular,
                f"d={d}, {chart}",
                reduced_fiber(d, 12, infinity),
            )
    run(singular, "d=1, infinity", reduced_fiber(1, 12, True))
    target = [f"a{j}^10" for j in range(1, 5)]
    run(singular, "d=1, affine radical", d_one_residual(12) + show_target_inclusion(target))
    run(
        singular,
        "d=1, exact initial forms",
        d_one_residual(16, rational=True) + TANGENT_CODE + show_target_inclusion(target),
    )


if __name__ == "__main__":
    main()
