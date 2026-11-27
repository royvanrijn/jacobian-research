#!/usr/bin/env python3
"""Certify a fixed start and survivor gap for all integer axis offsets."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from singular_runner import run_singular


SINGULAR_FORMS = r"""
option(noredefine);
ring rational=0,(z,t,a,c),dp;

// A cleared binary moment. c is M-a, and z is the higher/lower ratio.
proc moment(int m, int q, int gap, int tied)
{
    poly f=0;
    poly wt;
    int k,j;
    for(k=0;k<=m;k++){
        wt=1;
        for(j=1;j<=k;j++){ wt=wt*number(m-j+1)/j; }
        for(j=1;j<=gap*k;j++){ wt=wt*(m*a+j); }
        for(j=gap*k;j<gap*m;j++){ wt=wt*(m*(c-1/number(q))-j); }
        if(tied==1 && k==m){ wt=wt*(1+t^m); }
        f=f+wt*z^k;
    }
    if(tied==1){ f=subst(f,c,gap); }
    return(f);
}

// Every factor of allowed is nonzero on the indicated integer domain.
proc strip(poly f, poly allowed)
{
    if(f==0){return(0);}
    poly divisor=gcd(f,allowed);
    poly reducedF;
    while(deg(divisor)>0){
        reducedF=f/divisor;
        if(reducedF*divisor!=f){ERROR("nonexact resultant division");}
        f=reducedF;
        divisor=gcd(f,allowed);
    }
    return(cleardenom(f));
}

proc marked_resultants(int n, int gap, int tied)
{
    ideal result;
    int q,m,k;
    poly f,g,eliminant,allowed;
    ideal rule;
    for(q=n;q<=n+2;q++){
        f=0;
        for(m=n;m<=n+2;m++){
            if(m!=q){
                if(f==0){f=moment(m,q,gap,tied);}
                else{g=moment(m,q,gap,tied);}
            }
        }
        if(tied==1){
            rule=std(ideal(t^q+1));
            f=reduce(f,rule);
            g=reduce(g,rule);
        }
        eliminant=resultant(f,g,z);
        allowed=q*a+q*c-1;
        if(tied==1){
            eliminant=resultant(reduce(eliminant,rule),t^q+1,t);
            allowed=subst(allowed,c,gap);
            for(k=1;k<=q*gap;k++){ allowed=allowed*(a+k); }
        } else {
            for(k=1;k<=gap;k++){allowed=allowed*(a+k)*(q*c-1-q*(k-1));}
            allowed=allowed*(q*(c-a-gap)-1);
        }
        result[q-n+1]=strip(eliminant,allowed);
    }
    return(result);
}

proc print_coefficients(ideal source)
{
    int j;
    poly term;
    intvec indices;
    for(j=1;j<=size(source);j++){
        term=source[j];
        while(term!=0){
            indices=leadexp(term);
            if(indices[1]!=0 || indices[2]!=0){ERROR("incomplete elimination");}
            print("COEFFICIENT "+string(j)+" "+string(indices[3])+" "+
                string(indices[4])+" "+string(leadcoef(term)));
            term=term-lead(term);
        }
    }
}
"""


def verify_modular(start, gap, mode, depth):
    tied = int(mode == "tied")
    source = SINGULAR_FORMS + f"\nideal marks=marked_resultants({start},{gap},{tied});\n"
    if depth:
        source += "print_coefficients(marks);\n"
        output = run_singular(source, timeout=300)
        coefficients = [[], [], []]
        for line in output.splitlines():
            if line.startswith("COEFFICIENT "):
                order, a_power, c_power, coefficient = map(int, line.split()[1:])
                coefficients[order - 1].append((a_power, c_power, coefficient))
        if not any(coefficients):
            raise ArithmeticError("No nonzero resultant returned.")
        counts = surviving_residues(coefficients, depth)
        if counts[-1]:
            raise ArithmeticError(f"Dyadic residues remain through level {depth}.")
        print(f"Dyadic surviving residue counts {counts}")
    else:
        source += r"""
ring finite=65537,(z,t,a,c),dp;
ideal equations=imap(rational,marks);
ideal basis=std(equations);
if(reduce(1,basis)!=0){ERROR("marked resultants do not generate the unit ideal");}
print("UNIT_IDEAL_CERTIFIED");
"""
        run_singular(source, timeout=300, required_markers=("UNIT_IDEAL_CERTIFIED",))
    print(f"Certified start {start}, gap {gap}, {mode} maximum.", flush=True)


def specialize_at_a(coefficients, offset, modulus):
    max_a_power = max((i for poly in coefficients for i, _, _ in poly), default=0)
    powers = [pow(offset, i, modulus) for i in range(max_a_power + 1)]
    specialized = []
    for poly in coefficients:
        row = [0] * (1 + max((j for _, j, _ in poly), default=0))
        for i, j, coefficient in poly:
            row[j] += coefficient * powers[i]
        specialized.append([coefficient % modulus for coefficient in reversed(row)])
    return specialized


def evaluate(polynomial, value, modulus):
    result = 0
    for coefficient in polynomial:
        result = (result * value + coefficient) % modulus
    return result


def surviving_residues(coefficients, depth):
    points = {(0, 0)}
    counts = []
    for level in range(1, depth + 1):
        modulus = 2**level
        shift = modulus // 2
        reduced = [
            [(i, j, c % modulus) for i, j, c in polynomial if c % modulus]
            for polynomial in coefficients
        ]
        lifts = {}
        for a, c in points:
            for offset in (a, a + shift):
                lifts.setdefault(offset, []).extend((c, c + shift))
        points = set()
        for a, candidates in lifts.items():
            rows = specialize_at_a(reduced, a, modulus)
            points.update(
                (a, c)
                for c in candidates
                if all(evaluate(row, c, modulus) == 0 for row in rows)
            )
        counts.append(len(points))
        if not points:
            break
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", type=int)
    parser.add_argument("gap", type=int, nargs="+")
    parser.add_argument("--mode", choices=("unique", "tied", "both"), default="both")
    parser.add_argument("--dyadic-depth", type=int, default=0)
    args = parser.parse_args()
    if args.start < 3 or min(args.gap) < 1 or args.dyadic_depth < 0:
        parser.error("Use a start of at least three, positive gaps, and nonnegative depth.")
    if args.dyadic_depth and args.mode != "unique":
        parser.error("The two-coordinate dyadic replay is for --mode unique.")
    modes = ("unique", "tied") if args.mode == "both" else (args.mode,)
    for gap in args.gap:
        for mode in modes:
            verify_modular(args.start, gap, mode, args.dyadic_depth)


if __name__ == "__main__":
    main()
