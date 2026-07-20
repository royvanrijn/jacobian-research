#!/usr/bin/env python3
"""Exact chart-first pilot scan with obstruction certificates, degrees 5--20.

This exhausts only the explicitly generated sparse B families, not all Keller
maps of these degrees. Results are written as JSON for feedback into generators.
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.charts import cubic_chart, u, v, x, y
from jcsearch.laurent import linear_certificate, integer_matrix, rank_mod_prime

PRIMES = (1000003, 1000033, 1000037)

def polynomial_ansatz(degree):
    # A full total-degree ansatz is still linear for fixed B.
    coeffs, terms = [], []
    for total in range(degree+1):
        for iy in range(total+1):
            ix = total-iy
            symbol = sp.Symbol(f"a_{ix}_{iy}")
            coeffs.append(symbol); terms.append(symbol*x**ix*y**iy)
    return sp.Add(*terms), coeffs

def candidates(degree, collision_normalized=False):
    # Coordinate controls and sparse top Newton edges.
    if collision_normalized:
        # Every B below satisfies B(0,0)=B(1,0)=0 and dB(0,0)=(0,1).
        return {
            "coordinate_y": y,
            "one_vertex_edge": y+x**(degree-1)*(x-1),
            "mixed_edge": y+x**(degree-2)*y*(x-1),
            "two_vertex_edge": y+(x**(degree-1)+y**(degree-1))*(x-1),
            "near_edge": y+x**(degree-2)*(x+y)*(x-1),
            "split_edge": y+x**2*y**(degree-3)*(x-1),
        }
    return {
        "coordinate_x": x,
        "triangular": x+y**degree,
        "pure_power": x**degree,
        "mixed_edge": x**(degree-1)*y,
        "two_vertex_edge": x**degree+y**degree,
        "near_edge": x**degree+x*y,
    }

def solve_for_A(B, degree, collision_normalized=False):
    A, unknowns = polynomial_ansatz(degree)
    residual = sp.Poly(sp.expand(sp.diff(A,x)*sp.diff(B,y)-sp.diff(A,y)*sp.diff(B,x)-1), x, y)
    equations=list(residual.coeffs())
    if collision_normalized:
        # F(0)=F(1,0)=0 and DF(0)=I. B's half is checked separately.
        b_checks=[B.subs({x:0,y:0}),B.subs({x:1,y:0}),
                  sp.diff(B,x).subs({x:0,y:0}),sp.diff(B,y).subs({x:0,y:0})-1]
        assert all(sp.expand(e)==0 for e in b_checks), (B,b_checks)
        equations += [A.subs({x:0,y:0}),A.subs({x:1,y:0}),
                      sp.diff(A,x).subs({x:0,y:0})-1,
                      sp.diff(A,y).subs({x:0,y:0})]
    cert, matrix, rhs = linear_certificate(equations, unknowns)
    im, irhs, scale = integer_matrix(matrix, rhs)
    modular = {str(p): {"rank": rank_mod_prime(im,p),
                        "augmented_rank": rank_mod_prime(im.row_join(irhs),p)} for p in PRIMES}
    solution = None
    if cert.consistent:
        vector = sp.linsolve((matrix,rhs), unknowns)
        solution = sp.sstr(next(iter(vector)))
    return cert, modular, solution

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--min-degree",type=int,default=5)
    parser.add_argument("--max-degree",type=int,default=20)
    parser.add_argument("--output",default="results/scan_2d_5_20.json")
    parser.add_argument("--collision-normalized",action="store_true")
    args=parser.parse_args()
    chart=cubic_chart()
    records=[]
    for degree in range(args.min_degree,args.max_degree+1):
        family=candidates(degree,args.collision_normalized)
        for name,B in family.items():
            cert, modular, solution=solve_for_A(B,degree,args.collision_normalized)
            if cert.consistent:
                obstruction="birational_coordinate_family"
            else:
                obstruction="jacobian_linear_inconsistency"
            B_chart=sp.cancel(B.subs({x:2/v,y:u-v/2}, simultaneous=True))
            records.append({"degree":degree,"family":name,"B":sp.sstr(B),
                "B_in_chart":sp.sstr(B_chart),"chart_jacobian":sp.sstr(chart.jacobian),
                "status":"survivor" if cert.consistent else "rejected",
                "obstruction":obstruction,"certificate":cert.__dict__,
                "modular_ranks":modular,"solution":solution})
        print(f"degree {degree}: checked {len(family)} sparse edge families")
    path=Path(args.output); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps({"scope":"enumerated sparse families only","records":records},indent=2))
    survivors=[r for r in records if r["status"]=="survivor"]
    print(f"wrote {len(records)} exact certificates to {path}")
    if args.collision_normalized:
        print(f"survivors after collision normalization: {len(survivors)}")
    else:
        print(f"survivors: {len(survivors)}; all are coordinate/triangular and birational")

if __name__=="__main__": main()
