#!/usr/bin/env python3
"""Exact controls for the constant directional-derivative Keller exclusion."""
from pathlib import Path
import hashlib
import json
import subprocess
import sympy as sp

HERE=Path(__file__).resolve().parent
SOURCE="https://www.math.purdue.edu/~ttm/jacobian"

def main():
    x,y,z=sp.symbols("x y z")
    source=(x,y,z)
    u=1+x*y
    foundation=sp.Matrix([
        u**3*z+y**2*u*(4+3*x*y),
        y+3*x*u**2*z+3*x*y**2*(4+3*x*y),
        2*x-3*x*x*y-x**3*z,
    ])
    assert sp.expand(foundation.jacobian(source).det())==-2
    a={x:0,y:0,z:-sp.Rational(1,4)}
    b={x:1,y:-sp.Rational(3,2),z:sp.Rational(13,2)}
    assert foundation.subs(a)==foundation.subs(b)
    # Tangent-to-identity target normalization; affine equivalence preserves
    # the existence of a constant nonzero Jacobian pairing.
    G=sp.Matrix([foundation[2]/2,foundation[1],foundation[0]])
    J=G.jacobian(source)
    assert J.subs({x:0,y:0,z:0})==sp.eye(3)
    v=sp.Matrix(sp.symbols("v0:3"));w=sp.Matrix(sp.symbols("w0:3"))
    pairing=sp.expand((v.T*J*w)[0])
    equations=[c for monomial,c in sp.Poly(pairing,x,y,z).terms()
               if monomial!=(0,0,0)]
    equations.append((v.T*w)[0]-1)
    variables=(*v,*w)
    basis=sp.groebner(equations,*variables,order="grevlex")
    assert len(basis.polys)==1 and basis.polys[0].as_expr()==1
    script=("ring R=0,("+",".join(map(str,variables))+"),dp;\n"
            +"ideal I="+",".join(str(e).replace("**","^") for e in equations)+";\n"
            +"ideal B=std(I);\n"
            +'print("UNIT|"+string(reduce(1,B)==0));\nexit;\n')
    result=subprocess.run(["Singular","-q"],input=script,text=True,capture_output=True,check=True)
    assert "UNIT|1" in result.stdout
    (HERE/"direction-control.sing").write_text(script)
    (HERE/"direction-control.log").write_text(result.stdout+result.stderr)

    # An explicit nonlinear source shear has derivative one along z, and
    # substitution of its inverse makes the displayed coordinate a parameter.
    P=x+x*y+z*z
    Q=y+x*z+z*z
    R=z+x*x+x*y+y*y
    H=sp.Matrix([P,Q,R])
    assert sp.diff(H[2],z)==1
    t=sp.symbols("t")
    inverse_shift=t-x*x-x*y-y*y
    K=H.subs(z,inverse_shift)
    assert sp.expand(K[2])==t
    assert max(sp.Poly(sp.expand(entry),x,y,t).total_degree() for entry in K)<=4
    assert sp.expand(K.jacobian((x,y,t)).det()-
                     H.jacobian(source).det().subs(z,inverse_shift))==0

    # The two additional universal gates and the exceptional eight-term row.
    h=sp.Function("h")
    A,B,C,D=sp.symbols("A B C D")
    assert sp.expand(sp.det(sp.Matrix([[A,0,0],[B,C,D],[D,A,B]]))-A*(B*C-A*D))==0
    normal=1+sp.symbols("a")*x
    R0=sp.Matrix([[A,B,C],[0,normal,0],[D,A,B]])
    assert sp.expand(R0.det()-normal*(A*B-C*D))==0

    stable=sp.Matrix([*G,t])
    stable_J=stable.jacobian((x,y,z,t))
    assert sp.expand(stable_J.det())==1 and stable_J[3,3]==1
    assert stable.subs(a)==stable.subs(b)
    degrees=[sp.Poly(sp.expand(entry),x,y,z,t).total_degree() for entry in stable]
    assert max(degrees)==7
    receipt=dict(
        status="pass",source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        primary_plane_source=SOURCE,plane_degree_import=100,
        scalar_pairing_exclusion_max_three_variable_degree=10,
        normalized_foundation_degrees=degrees[:3],
        foundation_pairing_incidence="unit over Q in SymPy and Singular",
        foundation_incidence_equations=len(equations),
        nonlinear_shear_control="inverse substitution and Jacobian chain rule exact",
        four_variable_control="degree-seven determinant-one collision with fourth coordinate fixed",
        four_variable_degrees=degrees,
        sharp_dimension_scope="Minimum dimension four for degree-at-most-ten Keller collisions admitting a constant nonzero Jacobian pairing; examples exist already at degree seven.",
        proof_boundary="The all-coefficient affine-direction reduction and degree bound are written proofs using Moh's degree-100 plane theorem. Symbolic controls do not reprove Moh.",
    )
    (HERE/"coordinate-gates-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt),flush=True)
if __name__=="__main__":
    main()

