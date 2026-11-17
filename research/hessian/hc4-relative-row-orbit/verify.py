"""Exact row-orbit, gauge and singleton-fibre control calculations.

Orbit completeness, the singleton-zero reformulation and normal-surface
flatness are written proofs in PROOF.md, not consequences of a search.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
checks=[]


def check(name,condition):
    assert condition,name
    checks.append(name)


def zero(name,expression):
    values=list(expression) if isinstance(expression,S.MatrixBase) else [expression]
    check(name,all(S.cancel(value)==0 for value in values))


def degree(matrix,x,y):
    return max(int(S.Poly(value,x,y).total_degree()) for value in matrix)


def main():
    start=time.monotonic()
    x,y,v,h=S.symbols('x y v h')
    P=18*v**5+81*v**4+120*v**3+60*v*v-1
    Q=(v+1)*(v+2)
    aa,bb=-P,Q
    T=S.cancel((2-aa*bb)/v)
    c12=S.rem(S.invert(aa,2*bb)*T,2*bb,v)
    c21=S.cancel((T-aa*c12)/(2*bb))
    c11=S.cancel((1-v*bb*c21)/aa)
    c22=S.cancel((2-v*aa*c12)/bb)
    C0=S.Matrix([[c11,x*x*c12],[y*y*c21,c22]]).subs(v,x*y)
    gamma=S.Matrix([[-y*P.subs(v,x*y),x*Q.subs(v,x*y)]])
    ell=S.Matrix([[y,2*x]])
    zero('explicit polynomial normalizer row',gamma*C0-ell)
    zero('explicit normalizer determinant',C0.det()-1)
    check('normalizer degree ten',degree(C0,x,y)==10)
    check('inverse normalizer degree ten',degree(C0.adjugate(),x,y)==10)
    zero('local row Jacobian orientation',gamma.jacobian((x,y)).det().subs({x:0,y:0})+2)
    # Universal oriented-presentation correction, without an ideal decision.
    a11,a12,a21,a22,mu,nu=S.symbols('a11 a12 a21 a22 mu nu')
    A=S.Matrix([[a11,a12],[a21,a22]])
    g=ell*A
    kell=S.Matrix([-2*x,y])
    corrected=A+kell*S.Matrix([[nu,-mu]])
    zero('universal correction keeps the row',ell*corrected-g)
    zero('universal adjugate identity',A.adjugate()*kell-S.Matrix([-g[1],g[0]]))
    zero('universal determinant correction',corrected.det()-A.det()+nu*g[1]+mu*g[0])
    # A nontrivial literal ideal certificate and oriented lift.
    f=y+2*x**3
    gg=2*x+f*f
    target=S.Matrix([[f,gg]])
    naive=S.Matrix([[1,y+4*x**3],[x*x,1+2*x**5]])
    zero('example row decomposition',ell*naive-target)
    zero('example determinant defect is in target ideal',naive.det()-1+x*x*f)
    oriented=naive+kell*S.Matrix([[0,x*x]])
    zero('example determinant-corrected lift',oriented.det()-1)
    zero('example complete row lift',gamma*C0*oriented-target)
    zero('example x ideal-membership certificate',-f*f/2+gg/2-x)
    zero('example y ideal-membership certificate',(1+x*x*f)*f-x*x*gg-y)
    zero('example Jacobian remains Keller',target.jacobian((x,y)).det()+2)
    # Generic complete stabilizer formula; completeness is the written UFD proof.
    g1,g2=S.symbols('g1 g2')
    row=S.Matrix([[g1,g2]])
    kg=S.Matrix([-g2,g1])
    nilpotent=kg*row
    stabilizer=S.eye(2)+h*nilpotent
    zero('stabilizer fixes generic row',row*stabilizer-row)
    zero('stabilizer determinant one',stabilizer.det()-1)
    zero('stabilizer square-zero term',nilpotent*nilpotent)
    zero('stabilizer polynomial inverse',stabilizer*(S.eye(2)-h*nilpotent)-S.eye(2))
    kgamma=S.Matrix([-gamma[1],gamma[0]])
    zero('normalizer transports stabilizer kernel',C0*kell-kgamma)
    gauge=C0+h*kgamma*ell
    zero('all-polynomial gauge row identity',gamma*gauge-ell)
    zero('all-polynomial gauge determinant',gauge.det()-1)
    check('gauge coefficient matrix degree twelve',degree(kgamma*ell,x,y)==12)
    high=gauge.subs(h,x**100)
    check('degree112 correction still gives a linear row',degree(high,x,y)==112)
    normalized=S.Matrix([[1,0],[x**99,1]])
    C109=C0*normalized
    zero('normalized degree99 matrix determinant',normalized.det()-1)
    zero('normalized-box strict extension output',gamma*C109-S.Matrix([[y+2*x**100,2*x]]))
    check('normalized-box original matrix degree109',degree(C109,x,y)==109)
    # Actual cotangent determinant, with arbitrary upper Hessian block.
    h11,h12,h22,j11,j12,j21,j22=S.symbols('h11 h12 h22 j11 j12 j21 j22')
    J=S.Matrix([[j11,j12],[j21,j22]])
    top=S.Matrix([[h11,h12],[h12,h22]])
    H=top.row_join(2*J.T).col_join((2*J).row_join(S.zeros(2)))
    zero('full relative cotangent Hessian determinant',H.det()-16*J.det()**2)
    # Known dimension-three control, independently replayed for this boundary use.
    z,a,b,c,t,s,TT,SS=S.symbols('z a b c t s T S')
    u=1+x*y
    F=S.Matrix([u**3*z+y*y*u*(4+3*x*y),
                y+3*x*u*u*z+3*x*y*y*(4+3*x*y),
                2*x-3*x*x*y-x**3*z])
    zero('dimension-three control Jacobian',F.jacobian((x,y,z)).det()+2)
    cubic=c*TT**3-2*TT*TT*SS+b*TT*SS*SS-2*a*SS**3
    source_values={a:F[0],b:F[1],c:F[2],TT:u,SS:x}
    zero('control global binary-cubic identity',cubic.subs(source_values,simultaneous=True))
    zero('control homogeneous T derivative',S.diff(cubic,TT).subs(source_values,simultaneous=True)-2*x)
    zero('control homogeneous S derivative',S.diff(cubic,SS).subs(source_values,simultaneous=True)+2*u)
    zero('global projective coordinates generate the unit ideal',u-x*y-1)
    # Literal maximal-ideal certificate; no Groebner decision is used.
    aux=y*(1+3*x*y)+x*u*z
    row_y=S.Matrix([[-3*x+S.Rational(9,2)*x*x*aux,u-x*(F[1]+y)/2,0]])
    aux_z=x*(u*u+u+1)*z+y*u*(4+3*x*y)
    row_z=S.Matrix([[1,0,0]])-aux_z*row_y
    row_x=S.Matrix([[0,0,S.Rational(1,2)]])+S.Rational(3,2)*x*x*row_y+x**3*row_z/2
    for name,row,expected in [('x',row_x,x),('y',row_y,y),('z',row_z,z)]:
        zero('control literal maximal-ideal certificate '+name,(row*F)[0]-expected)
    # Both etale charts of the finite cubic normalization.
    D=3*c*t*t-4*t+b
    xr=2/D
    yr=t-D/2
    zr=(2*xr-3*xr*xr*yr-c)/xr**3
    target_a=(c*t**3-2*t*t+b*t)/2
    image=F.subs({x:xr,y:yr,z:zr},simultaneous=True)
    zero('simple finite-root chart inverse',image-S.Matrix([target_a,b,c]))
    zero('finite-root chart root recovered',(1+xr*yr)/xr-t)
    E=-2+2*b*s-6*a*s*s
    xs=-2*s/E
    ys=b-3*a*s
    zs=-a*E**3/8-ys*ys*E*E/4+3*ys*ys*E/2
    target_c=2*s-b*s*s+2*a*s**3
    image=F.subs({x:xs,y:ys,z:zs},simultaneous=True)
    zero('simple infinite-root chart inverse',image-S.Matrix([a,b,target_c]))
    zero('infinite-root chart root recovered',xs/(1+xs*ys)-s)
    zero('finite zero fibre has double-plus-simple support',cubic.subs({a:0,b:0,c:0})+2*TT*TT*SS)
    points=[S.Matrix([0,0,-t*t]),S.Matrix([-1/(2*t),3*t,26*t*t]),S.Matrix([1/(2*t),-3*t,26*t*t])]
    for i,point in enumerate(points):
        zero('complete escaping fibre point '+str(i),F.subs(dict(zip((x,y,z),point)),simultaneous=True)-S.Matrix([-t*t,0,0]))
    sources={}
    for row in json.loads((HERE/'sources/manifest.json').read_text()):
        actual=hashlib.sha256((HERE/'sources'/row['file']).read_bytes()).hexdigest()
        check('primary import source '+row['file'],actual==row['sha256'])
        sources[row['url']]=actual
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),checks=checks,
        normalizer=[list(map(str,C0.row(i))) for i in range(2)],normalizer_degree=10,
        orbit='All ordered generators of (x,y) with local Jacobian orientation -2',
        stabilizer='I+h*(-g2,g1)^T*(g1,g2)',
        unbounded_gauge_degree='degree(h)+12; unchanged linear output',
        normalized_sufficient_gate=99,strict_extension_original_degree_example=109,
        control_generic_degree=3,control_zero_source_fibre_length=1,
        control_zero_normalization_lengths=[1,2],
        control_literal_ideal_certificate=[list(map(str,row)) for row in [row_x,row_y,row_z]],
        general_missing_length='d-1; written finite-normal-surface flatness argument',
        unrestricted_singleton_zero_plane_problem='OPEN',
        formal_verification=False,external_review=False,uses_Tang=False,
        primary_source_sha256=sources,elapsed_seconds=round(time.monotonic()-start,3),
        runtime={'python':platform.python_version(),'sympy':S.__version__}),indent=2))


if __name__=='__main__':
    main()
