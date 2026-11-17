"""Exact controls for the written proper/Morse boundary and monic free basis."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item.doit()))==0 for item in values))

    x,y,s,t = S.symbols('x y s t')
    z = (x,y)
    w = S.Matrix([s,t])
    K = S.Matrix([x+y*y,y+(x+y*y)**2])
    J = K.jacobian(z)
    zero('the nonlinear positive Keller control has determinant one',J.det()-1)
    M = J.T.inv()
    zero('the inverse transpose is a polynomial unimodular matrix',M*J.T-S.eye(2))
    check('the chosen Keller control has two genuinely nonlinear outputs',
        [S.Poly(v,*z).total_degree() for v in K]==[2,4])
    D = -(K.dot(K))/2
    fM = (w.dot(K))+D
    relative = S.Matrix([S.diff(fM,a) for a in z])
    zero('the universal Morse relative equations have the exact invertible factor',relative-J.T*(w-K))
    zero('the full actual Morse-model Hessian is a unit',S.hessian(fM,(s,t,x,y)).det()-1)
    graph = {s:K[0],t:K[1]}
    vertical = S.hessian(fM,z).subs(graph,simultaneous=True)
    zero('every point on its relative critical graph has the stated vertical Hessian',vertical+J.T*J)
    zero('that vertical Hessian determinant is the squared Keller constant',vertical.det()-1)
    zero('the graph equation recovers the whole relative ideal',M*relative-(w-K))
    zero('the full critical potential value is zero on K=0',
        fM.subs({s:0,t:0,x:0,y:0},simultaneous=True))
    p,q = S.symbols('p q')
    Kinv = {y:q-p*p,x:p-(q-p*p)**2}
    zero('the positive control has an actual polynomial plane inverse',K.subs(Kinv,simultaneous=True)-S.Matrix([p,q]))
    # Full gradient inversion is separate from ambient source substitutions.
    eta1,eta2 = S.symbols('eta1 eta2')
    recover_w = K+M*S.Matrix([eta1,eta2])
    grad_z = S.Matrix([S.diff(fM,a) for a in z])
    zero('the cotangent part of the full gradient has a literal polynomial inverse',
        grad_z.subs({s:recover_w[0],t:recover_w[1]},simultaneous=True)-S.Matrix([eta1,eta2]))

    # Fermat construction for this same, nonconstant inverse-transpose matrix.
    B = J.T
    degreeB = max(S.Poly(v,*z).total_degree() for v in B if v!=0)
    N = 4
    check('the Fermat exponent strictly exceeds every inverse-matrix entry degree',N>degreeB)
    powers = S.Matrix([x**N,y**N])
    P = M*powers
    zero('the finite-free graph equations have the full polynomial inverse-matrix identity',B*P-powers)
    DN = -(x**(N+1)+y**(N+1))/S.Integer(N+1)
    fN = DN+w.dot(K)
    zero('the finite-free-model full actual Hessian remains a unit',S.hessian(fN,(s,t,x,y)).det()-1)
    relativeN = S.Matrix([S.diff(fN,a) for a in z])
    zero('the proper-model relative equations have the prescribed graph',relativeN-J.T*(w-P))
    graphN = {s:P[0],t:P[1]}
    HN = S.hessian(fN,z).subs(graphN,simultaneous=True)
    zero('its graph differential again gives the full vertical Hessian identity',HN+J.T*P.jacobian(z))
    zero('the finite-free projection ramifies at the origin',P.jacobian(z).subs({x:0,y:0}))
    zero('the vertical Hessian vanishes at that relative critical point',HN.subs({x:0,y:0}))
    zero('the full Hessian is still nondegenerate at that same point',S.hessian(fN,(s,t,x,y)).det().subs({s:0,t:0,x:0,y:0})-1)

    # Actual monic Groebner controls over the target polynomial coefficient ring.
    def monic_control(name,variables,Bmat,exponent):
        size = len(variables)
        targets = S.symbols('target0:'+str(size))
        vec = S.Matrix([v**exponent for v in variables])
        Minv = Bmat.inv()
        Pmap = Minv*vec
        equations = list(vec-Bmat*S.Matrix(targets))
        gb = S.groebner(equations,*variables,domain=S.QQ[targets],order='grlex')
        leads = {poly.LM(order=gb.order).exponents for poly in gb.polys}
        expected = {tuple(exponent if i==j else 0 for i in range(size)) for j in range(size)}
        check(name+' has exactly the pairwise coprime pure-power leading monomials',leads==expected)
        check(name+' is monic over the full target coefficient ring',
            all(poly.LC(order=gb.order)==1 for poly in gb.polys))
        basis = [S.prod(v**a for v,a in zip(variables,exponents))
            for exponents in itertools.product(range(exponent),repeat=size)]
        check(name+' has the stated unreduced standard monomial basis',
            len(basis)==exponent**size and all(gb.reduce(m)[1]==m for m in basis))
        zero(name+' graph equations vanish after the actual polynomial target substitution',
            S.Matrix(equations).subs(dict(zip(targets,Pmap)),simultaneous=True))
        zero(name+' zero fibre is exactly the pure-power complete intersection',Bmat*Pmap-vec)
        for i in range(size):
            for j in range(i):
                sp = variables[j]**exponent*equations[i]-variables[i]**exponent*equations[j]
                zero(name+' monic overlap '+str((i,j))+' reduces to zero',gb.reduce(sp)[1])
        return exponent**size

    dimensions = [monic_control('nonlinear Keller n=2,N=4',z,B,N)]
    a,b,c = S.symbols('a b c')
    B3 = S.Matrix([[1,a,0],[0,1,b],[0,0,1]])
    zero('the dimension-three algebra control is unimodular',B3.det()-1)
    dimensions.append(monic_control('unimodular n=3,N=2',(a,b,c),B3,2))
    # Equality in the degree bound really can destroy finiteness.
    Bbad = S.Matrix([[1,x*x],[0,1]])
    Pbad = Bbad.inv()*S.Matrix([x*x,y*y])
    zero('the equality-bound matrix is nevertheless unimodular',Bbad.det()-1)
    zero('its equality-bound map contracts a whole affine line in one fibre',
        Pbad.subs(y,1)-S.Matrix([0,1]))
    check('the strict bound has genuinely failed in that negative control',
        max(S.Poly(v,x,y).total_degree() for v in Bbad if v!=0)==2)
    # Proper and vertical-Morse together: an actual finite degree-one section.
    f0 = s*x+t*y-(x*x+y*y)/2
    zero('the finite-Morse positive control has unit full Hessian',S.hessian(f0,(s,t,x,y)).det()-1)
    zero('its relative critical projection is the identity section',
        S.Matrix([S.diff(f0,x),S.diff(f0,y)])-S.Matrix([s-x,t-y]))
    zero('its vertical critical Hessian is everywhere invertible',S.hessian(f0,z).det()-1)
    # No hidden independent ambient embedding is asserted by either model.
    check('the finite-free degrees belong to the relative projection, not the Keller control',dimensions==[16,8])

    pins = {}
    for receipt_name,name in [
        ('receipt-relative-critical-graph.json','RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md'),
        ('receipt-critical-components.json','CRITICAL_COMPONENTS.md'),
        ('receipt-global-value.json','GLOBAL_VALUE_ORBITS.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    retained = json.loads((HERE/'receipt-intrinsic-plane.json').read_text())
    for name,digest in retained['dependencies'].items():
        if 'tang-hc4-audit' in name:
            pins[name] = digest
    retained = json.loads((HERE/'receipt-stable-hypersurface.json').read_text())
    for name,digest in retained['dependencies'].items():
        if 'hc4-relative-row-orbit/PROOF.md' in name:
            pins[name] = digest
    record = json.loads((HERE/'sources/proper-morse-boundary-import.json').read_text())
    for item in record['sources']:
        pins[str((HERE/'sources'/item['retained_file']).relative_to(ROOT))] = item['sha256']
    for name,digest in pins.items():
        check('retained mathematical or primary source pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_proper_morse_boundary.py','PROPER_MORSE_CRITICAL_BOUNDARY.md','sources/proper-morse-boundary-import.json']},
        dependencies=pins,finite_free_control_ranks=dimensions,
        written_not_finite_exhaustion=[
            'one proper relative Morse component, or one section component, suffices in every Tang projection rank',
            'every plane Keller pair has an everywhere-Morse relative critical encoding whose projection is exactly that pair',
            'arbitrary-field, arbitrary-dimensional monic finite-free lemma for every unimodular polynomial matrix and every exponent above the inverse degree',
            'finite-free relative critical projection encoding of every Keller pair and the two exact separate JC2 equivalences'],
        finite_fibres_asserted_finite_morphism=False,
        properness_alone_asserted_to_force_gradient_inverse=False,
        relative_Morse_alone_asserted_to_force_gradient_inverse=False,
        properness_of_Fermat_projection_asserted_properness_of_Keller_pair=False,
        degree_N_squared_asserted_Keller_degree=False,
        universal_potentials_asserted_closed_ambient_embeddings=False,
        full_Riemann_existence_essential_surjectivity_used=False,
        flatness_theorem_used_for_finite_free_basis=False,
        nonlinear_source_Hessian_invariance_used=False,
        nonautomorphic_plane_Keller_pair_constructed=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        new_external_embedding_classification_imported=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-proper-morse-boundary.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
