"""Finite controls for the full filtered gate. Universal claims use the proof."""
from pathlib import Path
import hashlib,json
import sympy as s
import filtered_gaussian_full_safe_gate as G
import verify_semipositive_gaussian_reduction as R

p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
t,z,e=G.t,R.z,R.e


def main():
    x,y,w,c=s.symbols('x y w c');pair=s.Matrix([[0,1],[1,0]])
    cases=[('nonlinear safe',(z*z,),(z,),z*x**3,(x,y),(s.Rational(2,3),s.Rational(4,3)),pair,[]),
           ('pure unsafe',(1+z,),(z,),-(1+z)*x*y+z*x**3,(x,y),(s.Rational(2,3),s.Rational(4,3)),pair,[1]),
           ('constant and variance defects',(1,),(z,),c+e*x*x/2,(x,),(1,),s.eye(1),[e,c]),
           ('residual alone false positive',(e+z,),(z,),0,(),(),s.zeros(0),[1]),
           ('nilpotent polynomial lift',(1+e*z,),(z,),0,(),(),s.zeros(0),[e]),
           ('nilpotent nonpure',(1+e*z,),(z,),e*(x*x-1)/2,(x,),(1,),s.eye(1),[e]),
           ('nonzero critical action',(1,),(z,),x**3+y,(x,y),(s.Rational(2,3),s.Rational(4,3)),pair,[1]),
           ('empty zero block',(),(),e*(x*x-1)/2,(x,),(1,),s.eye(1),[e])]
    controls=[]
    for name,H,Z,F,X,weights,precision,expected in cases:
        data=G.full_gate(H,Z,F,X,weights,precision)
        generators=[s.sympify(v,locals={'e':e,'c':c}) for v in data['safe_ideal_generators']]
        assert s.groebner(generators,e,c,domain=s.QQ)==s.groebner(expected,e,c,domain=s.QQ),name
        controls.append(dict(name=name,gate=data,expected_ordinary_ideal=[str(v) for v in expected]))
    assert controls[3]['gate']['coefficient_blocks']['fixedpoint_residual'] != []
    assert '-1' in controls[3]['gate']['coefficient_blocks']['normal_Jacobian']
    # Ordinary polynomial Newton corrections modulo e^4: a genuine degree lift.
    trunc=lambda f:s.series(s.cancel(f),e,0,4).removeO().expand()
    q=t+e*t*t;newton=[]
    for i in range(2):
        error=trunc(q-t*(1+e*q));inverse=trunc(1/(1-e*t))
        q=trunc(q-inverse*error)
        newton.append(dict(iteration=i+1,polynomial=str(q),residual=str(trunc(q-t*(1+e*q)))))
    assert q==t+e*t*t+e**2*t**3+e**3*t**4
    # Direct original Wick controls with nonpure safe ordinary nilpotents.
    seed=1+e*t
    original=R.original(w*(1+e*z)+e*(x*x-1)/2,(z,w,x),s.diag(pair,s.eye(1)),
        [('unit',s.Integer(1),seed),('Z',z,t+2*e*t*t),('W',w,s.Integer(0)),
         ('T2',x*x,1+2*e*t),('ZW',z*w,1+2*e*t)],N=4,mod_nil=True)
    # A denominator-free safe two-coordinate branch exercises the vector gate.
    z2=s.Symbol('z2')
    vector=G.full_gate((1,z*z),(z,z2),0,(),(),s.zeros(0))
    assert vector['safe_ideal_generators']==[]
    assert vector['truncation']==[str(t),str(t**3)]
    # Every polynomial f*(a/b)^m has unbounded poles if b has a pole zero.
    valuation=[]
    for m in range(5):
        value=(1-t)**2*(t/(1-t))**m
        valuation.append(dict(m=m,period=str(s.factor(value)),pole_order=max(0,m-2)))
    assert s.denom(s.cancel((1-t)**2*(t/(1-t))**3))!=1
    rejects=[]
    for fn in [lambda:G.full_gate((1,),(z,),x**3,(x,),(1,),s.eye(1)),
               lambda:G.full_gate((t+z,),(z,),0,(),(),s.zeros(0)),
               lambda:G.full_gate((x+z,),(z,),x*x,(x,),(1,),s.eye(1)),
               lambda:G.full_gate((1,),(z,),z, (z,),(1,),s.eye(1))]:
        try:fn()
        except ValueError:rejects.append(True)
        else:raise AssertionError('Accepted invalid architecture')
    packet=dict(gate_controls=controls,vector_gate=vector,polynomial_Newton_controls=newton,
        original_nilpotent_Wick_controls=original,rational_pole_controls=valuation,rejections=rejects,
        universal_claims_use_written_proofs=True,purity_premise_required=False,
        general_unfiltered_safety_implies_purity_resolved=False,actual_image_kernel_computed=False)
    out=p/'filtered_gaussian_full_safe_gate_certificate.json';out.write_text(json.dumps(packet,indent=2)+'\n')
    imports=['verify_filtered_gaussian_critical_action.py','verify_semipositive_gaussian_reduction.py',
             'verify_algebraic_gaussian_lifts.py','FILTERED_GAUSSIAN_CRITICAL_ACTION_AND_EXACT_ALL_RING_GATES.md',
             'POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md','BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md']
    runtime=dict(python=__import__('sys').version,sympy=s.__version__,sealed_imports_sha256={n:sha(p/n) for n in imports})
    runtime_path=p/'filtered_gaussian_full_safe_gate_runtime.json';runtime_path.write_text(json.dumps(runtime,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'FILTERED_SAFETY_WITHOUT_A_PURITY_PREMISE.md'),
        producer_sha256=sha(p/'filtered_gaussian_full_safe_gate.py'),checker_sha256=sha(Path(__file__)),
        certificate_sha256=sha(out),runtime_sha256=sha(runtime_path),
        exact_controls=dict(gate_models=len(controls),vector_gates=1,Newton_steps=len(newton),
                            original_Wick_identities=25,rational_pole_controls=len(valuation),rejections=len(rejects)),
        independent_replay=False,formal_verification=False,external_review=False,
        universal_claims_use_written_proofs=True,general_unfiltered_safety_implies_purity_resolved=False)
    (p/'filtered_gaussian_full_safe_gate_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
