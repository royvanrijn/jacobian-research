"""Small exact controls for the geometric fiber-Hom proof;no general producer."""
from pathlib import Path
import hashlib,json,platform
import sympy as s


def shifted_degree(degree,shift):
    # H^i(C[k])=H^(i+k)(C).
    return degree-shift


def run():
    t,a,c,z,e,lam,xi,eta=s.symbols('t a c z e lambda xi eta')
    shifts=[]
    for d in range(5):
        # N=O_B[1] for M=O_(t,B).
        n_degree=-1
        dagger_degree=shifted_degree(n_degree,-d)
        final_degree=shifted_degree(dagger_degree,d-1)
        assert dagger_degree==d-1 and final_degree==0
        # N=O_B for M=delta_(t0) external O_B.
        point_dagger=shifted_degree(0,-d)
        point_ext=shifted_degree(point_dagger,d-1)
        assert point_ext==1
        shifts.append({'base_dimension':d,'constant_N_degree':n_degree,
                       'constant_dagger_degree':dagger_degree,'constant_Hom_degree':final_degree,
                       'time_point_Ext_degree':point_ext})

    # For exp(at),dual relative de Rham has cokernel delta_(a0).
    # t^j maps to (-1)^j partial_a^j delta. Verify the full-tower lowering
    # and surjectivity formulas on six indices;no truncated cokernel is used.
    point_rows=[]
    for j in range(6):
        v=z**j
        coordinate=-s.diff(v,z)
        preimage=-z**(j+1)/s.Integer(j+1)
        assert s.expand(-s.diff(preimage,z)-v)==0
        assert s.expand(-s.diff(z*v,z)-z*(-s.diff(v,z)))==-v
        de_rham=a*t**j-j*(t**(j-1) if j else 0)
        mapped=(-1)**j*coordinate-j*((-1)**(j-1)*z**(j-1) if j else 0)
        assert s.expand(mapped)==0
        point_rows.append({'tower_index':j,'coordinate_action':str(coordinate),
                           'coordinate_preimage':str(preimage),'de_Rham_relation':str(de_rham)})
    # Koszul coordinate kernel lies in degree-1;normalized codimension1
    # restriction shifts[-1],so the surviving fiber is degree0.
    assert shifted_degree(-1,-1)==0
    assert shifted_degree(shifted_degree(-1,-1),1-1)==0

    # The real Gaussian centered linear family has Hom0 generically,Hom1 at0.
    f=s.exp(a*a*t*t/2)
    assert s.diff(f,t)-a*a*t*f==0
    linear_rows=[]
    for n in range(4):
        coeffs=s.symbols('b:'+str(n+1));q=sum(coeffs[j]*t**j for j in range(n+1))
        top=s.Poly(s.diff(q,t)-a*a*t*q,t).nth(n+1)
        assert top==-a*a*coeffs[-1]
        assert s.expand((s.diff(q,t)-a*a*t*q).subs(a,0)-s.diff(q,t))==0
        linear_rows.append({'degree':n,'generic_top_equation':str(top),'zero_phase_polynomial_solution':'1'})

    # Graph c=t^2:conormal(-2t,1),with nonzero lambda.
    # Normalize lambda=1;zero time covector forces t0 and hence c0.
    gb=s.groebner([c-t*t,-2*t],t,c)
    eliminated=[g.as_expr() for g in gb.polys if t not in g.as_expr().free_symbols]
    assert eliminated==[c]
    assert s.Matrix([[-2*t,0],[1,1]]).det()==-2*t

    # Central Euler parameters are not an extra holonomic differential variable.
    probe=t**3*a*a+t*a
    def Euler(q):return t*s.diff(q,t)-a*q
    assert s.expand(s.diff(Euler(probe),a)-Euler(s.diff(probe,a)))==-probe
    central_rows=[]
    for n in range(5):
        assert s.expand(Euler(t**n).subs(a,n))==0
        central_rows.append({'central_parameter':n,'polynomial_solution':str(t**n)})

    # Exact genuine quartic input is retained elsewhere. These controls only
    # compare its known marker t/(1-et) with radicals and Artin thickness.
    artin=[]
    for length in range(1,5):
        mark=sum(e**(j-1)*t**j for j in range(1,length+1))
        assert s.Poly(mark,t).degree()==length
        assert s.rem(s.Poly(mark,t).nth(length),e**length,e)==e**(length-1)
        assert s.rem(e**length,e**length,e)==0
        # All positive powers have radical(e),but each next ordinary ideal
        # omits the preceding generator.
        assert s.Poly(e**length,e).sqf_part().as_expr()==e
        assert s.rem(e**length,e**(length+1),e)==e**length
        artin.append({'length':length,'last_marked_time_degree':length,
                      'cutoff_radical_generator':'e','ordinary_next_power_strict':True})

    return {'completed':True,'python_version':platform.python_version(),'sympy_version':s.__version__,
            'derived_degree_controls':shifts,'exponential_linear_point_tower_controls':point_rows,
            'derived_point_fiber_degree':0,'ordinary_point_fiber_zero_does_not_replace_derived':True,
            'centered_Gaussian_linear_Hom_controls':linear_rows,
            'graph_noncharacteristic_exception_ideal':['c'],'central_Euler_commutator':'-identity',
            'central_Euler_resonance_controls':central_rows,'quartic_Artin_and_radical_controls':artin,
            'Gaussian_field_uniformity_written_proof':True,
            'reduced_target_extension_written_proof':True,
            'universal_constructibility_inferred_from_samples':False,
            'general_holonomic_stratification_producer_executed':False,
            'general_cutoff_ideal_or_least_integer_producer_executed':False,
            'explicit_universal_numeric_cutoff_computed':False,
            'ordinary_cutoff_ideal_stabilization_claimed':False,
            'uniformity_over_all_nonreduced_algebras_claimed':False,
            'general_GVC_or_SIC_application_completed':False,
            'prior_quartic_all_order_input_reproved':False,
            'protected_ternary_frontier_changed':False,
            'independent_review':False,'formal_verification':False,'novelty_claim':False}


if __name__=='__main__':
    path=Path(__file__).resolve();result=run()
    result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('fiberwise_hom_uniformity_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
