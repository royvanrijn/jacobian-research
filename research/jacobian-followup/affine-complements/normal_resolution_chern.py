"""PF109 exact controls; no routine infers an actual end or constructs a Keller map."""
import sympy as S


def boundary_chern_gate(jacobian_coefficients, *, whole_plane_open):
    values=[S.sympify(x) for x in jacobian_coefficients]
    if not values or any(not x.is_Integer or x<0 for x in values):
        raise ValueError('nonempty nonnegative integral Jacobian vector required')
    return dict(boundary_components=len(values),nonzero_divisor=any(x!=0 for x in values),
                independent_boundary_classes_certified=whole_plane_open,
                rational_nonzero_chern_certified=whole_plane_open and any(x!=0 for x in values),
                proof_premises=['smooth rational algebraic resolution','entire marked plane open',
                                'no boundary unit relations','core retracts resolution'],
                necessary_only=True,realization_claimed=False)


def a1_wrong_open_model():
    v,t=S.symbols('v t')
    P=t*t*(1+v**3)
    Q=t**3*(1+S.Rational(3,2)*v**3+S.Rational(3,8)*v**6)
    J=S.factor(S.diff(P,v)*S.diff(Q,t)-S.diff(P,t)*S.diff(Q,v))
    return dict(v=v,t=t,P=P,Q=Q,Jacobian=J,canonical_vector=[8,4],
                principal_unit_relation=[2,1],relation_multiple=4,
                transverse_index=9,finite_degree=12,
                quotient_order=2,complement='C* times C',
                boundary_is_actual_A1_link=True,
                actual_Keller_source=False,rational_Chern_class_zero=True)


def ade_chern_comparison(kind,rank=None):
    import ade_end_exclusion as A
    data=A.ade_data(kind,rank)
    return dict(type=kind,rank=data['rank'],oriented_end=data['orientation'],
                unique_positive_tight=data['unique_positive_tight'],
                Stein_implies_minimal=True,
                minimal_filling_negative_definite=True,
                minimal_filling_rank=data['rank'],
                theta=data['canonical_theta'],
                every_Stein_filling_rational_chern_zero=True,
                singular_actual_Keller_normalization_excluded=True,
                actual_oriented_end_identification_still_required=True,
                classifications_imported=True)
