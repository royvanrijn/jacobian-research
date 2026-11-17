"""Exact local identities for the completed norm-target proof."""
def cofactor_root_cover(s,z):
    u=s*z;v=s*(s-z);w=-z*(s-z)
    D=s*z*(s-z);B=(2*z-s)*(2*s-z)*(s+z)
    return {'roots':(u,v,w),'D':D,'B':B,
      'h':D**2,'r_minus_c':-(u+v+w)}

def critical_weights(e):
    if e==2:return (1,1,1,1,1),(2,2)
    if e==3:return (2,2,3,3,3),(6,6)
    raise ValueError('The proof covers uniform ramification2 or3 only.')

def residue_exponent(weights,equation_degrees):
    return sum(weights)-1-sum(equation_degrees)
