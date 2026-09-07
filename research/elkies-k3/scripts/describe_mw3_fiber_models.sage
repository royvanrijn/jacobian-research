from sage.all import *

# Root lattice A10 + A2 + A1^2.
# In characteristic zero:
# A10 -> I11 only
# A2  -> I3 or IV
# A1  -> I2 or III
# K3 total Euler number = 24.
opts_A2=[("I3",3),("IV",4)]
opts_A1=[("I2",2),("III",3)]

print("MW3FIB|stage=start|ADE=A10+A2+A1+A1|torsion=trivial")
for a2,e2 in opts_A2:
    for a1,e1 in opts_A1:
        for b1,eb in opts_A1:
            used=11+e2+e1+eb
            residual=24-used
            print(f"MW3FIB|I11+{a2}+{a1}+{b1}|e_used={used}|residual_euler={residual}")
print("MW3FIB|semistable_preferred=I11+I3+I2+I2+6I1")
print("MW3FIB|short_weierstrass_normalization=I11@infinity,I3@0,I2@1,I2@lambda")
print("MW3FIB|Delta_degree=13|Delta=t^3*(t-1)^2*(t-lambda)^2*R6(t)")
print("MW3FIB|family_dimension_before_MW_jump=4")
print("MW3FIB|target_MW_rank=3|target_regulator=79/11")
