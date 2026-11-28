"""Hermite interpolation radial quadrature, trial implementation."""
from flint import arb,arb_mat
from check_three_cell_sharp_undamped_extension import upper_float


class HermiteQuadrature:
    def __init__(self, panels, radius, order):
        self.radius=arb(radius)
        self.width=self.radius/panels
        self.nodes=[j*self.width for j in range(panels+1)]
        self.order=order
        self.value_weights=arb_mat(panels+1,order+1)
        self.slope_weights=arb_mat(panels+1,order+1)
        self.panel_weights=arb_mat(panels,order+1)
        for i,(left,right) in enumerate(zip(self.nodes,self.nodes[1:])):
            left_power,right_power=arb(1),arb(1)
            moments=[]
            for n in range(order+4):
                left_power*=left
                right_power*=right
                moments.append((right_power-left_power)/(n+1))
            for n in range(order+1):
                z0,z1,z2,z3=moments[n:n+4]
                w=self.width
                s1=(z1-left*z0)/w
                s2=(z2-2*left*z1+left**2*z0)/w**2
                s3=(z3-3*left*z2+3*left**2*z1-left**3*z0)/w**3
                self.value_weights[i,n]+=2*s3-3*s2+z0
                self.value_weights[i+1,n]+=-2*s3+3*s2
                self.slope_weights[i,n]+=w*(s3-2*s2+s1)
                self.slope_weights[i+1,n]+=w*(s3-s2)
                self.panel_weights[i,n]=z0

    def radial_moments(self,d,kappa,beta,row):
        orders=(self.order,self.order//2,0)
        exponents=(beta,arb(1),1-kappa)
        assert all(exponent>0 for exponent in exponents)
        node_values=[[] for _ in exponents]
        node_slopes=[[] for _ in exponents]
        panel_errors=[[] for _ in exponents]
        for j,x in enumerate(self.nodes):
            q0,q1,q2,q3,q4=row[-1],arb(0),arb(0),arb(0),arb(0)
            for coefficient in reversed(row[:-1]):
                q4=q4*x+4*q3
                q3=q3*x+3*q2
                q2=q2*x+2*q1
                q1=q1*x+q0
                q0=q0*x+coefficient
            assert q0>0
            values=[8*d*(-exponent*x**2).exp()/q0 for exponent in exponents]
            slopes=[value*(-2*exponent*x-q1/q0) for exponent,value in zip(exponents,values)]
            if j:
                t1=q1/left_polynomial
                t2=q2/left_polynomial
                t3=q3/left_polynomial
                t4=q4/left_polynomial
                for errors,exponent,left_value in zip(panel_errors,exponents,left_values):
                    l1=2*exponent*x+t1
                    l2=2*exponent+t2+t1**2
                    l3=t3+3*t1*t2+2*t1**3
                    l4=t4+4*t1*t3+3*t2**2+12*t1**2*t2+6*t1**4
                    fourth_bound=left_value*(l1**4+6*l1**2*l2+3*l2**2+4*l1*l3+l4)
                    error=(fourth_bound*self.width**4/384).upper()
                    errors.append((-error).union(error))
            for output,value in zip(node_values,values):
                output.append(value)
            for output,value in zip(node_slopes,slopes):
                output.append(value)
            left_polynomial,left_values=q0,values

        radius=self.radius
        leading=max(range(len(row)),key=lambda j:upper_float(row[j]*radius**j))
        result=[]
        for values,slopes,errors,order,exponent in zip(node_values,node_slopes,panel_errors,orders,exponents):
            values_row=arb_mat(1,len(values),values)
            slope_row=arb_mat(1,len(slopes),slopes)
            errors_row=arb_mat(1,len(errors),errors)
            enclosure=values_row*self.value_weights+slope_row*self.slope_weights+errors_row*self.panel_weights
            moments=[]
            for n in range(order+1):
                power=n-leading-1
                decay=exponent-arb(max(0,power))/(2*radius**2)
                assert decay>0
                tail=4*d*radius**power*(-exponent*radius**2).exp()/(row[leading]*decay)
                moments.append(enclosure[0,n]+arb(0).union(tail))
            result.append(moments)
        tt,uu,(ww,)=result
        assert ww>0
        return tt,uu,ww
