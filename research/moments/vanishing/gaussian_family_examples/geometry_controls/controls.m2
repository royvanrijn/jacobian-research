load "/home/royvanrijn/src/jacobian-research/research/moments/vanishing/gaussian_family_core.m2";
W=QQ[x,y,dx,dy,WeylAlgebra=>{x=>dx,y=>dy}];
C=W^1/ideal(y-x^2,dx+2*x*dy);
O=W^1/ideal(dx,dy); B=QQ[x,y]; S=ideal(y-x^2);
r1=gGenericRank(C,S);r2=gGenericRank(C++C,S);r3=gGenericRank(C++O,S);
assert(r1#0==1 and r2#0==2 and r3#0==1);
W=QQ[t,c0,dt,dc0,WeylAlgebra=>{t=>dt,c0=>dc0}];
O=W^1/ideal(dt,dc0); Delta=W^1/ideal(dt,c0);
L1=localCohom({0,1},ideal(c0),O,Strategy=>OaTa);
L2=localCohom({0,1},ideal(c0),Delta,Strategy=>OaTa);
assert(gIsZero(L1#0) and not gIsZero(L1#1));
assert(not gIsZero(L2#0) and gIsZero(L2#1));
B=QQ[c0];
Graph=W^1/ideal(W_1-W_0^2,W_2+2*W_0*W_3);
bad=gSourceBad(Graph,ideal(0_B));
assert(#bad==1 and first bad==ideal(c0));
assert(not isHolonomic(W^1/ideal(W_0*W_2-W_1)));
W=QQ[x,dx,WeylAlgebra=>{x=>dx}];
K=Drestriction(ideal(x^2*dx+1),{1});
assert(gIsZero(K#0) and gIsZero(K#1));
H=new HashTable from {
  "nonlinear_conormal_ranks"=>{r1#0,r2#0,r3#0},
  "ordinary_restriction_uses_codimension_cohomology"=>true,
  "punctual_source_degree_zero_would_be_wrong"=>true,
  "vertical_graph_critical_parameter"=>gSerializeIdeal first bad,
  "central_Euler_is_not_jointly_holonomic"=>true,
  "irregular_origin_restriction_is_zero"=>true};
stdio << "HELPERS_COMPLETE|" << toJSON H << endl;
exit 0;
