needsPackage "PrimaryDecomposition";needsPackage "JSON";
R=QQ[c0,c1];I=ideal matrix(R,{{(1)*c1^4+(2)*c0^2*c1^2+(1)*c0^4}});
qs=primaryDecomposition I;
if intersect qs != I then error "primary intersection mismatch";
encode=f->apply(listForm f,u->{u#0,toString u#1});
rows=apply(qs,Q->(
 P:=radical Q;
 if not isPrime P then error "nonprime radical";
 U:=support ((independentSets(P,Limit=>1))_0);
 if #U != dim(R/P) then error "wrong independent-set cardinality";
 new HashTable from {"primary_ideal"=>apply(first entries gens Q,encode),
 "radical_ideal"=>apply(first entries gens P,encode),
 "independent_indices"=>select(toList(0..numgens R-1),i->member(R_i,U)),
 "dimension"=>dim(R/P),"prime_radical_checked"=>true}
 ));
stdio << "PRIMARY_JSON|" << toJSON rows << endl;
print("PRIMARY_COMPLETE");exit 0;
