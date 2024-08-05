(*============*)
(*  Geometry  *)
(*============*)

Comment@"Setup of manifold, metric, defining tensors";
DefManifold[M,4,IndexRange[{a,s}]]; 
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];
DefCovD[CDT[-a],Torsion->True, SymbolOfCovD->{"#","D"},FromMetric->metric];
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];
MetricInBasis[metric, -cartesian,{1,-1,-1,-1}];
MetricCompute[metric,cartesian, All];
bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];
