(*============*)
(*  Geometry  *)
(*============*)

Section@"Components of the fields, background, and perturbations";

DemonstrateComponents[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=ComponentArray;
	Expr//=ToValues;
	Expr//=ToValues;
	Expr//=MatrixForm;	
	(InputExpr->Expr)//DisplayExpression;
];

Comment@"The background metric tensor.";
DefManifold[M,4,IndexRange[{a,s}]]; 
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];
DefCovD[CDT[-a],Torsion->True, SymbolOfCovD->{"#","D"},FromMetric->metric];
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];
Block[{Print=Null},
	MetricInBasis[metric, -cartesian,{1,-1,-1,-1}];
	MetricCompute[metric,cartesian, All];
];
metric[-{a,cartesian},-{b,cartesian}]//DemonstrateComponents;
bgRules = SymmetricSpaceRules[CD,0];
Block[{Print=Null},
	SetOptions[ToBackground,BackgroundSolution->bgRules];
];
