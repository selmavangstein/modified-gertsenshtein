(*============*)
(*  Geometry  *)
(*============*)

DemonstrateComponents[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=ComponentArray;
	Expr//=ToValues;
	Expr//=ToValues;
	Expr//=MatrixForm;	
	(InputExpr->Expr)//DisplayExpression;
];

Comment@"Setup of manifold and metric.";
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
