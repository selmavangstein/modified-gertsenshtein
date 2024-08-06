(*==============*)
(*  Components  *)
(*==============*)

Comment@"The metric perturbation in the transverse-traceless gauge.";
(*\[ScriptH]:*)
zerovalues=Table[0,{i,0,3},{j,0,3}];
Block[{Print=Null},
	ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues];
];
DefScalarFunction[{h1,h2}];
PrintAs@h1^="\!\(\*SubscriptBox[\(\[ScriptH]\), \(+\)]\)";
PrintAs@h2^="\!\(\*SubscriptBox[\(\[ScriptH]\), \(\[Times]\)]\)";
Block[{Print=Null},
	ComponentValue[H[{1,-cartesian},{1,-cartesian}],h1[t[],z[]]];
	ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]];
	ComponentValue[H[{1,-cartesian},{2,-cartesian}],h2[t[],z[]]];
	ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]];
];
H[-{a,cartesian},-{b,-cartesian}]//DemonstrateComponents;

Comment@"The background Maxwell tensor.";
Block[{Print=Null},
	ComponentValue[ComponentArray@F[-{a,cartesian},-{b,-cartesian}],zerovalues];
	DefConstantSymbol[Bx,PrintAs->"\!\(\*SubscriptBox[\(B\), \(x\)]\)"];
	ComponentValue[F[{2,-cartesian},{3,-cartesian}],Bx];
	F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm;
	ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]];
];
F[-{a,cartesian},-{b,cartesian}]//DemonstrateComponents;

Comment@"The perturbative Maxwell tensor.";
Block[{Print=Null},
	AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues];
	DefScalarFunction[\[ScriptB]];
	ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],-\[ScriptB][t[],z[]]];
	DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}];
	PrintAs@\[ScriptCapitalE]x^="\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\(x\)\)]\)";
	PrintAs@\[ScriptCapitalE]y^="\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\(y\)\)]\)";
	PrintAs@\[ScriptCapitalE]z^="\!\(\*SubscriptBox[\(\[ScriptCapitalE]\), \(\(z\)\)]\)";
	ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}];
	ChangeComponents[pertF[{a,cartesian},{b,cartesian}],pertF[-{a,cartesian},-{b,cartesian}]];
];
pertF[-{a,cartesian},-{b,cartesian}]//DemonstrateComponents;

Comment@"The perturbative torsion vectors.";
pertTtoVec=MakeRule[{pertT[a,-b,-c],epsilonmetric[a,-b,-c,-d]pertQ[d]+1/2 (-delta[a,-c]  pertU[-b]+delta[a,-b]   pertU[-c])},MetricOn->All,ContractMetrics->True];
funcPertTtoVec[Expr_]:=Module[{expr=Expr},
	expr=expr/.pertTtoVec;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr
];
DefScalarFunction[{\[ScriptQ]0,\[ScriptQ]2,\[ScriptU]0,\[ScriptU]2}];
PrintAs@\[ScriptQ]0^="\!\(\*SubscriptBox[\(\[ScriptQ]\), \(0\)]\)";
PrintAs@\[ScriptQ]2^="\!\(\*SubscriptBox[\(\[ScriptQ]\), \(2\)]\)";
PrintAs@\[ScriptU]0^="\!\(\*SubscriptBox[\(\[ScriptU]\), \(0\)]\)";
PrintAs@\[ScriptU]2^="\!\(\*SubscriptBox[\(\[ScriptU]\), \(2\)]\)";
Block[{Print=Null},
	AllComponentValues[pertQ[-{a,cartesian}],{\[ScriptQ]0[t[],z[]],0,\[ScriptQ]2[t[],z[]],0}];
	ChangeComponents[pertQ[{a,cartesian}],pertQ[-{a,cartesian}]];
	AllComponentValues[pertU[-{a,cartesian}],{\[ScriptU]0[t[],z[]],0,\[ScriptU]2[t[],z[]],0}];
	ChangeComponents[pertU[{a,cartesian}],pertU[-{a,cartesian}]];
];
pertQ[-{a,cartesian}]//DemonstrateComponents;
pertU[-{a,cartesian}]//DemonstrateComponents;

Comment@"The background torsion vector.";
DefConstantSymbol[q0];
PrintAs@q0^="\!\(\*SubscriptBox[\(q\), \(0\)]\)";
Block[{Print=Null},
	AllComponentValues[Q[-{a,cartesian}],{q0,0,0,0}];
	ChangeComponents[Q[{a,cartesian}],Q[-{a,cartesian}]];
];
Q[-{a,cartesian}]//DemonstrateComponents;

Comment@"The unit vector field.";
funcChristCartZero[expr_]:=expr/.ChristoffelCDPDcartesian->Zero;
DefTensor[u[a],M];
u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True];
Block[{Print=Null},
	AllComponentValues[u[-{a,cartesian}],{1,0,0,0}];
];
u[-{a,cartesian}]//DemonstrateComponents;
