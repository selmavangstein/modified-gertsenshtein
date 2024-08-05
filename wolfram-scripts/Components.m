(*==============*)
(*  Components  *)
(*==============*)

Comment@"Setting chart and all components";
(*\[ScriptH]:*)
zerovalues=Table[0,{i,0,3},{j,0,3}];
ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues];
DefScalarFunction[{h1,h2}];
ComponentValue[H[{1,-cartesian},{1,-cartesian}],h1[t[],z[]]];
ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]];
ComponentValue[H[{1,-cartesian},{2,-cartesian}],h2[t[],z[]]];
ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]];
(*F:*)
ComponentValue[ComponentArray@F[-{a,cartesian},-{b,-cartesian}],zerovalues];
DefConstantSymbol[Bx];
ComponentValue[F[{2,-cartesian},{3,-cartesian}],Bx];
F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm;
ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]];
(*\[ScriptCapitalF]:*)
AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues];
DefScalarFunction[\[ScriptB]];
ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],-\[ScriptB][t[],z[]]];
DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}];
ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}];
ChangeComponents[pertF[{a,cartesian},{b,cartesian}],pertF[-{a,cartesian},-{b,cartesian}]];
(*\[ScriptCapitalT]:*)
pertTtoVec=MakeRule[{pertT[a,-b,-c],epsilonmetric[a,-b,-c,-d]pertQ[d]+1/2 (-delta[a,-c]  pertU[-b]+delta[a,-b]   pertU[-c])},MetricOn->All,ContractMetrics->True];
funcPertTtoVec[Expr_]:=Module[{expr=Expr},
	expr=expr/.pertTtoVec;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr
];
DefScalarFunction[{\[ScriptQ]0,\[ScriptQ]2,\[ScriptU]0,\[ScriptU]2}];
AllComponentValues[pertQ[-{a,cartesian}],{\[ScriptQ]0[t[],z[]],0,\[ScriptQ]2[t[],z[]],0}];
ChangeComponents[pertQ[{a,cartesian}],pertQ[-{a,cartesian}]];
AllComponentValues[pertU[-{a,cartesian}],{\[ScriptU]0[t[],z[]],0,\[ScriptU]2[t[],z[]],0}];
ChangeComponents[pertU[{a,cartesian}],pertU[-{a,cartesian}]];
(*T:*)
DefConstantSymbol[q0];
AllComponentValues[Q[-{a,cartesian}],{q0,0,0,0}];
ChangeComponents[Q[{a,cartesian}],Q[-{a,cartesian}]];
funcChristCartZero[expr_]:=expr/.ChristoffelCDPDcartesian->Zero;
DefTensor[u[a],M];
u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True];
AllComponentValues[u[-{a,cartesian}],{1,0,0,0}];
