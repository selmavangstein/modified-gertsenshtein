
$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
Off@ValidateSymbol::used;
	$ThisDirectory=If[NotebookDirectory[]==$Failed,
		Directory[],
		NotebookDirectory[],
		NotebookDirectory[]];
	<<xAct`xPlain`;
	<<xAct`xTensor`;
	<<xAct`xTras`;
	<<xAct`xCoba`;
ParallelNeeds["xAct`xPlain`"];
ParallelNeeds["xAct`xTensor`"];
ParallelNeeds["xAct`xTras`"];
ParallelNeeds["xAct`xCoba`"];
$DefInfoQ=False;
On@ValidateSymbol::used;
$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;
SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];

Title@"Angular momentum localisation";

DisplayRule~SetAttributes~HoldAll;
DisplayRule[InputExpr_,InputRule_]:=Module[{Expr=Evaluate@InputExpr,EqnLabelValue=ToString@Defer@InputRule},

	EqnLabelValue//=StringDelete[#,"Defer["]&;
	EqnLabelValue//=StringDelete[#,"]"]&;

	Expr=Expr/.Evaluate@InputRule;
	Expr//=ToCanonical;
	Expr//=ContractMetric;
	Expr//=ScreenDollarIndices;
	Expr//=CollectTensors;
	Expr//=ScreenDollarIndices;
	DisplayExpression[(InputExpr->Expr),EqnLabel->EqnLabelValue];
Expr];

DemonstrateComponents[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=ComponentArray;
	Expr//=ToValues;
	Expr//=ToValues;
	Expr//=MatrixForm;	
	(InputExpr->Expr)//DisplayExpression;
];

DefConstantSymbol[MyOmega,PrintAs->"\[Omega]"];
DefConstantSymbol[MyC,PrintAs->"\[ScriptC]"];
DefConstantSymbol[MyG,PrintAs->"\[ScriptCapitalG]"];
DefConstantSymbol[MyK,PrintAs->"\[ScriptK]"];
DefScalarFunction[MyF,PrintAs->"\[ScriptF]"];

Comment@"The Minkowski metric tensor.";
DefManifold[M,4,IndexRange[{a,z}]]; 
DefMetric[-1,G[-a,-b],CD,PrintAs->"\[Eta]",SymCovDQ->True];
DefChart[cartesian,M,{0,1,2,3},{ctc[],xc[],yc[],zc[]}];
Block[{Print=Null},
	MetricInBasis[G, -cartesian,{1,-1,-1,-1}];
	MetricCompute[G,cartesian, All];
];
G[-{a,cartesian},-{b,cartesian}]//DemonstrateComponents;
bgRules = SymmetricSpaceRules[CD,0];
Block[{Print=Null},
	SetOptions[ToBackground,BackgroundSolution->bgRules];
];

(*Yet to be implemented*)
(*
AllComponentValues[epsilonG[{a,cartesian},{b,cartesian},{c,cartesian},{d,cartesian}],-$epsilonSign*Normal@LeviCivitaTensor@4];
AllComponentValues[epsilonG[{-a,-cartesian},{-b,-cartesian},{-c,-cartesian},{-d,-cartesian}],$epsilonSign*Normal@LeviCivitaTensor@4];
*)

PrintAs@ctc^="\[ScriptC]\[ScriptT]";
PrintAs@xc^="\[ScriptX]";
PrintAs@yc^="\[ScriptY]";
PrintAs@zc^="\[ScriptZ]";

Comment@"The metric perturbation.";
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];
zerovalues=Table[0,{i,0,3},{j,0,3}];
Block[{Print=Null},
	ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues];
];
Block[{Print=Null},
	ComponentValue[H[{1,-cartesian},{1,-cartesian}],-MyF[xc[],yc[]]*Cos[MyK*zc[]-MyOmega*ctc[]/MyC]];

	ComponentValue[H[{1,-cartesian},{2,-cartesian}],-MyF[xc[],yc[]]*Sin[MyK*zc[]-MyOmega*ctc[]/MyC]];
	ComponentValue[H[{2,-cartesian},{1,-cartesian}],-MyF[xc[],yc[]]*Sin[MyK*zc[]-MyOmega*ctc[]/MyC]];

	ComponentValue[H[{2,-cartesian},{2,-cartesian}],MyF[xc[],yc[]]*Cos[MyK*zc[]-MyOmega*ctc[]/MyC]];

	ComponentValue[H[{1,-cartesian},{3,-cartesian}],-(1/MyK)*(D[MyF[xc[],yc[]],yc[]]*Cos[MyK*zc[]-MyOmega*ctc[]/MyC]-D[MyF[xc[],yc[]],xc[]]*Sin[MyK*zc[]-MyOmega*ctc[]/MyC])];
	ComponentValue[H[{3,-cartesian},{1,-cartesian}],-(1/MyK)*(D[MyF[xc[],yc[]],yc[]]*Cos[MyK*zc[]-MyOmega*ctc[]/MyC]-D[MyF[xc[],yc[]],xc[]]*Sin[MyK*zc[]-MyOmega*ctc[]/MyC])];

	ComponentValue[H[{2,-cartesian},{3,-cartesian}],-(1/MyK)*(D[MyF[xc[],yc[]],xc[]]*Cos[MyK*zc[]-MyOmega*ctc[]/MyC]+D[MyF[xc[],yc[]],yc[]]*Sin[MyK*zc[]-MyOmega*ctc[]/MyC])];
	ComponentValue[H[{3,-cartesian},{2,-cartesian}],-(1/MyK)*(D[MyF[xc[],yc[]],xc[]]*Cos[MyK*zc[]-MyOmega*ctc[]/MyC]+D[MyF[xc[],yc[]],yc[]]*Sin[MyK*zc[]-MyOmega*ctc[]/MyC])];

	ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]];
];
H[-{a,cartesian},-{b,cartesian}]//DemonstrateComponents;

Comment@"The unit vector along the beamline.";
DefTensor[V[a],M,PrintAs->"\[ScriptN]"];
Block[{Print=Null},
	ComponentValue[V[{0,-cartesian}],0];
	ComponentValue[V[{1,-cartesian}],0];
	ComponentValue[V[{2,-cartesian}],0];
	ComponentValue[V[{3,-cartesian}],1];
	ChangeComponents[V[{a,cartesian}],V[-{a,cartesian}]];
];
V[{a,cartesian}]//DemonstrateComponents;

FunctionCheck::culpret="Function `1` will now act on the expression `2`.";
FunctionCheck[InputFunc_][InputExpr_]:=(InputExpr//InputFunc);
funcChristCartZero[expr_]:=expr/.ChristoffelCDPDcartesian->Zero;
CycleAverage[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=Integrate[#,{ctc[],0,2*Pi*MyC/MyOmega}]&;
	Expr/=2*Pi*MyC/MyOmega;
	Expr=Expr/.{MyK->MyOmega/MyC};
	Expr//=FullSimplify;
Expr];
ComputeMyComponents[InputExpr_]:=Module[{Expr=InputExpr,LLC},
	Expr//ToCanonical;
	Expr//ContractMetric;
	Expr//ScreenDollarIndices;
	Expr//=ChristoffelToGradMetric;
	Expr//=SeparateMetric[metric];
	Expr//ToCanonical;
	Expr//=ToBasis[cartesian];
	Expr//=FunctionCheck@funcChristCartZero;
	Expr//=FunctionCheck@ToBasis[cartesian];
	Expr//=FunctionCheck@funcChristCartZero;
	Expr//=FunctionCheck@ToBasis[cartesian];
	Expr//=FunctionCheck@funcChristCartZero;
	Expr//=FunctionCheck@TraceBasisDummy;
	Expr//=FunctionCheck@TraceBasisDummy;
	Expr//=FunctionCheck@TraceBasisDummy;
	Expr//=FunctionCheck@ComponentArray;
	Expr//=FunctionCheck@ToValues;
	Expr//=FunctionCheck@ToValues;
	Expr//=FunctionCheck@ToValues;
	Expr//=FunctionCheck@ToCanonical;
	Expr//=(CycleAverage/@#)&;
	Expr//=MatrixForm;
	Expr//DisplayExpression;
];

Section@"The Landau-Lifshitz pseudotensor";
Expr=(MyC^4/(16*Pi*MyG))*(CD[-l]@H[i,k]*CD[-m]@H[l,m]-CD[-l]@H[i,l]*CD[-m]@H[k,m]+(1/2)*G[i,k]*CD[-p]@H[l,n]*CD[-n]@H[p,-l]-CD[-p]@H[k,n]*CD[i]@H[p,-n]-CD[-p]@H[i,n]*CD[k]@H[p,-n]+CD[-n]@H[i,l]*CD[n]@H[k,-l]+(1/2)*CD[i]@H[n,-q]*CD[k]@H[q,-n]-(1/4)*CD[i]@H[n,-n]*CD[k]@H[p,-p]-(1/4)*G[i,k]*CD[-l]@H[n,-q]*CD[l]@H[q,-n]+(1/8)*G[i,k]*CD[-l]@H[n,-n]*CD[l]@H[p,-p]);
DisplayExpression[Expr,EqnLabel->"LandauLifshtzPseudotensor"];
Expr//ComputeMyComponents;

Section@"The Einstein pseudotensor";
Expr=(MyC^4/(32*Pi*MyG))*(CD[-m]@H[-r,-s]*CD[-n]@H[r,s]-(1/2)*G[-m,-n]*CD[-a]@H[-r,-s]*CD[a]@H[r,s]);
DisplayExpression[Expr,EqnLabel->"LandauLifshtzPseudotensor"];
Expr//ComputeMyComponents;

Section@"The affine tensor of Butcher";
Expr=(MyC^4/(8*Pi*MyG))*((1/4)*CD[-p]@H[-a,-b]*CD[-q]@H[a,b]-(1/8)*CD[-p]@H[a,-a]*CD[-q]@H[b,-b]-(1/8)*G[-p,-q]*(CD[-c]@H[-a,-b]*CD[c]@H[a,b]-(1/2)*CD[-c]@H[a,-a]*CD[c]@H[b,-b]));
DisplayExpression[Expr,EqnLabel->"LandauLifshtzPseudotensor"];
Expr//ComputeMyComponents;

(*Yet to be implemented*)
(*
Expr=epsilonG[n,m,-u,c]*V[u]*(H[-b,-n]-(1/2)*G[-b,-n]*H[y,-y])*(CD[a]@(H[-m,b]-(1/2)*G[-m,b]*H[z,-z])-CD[b]@(H[-m,a]-(1/2)*G[-m,a]*H[z,-z]));
DisplayExpression[Expr,EqnLabel->"LandauLifshtzPseudotensor"];
Expr//ComputeMyComponents;
*)

Supercomment@"This is the end of the script.";

Quit[];
