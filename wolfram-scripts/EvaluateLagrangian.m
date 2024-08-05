(* ::Package:: *)

(*=====================*)
(*  EvaluateLagrangianBG  *)
(*=====================*)

EvaluateLagrangianBG[Lagrangian_,resultsFileName_]:=Module[{lagrangian=Lagrangian,linearizedAction,FtoA,AtoF,funcAtoF,toH,funcToH,
	topertA,funcToPertA,topertT,funcToPertT,pertFtoA,pertAtoF,funcPertAtoF,TtoVec,funcTtoVec,lorentz,commuteCD,
		funcLorentz,funcCD,torsionField,maxwellField,einsteinField,ToOrderH,ToOrderCDH,DeleteFirstOrderPart,zeroValues,
		pertTtoVec,funcPertTtoVec,funcChristCartZero,torsionC,einsteinC,maxwellC,maxwellCurl,
		 torsionExpr, einsteinExpr, maxwellExpr},
Comment@"Takes lagrangian, calculates field equations and evaluates the component equations with a constant torsion
	 background";
Print@"Takes lagrangian, calculates field equations and evaluates the component equations with a constant torsion
	 background";
Comment@"Everything is stashed together for now, will split up for readability later";

Get@FileNameJoin@{$ThisDirectory,"Parallelisation.m"};

Comment@"Defining tensor relationships - can this be a separate file?";
Comment@"Will - if I do this in a separate file, do I have to do something for these to be defined in this file as well?";
Print@"Defining tensor relationships";
(*\[Epsilon]:*)
Print@"epsilon 1";
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True];
Print@"epsilon 2";
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,-cartesian},{1,-cartesian},{2,-cartesian},{3,-cartesian}],1},MetricOn->All,ContractMetrics->True];

Print@"making rules";
FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True];
AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];
funcAtoF[Expr_]:=Module[{expr=Expr},
	expr=expr/.AtoF;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
Perturbationmetric[LI[n_],___]/;n>1:=0;
toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True];
funcToH[Expr_]:=Module[{expr=Expr},
	expr=expr/.toH;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M];
perturbationA[LI[n_],___]/;n>1:=0;
topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True];
funcToPertA[Expr_]:=Module[{expr=Expr},
	expr=expr/.topertA;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
DefTensorPerturbation[perturbationT[LI[order],a,-b,-c],TorsionCDT[a,-b,-c],M];
perturbationT[LI[n_],___]/;n>1:=0;
topertT= MakeRule[{perturbationT[LI[1],a,-b,-c],pertT[a,-b,-c]},MetricOn-> All, ContractMetrics-> True];
funcToPertT[Expr_]:=Module[{expr=Expr},
	expr=expr/.topertT;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];

pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True];
pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},
	MetricOn->All,ContractMetrics->True];
funcPertAtoF[Expr_]:=Module[{expr=Expr},
	expr=expr/.pertAtoF;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
DefTensor[Q[a],M];
TtoVec=MakeRule[{TorsionCDT[a,-b,-c],epsilonmetric[a,-b,-c,-d]Q[d]},MetricOn->All,ContractMetrics->True];
funcTtoVec[Expr_]:=Module[{expr=Expr},
	expr=expr/.TtoVec;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];



Comment@"Starting to manipulate the lagrangian";
Print@"Starting to manipulate the lagrangian";
lagrangian=lagrangian/.FtoA;
Print@"initial lagrangian";
Print@lagrangian;

Comment@"Writing the lagrangian in terms of the torsion tensor";
Print@"Writing the lagrangian in terms of the torsion tensor";
lagrangian=ChangeCovD[lagrangian,CDT,CD]//ToCanonical//ContractMetric//ScreenDollarIndices//CollectTensors;
Print@"ChangeCovD";
Print@lagrangian;
lagrangian=ChangeCurvature[lagrangian,CDT,CD];
lagrangian//=ChristoffelToGradMetric;
lagrangian//=ContractMetric;
lagrangian//=ToCanonical;
lagrangian//=ScreenDollarIndices;
Print["Written in terms of torsion:"];
Print[lagrangian];



Comment@"Expanding and linearizing the lagrangian";
Print@"Expanding and linearizing the lagrangian";
Comment@"This part is a little slow, but it breaks if I try to use ApplyParallel";
Print@"lilslow";
linearizedAction=PerturbBackground[lagrangian,2, BackgroundSolution->bgRules]//ExpandPerturbation//ToBackground//CollectTensors;
Comment@"Simplifying the linearized lagrangian";
Print@"Simplifying the linearized lagrangian";
linearizedAction//=Expand;
linearizedAction=ApplyParallel[linearizedAction,{funcToH,funcToPertA,funcAtoF,funcToPertT}];
linearizedAction = linearizedAction/.Sqrt[-Detmetric[]]->1;
linearizedAction//=ToCanonical;
linearizedAction//=ContractMetric;
linearizedAction//=ScreenDollarIndices;
Print["Expanded and linearized:"];
Print[lagrangian];



Comment@"Imposing background torsion";
Print@"Imposing background torsion";
linearizedAction//=Expand;
linearizedAction=ApplyParallel[linearizedAction,{funcTtoVec,ToCanonical,ContractMetric}];
Print["Imposed bg torsion:"];
Print[lagrangian];

Comment@"Adding in traceless \[ScriptH] and Lorentz gauge";
Print@"Adding in traceless \[ScriptH] and Lorentz gauge";
H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];
lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];
commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];
funcLorentz[Expr_]:=Module[{expr=Expr},
	expr=expr/.lorentz;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
funcCD[Expr_]:=Module[{expr=Expr},
	expr=expr/.commuteCD;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
linearizedAction=linearizedAction//funcLorentz//funcCD//BreakScalars;
linearizedAction//=Expand;
Print["Final linearized action:"];
Print[linearizedAction];
Quit[];









Comment@"Calculating torsion field equations";
Print@"Calculating torsion field equations";
torsionField=ApplyParallel[linearizedAction, {VarD[pertT[k,-l,-m],CDT]}];
torsionField//=Expand;
torsionField=ApplyParallel[linearizedAction, {ToCanonical,ContractMetric}];
torsionField=ChangeCovD[torsionField,CDT,CD]//ChristoffelToGradMetric//ToCanonical//ContractMetric//ScreenDollarIndices;
torsionField//=Expand;
torsionField=ApplyParallel[torsionField,{funcLorentz,funcCD,funcTtoVec}];

Comment@"Calculating einstein field equations";
Print@"Calculating einstein field equations";
einsteinField=ApplyParallel[linearizedAction, {VarD[H[r,s],CDT]}];
Print@"VarD done, toCanonical next";
einsteinField=ApplyParallel[linearizedAction, {ToCanonical,ContractMetric,ScreenDollarIndices}];
einsteinField//=Expand;
einsteinField=ApplyParallel[einsteinField,{funcPertAtoF}];
einsteinField=ChangeCovD[einsteinField,CDT,CD]//ChristoffelToGradMetric//ContractMetric;
einsteinField//=Expand;
einsteinField=ApplyParallel[einsteinField, {funcTtoVec,funcPertAtoF}];
einsteinField//=Expand;
einsteinField=ApplyParallel[einsteinField, {funcLorentz,funcCD}];

Comment@"Calculating maxwell field equations";
Print@"Calculating maxwell field equations";
maxwellField=ApplyParallel[linearizedAction,{VarD[pertA[k],CDT],ToCanonical,ContractMetric,ScreenDollarIndices}];
maxwellField//=Expand;
maxwellField=ApplyParallel[maxwellField,{funcPertAtoF}];
maxwellField=ChangeCovD[maxwellField,CDT,CD]//ChristoffelToGradMetric;
maxwellField//=Expand;
maxwellField=ApplyParallel[maxwellField,{funcTtoVec,funcLorentz,funcCD}];

Comment@"Simplifying field eqs";
Print@"Simplifying field eqs";
DefConstantSymbol[PerturbativeParameter, PrintAs->"\[Epsilon]"];
ToOrderH = MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},MetricOn->All,ContractMetrics->True];
ToOrderCDH=MakeRule[{CD[-c]@H[-a,-b],PerturbativeParameter*CD[-c]@H[-a,-b]},MetricOn->All,ContractMetrics->True];
DeleteFirstOrderPart[InputExpr_]:=Module[{OutputLagrangian=InputExpr,SecondOrderPart,FirstOrderPart,NullOrderPart},
	OutputLagrangian=OutputLagrangian/.ToOrderCDH/.ToOrderH;
	SecondOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,2}]&;
	SecondOrderPart//=Normal;
	SecondOrderPart=SecondOrderPart/.PerturbativeParameter->1;
	SecondOrderPart//=ToCanonical;
	SecondOrderPart//=ContractMetric;
	SecondOrderPart//=ScreenDollarIndices;
	SecondOrderPart//=CollectTensors;
	
	FirstOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,1}]&;
	FirstOrderPart//=Normal;
	FirstOrderPart=FirstOrderPart/.PerturbativeParameter->1;
	FirstOrderPart//=ToCanonical;
	FirstOrderPart//=ContractMetric;
	FirstOrderPart//=ScreenDollarIndices;
	FirstOrderPart//=CollectTensors;
	
	NullOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,0}]&;
	NullOrderPart//=Normal;
	NullOrderPart=NullOrderPart/.PerturbativeParameter->1;
	NullOrderPart//=ToCanonical;
	NullOrderPart//=ContractMetric;
	NullOrderPart//=ScreenDollarIndices;
	NullOrderPart//=CollectTensors;
	
	OutputLagrangian=SecondOrderPart-(FirstOrderPart-NullOrderPart);
	OutputLagrangian//=ToCanonical;
	OutputLagrangian//=ContractMetric;
	OutputLagrangian//=ScreenDollarIndices;
	OutputLagrangian//=CollectTensors;
OutputLagrangian];

einsteinField//=Expand;
maxwellField//=Expand;
torsionField//=Expand;
einsteinField=ApplyParallel[einsteinField,{DeleteFirstOrderPart}];
maxwellField=ApplyParallel[maxwellField,{DeleteFirstOrderPart}];
torsionField=ApplyParallel[torsionField,{DeleteFirstOrderPart}];

Print@"Will be doing coord stuff but quitting now";
Quit[];


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
DefTensor[pertQ[-a],M,PrintAs->"\[ScriptCapitalQ]"];
DefTensor[pertU[-a],M, PrintAs->"\[ScriptCapitalU]"];
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

Comment@"Evaluating torsion component eqs...";
funcChristCartZero[expr_]:=expr/.ChristoffelCDPDcartesian->Zero;
torsionField//=Expand;
torsionC=ApplyParallel[torsionField,{funcPertTtoVec,funcTtoVec,ToBasis[cartesian]}];
torsionC//Expand;
torsionExpr=ApplyParallel[torsionC,{funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}];

Comment["Evaluating Einstein component eqs..."];
einsteinField//=Expand;
einsteinC=ApplyParallel[einsteinField,{funcPertTtoVec,funcTtoVec,funcPertAtoF, SeparateMetric[metric],ToCanonical,ToBasis[cartesian]}];
einsteinC//=Expand;
einsteinExpr=ApplyParallel[einsteinC, {funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical}];

Comment@"Evaluating Maxwell component eqs";
DefTensor[u[a],M];
u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True];
AllComponentValues[u[-{a,cartesian}],{1,0,0,0}];
maxwellField//=Expand;
maxwellC=ApplyParallel[maxwellField,{funcPertTtoVec,funcTtoVec, ToCanonical,ToBasis[cartesian], ToBasis[cartesian]}];
maxwellC=maxwellC/.ChristoffelCDPDcartesian->Zero;
maxwellCurl=u[-{l,cartesian}]epsilonmetric[{l,cartesian},{i,cartesian},{f,cartesian},{k,cartesian}]CD[-{i,cartesian}][maxwellC];
maxwellCurl//=Expand;
maxwellExpr=ApplyParallel[maxwellCurl,{ContractMetric,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}];

Comment@"Saving results...";
DumpSave[FileNameJoin[{$ThisDirectory,"results",resultsFileName<>".mx"}],{\[ScriptCapitalL],maxwellCurl,maxwellC,maxwellField,maxwellExpr,einsteinField,einsteinExpr,torsionField,torsionExpr}];
Comment@"Goodbye and thanks for all the fish";
]
