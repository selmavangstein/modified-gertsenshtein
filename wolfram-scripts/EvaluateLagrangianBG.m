(* ::Package:: *)

(*========================*)
(*  EvaluateLagrangianBG  *)
(*========================*)

Get@FileNameJoin@{$ThisDirectory,"Parallelisation.m"};
Get@FileNameJoin@{$ThisDirectory,"Geometry.m"};
Get@FileNameJoin@{$ThisDirectory,"DefFields.m"};
Get@FileNameJoin@{$ThisDirectory,"DefRules.m"};
Get@FileNameJoin@{$ThisDirectory,"DeleteFirstOrderPart.m"};
Get@FileNameJoin@{$ThisDirectory,"Components.m"};
Get@FileNameJoin@{$ThisDirectory,"InitialExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"FurtherExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"SimplifyLagrangian.m"};
Get@FileNameJoin@{$ThisDirectory,"ImposeBackgroundTorsion.m"};
Get@FileNameJoin@{$ThisDirectory,"SuperChangeCovD.m"};

EvaluateLagrangianBG[Lagrangian_,resultsFileName_]:=Module[{lagrangian=Lagrangian,torsionC,einsteinC,maxwellC,maxwellCurl},

	lagrangian//=InitialExpand;
	lagrangian//=FurtherExpand;
	lagrangian//=SimplifyLagrangian;
	lagrangian//=ImposeBackgroundTorsion;

	Comment@"Calculating torsion field equations";
	torsionField=ApplyParallel[lagrangian, {VarD[pertT[k,-l,-m],CDT],ScreenDollarIndices}];
	torsionField//=Expand;
	torsionField=ApplyParallel[torsionField, {ToCanonical,ContractMetric,ScreenDollarIndices}];
	torsionField//=Expand;
	torsionField//=SuperChangeCovD;
	torsionField//=Expand;
	torsionField=ApplyParallel[torsionField,{funcLorentz,funcCD,funcTtoVec}];
	torsionField//=Expand;
	torsionField=ApplyParallel[torsionField,{DeleteFirstOrderPart}];

	Comment@"Calculating einstein field equations";
	einsteinField=ApplyParallel[lagrangian, {VarD[H[r,s],CDT]}];
	einsteinField//=Expand;
	einsteinField=ApplyParallel[einsteinField, {ToCanonical,ContractMetric,ScreenDollarIndices}];
	einsteinField//=Expand;
	einsteinField=ApplyParallel[einsteinField,{funcPertAtoF}];
	einsteinField//=Expand;
	einsteinField//=SuperChangeCovD;
	einsteinField//=Expand;
	einsteinField=ApplyParallel[einsteinField, {funcTtoVec,funcPertAtoF}];
	einsteinField//=Expand;
	einsteinField=ApplyParallel[einsteinField, {funcLorentz,funcCD}];
	einsteinField//=Expand;
	einsteinField=ApplyParallel[einsteinField,{DeleteFirstOrderPart}];

	Comment@"Calculating maxwell field equations";
	maxwellField=ApplyParallel[lagrangian,{VarD[pertA[k],CDT],ToCanonical,ContractMetric,ScreenDollarIndices}];
	maxwellField//=Expand;
	maxwellField=ApplyParallel[maxwellField,{funcPertAtoF}];
	maxwellField//=Expand;
	maxwellField//=SuperChangeCovD;
	maxwellField//=Expand;
	maxwellField=ApplyParallel[maxwellField,{funcTtoVec,funcLorentz,funcCD}];
	maxwellField//=Expand;
	maxwellField=ApplyParallel[maxwellField,{DeleteFirstOrderPart}];

	Comment@"Evaluating torsion component eqs...";
	torsionField//=Expand;
	torsionC=ApplyParallel[torsionField,{funcPertTtoVec,funcTtoVec,ToBasis[cartesian]}];
	torsionC//=Expand;
	torsionExpr=ApplyParallel[torsionC,{funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}];
	torsionExpr//DisplayExpression;

	Comment@"Evaluating Einstein component eqs...";
	einsteinField//=Expand;
	einsteinC=ApplyParallel[einsteinField,{funcPertTtoVec,funcTtoVec,funcPertAtoF, SeparateMetric[metric],ToCanonical,ToBasis[cartesian]}];
	einsteinC//=Expand;
	einsteinExpr=ApplyParallel[einsteinC, {funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical}];
	einsteinExpr//DisplayExpression;

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
	maxwellExpr//DisplayExpression;

	Comment@"Saving results...";
	DumpSave[FileNameJoin[{$ThisDirectory,"results",resultsFileName<>".mx"}],{lagrangian,maxwellField,maxwellExpr,einsteinField,einsteinExpr,torsionField,torsionExpr}];
	Comment@"Goodbye and thanks for all the fish";
{lagrangian,maxwellField,resultMaxwell,einsteinField,resultEinstein,torsionField,resultTorsion}];
