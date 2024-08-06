(* ::Package:: *)

(*========================*)
(*  EvaluateLagrangianBG  *)
(*========================*)

Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Parallelisation.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Geometry.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","DefFields.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","DefRules.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","DeleteFirstOrderPart.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Components.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","InitialExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","FurtherExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","SimplifyLagrangian.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","ImposeBackgroundTorsion.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","SuperChangeCovD.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","ParallelExpand.m"};

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
	einsteinField=ApplyParallel[einsteinField, {ToCanonical,ContractMetric,ScreenDollarIndices}];
	einsteinField//=Expand;
	einsteinField=ApplyParallel[einsteinField,{funcPertAtoF}];
	einsteinField//=Expand;
	einsteinField//=SuperChangeCovD;
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
	torsionField//=ParallelExpand;
	torsionC=ApplyParallel[torsionField,{funcPertTtoVec,funcTtoVec,ToBasis[cartesian]}];
	torsionField//=ParallelExpand;
	torsionExpr=ApplyParallel[torsionC,{funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}];
	torsionExpr//DisplayExpression;

	Comment@"Evaluating Einstein component eqs...";
	einsteinField//=ParallelExpand;
	einsteinC=ApplyParallel[einsteinField,{funcPertTtoVec,funcTtoVec,funcPertAtoF, SeparateMetric[metric],ToCanonical,ToBasis[cartesian]}];
	einsteinField//=ParallelExpand;
	einsteinExpr=ApplyParallel[einsteinC, {funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical}];
	einsteinExpr//DisplayExpression;

	Comment@"Evaluating Maxwell component eqs";
	maxwellField//=ParallelExpand;
	maxwellC=ApplyParallel[maxwellField,{funcPertTtoVec,funcTtoVec, ToCanonical,ToBasis[cartesian], ToBasis[cartesian]}];
	maxwellC=maxwellC/.ChristoffelCDPDcartesian->Zero;
	maxwellCurl=u[-{l,cartesian}]epsilonmetric[{l,cartesian},{i,cartesian},{f,cartesian},{k,cartesian}]CD[-{i,cartesian}][maxwellC];
	maxwellCurl=ApplyParallel[maxwellCurl,{Expand}];
	maxwellExpr=ApplyParallel[maxwellCurl,{ContractMetric,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}];
	maxwellExpr//DisplayExpression;

	Comment@"Saving results...";
	DumpSave[FileNameJoin[{$ThisDirectory,"results",resultsFileName<>".mx"}],{lagrangian,maxwellField,maxwellExpr,einsteinField,einsteinExpr,torsionField,torsionExpr}];
	Comment@"Goodbye and thanks for all the fish";
{lagrangian,maxwellField,resultMaxwell,einsteinField,resultEinstein,torsionField,resultTorsion}];
