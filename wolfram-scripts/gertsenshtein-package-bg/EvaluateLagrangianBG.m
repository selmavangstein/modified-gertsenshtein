(*========================*)
(*  EvaluateLagrangianBG  *)
(*========================*)

Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"Parallelisation.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"Geometry.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"DefFields.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"DefRules.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"DeleteFirstOrderPart.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"Components.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"InitialExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"FurtherExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"SimplifyLagrangian.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"ImposeBackgroundTorsion.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"SuperChangeCovD.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"ParallelExpand.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"ToCanonicalCheck.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"Lagrangian.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG",
	"Cleaning.m"};

EvaluateLagrangianBG[Lagrangian_,resultsFileName_]:=Module[{lagrangian=Lagrangian,maxwellCurl},

	lagrangian//=InitialExpand;
	lagrangian//=FurtherExpand;
	lagrangian//=SimplifyLagrangian;
	lagrangian//=ImposeBackgroundTorsion;

	Comment@"Calculating Cartan field equations.";
	torsionField=ApplyParallel[lagrangian, {VarD[pertT[k,-l,-m],CDT],ScreenDollarIndices}];
	torsionField//=Expand;
	torsionField=ApplyParallel[torsionField, {ToCanonical,ContractMetric,ScreenDollarIndices}];
	torsionField//=Expand;
	torsionField//=SuperChangeCovD;
	torsionField//=Expand;
	torsionField=ApplyParallel[torsionField,{funcLorentz,funcCD,funcTtoVec}];
	torsionField//=Expand;
	torsionField=ApplyParallel[torsionField,{DeleteFirstOrderPart}];

	Comment@"Calculating Einstein field equations.";
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

	Comment@"Calculating Maxwell field equations.";
	maxwellField=ApplyParallel[lagrangian,{VarD[pertA[k],CDT],ToCanonical,ContractMetric,ScreenDollarIndices}];
	maxwellField//=Expand;
	maxwellField=ApplyParallel[maxwellField,{funcPertAtoF}];
	maxwellField//=Expand;
	maxwellField//=SuperChangeCovD;
	maxwellField//=Expand;
	maxwellField=ApplyParallel[maxwellField,{funcTtoVec,funcLorentz,funcCD}];
	maxwellField//=Expand;
	maxwellField=ApplyParallel[maxwellField,{DeleteFirstOrderPart}];

	DumpSave[FileNameJoin[{$ThisDirectory,"results","fields"<>resultsFileName<>".mx"}],{lagrangian,maxwellField,einsteinField,torsionField}];
	(*DumpSave[FileNameJoin[{$ThisDirectory,"results","FO"<>"fields"<>resultsFileName<>".mx"}],{lagrangian,maxwellField,einsteinField,torsionField}];*)

	Get@FileNameJoin@{$ThisDirectory,"results","fields"<>resultsFileName<>".mx"};
	(*Get@FileNameJoin@{$ThisDirectory,"results","FO"<>"fields"<>resultsFileName<>".mx"};*)
	$ErrorFiles=FileNames["results/BadEvaluation*"];
	DeleteFile/@$ErrorFiles;

	Comment@"Calculating Cartan component equations.";
	Comment@"The first parallelisation.";
	torsionField//=ParallelExpand;
	torsionField//=funcPertAtoF;
	Comment@"The second parallelisation.";
	torsionField//=(ApplyParallel[#,{
		funcPertTtoVec,
		funcTtoVec
	}])&;
	torsionField//=ToCanonical;
	torsionField//=ContractMetric;
	torsionField//=ScreenDollarIndices;
	Comment@"The third parallelisation.";
	torsionField//=(ApplyParallel[#,{
		ToBasis[cartesian]
	}])&;
	Comment@"The fourth parallelisation.";
	torsionField//=ParallelExpandMinimal;
	Commment@"The fifth parallelisation.";
	torsionField//=(ApplyParallel[#,{
		MultipleStepsTorsion
	}])&;
	torsionField//DisplayExpression;

	Comment@"Evaluating Einstein component equations.";
	einsteinField//=ParallelExpand;
	einsteinField//=funcPertAtoF;
	einsteinField//=funcPertAtoF;
	einsteinField//=funcPertAtoF;
	einsteinField//=funcPertAtoF;
	Comment@"The first parallelisation.";
	einsteinField//=(ApplyParallel[#,{
		funcPertTtoVec,
		funcTtoVec,
		funcPertAtoF,
		SeparateMetric[metric],
		ToCanonical,
		ToBasis[cartesian]
	}])&;
	Comment@"The second parallelisation.";
	einsteinField//=ParallelExpandMinimal;
	Comment@"The third parallelisation.";
	einsteinField//=(ApplyParallel[#,{
		MultipleStepsEinstein
	}])&;
	einsteinField//DisplayExpression;

	Comment@"Evaluating Maxwell component equations.";
	maxwellField//=ParallelExpand;
	maxwellField//=(ApplyParallel[#,{
		funcPertTtoVec,
		funcTtoVec,
		ToCanonical,
		ToBasis[cartesian],
		ToBasis[cartesian]
	}])&;
	maxwellField=maxwellField/.ChristoffelCDPDcartesian->Zero;
	maxwellField//=(u[-{l,cartesian}]epsilonmetric[{l,cartesian},{i,cartesian},{f,cartesian},{k,cartesian}]CD[-{i,cartesian}][#])&;
	maxwellField//=ParallelExpandMinimal;
	maxwellField//=(ApplyParallel[#,{
		ContractMetric,
		TraceBasisDummy,
		TraceBasisDummy,
		ComponentArray,
		ToValues,
		ToValues,
		ToCanonical,
		SeparateMetric[metric],
		ToBasis[cartesian],
		ToBasis[cartesian],
		TraceBasisDummy,
		TraceBasisDummy,
		ToCanonical,
		ToValues
	}])&;
	maxwellField//DisplayExpression;

	Comment@"Saving results...";
	DumpSave[FileNameJoin[{$ThisDirectory,"results",resultsFileName<>".mx"}],
		{lagrangian,maxwellField,einsteinField,torsionField}];
	(*DumpSave[FileNameJoin[{$ThisDirectory,"results","FO"<>resultsFileName<>".mx"}],
		{lagrangian,maxwellField,einsteinField,torsionField}];*)
];

ShowComponents[InputComponents_,InputRules_]:=Module[{myComponents=InputComponents},
	myComponents=myComponents/.InputRules;
	myComponents//=Flatten;
	myComponents//=DeleteDuplicates;
	myComponents//=((#==0)&/@#)&;
	myComponents//=((#//FullSimplify)&/@#)&;
	myComponents//=DeleteDuplicates;
	myComponents//=DeleteCases[#,True]&;
	DisplayExpression@(myComponents~MyRaggedBlock~1);
];

StudySystem[InputRules_]:=Module[{rules=InputRules,myLagrangian,myComponents,tList,eList,mList,expr,finalExpr},

	Comment@"Here is the list of rules.";
	DisplayExpression@(rules~MyRaggedBlock~5);

	Comment@"Here is the non-linear Lagrangian.";
	myLagrangian=lagrangian/.rules;
	myLagrangian//ToCanonical;
	myLagrangian//ContractMetric;
	myLagrangian//ScreenDollarIndices;
	myLagrangian//DisplayExpression;

	Subsection@"Here are the zeroth-order equations.";
	Get@FileNameJoin@{$ThisDirectory,"results","FO"<>resultsFileName<>".mx"};
	Comment@"The Cartan components.";
	ShowComponents[torsionField,rules];
	Comment@"The Einstein components.";
	ShowComponents[einsteinField,rules];
	Comment@"The Maxwell components.";
	ShowComponents[maxwellField,rules];

	Subsection@"Here are the first-order equations.";
	Get@FileNameJoin@{$ThisDirectory,"results",resultsFileName<>".mx"};
	Comment@"The Cartan components.";
	ShowComponents[torsionField,rules];
	Comment@"The Einstein components.";
	ShowComponents[einsteinField,rules];
	Comment@"The Maxwell components.";
	ShowComponents[maxwellField,rules];

	Subsection@"Here is the reduced set of first-order equations.";

	tList=torsionField//Flatten;
	eList=einsteinField//Flatten;
	mList=maxwellField//Flatten;
	expr=Join[tList,eList,mList];

	finalExpr=fullySimplify[expr,rules];
	DisplayExpression@(finalExpr~MyRaggedBlock~1);
];
