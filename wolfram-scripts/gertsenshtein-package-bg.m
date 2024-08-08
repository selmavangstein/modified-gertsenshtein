(*============================*)
(*  gertsenshtein-package-bg  *)
(*============================*)

Off@ValidateSymbol::used;
Quiet@Block[{Print=Null},
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
];
$DefInfoQ=False;
On@ValidateSymbol::used;
Title@"The Gertsenshtein effect in gauge theories of gravity";

$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;
SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];

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

Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG.m"};
