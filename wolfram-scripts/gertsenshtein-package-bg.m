(* ::Package:: *)

(*============================*)
(*  gertsenshtein-package-bg  *)
(*============================*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
<<xAct`xPlain`;

Title@"Scripts for evaluating the field equations from a Lagrangian";

Comment@"Loading packages, if not loaded already...";
Off@ValidateSymbol::used;
<<xAct`xTensor`;
<<xAct`xTras`;
<<xAct`xCoba`;
On@ValidateSymbol::used;
Comment@"Using ParallelNeeds...";
ParallelNeeds["xAct`xPlain`"];
ParallelNeeds["xAct`xTensor`"];
ParallelNeeds["xAct`xTras`"];
ParallelNeeds["xAct`xCoba`"];
$DefInfoQ=False;

$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;
SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];

Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG.m"};
