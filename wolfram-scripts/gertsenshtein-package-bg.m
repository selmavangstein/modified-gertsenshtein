(* ::Package:: *)

(*==============================================*)
(*  From Lagrangian to component field equations  *)
(*==============================================*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
<<xAct`xPlain`;

Title@"Scripts for evaluating the field equations from a Lagrangian";

Comment@"Loading packages, if not loaded already...";
Print@"Loading packages, if not loaded already...";
Off@ValidateSymbol::used;
<<xAct`xTensor`;
<<xAct`xTras`;
<<xAct`xCoba`;
On@ValidateSymbol::used;
Comment@"Using ParallelNeeds...";
Print@"Using ParallelNeeds...b";
ParallelNeeds["xAct`xPlain`"];
ParallelNeeds["xAct`xTensor`"];
ParallelNeeds["xAct`xTras`"];
ParallelNeeds["xAct`xCoba`"];
$DefInfoQ=False;
(*Unprotect@AutomaticRules;
Options[AutomaticRules]={Verbose->False};
Protect@AutomaticRules;*)

$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;
SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];

Get@FileNameJoin@{$ThisDirectory,"EvaluateLagrangianBG.m"};
