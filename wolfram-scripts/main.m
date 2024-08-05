(* ::Package:: *)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package.m"};
Get@FileNameJoin@{$ThisDirectory, "EvaluateLagrangian.m"};


(*Put your lagrangian here*)
DefConstantSymbol[\[Kappa]];
DefConstantSymbol[\[Lambda]];
lagrangian=Sqrt[-Detmetric[]](\[Kappa] RicciScalarCDT[]+ \[Lambda] F[-a,-b]F[a,b]);
resultsFileName="testmain";

EvaluateLagrangianBG[lagrangian,resultsFileName];
(*EvaluateLagrangianNoBG[lagrangian,resultsFileName];*)(*-this will be a thing in the future*)
