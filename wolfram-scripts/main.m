(* ::Package:: *)

(*========*)
(*  main  *)
(*========*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg.m"};

DefConstantSymbol[\[Kappa]];
DefConstantSymbol[\[Kappa]4];
DefConstantSymbol[\[Lambda]];
lagrangian=Sqrt[-Detmetric[]](\[Kappa] RicciScalarCDT[]+ \[Lambda] F[-a,-b]F[a,b]);
(*lagrangian=(Sqrt[-Detmetric[]](RicciScalarCDT[]+\[Kappa]4 RiemannCDT[-a,-b,-c,-d]RiemannCDT[a,b,c,d]+\[Lambda] F[-a,-b]F[a,b]));*)
resultsFileName="rr4-bg";

resultList=EvaluateLagrangianBG[lagrangian,resultsFileName];
(*EvaluateLagrangianNoBG[lagrangian,resultsFileName];*)(*-this will be a thing in the future*)
Quit[];
