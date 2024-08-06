(* ::Package:: *)

(*========*)
(*  main  *)
(*========*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg.m"};

DefConstantSymbol[Alp0,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(0\)]\)"];
DefConstantSymbol[Alp1,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)"];
DefConstantSymbol[Alp2,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(2\)]\)"];
DefConstantSymbol[Alp3,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(3\)]\)"];
DefConstantSymbol[Alp4,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(4\)]\)"];
DefConstantSymbol[Alp5,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(5\)]\)"];
DefConstantSymbol[Alp6,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(6\)]\)"];
DefConstantSymbol[Bet1,PrintAs->"\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\)"];
DefConstantSymbol[Bet2,PrintAs->"\!\(\*SubscriptBox[\(\[Beta]\), \(2\)]\)"];
DefConstantSymbol[Bet3,PrintAs->"\!\(\*SubscriptBox[\(\[Beta]\), \(3\)]\)"];

Comment@"Conventions regarding the curvature";

Quit[];

lagrangian=Sqrt[-Detmetric[]](
	\[Kappa] RicciScalarCDT[]
	-(1/4)*F[-a,-b]F[a,b]);

(*lagrangian=Sqrt[-Detmetric[]](\[Kappa] RicciScalarCDT[]+ \[Lambda] F[-a,-b]F[a,b]);*)
lagrangian=(Sqrt[-Detmetric[]](RicciScalarCDT[]+\[Kappa]4 RiemannCDT[-a,-b,-c,-d]RiemannCDT[a,b,c,d]+\[Lambda] F[-a,-b]F[a,b]));
resultsFileName="rr4-bg";

resultList=EvaluateLagrangianBG[lagrangian,resultsFileName];
(*EvaluateLagrangianNoBG[lagrangian,resultsFileName];*)(*-this will be a thing in the future*)
Quit[];
