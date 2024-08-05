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

Get@FileNameJoin@{$ThisDirectory,"EvaluateLagrangian.m"};

Comment@"Setup of manifold, metric, defining tensors";
Print@"Setup of manifold, metric, defining tensors";
DefManifold[M,4,IndexRange[{a,s}]]; 
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];
DefCovD[CDT[-a],Torsion->True, SymbolOfCovD->{"#","D"},FromMetric->metric];

Print["Chart and basis on chart"];
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];
MetricInBasis[metric, -cartesian,{1,-1,-1,-1}];
MetricCompute[metric,cartesian, All];

bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];

DefTensor[A[-a], M];
DefTensor[F[-a,-b],M,Antisymmetric[{-a,-b}]];
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];
DefTensor[pertA[-a],M,PrintAs->"\[ScriptCapitalA]"];
DefTensor[pertF[-a,-b],M, Antisymmetric[{-a,-b}], PrintAs->"\[ScriptCapitalF]"];
DefTensor[pertT[a,-b,-c],M, Antisymmetric[{-b,-c}], PrintAs->"\[ScriptCapitalT]"];

Comment@"Ready for use";
Print@"Ready for use";
