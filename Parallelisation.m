(* ::Package:: *)

(*===================*)
(*  Parallelisation  *)
(*===================*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
<<xAct`xPlain`;

Title@"Scripts for parallel computations";

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
Unprotect@AutomaticRules;
Options[AutomaticRules]={Verbose->False};
Protect@AutomaticRules;

Get@FileNameJoin@{$ThisDirectory,"Parallelisation","PrepareFiles.m"};
Get@FileNameJoin@{$ThisDirectory,"Parallelisation","QuietParallelSubmit.m"};
Get@FileNameJoin@{$ThisDirectory,"Parallelisation","ProcessOperator.m"};
Get@FileNameJoin@{$ThisDirectory,"Parallelisation","ApplyParallel.m"};

DistributeDefinitions@$ThisDirectory;
Comment@"If there were no errors above, then the kernel should be ready for parallel computing.";

(*
Comment@"We have a very long expression consisting of many terms. We wish to break this expression into a list of terms, rather than a sum. Each term in the list will be exported to a binary file. Each file will then be reimported by a parallel kernel, processed, and then the results recombined. As a placeholder, our \"very long expression\" will just be the English alphabet.";
VeryLongExpression=Total@Alphabet[];
DisplayExpression[VeryLongExpression,EqnLabel->"VeryLongExpression"];

Comment@"Now we try to process this expression in parallel.";
VeryLongExpression//=ApplyParallel;
DisplayExpression[VeryLongExpression,EqnLabel->"VeryLongExpression"];
*)
