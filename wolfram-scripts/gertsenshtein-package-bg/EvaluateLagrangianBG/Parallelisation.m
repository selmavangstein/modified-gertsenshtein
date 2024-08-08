(*===================*)
(*  Parallelisation  *)
(*===================*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
<<xAct`xPlain`;
(*
Off@ValidateSymbol::used;
<<xAct`xTensor`;
<<xAct`xTras`;
<<xAct`xCoba`;
On@ValidateSymbol::used;
ParallelNeeds["xAct`xPlain`"];
ParallelNeeds["xAct`xTensor`"];
ParallelNeeds["xAct`xTras`"];
ParallelNeeds["xAct`xCoba`"];
*)
$DefInfoQ=False;
Unprotect@AutomaticRules;
Options[AutomaticRules]={Verbose->False};
Protect@AutomaticRules;

Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Parallelisation","PrepareFiles.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Parallelisation","QuietParallelSubmit.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Parallelisation","ProcessOperator.m"};
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Parallelisation","ApplyParallel.m"};

DistributeDefinitions@$ThisDirectory;
