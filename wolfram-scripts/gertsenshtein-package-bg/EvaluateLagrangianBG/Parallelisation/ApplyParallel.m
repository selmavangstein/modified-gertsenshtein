(* ::Package:: *)

(*=================*)
(*  ApplyParallel  *)
(*=================*)

ApplyParallel[VeryLongExpression_, commands_List]:=Module[{AllFileNames,Results},
	AllFileNames=PrepareFiles@VeryLongExpression;
	Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg","EvaluateLagrangianBG","Parallelisation","SaveBinaries.m"};
	Results=(QuietParallelSubmit@(	
		Off@(RuleDelayed::rhs);
		Get[FileNameJoin[{$ThisDirectory,#<>".mx"}],#]&/@{"xAct`xTensor`","xAct`xTensor`Private`","xAct`xCoba`","xAct`xCoba`Private`","TangentM`","Global`"};
		Get@#;
		On@(RuleDelayed::rhs);
	ProcessOperator[commands]))&~Map~AllFileNames;
	Results//=WaitAll;
	DeleteFile/@AllFileNames;
	Results//=Total;
Results];
