(* ::Package:: *)

(*=================*)
(*  ApplyParallel  *)
(*=================*)

ApplyParallel[VeryLongExpression_]:=Module[{AllFileNames,Results},
	AllFileNames=PrepareFiles@VeryLongExpression;
	Get@FileNameJoin@{$ThisDirectory,"Parallelisation","SaveBinaries.m"};
	Results=(QuietParallelSubmit@(	
		Off@(RuleDelayed::rhs);
		Get[FileNameJoin[{$ThisDirectory,#<>".mx"}],#]&/@{"xAct`xTensor`","xAct`xTensor`Private`","TangentM`","Global`"};
		Get@#;
		On@(RuleDelayed::rhs);
	ProcessOperator[]))&~Map~AllFileNames;
	Results//=WaitAll;
	DeleteFile/@AllFileNames;
	Results//=Total;
Results];
