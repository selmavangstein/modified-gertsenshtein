(* ::Package:: *)

(*================*)
(*  PrepareFiles  *)
(*================*)

PrepareFiles[VeryLongExpression_]:=Module[{Expr=VeryLongExpression,AllFileNames},
	toList[expr_]:=expr/. {Plus[a_,b__]:>List[a,b]};
	Expr=toList[Expr]//Flatten;	
	(*Expr=Expr/.Plus->List//Flatten*)
	AllFileNames=CreateFile[]~Table~{ii,1,Length@Expr};
	MapThread[(TheOperator=#2;#1~DumpSave~TheOperator)&,
	{AllFileNames,Expr}];
AllFileNames];

