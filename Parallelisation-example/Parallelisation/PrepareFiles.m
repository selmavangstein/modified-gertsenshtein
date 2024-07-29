(* ::Package:: *)

(*================*)
(*  PrepareFiles  *)
(*================*)

PrepareFiles[VeryLongExpression_]:=Module[{Expr=VeryLongExpression,AllFileNames},	
	Expr=List@@Expr;
	AllFileNames=CreateFile[]~Table~{ii,1,Length@Expr};
	MapThread[(TheOperator=#2;#1~DumpSave~TheOperator)&,
	{AllFileNames,Expr}];
AllFileNames];

