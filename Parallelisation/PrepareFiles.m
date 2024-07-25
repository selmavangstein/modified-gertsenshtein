(*================*)
(*  PrepareFiles  *)
(*================*)

PrepareFiles[VeryLongExpression_]:=Module[{Expr=VeryLongExpression,AllFileNames},	
	Expr=Expr/.{Plus->List};
	AllFileNames=Table["Operator"<>ToString@ii,{ii,1,Length@Expr}];
	MapThread[(TheOperator=#2;DumpSave[FileNameJoin[{$ThisDirectory,#1<>".mx"}],TheOperator])&,
	{AllFileNames,Expr}];
AllFileNames];

