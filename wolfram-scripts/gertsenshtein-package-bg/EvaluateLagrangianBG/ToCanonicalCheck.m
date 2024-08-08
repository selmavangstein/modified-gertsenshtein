(*====================*)
(*  ToCanonical  *)
(*====================*)

(*$MagicString="culpret";*)
$MagicString="etected";

FunctionCheck::culpret="Function `1` will now act on the expression `2`.";
FunctionCheck[InputFunc_][InputExpr_]:=(InputExpr//InputFunc);
(*FunctionCheck[InputFunc_][InputExpr_]:=(Print@InputFunc;Print@InputExpr;InputExpr//InputFunc);*)
(*FunctionCheck[InputFunc_][InputExpr_]:=(Message[FunctionCheck::culpret,InputFunc,InputExpr];
						InputExpr//InputFunc);*)
Discriminator[InputFileName_]:=Module[{FileContentData},
	FileContentData=Import[InputFileName,"Text"];
	FileContentData=StringDelete[FileContentData,"\n"];
	FileContentData//=StringReplace[#,{"  "->" "}]&;
	FileContentData//=StringReplace[#,{"  "->" "}]&;
	FileContentData//=StringReplace[#,{"  "->" "}]&;
	(FileContentData~StringContainsQ~$MagicString)~If~(
		(*WriteString[FileNameJoin@{$ThisDirectory,"results","BadEvaluation"<>ToString@Hash@FileContentData},FileContentData];*)
		FileContentString=FileContentData;
		DumpSave[FileNameJoin@{$ThisDirectory,"results","BadEvaluation"<>ToString@Hash@FileContentData<>".mx"},{FileContentString}];
	);
];

ReadErrors[InputFileName_]:=Module[{FileContentData},
	Comment@{"We read the error file \"",InputFileName,"\"..."};
	(*FileContentData=Import[InputFileName,"Text"];*)
	Get@InputFileName;
	FileContentData=FileContentString;
	FileContentData=StringSplit[FileContentData,"BoxForm`AbsoluteRawBoxes"];
	FileContentData=("BoxForm`AbsoluteRawBoxes"<>#)&/@FileContentData;
	FileContentData=(#~ToExpression~TraditionalForm)&/@FileContentData;
	Print/@FileContentData;
];

MultipleStepsTorsion[InputExpr_]:=Module[{Expr=InputExpr,NewStream=OpenWrite[],FileContentData},
	Block[{$Messages=Append[$Messages,NewStream]},
		Off@General::stop;
		Expr//=FunctionCheck@funcChristCartZero;
		Expr//=FunctionCheck@ToBasis[cartesian];
		Expr//=FunctionCheck@funcChristCartZero;
		Expr//=FunctionCheck@ToBasis[cartesian];
		Expr//=FunctionCheck@funcChristCartZero;
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@ComponentArray;
		Expr//=FunctionCheck@ToValues;
		Expr//=FunctionCheck@ToValues;
		Expr//=FunctionCheck@ToValues;
		Expr//=FunctionCheck@ToCanonical;
		Expr//=FunctionCheck@SeparateMetric[metric];
		Expr//=FunctionCheck@ToBasis[cartesian];
		Expr//=FunctionCheck@ToBasis[cartesian];
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@ToCanonical;
		Expr//=FunctionCheck@ToValues;
		On@General::stop;
	];
	Discriminator[NewStream[[1]]];
	Close@NewStream;
Expr];

MultipleStepsEinstein[InputExpr_]:=Module[{Expr=InputExpr,NewStream=OpenWrite[],FileContentData},
	Block[{$Messages=Append[$Messages,NewStream]},
		Off@General::stop;
		Expr//=FunctionCheck@funcChristCartZero;
		Expr//=FunctionCheck@ToBasis[cartesian];
		Expr//=FunctionCheck@funcChristCartZero;
		Expr//=FunctionCheck@ToBasis[cartesian];
		Expr//=FunctionCheck@funcChristCartZero;
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@TraceBasisDummy;
		Expr//=FunctionCheck@ComponentArray;
		Expr//=FunctionCheck@ToValues;
		Expr//=FunctionCheck@ToValues;
		Expr//=FunctionCheck@ToValues;
		Expr//=FunctionCheck@ToCanonical;
		On@General::stop;
	];
	Discriminator[NewStream[[1]]];
	Close@NewStream;
Expr];
