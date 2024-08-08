(*========*)
(*  Test  *)
(*========*)

<<xAct`xPlain`;

Comment@"We define some things...";
Code[
	MyNewFunc[x_]:=x/0;
	MyNewerFunc[x_]:=x[[3]];
	$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];

	$MagicString="is longer than depth of object";

	FunctionCheck::culpret="Function `1` will now act on the expression `2`.";
	FunctionCheck[InputFunc_][InputExpr_]:=(Message[FunctionCheck::culpret,InputFunc,InputExpr];
							InputExpr//InputFunc);
	Discriminator[InputFileName_]:=Module[{FileContentData},
		FileContentData=Import[InputFileName,"Text"];
		FileContentData=StringDelete[FileContentData,"\n"];
		FileContentData//=StringReplace[#,{"  "->" "}]&;
		FileContentData//=StringReplace[#,{"  "->" "}]&;
		FileContentData//=StringReplace[#,{"  "->" "}]&;
		(FileContentData~StringContainsQ~$MagicString)~If~(
			WriteString[FileNameJoin@{$ThisDirectory,"results","BadEvaluation"<>ToString@Hash@FileContentData},FileContentData];
		);
	];

	ReadErrors[InputFileName_]:=Module[{FileContentData},
		Comment@{"We read the error file \"",InputFileName,"\"..."};
		FileContentData=Import[InputFileName,"Text"];
		FileContentData=StringSplit[FileContentData,"BoxForm`AbsoluteRawBoxes"];
		FileContentData=("BoxForm`AbsoluteRawBoxes"<>#)&/@FileContentData;
		FileContentData=(#~ToExpression~TraditionalForm)&/@FileContentData;
		Print/@FileContentData;
	];

	MultipleSteps[InputExpr_]:=Module[{Expr=InputExpr,NewStream=OpenWrite[],FileContentData},
		Block[{$Messages=Append[$Messages,NewStream]},
		(*Block[{$Messages=NewStream,ApplyTo[x_,f_]:=(x//=FunctionCheck@f)},*)
			Expr//=FunctionCheck@MyNewerFunc;
			Expr=1;
			Expr//=FunctionCheck@MyNewFunc;
		];
		Discriminator[NewStream[[1]]];
		Close@NewStream;
	Expr];
];


Comment@"We test the function...";

Code[
	$ErrorFiles=FileNames["results/BadEvaluation*"];
	DeleteFile/@$ErrorFiles;
	MultipleSteps[1];
	$ErrorFiles=FileNames["results/BadEvaluation*"];
	ReadErrors/@$ErrorFiles;
];
