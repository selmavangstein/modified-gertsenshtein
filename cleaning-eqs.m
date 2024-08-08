(* ::Package:: *)

fullySimplify[inputExpr_List]:=Module[{currExpr=inputExpr, simplerExpr},
	currExpr=applyMaxwellsEquation[currExpr];
	While[
		True,
		simplerExpr=simplifyExpr[currExpr];
		Print["done with a simplification"];
		If[
			simplerExpr===currExpr,Break[];
		];
		currExpr=simplerExpr;
	];
	Print["All done, hope you had a good experience"];
	currExpr
];

applyMaxwellsEquation[inputExpr_]:=Module[{outputExpr=inputExpr,maxwellRules},
	(*These are for sure not perfect, but they work until they don't*)
	(*Also, if Will does his job, this won't be necessarry I think*)
	maxwellRules={
	\!\(\*SuperscriptBox[\(\[ScriptCapitalE]x\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],x[],y[],z[]]->Derivative[2,0][\[ScriptB]][t[],z[]]+\!\(\*SuperscriptBox[\(\[ScriptCapitalE]z\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],x[],y[],z[]],
	\!\(\*SuperscriptBox[\(\[ScriptCapitalE]x\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],x[],y[],z[]]->\!\(\*SuperscriptBox[\(\[ScriptCapitalE]y\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],x[],y[],z[]],
	\!\(\*SuperscriptBox[\(\[ScriptCapitalE]z\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],x[],y[],z[]]->\!\(\*SuperscriptBox[\(\[ScriptCapitalE]y\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],x[],y[],z[]]};
	outputExpr=Expand[Simplify[outputExpr/.maxwellRules]];
	outputExpr
];

simplifyExpr[inputExpr_]:=Module[{outputExpr=inputExpr,eqs,allParts,variables,sols,zeroSolutions,extendedZeroSolutions},
	outputExpr//=DeleteDuplicates;
	outputExpr//=Expand;
	eqs=Thread[outputExpr==0];
	allParts=DeleteDuplicates@Flatten[outputExpr/. Plus|Times->List];
	variables=Select[allParts,!FreeQ[#,_Symbol]||!FreeQ[#,_[__]]&];
	variables=DeleteDuplicates[variables/.{x_^n_:>x}];
	constants=Select[variables, FreeQ[#,_[__]]&];
	variables=Select[variables,FreeQ[#,Derivative]&];
	Print["Looking for trivial solutions..."];
	outputExpr=applyTrivialSolutions[outputExpr,constants];
	Print["trivial solutions applied"];
	Print["solving the full thing..."];
	sols=Assuming[Thread[constants!=0],Solve[eqs,variables]];
	(*I have only seen this produce one solution, but if not it just picks one*)
	If[sols=!={},sols=sols[[1]]];
	Print["done solving, simplifying..."];
	zeroSolutions=Select[sols,#[[2]]===0&];
	extendedZeroSolutions=Flatten[generateDerivativeRules/@zeroSolutions]//DeleteDuplicates;
	outputExpr=Expand[Simplify[outputExpr/.extendedZeroSolutions]];
	outputExpr=DeleteCases[outputExpr,0];
	outputExpr
];

applyTrivialSolutions[inputExpr_,constants_List]:=Module[{currExpr=inputExpr,i,simplerExpr,trivialEqs,trivialSols,trivialExpr, trivialVariables,prodElements,zeroFuncs,zeroFuncsRule,dRules},
	simplerExpr=0;
	i=0;
	While[
		True,
		(*We find trivial solutions of from stuff=0*)
		trivialEqs=Thread[Select[currExpr, Head[#]=!=Plus&]==0];
		trivialSols=Assuming[Thread[constants!=0],Solve[trivialEqs]];
		(*I have only seen this produce one solution, but if not it just picks one*)
		trivialSols=trivialSols[[1]];
		(*Before applying the rules, we want to set the derivatives of the identified fcns to zero as well, so we don't lose that information*)
		dRules=Flatten[generateDerivativeRules/@trivialSols];
		(*We can now apply the solutions we have found*)
		simplerExpr=DeleteDuplicates[currExpr/.trivialSols/.dRules];
		If[currExpr===simplerExpr,Break[]];
		currExpr=simplerExpr;
		i+=1;
	];
	simplerExpr
]

generateDerivativeRules[rule_]:=Module[{var,args,derivatives,derivTuples},
	var=rule[[1]];
	If[
		Head[var]===Symbol,
		(*If the variable is a symbol,return nothing*)
		{},
		args=List@@var;
		(*Generate derivatives of the function*)
		If[
			Length[args]==2,
			derivTuples=Flatten[Table[{i,j},{i,0,4},{j,0,4}],1];
			derivTuples=Select[derivTuples,#[[1]]+#[[2]]<=4&];
		];
		(*I don't think we ever have E-stuff up to 4th derivative, so a little redundant, but its pretty quick*)
		If[
			Length[args]==4,
			derivTuples=Flatten[Table[{i,j,k,l},{i,0,4},{j,0,4},{k,0,4},{l,0,4}],3];
			derivTuples=Select[derivTuples,#[[1]]+#[[2]]+#[[3]]+#[[4]]<=4&];
		];
		derivatives=Table[Derivative[Sequence@@derivSeq][var[[0]]][Sequence@@args],{derivSeq,derivTuples}];
		(*Create rules to set these derivatives to zero*)
		Rule[#,0]&/@derivatives
	]
];


$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
resultsFileName = "cleanEquations";

Get@FileNameJoin[{$ThisDirectory,"wolfram-scripts/results/rr4-bg.mx"}]
tList=torsionExpr//Flatten;
eList=einsteinExpr//Flatten;
mList=maxwellExpr//Flatten;
expr=Join[tList,eList,mList];

finalExpr=fullySimplify[expr];
finalEqs=Thread[finalExpr==0];

DumpSave[FileNameJoin[{$ThisDirectory,resultsFileName<>".mx"}],finalEqs];
Print["Results saved in ", resultsFileName]
