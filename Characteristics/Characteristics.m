<<xAct`xPlain`

Comment@"Here is the beginning of a new tool."


HarmonicCharacteristics[InputSystem_]:=Module[{Expr=InputSystem},
	Expr=Flatten@(Variables/@Expr);
	Expr//=DeleteDuplicates;
	AllVariables=Expr;
	Block[{Derivative},Derivative[Orders___][Target___][Vars___]:=Target[Vars];Expr//=FullSimplify];
	Expr//=DeleteDuplicates;
	Couplings=Expr//Cases[_Symbol];
	Couplings//=Sort;
	Comment@"Here are the couplings in your system:";
	Couplings//DisplayExpression;
	Functions=Expr~Complement~Couplings;
	Functions//=Sort;
	Comment@"Here are the functions in your system:";
	Functions//DisplayExpression;

	Comment@"Here is the harmonic ansatz:";
	Format@V^=ToExpression@"\[ScriptV]";
	HarmonicAnsatz=(#->-D[#,t,z]*V^2)&/@Functions;

	AllDerivatives=AllVariables~Complement~(Functions~Join~Couplings);
	AllDerivatives//=Sort;
	AllDerivatives//DisplayExpression;
	RandomStrings=Table[Symbol@(ResourceFunction["RandomString"][10]),{Length@AllDerivatives}];
	ToInertDerivative=Thread[AllDerivatives->RandomStrings];
	FromInertDerivative=Thread[RandomStrings->AllDerivatives];

	Comment@"Here are the derivatives of the functions in your system:";
	FunctionHeads=Head/@Functions;
	Derivatives={
			Symbol@(ToString@#<>"\[ScriptT]"),
			Symbol@(ToString@#<>"\[ScriptZ]")
			}&/@FunctionHeads;
	Derivatives//=Flatten;
	Derivatives=(#@@{t,z})&/@Derivatives;
	Derivatives//DisplayExpression;
	SecondDerivatives={D[#,{t,1}],D[#,{z,1}]}&/@Derivatives;
	SecondDerivatives//=Flatten;
	DerivativeRules={
			D[#,{t,1}],
			D[#,{z,1}]
		}&/@Functions;
	DerivativeRules//=Flatten;
	SecondDerivativeRules={D[#,{t,1}],D[#,{z,1}]}&/@DerivativeRules;
	SecondDerivativeRules//=Flatten;
	DerivativeRules=Thread[SecondDerivatives->SecondDerivativeRules];

	TemporalMatrix=Table[Symbol@("T"<>ToString@i<>ToString@j),{i,Length@Derivatives},{j,Length@Derivatives}];
	TemporalMatrix=IdentityMatrix[Length@Derivatives];
	SpatialMatrix=Table[Symbol@("Z"<>ToString@i<>ToString@j),{i,Length@Derivatives},{j,Length@Derivatives}];
	MatrixVariables=Flatten@TemporalMatrix~Join~Flatten@SpatialMatrix;
	MatrixVariables=Flatten@SpatialMatrix;
	MatrixVariables//DisplayExpression;

	

	Expr=InputSystem/.ToInertDerivative;
	Expr=Expr/.HarmonicAnsatz;
	Expr=Expr/.FromInertDerivative;
	Expr//DisplayExpression;

	FormalSystem=TemporalMatrix.D[Derivatives,t]+SpatialMatrix.D[Derivatives,z];
	FormalSystem=FormalSystem/.DerivativeRules;
	FormalSystem//=FullSimplify;
	FormalSystem//MatrixForm//DisplayExpression;
	FormalSystem=MapThread[Equal,{FormalSystem,Flatten@({#,0}&/@Expr)}];
	FormalSystem//=FullSimplify;
	FormalSystem//DisplayExpression;
	Coeffs=CoefficientArrays[FormalSystem,DeleteDuplicates@SecondDerivativeRules][[2]];
	Coeffs//=Normal;
	Coeffs//=Flatten;
	Coeffs=(#==0)&/@Coeffs;
	Solution=First@Quiet@Solve[Coeffs,MatrixVariables];
	Solution//DisplayExpression;

	TemporalMatrix=TemporalMatrix/.Solution;
	SpatialMatrix=SpatialMatrix/.Solution;
	TemporalMatrix//MatrixForm//DisplayExpression;
	SpatialMatrix//MatrixForm//DisplayExpression;
	Format@Lamb^=ToExpression@"\[Lambda]";
	MyDet=Det@(TemporalMatrix-Lamb*SpatialMatrix);
	MyDet//=Normal;
	MyDet//=Numerator;
	MyDet//=FullSimplify;
	AllAssumptions=((#>0)&/@MatrixVariables)~Join~(Element[#,Reals]&/@(Couplings));
	AllAssumptions=((#>0)&/@(MatrixVariables~Join~Couplings));
	AllAssumptions=AllAssumptions~Join~{V>0};
	MyDet=Assuming[AllAssumptions,Solve[MyDet==0,Lamb]];
	MyDet=First/@MyDet;
	MyDet//=DeleteDuplicates;
	FullSols=MyDet;
	MyDet=(Lamb/.#)&/@MyDet;
	MyDet//DisplayExpression;
	MyDet1=Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]]&/@MyDet;
	MyDet1//DisplayExpression;
	MyDet0=Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,0}]]&/@MyDet;
	MyDet0//DisplayExpression;

	Reality=(Assuming[AllAssumptions,FullSimplify@Im[#]]===0)&/@MyDet1;
	Reality//=Position[#,False]&;
	Reality//DisplayExpression;
	(*
	MyDet//=Delete[#,Reality]&;
	FullSols//=Delete[#,Reality]&;
	MyDet1//=Delete[#,Reality]&;
	MyDet0//=Delete[#,Reality]&;
	*)

	Comment@"Here are the solutions for the characteristic equation:";
	FindNullVector[LambSol_,Direct_]:=Module[{NullVector},
		Block[{Mag=1,V=-0.1},
			Expr=(TemporalMatrix-Lamb*SpatialMatrix)/.LambSol;
			(*Expr=Assuming[AllAssumptions,FullSimplify@Expr];*)
			(*Expr=Map[(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,2}]])&,Expr,{2}];*)
			(*
			Expr//=Eigensystem[#,Quartics->True]&;
			EigenValues=First@Expr;
			EigenVectors=Last@Expr;
			*)
			(*EigenValues=Map[(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]])&,EigenValues];
			EigenVectors=Map[(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]])&,EigenVectors,{2}];*)
			(*
			EigenValues//DisplayExpression;
			EigenVectors//DisplayExpression;
			*)
			(**)
			Expr=Assuming[AllAssumptions,NullSpace[Transpose@Expr]];
			(*Expr=(Assuming[AllAssumptions,Sqrt[2]*Normalize[#]])&/@Expr;*)
			Expr=(Assuming[AllAssumptions,Normalize[#]])&/@Expr;
			(*Expr=Map[(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]])&,Expr,{2}];*)
			Comment@"Solution array:";
			NewVar=(Transpose[#].TemporalMatrix)&/@Expr;
			(*NewVar=Map[(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]])&,NewVar,{2}];*)
			NewVar//DisplayExpression;
			Comment@"End of sol array.";
			firVar=Table[If[i==j,{Direct*Lamb,Direct},{0,0}],{i,Length@Functions},{j,Length@Functions}];
			firVar=Flatten/@firVar;
			firVar=(#/.LambSol)&/@firVar;
			(*
			firVar=(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]]&/@#)&/@firVar;
			*)
			(*firVar=(Assuming[AllAssumptions,Sqrt[2]*Normalize[#]])&/@firVar;*)
			firVar=(Assuming[AllAssumptions,Normalize[#]])&/@firVar;
			(*firVar=Map[(Assuming[AllAssumptions,Normal@FullSimplify@Series[#,{V,0,1}]])&,firVar,{2}];*)
		];
		(**)
	{NewVar,firVar}];	
	FinalSols=MapThread[FindNullVector,{FullSols,MyDet0}];	
	(*FinalSols=FinalSols~Join~(FinalSols/.{V->-V});*)
	Comment@"Solution array:";
	DisplayExpression/@FinalSols;
	Comment@"Now for pairing.";
];


Format@Phi^=ToExpression@"\[Phi]";
Format@Psi^=ToExpression@"\[Psi]";
Format@M^=ToExpression@"\[ScriptM]";
InputSystem={D[Phi[t,z],{t,2}]-D[Phi[t,z],{z,2}]+M^2*Phi[t,z],D[Psi[t,z],{t,2}]-D[Psi[t,z],{z,2}]+M^2*Psi[t,z]};

Format@H^=ToExpression@"\[\ScriptH]";
Format@B^=ToExpression@"\[ScriptB]";
Format@Mag^=ToExpression@"\[ScriptCapitalB]";
InputSystem={D[B[t,z],{t,2}]-D[B[t,z],{z,2}]-Mag*D[H[t,z],{z,2}],D[H[t,z],{t,2}]-D[H[t,z],{z,2}]+4*Mag*B[t,z]};

DisplayExpression/@InputSystem;

InputSystem//HarmonicCharacteristics;

Supercomment@"This is the end of the script."

Quit[];

Matr={
{1,l,0,B*l},
{l,1,0,0},
{-y*l,x,1,l},
{0,0,l,1}};

(*
Matr={
{1,0,-1,0},
{0,1,0,-1},
{1,0,1,0},
{0,1,0,1}};
*)

Matr//=Det;
Matr//=FullSimplify;
Matr//DisplayExpression;
Matr=Matr/.{x->4B/w^2-y};
Matr//=FullSimplify;
Sols=Solve[Matr==0,l];
Sols//DisplayExpression;
Sols=(l/.#)&/@Sols;
Sols//DisplayExpression;
Sols=Assuming[w>0,FullSimplify@Series[#,{B,0,2}]]&/@Sols;
Sols//DisplayExpression;


Quit[]
