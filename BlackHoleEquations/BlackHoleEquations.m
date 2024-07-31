(*=====================*)
(*  xCobaCalculations  *)
(*=====================*)

<<xAct`xPlain`;

Title@"Stefan's modified black hole equations";

PartIIIProject@"The purpose of this script is to give a concrete set of tools for extracting the reduced field equations in a static, spherical spacetime. With these tools, it is hoped that the remainder of the project time can be spent on analysis of the physics.";

SetOptions[$FrontEndSession,EvaluationCompletionAction->"ScrollToOutput"];

<<xAct`xTensor`;
<<xAct`xPerm`;
<<xAct`xTras`;
<<xAct`xCoba`;

(*-------------------------------------*)
(*  Definition of Riemannian manifold  *)
(*-------------------------------------*)

Get@FileNameJoin@{NotebookDirectory[],"DefManifold.m"};

(*------------------*)
(*  SchwarzschildLike spacetime  *)
(*------------------*)

Comment@"Define the Schwarzschild functions.";

DefScalarFunction[B1,PrintAs->"\[CapitalTheta]"];
DefScalarFunction[B2,PrintAs->"\[CapitalOmega]"];

Comment@"Define the chart for spherical polar coordinates.";

DefChart[SphericalPolar,M4,{0,1,2,3},{ct[],cr[],ctheta[],cphi[]},ChartColor->Purple];

PrintAs@ct^="\[ScriptT]";
PrintAs@cr^="\[ScriptR]";
PrintAs@ctheta^="\[Theta]";
PrintAs@cphi^="\[Phi]";

Comment@"Set the components of the metric to those of Schwarzschild.";

MatrixForm[MatrixSphericalPolar=DiagonalMatrix[{B1[cr[]]^2,-B2[cr[]]^2,-cr[]^2,-cr[]^2*Sin[ctheta[]]^2}]];
MatrixForm@MetricInBasis[G,-SphericalPolar,MatrixSphericalPolar];
MetricCompute[G,SphericalPolar,All,Verbose->False];

Comment@"Just to check the line element which we built here, let's have a look!";

Expr=G[{-a,-SphericalPolar},{-b,-SphericalPolar}];
DisplayExpression@Expr;
Expr//=ComponentArray;
DisplayExpression@MatrixForm@Expr;
Expr//=ToValues;
DisplayExpression@MatrixForm@Expr;

Comment@"Now how about the inverse metric?";

Expr=G[{a,SphericalPolar},{b,SphericalPolar}];
DisplayExpression@Expr;
Expr//=ComponentArray;
DisplayExpression@MatrixForm@Expr;
Expr//=ToValues;
DisplayExpression@MatrixForm@Expr;

(*-----------------------------*)
(*  Minisuperspace Lagrangian  *)
(*-----------------------------*)

Comment@"Define scalar functions of the radial coordinate scalar, and functions which will represent the ADM quantities.";

DefTensor[Phi[],M4,PrintAs->"\[Phi]"];
DefTensor[Psi[],M4,PrintAs->"\[Psi]"];

DefScalarFunction[Psis,PrintAs->"\[Psi]"];
DefScalarFunction[Phis,PrintAs->"\[Phi]"];

(*----------------------------------------------------------*)
(*  Some routines to simplify transfer to SchwarzschildLike coordinates  *)
(*----------------------------------------------------------*)

RawSchwarzschildLike[Expr_]:=Module[{SchwarzschildLikeExpr=Expr},
	SchwarzschildLikeExpr//=NoScalar;	
	SchwarzschildLikeExpr=SchwarzschildLikeExpr/.{
		Phi[]->Phis[cr[]],
		Psi[]->Psis[cr[]]
		};	
	SchwarzschildLikeExpr//=ToBasis[SphericalPolar];
	Print@SchwarzschildLikeExpr;
	SchwarzschildLikeExpr//=TraceBasisDummy;
	Print@SchwarzschildLikeExpr;
	SchwarzschildLikeExpr//=ToValues;
	Print@SchwarzschildLikeExpr;
	SchwarzschildLikeExpr//=ToValues;
	Print@SchwarzschildLikeExpr;

	(*----------------------------------------*)
	(*  Evaluate at the equator of the chart  *)
	(*----------------------------------------*)

	SchwarzschildLikeExpr=SchwarzschildLikeExpr/.{ctheta[]->Pi/2};
	Print@SchwarzschildLikeExpr;

	(*--------------------------------------------------------------------*)
	(*  Now that variations are done, we reveal the square root function  *)
	(*--------------------------------------------------------------------*)

	SchwarzschildLikeExpr=SchwarzschildLikeExpr/.{Sq->Sqrt};
	Print@SchwarzschildLikeExpr;
	SchwarzschildLikeExpr];

(*------------------------------------------------------------*)
(*  A wrapper for the above, which looks inside Scalar heads  *)
(*------------------------------------------------------------*)

ToSchwarzschildLike[Expr_]:=Module[{SchwarzschildLikeExpr=Expr},

	(*-----------------------------------------------------*)
	(*  Even inside scalars, impose the coordinate ansatz  *)
	(*-----------------------------------------------------*)

	SchwarzschildLikeExpr=SchwarzschildLikeExpr/.{xAct`xTensor`Scalar->RawSchwarzschildLike};
	Print@SchwarzschildLikeExpr;

	(*---------------------------------------*)
	(*  Then impose the ansatz on the whole  *)
	(*---------------------------------------*)

	SchwarzschildLikeExpr//=RawSchwarzschildLike;	
	Print@SchwarzschildLikeExpr;

	(*------------------------------------------------------------------*)
	(*  Assume we evaluate at a point in the chart with nonzero radius  *)
	(*------------------------------------------------------------------*)
(*
	SchwarzschildLikeExpr=SchwarzschildLikeExpr/(MPl*2*3*5*7*11*Sf[ct[]]*cr[])^10;
*)
	SchwarzschildLikeExpr//=Simplify;	
	Print@SchwarzschildLikeExpr;
(*
	SchwarzschildLikeExpr//=Normal;
	SchwarzschildLikeExpr//=Numerator;

	SchwarzschildLikeExpr//=Simplify;	
*)
	SchwarzschildLikeExpr];

PartIIIProject@"Something I omitted in the video tutorial above was a discussion of how to set up reduced radial-function components of a tensor field which is not the metric. There were no such tensors in the torsionful effective theory I gave to Oliver, but of course in MoND we have the unit-timelike vector field. How to work with this?";

Comment@"We wish to define a unit-timelike one-form.";

DefTensor[A[-a],M4,PrintAs->"\[ScriptCapitalA]"];

Expr=A[{-a,-SphericalPolar}];
DisplayExpression@Expr;
Expr//=ComponentArray;
Expr//=ToValues;
DisplayExpression@Expr;

Comment@"We know that this one-form has the property of being unit-timelike, and that the Killing algebra is spherical. So, we define a couple of scalar functions.";

DefScalarFunction[A1,PrintAs->"\[CapitalPhi]"];
DefScalarFunction[A2,PrintAs->"\[CapitalPsi]"];

AllComponentValues[A[{-a,-SphericalPolar}],{A1[cr[]],B2[cr[]]*Sqrt[A1[cr[]]^2/B1[cr[]]^2-1],0,0}];

Expr=A[{-a,-SphericalPolar}];
DisplayExpression@Expr;
Expr//=ComponentArray;
DisplayExpression@Expr;
Expr//=ToValues;
DisplayExpression@Expr;

Comment@"Now check that a simple expression is carried properly to the component form.";

Expr=CD[-a]@Phi[]*CD[-b]@Psi[]*G[a,b];
DisplayExpression@Expr;
Expr//=ToSchwarzschildLike;
DisplayExpression@Expr;

Comment@"Next something more advanced: we want to show that our one-form really does the job.";

Expr=A[-a]*A[a];
DisplayExpression@Expr;
Expr//=SeparateMetric[G];
Expr//=ScreenDollarIndices;
DisplayExpression@Expr;
Expr//=ToSchwarzschildLike;
DisplayExpression@Expr;

Comment@"Great. Now let's check that this can be used for more complex expressions.";

Expr=A[-a]*A[-b]*RicciCD[a,b];
DisplayExpression@Expr;
Expr//=SeparateMetric[G];
Expr//=ScreenDollarIndices;
DisplayExpression@Expr;
Expr//=ToSchwarzschildLike;
DisplayExpression@Expr;

PartIIIProject@"Brilliant, that computation took a fraction of a second. So with these tools all the field equations can be reduced to ODEs in the radius within no more than a couple of hours' tinkering: it is basically just data entry.";

Quit[];
