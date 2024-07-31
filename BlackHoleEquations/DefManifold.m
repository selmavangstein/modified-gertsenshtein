(*===============*)
(*  DefManifold  *)
(*===============*)

DefManifold[M4,4,IndexRange[{a,z}]];
AddIndices[TangentM4,{a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,n1,m1,o1,p1,q1,r1,s1,t1,u1,v1,w1,x1,y1,z1}];

GSymb="\!\(\*OverscriptBox[\(\[ScriptG]\), \(\[Degree]\)]\)";
DefMetric[-1,G[-a,-b],CD,{";","\!\(\*OverscriptBox[\(\[Del]\), \(\[SmallCircle]\)]\)"},PrintAs->GSymb,SymCovDQ->True];

(*-----------------------------------------------------------------------*)
(*  Relabeling of the indices so that we can type Roman and look Greek!  *)
(*-----------------------------------------------------------------------*)

StandardIndices=ToString/@Alphabet[];
ExtendedIndices=ToExpression@(ToString@#<>"1")&/@StandardIndices;

GeoStandardIndicesSymb=(ToString@#)&/@Evaluate@((#[[2]])&/@{
	{a,"\[Alpha]"},
	{b,"\[Beta]"},
	{c,"\[Chi]"},
	{d,"\[Delta]"},
	{e,"\[Epsilon]"},
	{f,"\[Phi]"},
	{g,"\[Gamma]"},
	{h,"\[Eta]"},
	{i,"\[Iota]"},
	{j,"\[Theta]"},
	{k,"\[Kappa]"},
	{l,"\[Lambda]"},
	{m,"\[Mu]"},
	{n,"\[Nu]"},
	{o,"\[Omicron]"},
	{p,"\[Pi]"},
	{q,"\[Omega]"},
	{r,"\[Rho]"},
	{s,"\[Sigma]"},
	{t,"\[Tau]"},
	{u,"\[Upsilon]"},
	{v,"\[Psi]"},
	{w,"\[Omega]"},
	{x,"\[Xi]"},
	{y,"\[CurlyPhi]"},
	{z,"\[Zeta]"}});
GeoExtendedIndicesSymb=ToString@ToExpression@(ToString@#<>"'")&/@GeoStandardIndicesSymb;

(PrintAs@Evaluate@#1^=Evaluate@#2)&~MapThread~{ToExpression/@StandardIndices,GeoStandardIndicesSymb};
(PrintAs@Evaluate@#1^=Evaluate@#2)&~MapThread~{ToExpression/@ExtendedIndices,GeoExtendedIndicesSymb};

(*----------------------------------------------------*)
(*  Coupling constants used in the metrical analogue  *)
(*----------------------------------------------------*)

Comment@"Define a Planck mass.";

DefConstantSymbol[MPl,PrintAs->"\(\*SubscriptBox[\(\[ScriptCapitalM]\), \(Pl\)]\)"];
