(*==================*)
(*  ParallelExpand  *)
(*==================*)

ParallelExpand[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=ApplyParallel[#,{Expand,ToCanonical,ContractMetric,ScreenDollarIndices}]&;
	Expr//=ToCanonical;
	Expr//=ContractMetric;
	Expr//=ScreenDollarIndices;
Expr];
