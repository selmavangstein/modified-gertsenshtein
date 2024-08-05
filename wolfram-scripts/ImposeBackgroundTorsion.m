(*===========================*)
(*  ImposeBackgroundTorsion  *)
(*===========================*)

ImposeBackgroundTorsion[InputExpr_]:=Module[{linearizedAction=InputExpr},
	linearizedAction//=Expand;
	linearizedAction=ApplyParallel[linearizedAction,{funcTtoVec,ToCanonical,ContractMetric}];
	linearizedAction//DisplayExpression;
	linearizedAction//=funcLorentz;
	linearizedAction//=funcCD;
	linearizedAction//=BreakScalars;
	linearizedAction//=Expand;
linearizedAction];
