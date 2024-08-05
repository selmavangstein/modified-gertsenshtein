(*======================*)
(*  SimplifyLagrangian  *)
(*======================*)

SimplifyLagrangian[InputExpr_]:=Module[{linearizedAction=InputExpr},
	linearizedAction//=Expand;
	linearizedAction=ApplyParallel[linearizedAction,{funcToH,funcToPertA,funcAtoF,funcToPertT}];
	linearizedAction = linearizedAction/.Sqrt[-Detmetric[]]->1;
	linearizedAction//=ToCanonical;
	linearizedAction//=ContractMetric;
	linearizedAction//=ScreenDollarIndices;
linearizedAction];
