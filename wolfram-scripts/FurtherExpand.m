(*=================*)
(*  FurtherExpand  *)
(*=================*)

FurtherExpand[InputExpr_]:=Module[{lagrangian=InputExpr},
	Comment@"FurtherExpand";
	lagrangian//=NoScalar;
	lagrangian//=ContractMetric;
	lagrangian//=ToCanonical;
	lagrangian//=ScreenDollarIndices;
	lagrangian//=NoScalar;
	lagrangian//=ContractMetric;
	lagrangian//=ToCanonical;
	lagrangian//=ScreenDollarIndices;
	linearizedAction=PerturbBackground[lagrangian,2, BackgroundSolution->bgRules];
	linearizedAction//=ExpandPerturbation;
	linearizedAction//=ToBackground;
	linearizedAction//=CollectTensors;
linearizedAction];
