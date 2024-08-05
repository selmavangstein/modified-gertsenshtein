(*=================*)
(*  InitialExpand  *)
(*=================*)

InitialExpand[InputExpr_]:=Module[{lagrangian=InputExpr},
	lagrangian=lagrangian/.FtoA;
	lagrangian=ChangeCovD[lagrangian,CDT,CD];
	lagrangian//=ToCanonical;
	lagrangian//=ContractMetric;
	lagrangian//=ScreenDollarIndices;
	lagrangian//=CollectTensors;
	lagrangian=ChangeCurvature[lagrangian,CDT,CD];
	lagrangian//=ScreenDollarIndices;
	lagrangian//=ChristoffelToGradMetric;
	lagrangian//=ContractMetric;
	lagrangian//=ToCanonical;
	lagrangian//=ScreenDollarIndices;
lagrangian];
