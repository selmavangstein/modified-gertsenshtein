(*===================*)
(*  SuperChangeCovD  *)
(*===================*)

SuperChangeCovD[InputExpr_]:=ApplyParallel[InputExpr,{ChangeCovD[#,CDT,CD]&,ChristoffelToGradMetric,ToCanonical,ContractMetric,ScreenDollarIndices}];
