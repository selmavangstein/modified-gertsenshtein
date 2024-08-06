(*========================*)
(*  DeleteFirstOrderPart  *)
(*========================*)

DefConstantSymbol[PerturbativeParameter, PrintAs->"\[Epsilon]"];
ToOrderH = MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},MetricOn->All,ContractMetrics->True];
ToOrderCDH=MakeRule[{CD[-c]@H[-a,-b],PerturbativeParameter*CD[-c]@H[-a,-b]},MetricOn->All,ContractMetrics->True];
DeleteFirstOrderPart[InputExpr_]:=Module[{OutputLagrangian=InputExpr,SecondOrderPart,FirstOrderPart,NullOrderPart},
	OutputLagrangian=OutputLagrangian/.ToOrderCDH/.ToOrderH;
	SecondOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,2}]&;
	SecondOrderPart//=Normal;
	SecondOrderPart=SecondOrderPart/.PerturbativeParameter->1;
	SecondOrderPart//=ToCanonical;
	SecondOrderPart//=ContractMetric;
	SecondOrderPart//=ScreenDollarIndices;
	SecondOrderPart//=CollectTensors;
	
	FirstOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,1}]&;
	FirstOrderPart//=Normal;
	FirstOrderPart=FirstOrderPart/.PerturbativeParameter->1;
	FirstOrderPart//=ToCanonical;
	FirstOrderPart//=ContractMetric;
	FirstOrderPart//=ScreenDollarIndices;
	FirstOrderPart//=CollectTensors;
	
	NullOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,0}]&;
	NullOrderPart//=Normal;
	NullOrderPart=NullOrderPart/.PerturbativeParameter->1;
	NullOrderPart//=ToCanonical;
	NullOrderPart//=ContractMetric;
	NullOrderPart//=ScreenDollarIndices;
	NullOrderPart//=CollectTensors;
	
	OutputLagrangian=SecondOrderPart-(FirstOrderPart-NullOrderPart);
	OutputLagrangian//=ToCanonical;
	OutputLagrangian//=ContractMetric;
	OutputLagrangian//=ScreenDollarIndices;
	OutputLagrangian//=CollectTensors;
OutputLagrangian];
