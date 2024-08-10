(*===============*)
(*  ExampleCode  *)
(*===============*)

(*How to do it*)

(*
	Tasks....

- Write down your Lagrangian (or tweak it for a new model...)
- Expand F
- Use xPert to get everything to 2nd order (i.e. both metric perturbation and electromagnetic gauge potential)
- Use your rules to replace triangle symbols with actual tensors
- Try to get rid of explicit background electromagnetic gauge potential, replacing with background Faraday (!)
- VarD for two sets of fields equations
- xCoba (good luck!)

*)

MakeRule[{CD[-a]@A[-b],
	(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},
	MetricOn->All,ContractMetrics->True];

(*
FieldEqOfMetric(H,PertA,F)=0
FieldEqOfA(H,PertA,F)=0
*)

DefConstantSymbol[PerturbativeParameter,PrintAs->"\[Epsilon]"];
ToOrder=MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},...]~Join~MakeRule[{PertA[-a],...},...];

ExtractQuadraticPart[InputExpr_]:=Module[{LinearLagrangian=InputExpr,SecondOrderPart,FirstOrderPart},

	LinearLagrangian=LinearLagrangian/.ToOrder;

	SecondOrderPart=LinearLagrangian//Series[#,{PerturbativeParameter,0,2}]&;
	SecondOrderPart//=Normal;
	SecondOrderPart=SecondOrderPart/.PerturbativeParameter->1;
	SecondOrderPart//=xAct`PSALTer`Private`ToNewCanonical;

	FirstOrderPart=LinearLagrangian//Series[#,{PerturbativeParameter,0,1}]&;
	FirstOrderPart//=Normal;
	FirstOrderPart=FirstOrderPart/.PerturbativeParameter->1;
	FirstOrderPart//=xAct`PSALTer`Private`ToNewCanonical;

	LinearLagrangian=SecondOrderPart-FirstOrderPart;
	LinearLagrangian//=ToCanonical;
	LinearLagrangian//=ContractMetric;
	LinearLagrangian//=ScreenDollarIndices;
	LinearLagrangian//=CollectTensors;
LinearLagrangian];
