(*===============*)
(*  ExampleCode  *)
(*===============*)

(*How to do it*)

MakeRule[{CD[-a]@A[-b],
	(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},
	MetricOn->All,ContractMetrics->True];

FieldEqOfMetric(H,PertA,F)=0
FieldEqOfA(H,PertA,F)=0
