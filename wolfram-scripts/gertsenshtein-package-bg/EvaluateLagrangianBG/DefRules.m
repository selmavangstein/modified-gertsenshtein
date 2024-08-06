(*============*)
(*  DefRules  *)
(*============*)

epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True];
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,-cartesian},{1,-cartesian},{2,-cartesian},{3,-cartesian}],1},MetricOn->All,ContractMetrics->True];

FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True];
AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];
funcAtoF[Expr_]:=Module[{expr=Expr},
	expr=expr/.AtoF;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
Perturbationmetric[LI[n_],___]/;n>1:=0;
toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True];
funcToH[Expr_]:=Module[{expr=Expr},
	expr=expr/.toH;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M];
perturbationA[LI[n_],___]/;n>1:=0;
topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True];
funcToPertA[Expr_]:=Module[{expr=Expr},
	expr=expr/.topertA;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
DefTensorPerturbation[perturbationT[LI[order],a,-b,-c],TorsionCDT[a,-b,-c],M];
perturbationT[LI[n_],___]/;n>1:=0;
topertT= MakeRule[{perturbationT[LI[1],a,-b,-c],pertT[a,-b,-c]},MetricOn-> All, ContractMetrics-> True];
funcToPertT[Expr_]:=Module[{expr=Expr},
	expr=expr/.topertT;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True];
pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},
	MetricOn->All,ContractMetrics->True];
funcPertAtoF[Expr_]:=Module[{expr=Expr},
	expr=expr/.pertAtoF;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
TtoVec=MakeRule[{TorsionCDT[a,-b,-c],epsilonmetric[a,-b,-c,-d]Q[d]},MetricOn->All,ContractMetrics->True];
funcTtoVec[Expr_]:=Module[{expr=Expr},
	expr=expr/.TtoVec;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];
lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];
commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];
funcLorentz[Expr_]:=Module[{expr=Expr},
	expr=expr/.lorentz;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
funcCD[Expr_]:=Module[{expr=Expr},
	expr=expr/.commuteCD;
	expr//=ToCanonical;
	expr//=ContractMetric;
	expr//=ScreenDollarIndices;
expr];
