(* ::Package:: *)

(* ::Title:: *)
(*Modified Gertsenshtein Effect*)


(* ::Text:: *)
(*First Riemann permutation terms in lagrangian, no background torsion*)


(* ::Section::Closed:: *)
(*Setup*)


(* ::Subsection:: *)
(*Loading packages*)


(* ::Text:: *)
(*xTras should load all the packages we need*)


(* ::Input::Initialization:: *)
<<xAct`xTras`


(* ::Text:: *)
(*Some convenient settings*)


(* ::Input::Initialization:: *)
$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;

SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];


(* ::Subsection:: *)
(*Manifold, basis, metric*)


(* ::Text:: *)
(*Defining a manifold, a metric, a small perturbation h, and making the background flat with SymmetricSpaceRules*)


(* ::Input::Initialization:: *)
DefManifold[M,4,IndexRange[{a,s}]]; (*might fix formatting: type roman, look greek. See DefManifold notebook*)
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];
DefCovD[CDT[-a],Torsion->True, SymbolOfCovD->{"#","D"},FromMetric->metric];


(* ::Input::Initialization:: *)
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];


(* ::Text:: *)
(*Place metric on that chart:*)


(* ::Input::Initialization:: *)
MetricInBasis[metric, -cartesian,{1,-1,-1,-1}]


(* ::Input::Initialization:: *)
MetricCompute[metric,cartesian, All]


(* ::Text:: *)
(*Setting a flat background*)


(* ::Input::Initialization:: *)
bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];


(* ::Subsection::Initialization:: *)
(*(*Defining our tensors*)*)


(* ::Input::Initialization:: *)
DefTensor[A[-a], M]


(* ::Input::Initialization:: *)
DefTensor[F[-a,-b],M,Antisymmetric[{-a,-b}]]


(* ::Text::Initialization:: *)
(*(*Perturbations:*)*)


(* ::Input::Initialization:: *)
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];


(* ::Input::Initialization:: *)
DefTensor[pertA[-a],M,PrintAs->"\[ScriptCapitalA]"]


(* ::Input::Initialization:: *)
DefTensor[pertF[-a,-b],M, Antisymmetric[{-a,-b}], PrintAs->"\[ScriptCapitalF]"]


(* ::Input::Initialization:: *)
DefTensor[pertT[a,-b,-c],M, Antisymmetric[{-b,-c}], PrintAs->"\[ScriptCapitalT]"]


(* ::Section::Closed:: *)
(*Rules*)


(* ::Subsection::Initialization:: *)
(*(*Going between F and A*)*)


(* ::Input::Initialization:: *)
FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];


(* ::Subsection::Initialization:: *)
(*(*Perturbation rules*)*)


(* ::Text::Initialization:: *)
(*(*Setting h as perturbation on the metric*)*)


(* ::Input::Initialization:: *)
Perturbationmetric[LI[n_],___]/;n>1:=0
toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True]


(* ::Text::Initialization:: *)
(*(*Defining a perturbation on A, and setting it to pertA*)*)


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M];
perturbationA[LI[n_],___]/;n>1:=0;
topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True]


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationT[LI[order],a,-b,-c],TorsionCDT[a,-b,-c],M];
perturbationT[LI[n_],___]/;n>1:=0;
topertT= MakeRule[{perturbationT[LI[1],a,-b,-c],pertT[a,-b,-c]},MetricOn-> All, ContractMetrics-> True];


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationR[LI[order],a,-b,-c],RicciScalarCDT[],M];
perturbationR[LI[n_],___]/;n>1:=0;


(* ::Subsection::Initialization:: *)
(*(*Between \[ScriptCapitalA] and \[ScriptCapitalF]*)*)


(* ::Text::Initialization:: *)
(*(*Now we connect the two perts:*)*)


(* ::Input::Initialization:: *)
pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},MetricOn->All,ContractMetrics->True];


(* ::Section::Initialization:: *)
(*(*Defining and expanding Lagrangian*)*)


(* ::Subsection::Initialization:: *)
(*(*Defining Lagrangian*)*)


(* ::Input::Initialization:: *)
DefConstantSymbol[\[Kappa]]
DefConstantSymbol[\[Lambda]]


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=(Sqrt[-Detmetric[]](RicciScalarCDT[]+\[Lambda] RiemannCDT[-a,-b,-c,-d]RiemannCDT[a,b,c,d]+\[Kappa] F[a,b]F[-a,-b]))/.FtoA


(* ::Subsection::Initialization:: *)
(*(*Writing it in terms of torsion*)*)


(* ::Text::Initialization:: *)
(*(*I have these two steps to get T and R specifically out there, and can still work easily with CD and the metric*)*)


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=ChangeCovD[\[ScriptCapitalL],CDT,CD](*redundant step - no CDT's in lagrangian*)


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=ChangeCurvature[\[ScriptCapitalL],CDT,CD]//ChristoffelToGradMetric//ContractMetric//ToCanonical


(* ::Subsection::Initialization:: *)
(*(*Expanding Lagrangian*)*)


(* ::Input::Initialization:: *)
pert\[ScriptCapitalL]=PerturbBackground[  \[ScriptCapitalL],2,
BackgroundSolution->bgRules]//ExpandPerturbation//ToBackground//CollectTensors;


(* ::Input::Initialization:: *)
pert\[ScriptCapitalL]=pert\[ScriptCapitalL]/.toH/.topertA/.topertT/.AtoF//ToCanonical;
linearizedAction = pert\[ScriptCapitalL]/.Sqrt[-Detmetric[]]->1


(* ::Subsection::Initialization:: *)
(*(*Defining traceless \[ScriptH] and Lorentz gauge*)*)


(* ::Input::Initialization:: *)
H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];
commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
linearizedAction=linearizedAction/.commuteCD/.lorentz//ToCanonical//CollectTensors


(* ::Section::Initialization::Closed:: *)
(*(*Field Equations*)*)


(* ::Subsection::Initialization:: *)
(*(*With respect to H:*)*)


(* ::Input::Initialization:: *)
einsteinField=VarD[H[a,b],CDT][linearizedAction]//ContractMetric//ToCanonical


(* ::Text::Initialization:: *)
(*(*Now we convert any CDT's in there, which will give new torsion tensors, and then we simplify further*)*)


(* ::Input::Initialization:: *)
einsteinField=ChangeCovD[einsteinField,CDT,CD]//ChristoffelToGradMetric;
einsteinField=einsteinField/.TorsionCDT->Zero/.pertAtoF//ToCanonical//ContractMetric;
einsteinField=einsteinField/.lorentz/.commuteCD


(* ::Subsection::Initialization::Closed:: *)
(*(*With respect to pertA:*)*)


(* ::Input::Initialization:: *)
maxwellField=VarD[pertA[a],CDT][linearizedAction]/.TorsionCDT->Zero//ContractMetric//ToCanonical


(* ::Text::Initialization:: *)
(*(*We do a similar simplification to above:*)*)


(* ::Input::Initialization:: *)
maxwellField=ChangeCovD[maxwellField,CDT,CD]//ChristoffelToGradMetric;
maxwellField=maxwellField/.TorsionCDT->Zero/.pertAtoF//ToCanonical//ContractMetric;
maxwellField=maxwellField/.lorentz/.commuteCD


(* ::Subsection::Initialization::Closed:: *)
(*(*With respect to torsion tensor:*)*)


(* ::Input::Initialization:: *)
torsionField=VarD[pertT[a,-b,-c],CDT][linearizedAction]/.TorsionCDT->Zero//ContractMetric//ToCanonical;
torsionField=ChangeCovD[torsionField,CDT,CD]//ChristoffelToGradMetric;
torsionField=torsionField/.TorsionCDT->Zero/.lorentz/.commuteCD//ToCanonical//ContractMetric


(* ::Subsection::Initialization:: *)
(*(*Removing terms with \[ScriptH] but keeping CD[\[ScriptH]]*)*)


(* ::Text::Initialization:: *)
(*(*Note that this only works with equations linear in h*)*)


(* ::Input::Initialization:: *)
DefConstantSymbol[PerturbativeParameter, PrintAs->"\[Epsilon]"]


(* ::Input::Initialization:: *)
ToOrderH = MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},MetricOn->All,ContractMetrics->True]


(* ::Text::Initialization:: *)
(*(*Since we have converted all derivatives from CDT to CD, this should still work:*)*)


(* ::Input::Initialization:: *)
ToOrderCDH=MakeRule[{CD[-c]@H[-a,-b],PerturbativeParameter*CD[-c]@H[-a,-b]},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
DeleteFirstOrderPart[InputExpr_]:=Module[{OutputLagrangian=InputExpr,SecondOrderPart,FirstOrderPart,NullOrderPart},OutputLagrangian=OutputLagrangian/.ToOrderCDH/.ToOrderH;
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


(* ::Input::Initialization:: *)
einsteinField=DeleteFirstOrderPart[einsteinField]


(* ::Input::Initialization:: *)
maxwellField=DeleteFirstOrderPart[maxwellField]


(* ::Input::Initialization:: *)
torsioneqField=DeleteFirstOrderPart[torsioneqField]


(* ::Section::Initialization:: *)
(*(*Components - enter xCoba*)*)


(* ::Text::Initialization:: *)
(*(*We now want to define the components of our tensors.*)*)


(* ::Subsection::Initialization::Closed:: *)
(*(*\[Epsilon]:*)*)


(* ::Input::Initialization:: *)
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,-cartesian},{1,-cartesian},{2,-cartesian},{3,-cartesian}],1},MetricOn->All,ContractMetrics->True]


(* ::Subsection::Initialization::Closed:: *)
(*(*\[ScriptH]:*)*)


(* ::Text::Initialization:: *)
(*(*Setting them all to zero first:*)*)


(* ::Input::Initialization:: *)
zerovalues=Table[0,{i,0,3},{j,0,3}]


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues]


(* ::Text::Initialization:: *)
(*(*Setting the non-zero components*)*)


(* ::Input::Initialization:: *)
DefScalarFunction[{h1,h2}]


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{1,-cartesian}],h1[t[],z[]]]


(* ::Input::Initialization:: *)
ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]]


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{2,-cartesian}],h2[t[],z[]]]


(* ::Text::Initialization:: *)
(*(*With indices up*)*)


(* ::Input::Initialization:: *)
ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]];


(* ::Subsection::Initialization::Closed:: *)
(*(*F:*)*)


(* ::Text::Initialization:: *)
(*(*Setting all to zero:*)*)


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray@F[-{a,cartesian},-{b,-cartesian}],zerovalues]


(* ::Text::Initialization:: *)
(*(*Giving it a constant B-field in the x-direction:*)*)


(* ::Input::Initialization:: *)
DefConstantSymbol[Bx]


(* ::Input::Initialization:: *)
ComponentValue[F[{2,-cartesian},{3,-cartesian}],Bx]


(* ::Input::Initialization:: *)
F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm


(* ::Input::Initialization:: *)
ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]];


(* ::Subsection::Initialization::Closed:: *)
(*(*\[ScriptCapitalF]:*)*)


(* ::Input::Initialization:: *)
AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues]


(* ::Input::Initialization:: *)
DefScalarFunction[\[ScriptB]]


(* ::Input::Initialization:: *)
ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],-\[ScriptB][t[],z[]]]


(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}]


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}]


(* ::Input::Initialization:: *)
ChangeComponents[pertF[{a,cartesian},{b,cartesian}],pertF[-{a,cartesian},-{b,cartesian}]];


(* ::Subsection::Initialization::Closed:: *)
(*(*\[ScriptCapitalT]:*)*)


(* ::Text::Initialization:: *)
(*(*We define two 4-vectors that \[ScriptCapitalT] is built from:*)*)


(* ::Input::Initialization:: *)
DefTensor[pertQ[-a],M,PrintAs->"\[ScriptCapitalQ]"];
DefTensor[pertU[-a],M, PrintAs->"\[ScriptCapitalU]"];


(* ::Input::Initialization:: *)
(*for now I have *)


(* ::Input::Initialization:: *)
pertTtoVec=MakeRule[{pertT[a,-b,-c],epsilonmetric[a,-b,-c,-d]pertQ[d]+1/2 (-delta[a,-c]  pertU[-b]+delta[a,-b]   pertU[-c])},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptQ]0,\[ScriptQ]2,\[ScriptU]0,\[ScriptU]2}]


(* ::Input::Initialization:: *)
AllComponentValues[pertQ[-{a,cartesian}],{\[ScriptQ]0[t[],z[]],0,\[ScriptQ]2[t[],z[]],0}]


(* ::Input::Initialization:: *)
ChangeComponents[pertQ[{a,cartesian}],pertQ[-{a,cartesian}]]


(* ::Input::Initialization:: *)
AllComponentValues[pertU[-{a,cartesian}],{\[ScriptU]0[t[],z[]],0,\[ScriptU]2[t[],z[]],0}]


(* ::Input::Initialization:: *)
ChangeComponents[pertU[{a,cartesian}],pertU[-{a,cartesian}]]


(* ::Section:: *)
(*Evaluating our field equations*)


(* ::Subsection::Initialization:: *)
(*(*Evaluating Torsion field equation*)*)


(* ::Text::Initialization:: *)
(*(*We run separatemetric to deal with epsilonmetrics after everything else is evaluated, to limit computation time*)*)


(* ::Input::Initialization:: *)
torsionC=torsionField/.pertTtoVec//ToCanonical//ContractMetric//ToBasis[cartesian];
torsionC=torsionC/.ChristoffelCDPDcartesian->Zero//ToBasis[cartesian]/.ChristoffelCDPDcartesian->Zero//TraceBasisDummy//ComponentArray//ToValues//ToValues//ToCanonical//SeparateMetric[metric]//ToBasis[cartesian];
torsionC=torsionC/.ChristoffelCDPDcartesian->Zero//TraceBasisDummy//ToCanonical//ToValues;


(* ::Input::Initialization:: *)
torsionC//MatrixForm


(* ::Subsection::Initialization:: *)
(*(*Evaluating Einstein*)*)


(* ::Text::Initialization:: *)
(*(*Might want to think about when we need/should separate and contract the metric*)*)


(* ::Input::Initialization:: *)
einsteinC=einsteinC/.pertTtoVec//ToCanonical//ToBasis[cartesian]//ToBasis[cartesian];
einsteinC=einsteinC/.ChristoffelCDPDcartesian->Zero//ContractMetric//TraceBasisDummy//ComponentArray//ToValues//ToValues//ToCanonical;


(* ::Input::Initialization:: *)
einsteinC//MatrixForm


(* ::Text::Initialization:: *)
(*(*No point in doing Maxwell - same field eqs as before.*)*)
