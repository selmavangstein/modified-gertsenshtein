(* ::Package:: *)

(* ::Title::Initialization:: *)
(*Modified Gertsenshtein Effect*)


(* ::Text::Initialization:: *)
(*ERF-term in lagrangian, no background torsion*)


Print["Running"]


(* ::Section::Initialization:: *)
(*Setup*)


(* ::Subsection::Initialization:: *)
(*Loading packages*)


(* ::Text::Initialization:: *)
(*xTras should load all the packages we need*)


(* ::Input::Initialization:: *)
<<xAct`xTras`
<<xAct`xPlain`


(* ::Text::Initialization:: *)
(*Some convenient settings*)


(* ::Input::Initialization:: *)
$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;

SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];


(* ::Input::Initialization:: *)
Print["getting parallelization"]
directory=DirectoryName[$InputFileName];
Get@FileNameJoin@{directory,"Parallelisation.m"};


(* ::Subsection::Initialization:: *)
(*Manifold, basis, metric*)


(* ::Text::Initialization:: *)
(*Defining a manifold, a metric, a small perturbation h, and making the background flat with SymmetricSpaceRules*)


(* ::Input::Initialization:: *)
Print["Setup manifold, metric, chart, defining tensors and their relationships"];
DefManifold[M,4,IndexRange[{a,s}]]; (*might fix formatting: type roman, look greek. See DefManifold notebook*)


(* ::Input::Initialization:: *)
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];


(* ::Input::Initialization:: *)
DefCovD[CDT[-a],Torsion->True, SymbolOfCovD->{"#","D"},FromMetric->metric];


(* ::Text::Initialization:: *)
(*Chart first:*)


(* ::Input::Initialization:: *)
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];


(* ::Text::Initialization:: *)
(*Place metric on that chart:*)


(* ::Input::Initialization:: *)
MatrixForm[MetricInBasis[metric, -cartesian,{1,-1,-1,-1}]];


(* ::Input::Initialization:: *)
MetricCompute[metric,cartesian, All];


(* ::Text::Initialization:: *)
(*Setting a flat background*)


(* ::Input::Initialization:: *)
bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];


(* ::Subsection::Initialization:: *)
(*Defining our tensors*)


(* ::Input::Initialization:: *)
DefTensor[A[-a], M];


(* ::Input::Initialization:: *)
DefTensor[F[-a,-b],M,Antisymmetric[{-a,-b}]];


(* ::Text::Initialization:: *)
(*Perturbations:*)


(* ::Input::Initialization:: *)
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];


(* ::Input::Initialization:: *)
DefTensor[pertA[-a],M,PrintAs->"\[ScriptCapitalA]"];


(* ::Input::Initialization:: *)
DefTensor[pertF[-a,-b],M, Antisymmetric[{-a,-b}], PrintAs->"\[ScriptCapitalF]"];


(* ::Input::Initialization:: *)
DefTensor[pertT[a,-b,-c],M, Antisymmetric[{-b,-c}], PrintAs->"\[ScriptCapitalT]"];


(* ::Section::Initialization:: *)
(*Rules*)


(* ::Subsection::Initialization:: *)
(*Going between F and A*)


(* ::Input::Initialization:: *)
FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];


(* ::Subsection::Initialization:: *)
(*Perturbation rules*)


(* ::Text::Initialization:: *)
(*Setting h as perturbation on the metric*)


(* ::Input::Initialization:: *)
Perturbationmetric[LI[n_],___]/;n>1:=0
toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True];


(* ::Text::Initialization:: *)
(*Defining a perturbation on A, and setting it to pertA*)


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M];


(* ::Input::Initialization:: *)
perturbationA[LI[n_],___]/;n>1:=0;


(* ::Input::Initialization:: *)
topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True];


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationT[LI[order],a,-b,-c],TorsionCDT[a,-b,-c],M];
perturbationT[LI[n_],___]/;n>1:=0;
topertT= MakeRule[{perturbationT[LI[1],a,-b,-c],pertT[a,-b,-c]},MetricOn-> All, ContractMetrics-> True];


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationR[LI[order],a,-b,-c],RicciScalarCDT[],M];


(* ::Input::Initialization:: *)
perturbationR[LI[n_],___]/;n>1:=0;


(* ::Subsection::Initialization:: *)
(*Between \[ScriptCapitalA] and \[ScriptCapitalF]*)


(* ::Text::Initialization:: *)
(*Now we connect the two perts:*)


(* ::Input::Initialization:: *)
pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},MetricOn->All,ContractMetrics->True];


(* ::Section::Initialization:: *)
(*Defining and expanding Lagrangian*)


(* ::Subsection::Initialization:: *)
(*Defining Lagrangian*)


Comment@"Defining and expanding lagrangian";
Print@"Defining and expanding lagrangian";


(* ::Input::Initialization:: *)
DefConstantSymbol[\[Kappa]];
DefConstantSymbol[\[Lambda]];


(* ::Text::Initialization:: *)
(*This is our modified lagrangian:*)


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=Sqrt[-Detmetric[]](RicciScalarCDT[]+epsilonmetric[a,b,-c,-d]\[Lambda] RicciCDT[c,d]F[-a,-b]+\[Kappa] F[a,b]F[-a,-b]);


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=\[ScriptCapitalL]/.FtoA;


(* ::Subsection::Initialization:: *)
(*Writing it in terms of torsion*)


(* ::Text::Initialization:: *)
(*I have these two steps to get T and R specifically out there, and can still work easily with CD and the metric*)


(* ::Input::Initialization:: *)
ChangeCovD[%,CDT,CD](*redundant step - no CDT's in lagrangian*)


(* ::Input::Initialization:: *)
ChangeCurvature[%,CDT,CD]


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=%//ChristoffelToGradMetric//ContractMetric//ToCanonical


(* ::Subsection::Initialization:: *)
(*Expanding Lagrangian*)


(* ::Input::Initialization:: *)
PerturbBackground[  \[ScriptCapitalL],2,
BackgroundSolution->bgRules];
%//ExpandPerturbation//ToBackground//CollectTensors;


(* ::Input::Initialization:: *)
%/.toH/.topertA/.topertT;
%/.AtoF//ToCanonical;
linearizedAction = %/.Sqrt[-Detmetric[]]->1


(* ::Subsection::Initialization:: *)
(*Defining traceless \[ScriptH] and Lorentz gauge*)


(* ::Text::Initialization:: *)
(*Okay, weird observation here. If I use H~automaticRules on the two containing CD I get a bug with CD being incompatible with the metric at the very end! ToCanonical[] seems to still struggle with metric, but we just don't use it and we are good-ish*)


(* ::Input::Initialization:: *)
H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];
commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
linearizedAction=linearizedAction/.commuteCD/.lorentz//ToCanonical//CollectTensors


(* ::Section::Initialization:: *)
(*Field Equations*)


(* ::Subsection::Initialization:: *)
(*With respect to H:*)


(* ::Input::Initialization:: *)
VarD[H[a,b],CDT][linearizedAction]


(* ::Input::Initialization:: *)
einstein =%/.TorsionCDT->Zero//ContractMetric//ToCanonical


(* ::Text::Initialization:: *)
(*Now we convert any CDT's in there and set background torsion to zero*)


(* ::Input::Initialization:: *)
ChangeCovD[einstein,CDT,CD]//ChristoffelToGradMetric;


(* ::Input::Initialization:: *)
%/.TorsionCDT->Zero/.pertAtoF//ToCanonical//ContractMetric


(* ::Text::Initialization:: *)
(*We can now apply the lorentz gauge as well:*)


(* ::Input::Initialization:: *)
einstein=%/.lorentz/.commuteCD


(* ::Subsection::Initialization:: *)
(*With respect to pertA:*)


(* ::Input::Initialization:: *)
VarD[pertA[a],CDT][linearizedAction];


(* ::Input::Initialization:: *)
maxwell=%/.TorsionCDT->Zero//ContractMetric//ToCanonical


(* ::Text::Initialization:: *)
(*We do a similar simplification to above:*)


(* ::Input::Initialization:: *)
ChangeCovD[maxwell,CDT,CD]//ChristoffelToGradMetric;


(* ::Input::Initialization:: *)
%/.TorsionCDT->Zero/.pertAtoF//ToCanonical//ContractMetric


(* ::Input::Initialization:: *)
maxwell=%/.lorentz/.commuteCD


(* ::Subsection::Initialization:: *)
(*With respect to torsion tensor:*)


(* ::Input::Initialization:: *)
VarD[pertT[a,-b,-c],CDT][linearizedAction]


(* ::Input::Initialization:: *)
torsioneq=%/.TorsionCDT->Zero//ContractMetric//ToCanonical


(* ::Input::Initialization:: *)
ChangeCovD[torsioneq,CDT,CD]//ChristoffelToGradMetric;


(* ::Input::Initialization:: *)
torsioneq=%/.TorsionCDT->Zero/.lorentz/.commuteCD/.pertAtoF//ToCanonical//ContractMetric


(* ::Subsection::Initialization:: *)
(*Removing terms with \[ScriptH] but keeping CD[\[ScriptH]]*)


(* ::Text::Initialization:: *)
(*Note that this only works with equations linear in h*)


(* ::Input::Initialization:: *)
DefConstantSymbol[PerturbativeParameter, PrintAs->"\[Epsilon]"]


(* ::Input::Initialization:: *)
ToOrderH = MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},MetricOn->All,ContractMetrics->True]


(* ::Text::Initialization:: *)
(*Since we have converted all derivatives from CDT to CD, this should still work:*)


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
einstein=DeleteFirstOrderPart[einstein]


(* ::Input::Initialization:: *)
maxwell=DeleteFirstOrderPart[maxwell]


(* ::Input::Initialization:: *)
torsioneq=DeleteFirstOrderPart[torsioneq]


(* ::Section::Initialization:: *)
(*Components - enter xCoba*)


(* ::Text::Initialization:: *)
(*We now want to define the components of our tensors.*)


(* ::Subsection::Initialization:: *)
(*\[Epsilon]:*)


(* ::Input::Initialization:: *)
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,-cartesian},{1,-cartesian},{2,-cartesian},{3,-cartesian}],1},MetricOn->All,ContractMetrics->True]


(* ::Subsection::Initialization:: *)
(*\[ScriptH]:*)


(* ::Text::Initialization:: *)
(*Setting them all to zero first:*)


(* ::Input::Initialization:: *)
zerovalues=Table[0,{i,0,3},{j,0,3}]


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues]


(* ::Text::Initialization:: *)
(*Setting the non-zero components*)


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptH]1,\[ScriptH]2}]


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{1,-cartesian}],\[ScriptH]1[t[],z[]]]


(* ::Input::Initialization:: *)
ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]]


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{2,-cartesian}],\[ScriptH]2[t[],z[]]]


(* ::Text::Initialization:: *)
(*With indices up*)


(* ::Input::Initialization:: *)
ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]]


(* ::Input::Initialization:: *)
ComponentArray[H[{a,cartesian},{b,cartesian}]]//ToValues


(* ::Input::Initialization:: *)
%//ToValues//ToValues//MatrixForm


(* ::Subsection::Initialization:: *)
(*F:*)


(* ::Text::Initialization:: *)
(*Setting all to zero:*)


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray@F[-{a,cartesian},-{b,-cartesian}],zerovalues]


(* ::Text::Initialization:: *)
(*Giving it a constant B-field in the x-direction:*)


(* ::Input::Initialization:: *)
DefConstantSymbol[Bx]


(* ::Input::Initialization:: *)
ComponentValue[F[{2,-cartesian},{3,-cartesian}],Bx]


(* ::Input::Initialization:: *)
F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm


(* ::Input::Initialization:: *)
ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]]


(* ::Subsection::Initialization:: *)
(*\[ScriptCapitalF]:*)


(* ::Input::Initialization:: *)
AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues]


(* ::Input::Initialization:: *)
DefScalarFunction[\[ScriptB]]


(* ::Input::Initialization:: *)
ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],-\[ScriptB][t[],z[]]]


(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}]


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}](*not completely sure if it is a fcn of all three, or just t,z*)


(* ::Input::Initialization:: *)
ToValues[ToValues[ComponentArray[pertF[-{a,cartesian},-{b,cartesian}]]]//Simplification]//MatrixForm


(* ::Input::Initialization:: *)
ChangeComponents[pertF[{a,cartesian},{b,cartesian}],pertF[-{a,cartesian},-{b,cartesian}]]


(* ::Subsection::Initialization:: *)
(*\[ScriptCapitalT]:*)


(* ::Text::Initialization:: *)
(*We define two 4-vectors that \[ScriptCapitalT] is built from:*)


(* ::Input::Initialization:: *)
DefTensor[pertQ[-a],M,PrintAs->"\[ScriptCapitalQ]"];
DefTensor[pertU[-a],M, PrintAs->"\[ScriptCapitalU]"];


(* ::Input::Initialization:: *)
(*for now I have *)


(* ::Input::Initialization:: *)
pertTtoVec=MakeRule[{pertT[a,-b,-c],epsilonmetric[a,-b,-c,-d]pertQ[d]+1/2 (-delta[a,-c]  pertU[-b]+delta[a,-b]   pertU[-c])},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)



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


(* ::Text::Initialization:: *)
(*We look at what we end up with for our perturbed torsion tensor:*)


(* ::Input::Initialization:: *)
pertT[a,-b,-c]/.pertTtoVec//SeparateMetric[metric]//ToBasis[cartesian]


(* ::Input::Initialization:: *)
%//TraceBasisDummy//ComponentArray


(* ::Input::Initialization:: *)
%//ToCanonical


(* ::Input::Initialization:: *)
%//ToValues


(* ::Input::Initialization:: *)
%//MatrixForm


(* ::Text::Initialization:: *)
(*This is now our torsion vector. Note that it took quite a lot of work to get just here, so it seems probable that we need to limit the x,z-components like we have for our perturbed maxwell field*)


(* ::Section::Initialization:: *)
(*Evaluating our field equations*)


(* ::Subsection::Initialization:: *)
(*Evaluating Torsion field equation*)


(* ::Input::Initialization:: *)
torsioneqCart=torsioneq/.pertTtoVec//ToCanonical//ContractMetric//ToBasis[cartesian]//ToBasis[cartesian]


(* ::Input::Initialization:: *)
torsioneqCart1=torsioneqCart/.ChristoffelCDPDcartesian->Zero//TraceBasisDummy//ComponentArray


(* ::Input::Initialization:: *)
torsioneqCart2=torsioneqCart1//ToValues//ToValues


(* ::Input::Initialization:: *)
%//SeparateMetric[metric]


(* ::Input::Initialization:: *)
%//ToBasis[cartesian]


(* ::Input::Initialization:: *)
%/.ChristoffelCDPDcartesian->Zero


(* ::Input::Initialization:: *)
%//TraceBasisDummy


(* ::Input::Initialization:: *)
torsioneqCart3=%//ComponentArray


(* ::Input::Initialization:: *)
%//ToCanonical


(* ::Input::Initialization:: *)
%//ToValues


(* ::Input::Initialization:: *)
torsionExpr = %//MatrixForm


(* ::Subsection::Initialization:: *)
(*Evaluating Einstein*)


(* ::Text::Initialization:: *)
(*We are now ready to evaluate it:*)


(* ::Input::Initialization:: *)
einstein/.pertTtoVec//ToCanonical//ToBasis[cartesian]//ToBasis[cartesian]


(* ::Input::Initialization:: *)
einsteinCart=%/.ChristoffelCDPDcartesian->Zero


(* ::Input::Initialization:: *)
%//ContractMetric


(* ::Input::Initialization:: *)
einsteinCart1=%//TraceBasisDummy//ComponentArray


(* ::Input::Initialization:: *)
einsteinCart2=einsteinCart1//ToValues//ToValues


(* ::Input::Initialization:: *)
%//ToCanonical


(* ::Input::Initialization:: *)
einsteinExpr=%//MatrixForm


(* ::Text::Initialization:: *)
(*Eqs with h-b-coupling: normal*)


(* ::Subsection::Initialization:: *)
(*Curl of Maxwell*)


(* ::Input::Initialization:: *)
DefTensor[u[a],M]


(* ::Input::Initialization:: *)
u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
AllComponentValues[u[-{a,cartesian}],{1,0,0,0}]


(* ::Input::Initialization:: *)
maxwell


(* ::Input::Initialization:: *)
maxwellCart=maxwell/.pertTtoVec//ToCanonical//ToBasis[cartesian]//ToBasis[cartesian]


(* ::Input::Initialization:: *)
%/.ChristoffelCDPDcartesian->Zero


(* ::Input::Initialization:: *)
maxwellCurl=u[-{l,cartesian}]epsilonmetric[{l,cartesian},{i,cartesian},{f,cartesian},{a,cartesian}]CD[-{i,cartesian}][%]


(* ::Input::Initialization:: *)
maxwellCurl1=maxwellCurl//ContractMetric//TraceBasisDummy


(* ::Input::Initialization:: *)
maxwellCurl2=maxwellCurl1//ComponentArray


(* ::Input::Initialization:: *)
maxwellCurl3=maxwellCurl2//ToValues


(* ::Input::Initialization:: *)
maxwellCurl4=maxwellCurl3//ToValues


(* ::Input::Initialization:: *)
maxwellExpr=maxwellCurl4//ToCanonical


(* ::Input::Initialization:: *)
%//MatrixForm


(* ::Input::Initialization:: *)
%//SeparateMetric[metric]


(* ::Input::Initialization:: *)
%//ToBasis[cartesian]


(* ::Input::Initialization:: *)
%//TraceBasisDummy


(* ::Input::Initialization:: *)
%//ToCanonical


(* ::Input::Initialization:: *)
%//ToValues


(* ::Input::Initialization:: *)
maxwellExpr=%


(* ::Input::Initialization:: *)
%//MatrixForm
