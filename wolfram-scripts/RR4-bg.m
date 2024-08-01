(* ::Package:: *)

(* ::Title:: *)
(*The Gertsenshtein Effect*)


(* ::Subtitle:: *)
(*Following Palessandro and Rothman 2023*)


(* ::Section::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*Setup*)*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*Loading packages*)*)*)*)*)*)*)*)*)


(* ::Input:: *)
(**)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*xTras should load all the packages we need*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
<<xAct`xTras`
<<xAct`xPlain`


Title@"Welcome to the SANDBOX"


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*Some convenient settings*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
$PrePrint = ScreenDollarIndices;
$CovDFormat = "Prefix";
$CommuteCovDsOnScalars=True;
$CVVerbose=False;

SetOptions[ToCanonical,UseMetricOnVBundle->All];
SetOptions[ContractMetric,AllowUpperDerivatives->True];


(* ::Input::Initialization:: *)
Comment@"getting parallelization";
directory=DirectoryName[$InputFileName];
Get@FileNameJoin@{directory,"Parallelisation.m"};


(* ::Subsection::Closed:: *)
(*Manifold, basis, metric*)


(* ::Text:: *)
(*Defining a manifold, a metric, a small perturbation h, and making the background flat with SymmetricSpaceRules*)


(* ::Input::Initialization:: *)
Comment@"Setup manifold, metric, chart, defining tensors and their relationships"
DefManifold[M,4,IndexRange[{a,s}]]; (*might fix formatting: type roman, look greek. See DefManifold notebook*)
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];
DefCovD[CDT[-a],Torsion->True, SymbolOfCovD->{"#","D"},FromMetric->metric];


(* ::Input::Initialization:: *)
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];


(* ::Text:: *)
(*Place metric on that chart:*)


(* ::Input::Initialization:: *)
MetricInBasis[metric, -cartesian,{1,-1,-1,-1}];


(* ::Input::Initialization:: *)
MetricCompute[metric,cartesian, All];


(* ::Text:: *)
(*Setting a flat background*)


(* ::Input::Initialization:: *)
bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*Defining our tensors*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefTensor[A[-a], M]


(* ::Input::Initialization:: *)
DefTensor[F[-a,-b],M,Antisymmetric[{-a,-b}]]


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*Perturbations:*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];


(* ::Input::Initialization:: *)
DefTensor[pertA[-a],M,PrintAs->"\[ScriptCapitalA]"]


(* ::Input::Initialization:: *)
DefTensor[pertF[-a,-b],M, Antisymmetric[{-a,-b}], PrintAs->"\[ScriptCapitalF]"]


DefTensor[pertT[a,-b,-c],M, Antisymmetric[{-b,-c}], PrintAs->"\[ScriptCapitalT]"];


(* ::Section::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*Rules*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*Going between F and A*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];
funcAtoF[expr_]:=expr/.AtoF;


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*Perturbation rules*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*Setting h as perturbation on the metric*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
Perturbationmetric[LI[n_],___]/;n>1:=0
toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True];
funcToH[expr_]:=expr/.toH;


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*Defining a perturbation on A, and setting it to pertA*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M]


(* ::Input::Initialization:: *)
perturbationA[LI[n_],___]/;n>1:=0


(* ::Input::Initialization:: *)
topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True];
funcToPertA[expr_]:=expr/.topertA;


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationT[LI[order],a,-b,-c],TorsionCDT[a,-b,-c],M];
perturbationT[LI[n_],___]/;n>1:=0;
topertT= MakeRule[{perturbationT[LI[1],a,-b,-c],pertT[a,-b,-c]},MetricOn-> All, ContractMetrics-> True];
funcToPertT[expr_]:=expr/.topertT;


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationR[LI[order],a,-b,-c],RicciScalarCDT[],M];
perturbationR[LI[n_],___]/;n>1:=0;


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*Note that both of these get rid of perturbations of second order and higher. They do that in the paper as well.*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*Between \[ScriptCapitalA] and \[ScriptCapitalF]*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*Now we connect the two perts:*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},MetricOn->All,ContractMetrics->True];
funcPertAtoF[expr_]:=expr/.pertAtoF;


(* ::Section::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*Defining and expanding Lagrangian*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*Defining Lagrangian*)*)*)*)*)*)*)*)


Comment@"Defining and expanding lagrangian"


(* ::Input::Initialization:: *)
DefConstantSymbol[{\[Alpha],\[Kappa]4,\[Lambda]}];


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=(Sqrt[-Detmetric[]](RicciScalarCDT[]+\[Kappa]4 RiemannCDT[-a,-b,-c,-d]RiemannCDT[a,b,c,d]+\[Lambda] F[-a,-b]F[a,b]))/.FtoA


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*Writing it in terms of torsion*)*)*)*)*)*)*)*)


Comment["Writing in terms of torsion"]


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=ChangeCovD[\[ScriptCapitalL],CDT,CD];(*redundant step - no CDT's in lagrangian*)


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=ChangeCurvature[\[ScriptCapitalL],CDT,CD]//ChristoffelToGradMetric//ContractMetric//ToCanonical


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*Expanding Lagrangian*)*)*)*)*)*)*)*)


Comment["expanding lagrangian"]


Comment["This part is a little slow, but it breaks if I try to use ApplyParallel"]


(* ::Input::Initialization:: *)
linearizedAction=PerturbBackground[ \[ScriptCapitalL] ,2,
BackgroundSolution->bgRules]//ExpandPerturbation//ToBackground//CollectTensors;
(*linearizedAction=linearizedAction//ExpandPerturbation//ToBackground//CollectTensors;
linearizedAction=ApplyParallel[linearizedAction,{ExpandPerturbation,ToBackground,CollectTensors}];*)
(*I am struggling to get ApplyParallel to work with ExpandPertrubation for some reason...*)
(*linearizedAction=ApplyParallel[linearizedAction,{ToBackground,CollectTensors}]*)


Comment["Simplifying linearized action"]


(* ::Input::Initialization:: *)
linearizedAction=ApplyParallel[linearizedAction,{funcToH,funcToPertA, funcAtoF,funcToPertT,ToCanonical}]
(*linearizedAction=linearizedAction/.toH/.topertA/.AtoF//ToCanonical;*)


(* ::Input::Initialization:: *)
linearizedAction = linearizedAction/.Sqrt[-Detmetric[]]->1;


(* ::Subsection:: *)
(*Imposing background torsion*)


Comment@"Imposing background torsion";


(* ::Input::Initialization:: *)
DefTensor[Q[a],M];


(* ::Input::Initialization:: *)
TtoVec=MakeRule[{TorsionCDT[a,-b,-c],epsilonmetric[a,-b,-c,-d]Q[d]},MetricOn->All,ContractMetrics->True];
funcTtoVec[expr_]:=expr/.TtoVec;


(* ::Input::Initialization:: *)
linearizedAction=ApplyParallel[linearizedAction,{funcTtoVec,ToCanonical,ContractMetric}]


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*Adding in traceless \[ScriptH] and Lorentz gauge*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];(*H~AutomaticRules~MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];*)
(*CD~AutomaticRules~MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];*)


(* ::Input::Initialization:: *)
lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];
commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];
funcLorentz[expr_]:=expr/.lorentz;
funcCD[expr_]:=expr/.commuteCD;


(* ::Input::Initialization:: *)
linearizedAction=linearizedAction/.commuteCD/.lorentz//ToCanonical//CollectTensors;


(* ::Section::Initialization:: *)
(*(*(*(*(*(*(*(*Field Equations*)*)*)*)*)*)*)*)


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*With respect to torsion tensor:*)*)*)*)*)*)*)*)


Comment@"Calculating torsion field equations"


(* ::Input::Initialization:: *)
torsionField=ApplyParallel[linearizedAction, {VarD[pertT[k,-l,-m],CDT]}];
torsionField=ApplyParallel[torsionField, {ToCanonical,ContractMetric}];


(* ::Input::Initialization:: *)
torsionField=ChangeCovD[torsionField,CDT,CD]//ChristoffelToGradMetric;


(* ::Input::Initialization:: *)
torsionField=ApplyParallel[torsionField,{funcLorentz,funcCD,funcTtoVec,ToCanonical,ContractMetric}];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*With respect to H:*)*)*)*)*)*)*)*)


Comment@"Calculating einstein field equations"


(* ::Input::Initialization:: *)
einsteinField=ApplyParallel[linearizedAction, {VarD[H[k,l],CDT]}];


(* ::Input::Initialization:: *)
einsteinField=ApplyParallel[einsteinField,{funcPertAtoF,ToCanonical,ContractMetric}];


(* ::Input::Initialization:: *)
einsteinField=ChangeCovD[einsteinField,CDT,CD]//ChristoffelToGradMetric//ContractMetric;


(* ::Input::Initialization:: *)
einsteinField=ApplyParallel[einsteinField, {funcTtoVec,funcPertAtoF}];
einsteinField=ApplyParallel[einsteinField, {funcLorentz,funcCD,ToCanonical,ContractMetric}];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*With respect to pertA:*)*)*)*)*)*)*)*)


Comment@"Calculating maxwell field equations"


(* ::Input::Initialization:: *)
maxwellField=ApplyParallel[linearizedAction,{VarD[pertA[k],CDT]}];


(* ::Input::Initialization:: *)
maxwellField=ApplyParallel[maxwellField,{funcPertAtoF, ToCanonical,ContractMetric}];


(* ::Input::Initialization:: *)
maxwellField=ChangeCovD[maxwellField,CDT,CD]//ChristoffelToGradMetric;


(* ::Input::Initialization:: *)
maxwellField=ApplyParallel[maxwellField,{funcTtoVec,funcLorentz,funcCD,ToCanonical,ContractMetric}];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*Removing terms with \[ScriptH] but keeping CD[\[ScriptH]]*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*Note that this only works with equations linear in h*)*)*)*)*)*)*)*)


Comment@"Simplifying field eqs"


(* ::Input::Initialization:: *)
DefConstantSymbol[PerturbativeParameter, PrintAs->"\[Epsilon]"]


(* ::Input::Initialization:: *)
ToOrderH = MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},MetricOn->All,ContractMetrics->True]


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
einsteinField=ApplyParallel[einsteinField,{DeleteFirstOrderPart}];


(* ::Input::Initialization:: *)
maxwellField=ApplyParallel[maxwellField,{DeleteFirstOrderPart}];


(* ::Input::Initialization:: *)
torsionField=ApplyParallel[torsionField,{DeleteFirstOrderPart}];


(* ::Section::Initialization:: *)
(*(*(*(*(*(*(*(*Components - enter xCoba*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*We now want to define the components of our tensors.*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


Comment@"Setting all components"


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*\[Epsilon]:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,-cartesian},{1,-cartesian},{2,-cartesian},{3,-cartesian}],1},MetricOn->All,ContractMetrics->True];


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*\[ScriptH]:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*Setting them all to zero first:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
zerovalues=Table[0,{i,0,3},{j,0,3}];


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues];


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*Setting the non-zero components*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefScalarFunction[{h1,h2}];


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{1,-cartesian}],h1[t[],z[]]];


(* ::Input::Initialization:: *)
ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]];


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{2,-cartesian}],h2[t[],z[]]];


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*With indices up*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]];


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*F:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*Setting all to zero:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray@F[-{a,cartesian},-{b,-cartesian}],zerovalues];


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*Giving it a constant B-field in the x-direction:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefConstantSymbol[Bx];


(* ::Input::Initialization:: *)
ComponentValue[F[{2,-cartesian},{3,-cartesian}],Bx];


(* ::Input::Initialization:: *)
F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm;


(* ::Input::Initialization:: *)
ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]];


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*\[ScriptCapitalF]:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues];


(* ::Input::Initialization:: *)
DefScalarFunction[\[ScriptB]];


(* ::Input::Initialization:: *)
ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],-\[ScriptB][t[],z[]]];


(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}];


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}];


(* ::Input::Initialization:: *)
ChangeComponents[pertF[{a,cartesian},{b,cartesian}],pertF[-{a,cartesian},-{b,cartesian}]];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*\[ScriptCapitalT]:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*We define two 4-vectors that \[ScriptCapitalT] is built from:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefTensor[pertQ[-a],M,PrintAs->"\[ScriptCapitalQ]"];
DefTensor[pertU[-a],M, PrintAs->"\[ScriptCapitalU]"];


(* ::Input::Initialization:: *)
pertTtoVec=MakeRule[{pertT[a,-b,-c],epsilonmetric[a,-b,-c,-d]pertQ[d]+1/2 (-delta[a,-c]  pertU[-b]+delta[a,-b]   pertU[-c])},MetricOn->All,ContractMetrics->True];
funcPertTtoVec[expr_]:=expr/.pertTtoVec;


(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptQ]0,\[ScriptQ]2,\[ScriptU]0,\[ScriptU]2}];


(* ::Input::Initialization:: *)
AllComponentValues[pertQ[-{a,cartesian}],{\[ScriptQ]0[t[],z[]],0,\[ScriptQ]2[t[],z[]],0}];


(* ::Input::Initialization:: *)
ChangeComponents[pertQ[{a,cartesian}],pertQ[-{a,cartesian}]];


(* ::Input::Initialization:: *)
AllComponentValues[pertU[-{a,cartesian}],{\[ScriptU]0[t[],z[]],0,\[ScriptU]2[t[],z[]],0}];


(* ::Input::Initialization:: *)
ChangeComponents[pertU[{a,cartesian}],pertU[-{a,cartesian}]];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*(*T:*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
DefConstantSymbol[q0];


(* ::Input::Initialization:: *)
AllComponentValues[Q[-{a,cartesian}],{q0,0,0,0}];


(* ::Input::Initialization:: *)
ChangeComponents[Q[{a,cartesian}],Q[-{a,cartesian}]];


(* ::Section:: *)
(*Evaluating our field equations*)


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*Evaluating Torsion field equation*)*)*)*)*)*)*)*)


Comment@"Evaluating torsion component eqs..."


(* ::Input::Initialization:: *)
funcChristCartZero[expr_]:=expr/.ChristoffelCDPDcartesian->Zero;


(* ::Input::Initialization:: *)
torsionC=ApplyParallel[torsionField,{funcPertTtoVec,funcTtoVec,ToCanonical,ToBasis[cartesian]}]


(* ::Input::Initialization:: *)
torsionExpr=ApplyParallel[torsionC,{funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*Evaluating Einstein*)*)*)*)*)*)*)*)


Comment["Evaluating Einstein component eqs..."];


(* ::Input::Initialization:: *)
einsteinC=ApplyParallel[einsteinField,{funcPertTtoVec,funcTtoVec,funcPertAtoF, SeparateMetric[metric] ,ToCanonical,ToBasis[cartesian]}];
einsteinExpr=ApplyParallel[einsteinC, {funcChristCartZero,ToBasis[cartesian],funcChristCartZero,ToBasis[cartesian],funcChristCartZero,TraceBasisDummy,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToValues,ToCanonical}]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*(*(*(*(*(*(*(*Curl of Maxwell*)*)*)*)*)*)*)*)*)*)*)*)*)*)


Comment@"Evaluating Maxwell component eqs"


(* ::Input::Initialization:: *)
DefTensor[u[a],M];
u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True];
AllComponentValues[u[-{a,cartesian}],{1,0,0,0}];


(* ::Input::Initialization:: *)
maxwellC=ApplyParallel[maxwellField,{funcPertTtoVec,funcTtoVec, ToCanonical,ToBasis[cartesian], ToBasis[cartesian]}];
maxwellC=maxwellC/.ChristoffelCDPDcartesian->Zero;
maxwellCurl=u[-{l,cartesian}]epsilonmetric[{l,cartesian},{i,cartesian},{f,cartesian},{k,cartesian}]CD[-{i,cartesian}][maxwellC];


(* ::Input::Initialization:: *)
maxwellExpr=ApplyParallel[maxwellCurl,{ContractMetric,TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,TraceBasisDummy,ToCanonical,ToValues}];


(* ::Subsection:: *)
(*Saving and goodbye*)


Comment@"Saving results..."


DumpSave[FileNameJoin[{directory,"results/rr4-bg-results.mx"}],{\[ScriptCapitalL],maxwellCurl,maxwellC,maxwellField,maxwellExpr,einsteinField,einsteinExpr,torsionField,torsionExpr}]


Comment@"Goodbye and thanks for all the fish"


Quit[]
