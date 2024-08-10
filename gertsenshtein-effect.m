(* ::Package:: *)

(* ::Title:: *)
(*The Gertsenshtein Effect*)


(* ::Subtitle:: *)
(*Following Palessandro and Rothman 2023*)


(* ::Text:: *)
(*In this file, we:*)
(*-  take a Lagrangian, expand and linearize (currently a Einstein-Maxwell Lagrangian)*)
(*-  setting a Lorentz gauge with a traceless \[ScriptH]*)
(*-  vary it with respect to \[ScriptH] (a perturbation on g) to get the Einstein field equations*)
(*-  vary it  with respect to \[ScriptCapitalA] (a perturbation on A) to get the Maxwell field equations*)
(*-  Set background field to be some constant B=Bx, and the perturbed field to be b=by(t,z), and \[ScriptH]=\[ScriptH](z,t)*)
(*- Evaluate the resulting expressions using specific components*)
(*- Take the curl of the Maxwell field equations*)
(*- Result: two PDE's connecting h and b*)


(* ::Section::Initialization::Closed:: *)
(*(*Setup*)*)


(* ::Subsection::Initialization:: *)
(*(*Loading packages*)*)


(* ::Input:: *)
(**)


(* ::Text::Initialization:: *)
(*(*xTras should load all the packages we need*)*)


(* ::Input::Initialization:: *)
<<xAct`xTras`
<<xAct`xPlain`


(* ::Text::Initialization:: *)
(*(*Some convenient settings*)*)


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
Print["The directory of this script is: ",directory];
Get@FileNameJoin@{directory,"Parallelisation.m"};


(* ::Subsection::Initialization::Closed:: *)
(*(*Manifold, basis, metric*)*)


(* ::Text::Initialization:: *)
(*(*Defining a manifold, a metric, a small perturbation h, and making the background flat with SymmetricSpaceRules*)*)


(* ::Input::Initialization:: *)
Comment@"Setup manifold, metric, chart, defining tensors and their relationships"
DefManifold[M,4,IndexRange[{a,l}]]; (*might fix formatting: type roman, look greek. See DefManifold notebook*)
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];



(* ::Text::Initialization:: *)
(*(*Chart first:*)*)


(* ::Input::Initialization:: *)
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];


(* ::Text::Initialization:: *)
(*(*Metric:*)*)


(* ::Input::Initialization:: *)
MatrixForm[MetricInBasis[metric, -cartesian,{1,-1,-1,-1}]]


(* ::Input::Initialization:: *)
MetricCompute[metric,cartesian, All]


(* ::Text::Initialization:: *)
(*(*Setting a flat background*)*)


(* ::Input::Initialization:: *)
bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];


(* ::Subsection::Initialization::Closed:: *)
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


(* ::Section::Initialization::Closed:: *)
(*Rules*)


(* ::Subsection::Initialization:: *)
(*Going between F and A*)


(* ::Input::Initialization:: *)
FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];


(* ::Subsection::Initialization:: *)
(*Perturbation rules*)


(* ::Text::Initialization:: *)
(*Setting h as perturbation on the metric*)


(* ::Input::Initialization:: *)
Perturbationmetric[LI[n_],___]/;n>1:=0
toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True]


(* ::Text::Initialization:: *)
(*Defining a perturbation on A, and setting it to pertA*)


(* ::Input::Initialization:: *)
DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M]


(* ::Input::Initialization:: *)
perturbationA[LI[n_],___]/;n>1:=0


(* ::Input::Initialization:: *)
topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True]


(* ::Text::Initialization:: *)
(*Note that both of these get rid of perturbations of second order and higher. They do that in the paper as well.*)


(* ::Subsection::Initialization:: *)
(*Between \[ScriptCapitalA] and \[ScriptCapitalF]*)


(* ::Text::Initialization:: *)
(*Now we connect the two perts:*)


(* ::Input::Initialization:: *)
pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True]


(* ::Input::Initialization:: *)
pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},MetricOn->All,ContractMetrics->True];


(* ::Section::Initialization:: *)
(*Defining and expanding Lagrangian*)


(* ::Subsection::Initialization:: *)
(*Defining Lagrangian*)


(* ::Input::Initialization:: *)
Comment@"Defining and expanding lagrangian"
EHLagrangian =Sqrt[-Detmetric[]]RicciScalarCD[];


(* ::Input::Initialization:: *)
maxwellLagrangian = Sqrt[-Detmetric[]]F[a,b]F[-a,-b];


(* ::Input::Initialization:: *)
DefConstantSymbol[\[Kappa]];


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=EHLagrangian+\[Kappa] maxwellLagrangian;


(* ::Input::Initialization:: *)
Comment@"We look at the Lagrangian"
Print["We look at the Lagrangian"]
Print[\[ScriptCapitalL]]


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=\[ScriptCapitalL]/.FtoA


(* ::Input::Initialization:: *)
Comment@"We did FtoA"
Print[\[ScriptCapitalL]]


(* ::Input::Initialization:: *)
Comment@"We did contractmetric"
Print["We did contractmetric"]
Print[\[ScriptCapitalL]]


(* ::Input::Initialization:: *)
\[ScriptCapitalL]//ContractMetric


(* ::Subsection::Initialization:: *)
(*Expanding Lagrangian*)


(* ::Input::Initialization:: *)
linearizedAction=PerturbBackground[ \[ScriptCapitalL] ,2,
BackgroundSolution->bgRules]//ExpandPerturbation;
(*linearizedAction=linearizedAction//ExpandPerturbation//ToBackground//CollectTensors;
linearizedAction=ApplyParallel[linearizedAction,{ExpandPerturbation,ToBackground,CollectTensors}];*)
(*I am struggling to get ApplyParallel to work with ExpandPertrubation for some reason...*)
linearizedAction=ApplyParallel[linearizedAction,{ToBackground,CollectTensors}];


(* ::Input::Initialization:: *)
funcToH[expr_]:=expr/.toH;
funcToPertA[expr_]:=expr/.topertA;
funcAtoF[expr_]:=expr/.AtoF;


(* ::Input::Initialization:: *)
linearizedAction=ApplyParallel[linearizedAction,{funcToH,funcToPertA, funcAtoF,ToCanonical}];
(*linearizedAction=linearizedAction/.toH/.topertA/.AtoF//ToCanonical;*)


(* ::Input::Initialization:: *)
linearizedAction = linearizedAction/.Sqrt[-Detmetric[]]->1;


(* ::Input::Initialization:: *)
Print["We ended up with"]
Print[linearizedAction]


(* ::Section::Initialization::Closed:: *)
(*Adding in traceless \[ScriptH] and Lorentz gauge*)


(* ::Text::Initialization:: *)
(*Okay, weird observation here. If I use H~automaticRules on the two containing CD I get a bug with CD being incompatible with the metric at the very end! ToCanonical[] seems to still struggle with metric, but we just don't use it and we are good-ish*)


(* ::Input::Initialization:: *)
H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];(*H~AutomaticRules~MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];*)
(*CD~AutomaticRules~MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];*)


(* ::Input::Initialization:: *)
lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];
commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
linearizedAction=linearizedAction/.commuteCD/.lorentz//ToCanonical//CollectTensors
Print["We now have"]
Print[linearizedAction//ScreenDollarIndices]


(* ::Section::Initialization::Closed:: *)
(*Field Equations*)


(* ::Subsection::Initialization:: *)
(*With respect to H:*)


(* ::Input::Initialization:: *)
einstein=ApplyParallel[linearizedAction, {VarD[H[k,l],CD]}]


funcCD[expr_]:=expr/.commuteCD;


einstein=ApplyParallel[einstein,{funcCD,ContractMetric,ToCanonical}]


Print["Einstein is now:"]
Print[einstein//ScreenDollarIndices]


(* ::Subsection::Initialization:: *)
(*With respect to pertA:*)


maxwell=ApplyParallel[linearizedAction,{VarD[pertA[k],CD]}]


(* ::Input::Initialization:: *)
maxwell=maxwell/.lorentz//ContractMetric//ToCanonical


(* ::Text::Initialization:: *)
(*In the paper they removed all terms that had h but not its derivative, because the perturbation is small.*)


(* ::Subsection::Initialization:: *)
(*Removing terms with \[ScriptH] but keeping CD[\[ScriptH]]*)


(* ::Text::Initialization:: *)
(*Note that this only works with equations linear in h*)


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
einstein=ApplyParallel[einstein,{DeleteFirstOrderPart}]


(* ::Input::Initialization:: *)
maxwell=DeleteFirstOrderPart[maxwell]


(* ::Section::Initialization::Closed:: *)
(*Components - enter xCoba*)


(* ::Text::Initialization:: *)
(*We now want to define the components of our tensors. We will create a cartesian basis, put a flat metric on that basis, set the components of \[ScriptH] and F, define a new tensor \[ScriptCapitalF] and set its components and its relationship to \[ScriptCapitalA]. For now we set all components containing a 0-index to zero, because we don't care about them.*)


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
DefScalarFunction[{h1,h2}]


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{1,-cartesian}],h1[t[],z[]]]


(* ::Input::Initialization:: *)
ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]]


(* ::Input::Initialization:: *)
ComponentValue[H[{1,-cartesian},{2,-cartesian}],h2[t[],z[]]]


(* ::Text::Initialization:: *)
(*With indices up*)


(* ::Input::Initialization:: *)
ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]]


(* ::Input::Initialization:: *)
ComponentArray[H[{a,cartesian},{b,cartesian}]]//ToValues//ToValues//MatrixForm


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
ComponentValue[F[{2,-cartesian},{3,-cartesian}],-Bx]


(* ::Input::Initialization:: *)
F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm


(* ::Input::Initialization:: *)
ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]]


(* ::Subsection::Initialization::Closed:: *)
(*\[ScriptCapitalF]:*)


(* ::Input::Initialization:: *)
AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues]


(* ::Input::Initialization:: *)
DefScalarFunction[\[ScriptB]]


(* ::Input::Initialization:: *)
ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],\[ScriptB][t[],z[]]]


(* ::Input::Initialization:: *)
DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}]


(* ::Input::Initialization:: *)
ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}](*not completely sure if it is a fcn of all three, or just t,z*)


(* ::Input::Initialization:: *)
ToValues[ToValues[ComponentArray[pertF[-{a,cartesian},-{b,cartesian}]]]//Simplification]//MatrixForm


(* ::Input::Initialization:: *)
ChangeComponents[pertF[{a, cartesian}, {b, cartesian}], pertF[-{a, cartesian}, -{b, cartesian}]]


(* ::Section:: *)
(*Evaluating our field equations*)


(* ::Subsection::Initialization::Closed:: *)
(*Evaluating Einstein*)


Print["starting to evaluate einstein"];


funcPertAtoF[expr_]:=expr/.pertAtoF;


einsteinExpr=ApplyParallel[einstein,{funcPertAtoF,ToCanonical,SeparateMetric[metric],ToBasis[cartesian],ToBasis[cartesian],TraceBasisDummy,ComponentArray,ToValues,ToValues}]


Print["Done evaluating einstein"];


(* ::Text::Initialization:: *)
(*The interesting entries here are those at (1,2) and (2,1) - they are identical, as expected*)


(* ::Input::Initialization:: *)
-4 Bx \[Kappa] \[ScriptB][t[],z[]]-\!\(\*SuperscriptBox[
InterpretationBox[
StyleBox["h2",
ShowAutoStyles->False,
AutoSpacing->False],
$CellContext`h2,
Editable->False], 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],z[]]+\!\(\*SuperscriptBox[
InterpretationBox[
StyleBox["h2",
ShowAutoStyles->False,
AutoSpacing->False],
$CellContext`h2,
Editable->False], 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[t[],z[]]


DumpSave[FileNameJoin[{directory,"gertEinsteinExpr.mx"}],einsteinExpr];


Print["Just saved einstein to file"];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*Curl of Maxwell*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
Print["moving on to Maxwell"]
DefTensor[u[a],M];
u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True];
AllComponentValues[u[-{a,cartesian}],{1,0,0,0}];


epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True];


(* ::Input::Initialization:: *)
Print["calculating maxwell"]
maxwellC=ApplyParallel[maxwell,{funcPertAtoF,ToCanonical,ToBasis[cartesian], ToBasis[cartesian]}];


(* ::Input::Initialization:: *)
maxwellCa=u[-{l,cartesian}]epsilonmetric[{l,cartesian},{i,cartesian},{f,cartesian},{k,cartesian}]CD[-{i,cartesian}][maxwellC];
maxwellExpr=ApplyParallel[maxwellCa,{TraceBasisDummy,TraceBasisDummy,ComponentArray,ToValues,ToValues,ToCanonical}]


DumpSave[FileNameJoin[{directory,"gertMaxwellExpr.mx"}],{maxwellExpr,maxwell,maxwellC,maxwellCa}];


Comment@"Comment1"
Print@"Comment2"


Quit[]
