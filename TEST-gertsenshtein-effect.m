(* ::Package:: *)

(* ::Title:: *)
(*The Gertsenshtein Effect*)


(* ::Subtitle:: *)
(*Following Palessandro and Rothman 2023*)


(* ::Input::Initialization:: *)
Title@"Gertsenshtein effect script"
Comment@"Starting the script"


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


(* ::Section::Initialization:: *)
(*Setup*)


(* ::Subsection::Initialization:: *)
(*Loading packages*)


(* ::Text::Initialization:: *)
(*xTras should load all the packages we need*)


(* ::Input::Initialization:: *)
<<xAct`xTras`


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
Get@FileNameJoin@{NotebookDirectory[],"Parallelisation.m"};


(* ::Subsection::Initialization:: *)
(*Manifold, basis, metric*)


(* ::Text::Initialization:: *)
(*Defining a manifold, a metric, a small perturbation h, and making the background flat with SymmetricSpaceRules*)


(* ::Input::Initialization:: *)
DefManifold[M,4,IndexRange[{a,l}]]; (*might fix formatting: type roman, look greek. See DefManifold notebook*)
DefMetric[-1,metric[-a,-b],CD,PrintAs->"g",SymCovDQ->True];



(* ::Text::Initialization:: *)
(*Chart first:*)


(* ::Input::Initialization:: *)
DefChart[cartesian,M,{0,1,2,3},{t[],x[],y[],z[]}];


(* ::Text::Initialization:: *)
(*Metric:*)


(* ::Input::Initialization:: *)
MatrixForm[MetricInBasis[metric, -cartesian,{1,-1,-1,-1}]]


(* ::Input::Initialization:: *)
MetricCompute[metric,cartesian, All]


(* ::Text::Initialization:: *)
(*Setting a flat background*)


(* ::Input::Initialization:: *)
bgRules = SymmetricSpaceRules[CD,0];
SetOptions[ToBackground,BackgroundSolution->bgRules];


(* ::Subsection::Initialization:: *)
(*Defining our tensors*)


(* ::Input::Initialization:: *)
DefTensor[A[-a], M]


(* ::Input::Initialization:: *)
DefTensor[F[-a,-b],M,Antisymmetric[{-a,-b}]]


(* ::Text::Initialization:: *)
(*Perturbations:*)


(* ::Input::Initialization:: *)
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];


(* ::Input::Initialization:: *)
DefTensor[pertA[-a],M,PrintAs->"\[ScriptCapitalA]"]


(* ::Input::Initialization:: *)
DefTensor[pertF[-a,-b],M, Antisymmetric[{-a,-b}], PrintAs->"\[ScriptCapitalF]"]


(* ::Section:: *)
(*Rules*)


(* ::Subsection:: *)
(*Going between F and A*)


(* ::Input:: *)
(*FtoA = MakeRule[{F[-a,-b],CD[-a][A[-b]]-CD[-b][A[-a]]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*AtoF = MakeRule[{CD[-a]@A[-b],(1/2)*F[-a,-b]+(1/2)*(CD[-a]@A[-b]+CD[-b]@A[-a])},MetricOn->All,ContractMetrics->True];*)


(* ::Subsection:: *)
(*Perturbation rules*)


(* ::Text:: *)
(*Setting h as perturbation on the metric*)


(* ::Input:: *)
(*Perturbationmetric[LI[n_],___]/;n>1:=0*)
(*toH =MakeRule[{Perturbationmetric[LI[1],-a,-b],H[-a,-b]},MetricOn-> All, ContractMetrics-> True]*)


(* ::Text:: *)
(*Defining a perturbation on A, and setting it to pertA*)


(* ::Input:: *)
(*DefTensorPerturbation[perturbationA[LI[order],-a],A[-a],M]*)


(* ::Input:: *)
(*perturbationA[LI[n_],___]/;n>1:=0*)


(* ::Input:: *)
(*topertA = MakeRule[{perturbationA[LI[1],-a],pertA[-a]},MetricOn-> All, ContractMetrics-> True]*)


(* ::Text:: *)
(*Note that both of these get rid of perturbations of second order and higher. They do that in the paper as well.*)


(* ::Subsection:: *)
(*Between \[ScriptCapitalA] and \[ScriptCapitalF]*)


(* ::Text:: *)
(*Now we connect the two perts:*)


(* ::Input:: *)
(*pertFtoA = MakeRule[{pertF[-a,-b],CD[-a][pertA[-b]]-CD[-b][pertA[-a]]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*pertAtoF = MakeRule[{CD[-a]@pertA[-b],(1/2)*pertF[-a,-b]+(1/2)*(CD[-a]@pertA[-b]+CD[-b]@pertA[-a])},MetricOn->All,ContractMetrics->True];*)


(* ::Section::Initialization:: *)
(*Defining and expanding Lagrangian*)


(* ::Subsection::Initialization:: *)
(*Defining Lagrangian*)


(* ::Input::Initialization:: *)
EHLagrangian =Sqrt[-Detmetric[]]RicciScalarCD[]


(* ::Input::Initialization:: *)
maxwellLagrangian = Sqrt[-Detmetric[]]F[a,b]F[-a,-b]


(* ::Input::Initialization:: *)
DefConstantSymbol[\[Kappa]]


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=EHLagrangian+\[Kappa] maxwellLagrangian


(* ::Input::Initialization:: *)
\[ScriptCapitalL]=\[ScriptCapitalL]/.FtoA


(* ::Input::Initialization:: *)
Comment@"The adapted lagrangian is"
Print[\[ScriptCapitalL]]


(* ::Subsection::Initialization:: *)
(*Expanding Lagrangian*)


(* ::Input::Initialization:: *)
PerturbBackground[ \[ScriptCapitalL] ,2,
BackgroundSolution->bgRules]
%//ExpandPerturbation//ToBackground//CollectTensors


(* ::Input::Initialization:: *)
%/.toH/.topertA;
%/.AtoF//ToCanonical;
linearizedAction = %/.Sqrt[-Detmetric[]]->1


(* ::Section::Closed:: *)
(*Adding in traceless \[ScriptH] and Lorentz gauge*)


(* ::Text:: *)
(*Okay, weird observation here. If I use H~automaticRules on the two containing CD I get a bug with CD being incompatible with the metric at the very end! ToCanonical[] seems to still struggle with metric, but we just don't use it and we are good-ish*)


(* ::Input:: *)
(*H~AutomaticRules~ MakeRule[{H[-a,a],0},MetricOn->All,ContractMetrics->True];(*H~AutomaticRules~MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];*)*)
(*(*CD~AutomaticRules~MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];*)*)


(* ::Input:: *)
(*lorentz = MakeRule[{CD[-a][H[a,b]],0},MetricOn->All,ContractMetrics->True];*)
(*commuteCD = MakeRule[{CD[-a]@CD[c]@H[a,b],0},MetricOn->All,ContractMetrics->True];*)


(* ::Input:: *)
(*linearizedAction=linearizedAction/.commuteCD/.lorentz//ToCanonical//CollectTensors*)


(* ::Section::Closed:: *)
(*Field Equations*)


(* ::Subsection::Closed:: *)
(*With respect to H:*)


(* ::Input:: *)
(*VarD[H[a,b],CD][linearizedAction]*)


(* ::Input:: *)
(*einstein = %/.commuteCD//ContractMetric//ToCanonical*)


(* ::Subsection::Closed:: *)
(*With respect to pertA:*)


(* ::Input:: *)
(*VarD[pertA[a],CD][linearizedAction]*)


(* ::Input:: *)
(*maxwell=%/.lorentz//ContractMetric//ToCanonical*)


(* ::Text:: *)
(*In the paper they removed all terms that had h but not its derivative, because the perturbation is small.*)


(* ::Subsection:: *)
(*Removing terms with \[ScriptH] but keeping CD[\[ScriptH]]*)


(* ::Text:: *)
(*Note that this only works with equations linear in h*)


(* ::Input:: *)
(*DefConstantSymbol[PerturbativeParameter, PrintAs->"\[Epsilon]"]*)


(* ::Input:: *)
(*ToOrderH = MakeRule[{H[-a,-b],PerturbativeParameter*H[-a,-b]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*ToOrderCDH=MakeRule[{CD[-c]@H[-a,-b],PerturbativeParameter*CD[-c]@H[-a,-b]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*DeleteFirstOrderPart[InputExpr_]:=Module[{OutputLagrangian=InputExpr,SecondOrderPart,FirstOrderPart,NullOrderPart},OutputLagrangian=OutputLagrangian/.ToOrderCDH/.ToOrderH;*)
(*SecondOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,2}]&;*)
(*SecondOrderPart//=Normal;*)
(*SecondOrderPart=SecondOrderPart/.PerturbativeParameter->1;*)
(*SecondOrderPart//=ToCanonical;*)
(*SecondOrderPart//=ContractMetric;*)
(*SecondOrderPart//=ScreenDollarIndices;*)
(*SecondOrderPart//=CollectTensors;*)
(**)
(*FirstOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,1}]&;*)
(*FirstOrderPart//=Normal;*)
(*FirstOrderPart=FirstOrderPart/.PerturbativeParameter->1;*)
(*FirstOrderPart//=ToCanonical;*)
(*FirstOrderPart//=ContractMetric;*)
(*FirstOrderPart//=ScreenDollarIndices;*)
(*FirstOrderPart//=CollectTensors;*)
(**)
(*NullOrderPart=OutputLagrangian//Series[#,{PerturbativeParameter,0,0}]&;*)
(*NullOrderPart//=Normal;*)
(*NullOrderPart=NullOrderPart/.PerturbativeParameter->1;*)
(*NullOrderPart//=ToCanonical;*)
(*NullOrderPart//=ContractMetric;*)
(*NullOrderPart//=ScreenDollarIndices;*)
(*NullOrderPart//=CollectTensors;*)
(**)
(*OutputLagrangian=SecondOrderPart-(FirstOrderPart-NullOrderPart);*)
(*OutputLagrangian//=ToCanonical;*)
(*OutputLagrangian//=ContractMetric;*)
(*OutputLagrangian//=ScreenDollarIndices;*)
(*OutputLagrangian//=CollectTensors;*)
(*OutputLagrangian];*)


(* ::Input:: *)
(*einstein=DeleteFirstOrderPart[einstein]*)


(* ::Input:: *)
(*maxwell=DeleteFirstOrderPart[maxwell]*)


(* ::Section:: *)
(*Components - enter xCoba*)


(* ::Text:: *)
(*We now want to define the components of our tensors. We will create a cartesian basis, put a flat metric on that basis, set the components of \[ScriptH] and F, define a new tensor \[ScriptCapitalF] and set its components and its relationship to \[ScriptCapitalA]. For now we set all components containing a 0-index to zero, because we don't care about them.*)


(* ::Subsection::Closed:: *)
(*\[ScriptH]:*)


(* ::Text:: *)
(*Setting them all to zero first:*)


(* ::Input:: *)
(*zerovalues=Table[0,{i,0,3},{j,0,3}]*)


(* ::Input:: *)
(*ComponentValue[ComponentArray@H[-{a,cartesian},-{b,-cartesian}],zerovalues]*)


(* ::Text:: *)
(*Setting the non-zero components*)


(* ::Input:: *)
(*DefScalarFunction[{h1,h2}]*)


(* ::Input:: *)
(*ComponentValue[H[{1,-cartesian},{1,-cartesian}],h1[t[],z[]]]*)


(* ::Input:: *)
(*ComponentValue[H[{2,-cartesian},{2,-cartesian}],-H[{1,-cartesian},{1,-cartesian}]]*)


(* ::Input:: *)
(*ComponentValue[H[{1,-cartesian},{2,-cartesian}],h2[t[],z[]]]*)


(* ::Text:: *)
(*With indices up*)


(* ::Input:: *)
(*ChangeComponents[H[{a,cartesian},{b,cartesian}],H[-{a,cartesian},-{b,cartesian}]]*)


(* ::Input:: *)
(*ComponentArray[H[{a,cartesian},{b,cartesian}]]//ToValues*)


(* ::Input:: *)
(*%//ToValues//MatrixForm*)


(* ::Subsection:: *)
(*F:*)


(* ::Text:: *)
(*Setting all to zero:*)


(* ::Input:: *)
(*ComponentValue[ComponentArray@F[-{a,cartesian},-{b,-cartesian}],zerovalues]*)


(* ::Text:: *)
(*Giving it a constant B-field in the x-direction:*)


(* ::Input:: *)
(*DefConstantSymbol[Bx]*)


(* ::Input:: *)
(*ComponentValue[F[{2,-cartesian},{3,-cartesian}],-Bx]*)


(* ::Input:: *)
(*F[-{a,cartesian},-{b,cartesian}]//ComponentArray//ToValues//MatrixForm*)


(* ::Input:: *)
(*ChangeComponents[F[{a,cartesian},{b,cartesian}],F[-{a,cartesian},-{b,cartesian}]]*)


(* ::Subsection:: *)
(*\[ScriptCapitalF]:*)


(* ::Input:: *)
(*AllComponentValues[pertF[-{a,cartesian},-{b,cartesian}],zerovalues]*)


(* ::Input:: *)
(*DefScalarFunction[\[ScriptB]]*)


(* ::Input:: *)
(*ComponentValue[pertF[{1,-cartesian},{3,-cartesian}],\[ScriptB][t[],z[]]]*)


(* ::Input:: *)
(*DefScalarFunction[{\[ScriptCapitalE]x,\[ScriptCapitalE]y,\[ScriptCapitalE]z}]*)


(* ::Input:: *)
(*ComponentValue[ComponentArray[pertF[{0,-cartesian},-{a,cartesian}]],{0,\[ScriptCapitalE]x[t[],x[],y[],z[]],\[ScriptCapitalE]y[t[],x[],y[],z[]],\[ScriptCapitalE]z[t[],x[],y[],z[]]}](*not completely sure if it is a fcn of all three, or just t,z*)*)


(* ::Input:: *)
(*ToValues[ToValues[ComponentArray[pertF[-{a,cartesian},-{b,cartesian}]]]//Simplification]//MatrixForm*)


(* ::Input:: *)
(*ChangeComponents[pertF[{a,cartesian},{b,cartesian}],pertF[-{a,cartesian},-{b,cartesian}]]*)


(* ::Section:: *)
(*Evaluating our field equations*)


(* ::Subsection:: *)
(*Evaluating Einstein*)


(* ::Input:: *)
(*einstein/.pertAtoF//ToCanonical*)


(* ::Input:: *)
(*einsteinExpr=einstein/.pertAtoF//ToCanonical//SeparateMetric[metric] //ToBasis[cartesian]//ToBasis[cartesian]//TraceBasisDummy//ComponentArray//ToValues//ToValues//MatrixForm*)


(* ::Text:: *)
(*The interesting entries here are those at (1,2) and (2,1) - they are identical, as expected*)


(* ::Input:: *)
(*-4 Bx \[Kappa] \[ScriptB][t[],z[]]-\!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["h2",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`h2,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "2"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],z[]]+\!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["h2",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`h2,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"2", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],z[]]*)


(* ::Subsection:: *)
(*Curl of Maxwell*)


(* ::Input:: *)
(*DefTensor[u[a],M]*)


(* ::Input:: *)
(*u~AutomaticRules ~MakeRule[{u[a]u[-a],1},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*AllComponentValues[u[-{a,cartesian}],{1,0,0,0}]*)


(* ::Input:: *)
(*epsilonmetric~AutomaticRules~MakeRule[{epsilonmetric[{0,cartesian},{1,cartesian},{2,cartesian},{3,cartesian}],1},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*maxwell/.pertAtoF//ToCanonical*)


(* ::Input:: *)
(*curlfunc[expr_,a_,uvec_,coords_]:=Module[{maxwell=expr,maxcurl,l,i,f},*)
(*maxwell=maxwell/.pertAtoF//ToCanonical;*)
(*maxwell=maxwell//ToBasis[cartesian]//ToBasis[cartesian];*)
(*maxcurl=uvec[-{l,coords}]epsilonmetric[{l,coords},{i,coords},{f,coords},{a,coords}]CD[-{i,coords}][maxwell];*)
(*maxcurl=maxcurl//TraceBasisDummy//TraceBasisDummy//ComponentArray//ToValues//ToValues;*)
(*maxcurl=maxcurl//ToCanonical//MatrixForm;*)
(*maxcurl];*)


(* ::Input:: *)
(*maxwellExpr=curlfunc[maxwell,a,u,cartesian]*)


(* ::Text:: *)
(*So we need Maxwell again: 8 \[Kappa] \!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptCapitalE]x",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptCapitalE]x,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],x[],y[],z[]]-8 \[Kappa] \!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptCapitalE]z",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptCapitalE]z,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],x[],y[],z[]] is -b^(2,0), the term we were missing. How can I apply it?*)


(* ::Input:: *)
(*cheatEtoB=MakeRule[{\!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptCapitalE]x",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptCapitalE]x,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[], x[], y[], z[]],-\!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptB]",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptB],*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"2", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[], z[]]+\!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptCapitalE]z",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptCapitalE]z,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[], x[], y[], z[]]},MetricOn->All,ContractMetrics->True]*)


(* ::Input:: *)
(*-8 Bx \[Kappa] \!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["h2",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`h2,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "2"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],z[]]+8 \[Kappa] \!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptB]",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptB],*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "2"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],z[]]+8 \[Kappa] \!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptCapitalE]x",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptCapitalE]x,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],x[],y[],z[]]-8 \[Kappa] \!\(\*SuperscriptBox[*)
(*InterpretationBox[*)
(*StyleBox["\[ScriptCapitalE]z",*)
(*ShowAutoStyles->False,*)
(*AutoSpacing->False],*)
(*$CellContext`\[ScriptCapitalE]z,*)
(*Editable->False], *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[t[],x[],y[],z[]]/.cheatEtoB//ToCanonical*)


(* ::Text:: *)
(*I mean it looks right at least. This is not the way I want to be doing it though.*)


(* ::Text:: *)
(*I could maybe just use makerule to do the em field equations for both \[ScriptCapitalF] and F - I am not sure how well makerule does with components though, thus if that would do everything I wanted to.*)


(* ::Input::Initialization:: *)
Quit[]
