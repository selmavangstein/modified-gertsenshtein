(*==============*)
(*  Lagrangian  *)
(*==============*)

MyRaggedBlock[SomeList_,SomeWidth_]:=Module[{DataArray,PaddedArray,FrameRules},
	DataArray=SomeList~Partition~(UpTo@SomeWidth);
	FrameRules=DeleteCases[Flatten@Map[((First@Evaluate@Position[Evaluate@DataArray,#])->True)&,DataArray,{2}],Null];
	TextGrid[DataArray,Frame->{None,None,FrameRules}]];

Section@"Making contact with Mike's conventions";

Comment@"Since my work from the pandemic era is based on Mike's conventions, and Selma's work is based on the xAct conventions for Riemann-Cartan geometry, we need to make contact between these in order to correctly understand the constant-torsion background.";

Comment@"The first step is to define Mike's Riemann tensor. Obviously, xAct won't know that this is the Riemann tensor.";
DefTensor[RiemannMike[a,b,-c,-d],M,{Antisymmetric[{a,b}],Antisymmetric[{-c,-d}]},PrintAs->"\[ScriptCapitalR]"];
RiemannMike[a,b,-c,-d]//DisplayExpression;

Comment@"Now we define the Ricci tensor in Mike's conventions.";
DefTensor[RicciMike[a,-b],M,PrintAs->"\[ScriptCapitalR]"];
RicciMikeToRiemannMike=MakeRule[{RicciMike[a,-c],RiemannMike[a,b,-c,-b]},
	MetricOn->All,ContractMetrics->True];
RicciMike[-a,-b]~DisplayRule~RicciMikeToRiemannMike;

Comment@"Now we define the Ricci scalar in Mike's conventions.";
DefTensor[RicciScalarMike[],M,PrintAs->"\[ScriptCapitalR]"];
RicciScalarMikeToRicciMike=MakeRule[{RicciScalarMike[],RicciMike[a,-a]},
	MetricOn->All,ContractMetrics->True];
RicciScalarMike[]~DisplayRule~RicciScalarMikeToRicciMike;

Comment@"Now we define the torsion tensor in Mike's conventions.";
DefTensor[TorsionMike[a,-b,-c],M,Antisymmetric[{-b,-c}],PrintAs->"\[ScriptCapitalT]"];
TorsionMike[a,-b,-c]//DisplayExpression;

Comment@"We can probably identify the torsion tensor in Mike's conventions with the torsion tensor in xAct without worrying too much about numerical factors (we should check this eventually).";
TorsionMikeToTorsionCDT=MakeRule[{TorsionMike[a,-b,-c],TorsionCDT[a,-b,-c]},
	MetricOn->All,ContractMetrics->True];
TorsionMike[a,-b,-c]~DisplayRule~TorsionMikeToTorsionCDT;

Comment@"Now we define the contraction of the torsion tensor in Mike's conventions.";
DefTensor[TorsionContractionMike[a],M,PrintAs->"\[ScriptCapitalT]"];
TorsionContractionMikeToTorsionMike=MakeRule[{TorsionContractionMike[-a],TorsionMike[b,-b,-a]},
	MetricOn->All,ContractMetrics->True];
TorsionContractionMike[-a]~DisplayRule~TorsionContractionMikeToTorsionMike;

Comment@"Now we define the contorsion tensor in Mike's conventions, as defined below Eq. (47) on page 12 of arXiv:1510.06699.";
DefTensor[ContorsionMike[-a,-b,-c],M,Antisymmetric[{-a,-b}],PrintAs->"\[ScriptCapitalK]"];
ContorsionMikeToTorsionMike=MakeRule[{ContorsionMike[-a,-b,-c],-(1/2)*(TorsionMike[-a,-b,-c]+TorsionMike[-b,-c,-a]-TorsionMike[-c,-a,-b])},
	MetricOn->All,ContractMetrics->True];
ContorsionMike[-a,-b,-c]~DisplayRule~ContorsionMikeToTorsionMike;

Comment@"Now we use Eq. (110a) on page 20 of arXiv:1510.06699 to define the (flat) post-Riemannian expansion of Mike's curvature tensor.";
RiemannMikeExpand=MakeRule[{RiemannMike[a,b,-c,-d],CD[-c]@ContorsionMike[a,b,-d]-CD[-d]@ContorsionMike[a,b,-c]+ContorsionMike[a,-e,-c]*ContorsionMike[e,b,-d]-ContorsionMike[a,-e,-d]*ContorsionMike[e,b,-c]},
	MetricOn->All,ContractMetrics->True];
RiemannMike[a,b,-c,-d]~DisplayRule~RiemannMikeExpand;

Comment@{"By expanding the contorsion using",Cref@{"RiemannMikeExpand","ContorsionMikeToTorsionMike","TorsionMikeToTorsionCDT"}," we obtain the following."};
Expr=RiemannMike[a,b,-c,-d];
Expr//DisplayExpression;
Expr=Expr/.RiemannMikeExpand;
Expr//=ToCanonical;
Expr//=ContractMetric;
Expr//=ScreenDollarIndices;
Expr=Expr/.ContorsionMikeToTorsionMike;
Expr//=ToCanonical;
Expr//=ContractMetric;
Expr//=ScreenDollarIndices;
Expr=Expr/.TorsionMikeToTorsionCDT;
Expr//=ToCanonical;
Expr//=ContractMetric;
Expr//=ScreenDollarIndices;
DisplayExpression[Expr,EqnLabel->"MikeVersion"];

Comment@"Now we try the same expansion using the xAct conventions (note the different order of indices, which is the point of this whole exercise).";
ExprxAct=RiemannCDT[-c,-d,a,b];
ExprxAct//DisplayExpression;
ExprxAct//DisplayExprxActession;
ExprxAct//=ChangeCurvature[#,CDT,CD]&;
ExprxAct//=ChristoffelToGradMetric;
ExprxAct//=ToCanonical;
ExprxAct//=ContractMetric;
ExprxAct//=ScreenDollarIndices;
ExprxAct=ExprxAct/.{RiemannCD->Zero};
DisplayExpression[ExprxAct,EqnLabel->"xActVersion"];

Comment@{"Now we subtract the two expressions in",Cref@{"MikeVersion","xActVersion"}," to see if they are equal."};
Expr-=ExprxAct;
Expr//=ToCanonical;
Expr//=ContractMetric;
Expr//=ScreenDollarIndices;
Expr//DisplayExpression;

Comment@"Good: we are now able to claim the following index matchings when going from Mike's conventions to those of xAct.";
RiemannMikeToRiemannCDT=MakeRule[{RiemannMike[a,b,-c,-d],RiemannCDT[-c,-d,a,b]},
	MetricOn->All,ContractMetrics->True];
RiemannMike[a,b,-c,-d]~DisplayRule~RiemannMikeToRiemannCDT;

Comment@"As a final definition, we define the Weyl field strength tensor in Mike's conventions.";
DefTensor[HMike[-a,-b],M,Antisymmetric[{-a,-b}],PrintAs->"\[ScriptCapitalH]"];
HMikeToCDT=MakeRule[{HMike[-a,-b],Evaluate@((CD[-a]@TorsionContractionMike[-b]-CD[-b]@TorsionContractionMike[-a])/.TorsionContractionMikeToTorsionMike/.TorsionMikeToTorsionCDT//ToCanonical//ContractMetric//ScreenDollarIndices)},
	MetricOn->All,ContractMetrics->True];
HMike[-a,-b]~DisplayRule~HMikeToCDT;

Supercomment@{"We think of",Cref@{"TorsionMikeToTorsionCDT","RiemannMikeToRiemannCDT","HMikeToCDT"}," as being on the same footing."};

Section@"Defining the Lagrangian";

DefConstantSymbol[Alp0,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(0\)]\)"];
DefConstantSymbol[Alp1,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)"];
DefConstantSymbol[Alp2,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(2\)]\)"];
DefConstantSymbol[Alp3,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(3\)]\)"];
DefConstantSymbol[Alp4,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(4\)]\)"];
DefConstantSymbol[Alp5,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(5\)]\)"];
DefConstantSymbol[Alp6,PrintAs->"\!\(\*SubscriptBox[\(\[Alpha]\), \(6\)]\)"];
DefConstantSymbol[Bet1,PrintAs->"\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\)"];
DefConstantSymbol[Bet2,PrintAs->"\!\(\*SubscriptBox[\(\[Beta]\), \(2\)]\)"];
DefConstantSymbol[Bet3,PrintAs->"\!\(\*SubscriptBox[\(\[Beta]\), \(3\)]\)"];
DefConstantSymbol[Zet,PrintAs->"\[Zeta]"];
DefConstantSymbol[Thet,PrintAs->"\[Theta]"];
DefConstantSymbol[ChiCon,PrintAs->"\[Chi]"];
DefConstantSymbol[XiCon,PrintAs->"\[Xi]"];

Comment@"First we define the fully covariant Lagrangian as given in Eqs. (2) and (3) on page 3 of arXiv:1510.06699.";
lagrangian=Sqrt[-Detmetric[]](
	Alp0*RicciScalarMike[]
	+Alp1*RicciScalarMike[]*RicciScalarMike[]
	+Alp2*RicciMike[-a,-b]*RicciMike[a,b]
	+Alp3*RicciMike[-a,-b]*RicciMike[b,a]
	+Alp4*RiemannMike[-a,-b,-c,-d]*RiemannMike[a,b,c,d]
	+Alp5*RiemannMike[-a,-b,-c,-d]*RiemannMike[a,c,b,d]
	+Alp6*RiemannMike[-a,-b,-c,-d]*RiemannMike[c,d,a,b]
	+Bet1*TorsionMike[-a,-b,-c]*TorsionMike[a,b,c]
	+Bet2*TorsionMike[-a,-b,-c]*TorsionMike[b,a,c]
	+Bet3*TorsionContractionMike[-a]*TorsionContractionMike[a]
	+Zet*F[-a,-b]*RicciMike[a,b]
	+Thet*F[-a,-b]*HMike[a,b]
	+(ChiCon/3)*RicciMike[-m,-n]*HMike[m,n]
	-(XiCon/36)*HMike[-m,-n]*HMike[m,n]
	-(1/4)*F[-a,-b]F[a,b]);
lagrangian//=ToCanonical;
lagrangian//=ContractMetric;
lagrangian//=ScreenDollarIndices;
lagrangian//=Simplify;
DisplayExpression[lagrangian,EqnLabel->"SuperLagrangian"];

Supercomment@{"We use",Cref@{"TorsionMikeToTorsionCDT","RiemannMikeToRiemannCDT","HMikeToCDT"}," to move",Cref@"SuperLagrangian"," into Selma's xAct conventions."};
lagrangian=lagrangian/.RicciScalarMikeToRicciMike/.RicciMikeToRiemannMike/.RiemannMikeToRiemannCDT/.TorsionMikeToTorsionCDT/.TorsionContractionMikeToTorsionMike/.HMikeToCDT;
lagrangian=lagrangian/.RicciScalarMikeToRicciMike/.RicciMikeToRiemannMike/.RiemannMikeToRiemannCDT/.TorsionMikeToTorsionCDT/.TorsionContractionMikeToTorsionMike/.HMikeToCDT;
lagrangian=lagrangian/.RicciScalarMikeToRicciMike/.RicciMikeToRiemannMike/.RiemannMikeToRiemannCDT/.TorsionMikeToTorsionCDT/.TorsionContractionMikeToTorsionMike/.HMikeToCDT;
lagrangian=lagrangian/.RicciScalarMikeToRicciMike/.RicciMikeToRiemannMike/.RiemannMikeToRiemannCDT/.TorsionMikeToTorsionCDT/.TorsionContractionMikeToTorsionMike/.HMikeToCDT;
lagrangian//=ToCanonical;
lagrangian//=ContractMetric;
lagrangian//=ScreenDollarIndices;
lagrangian//=Simplify;
lagrangian//DisplayExpression;

Section@"Special theories of interest";

Comment@"Now we define the Karananas couplings, and give them in terms of Mike's couplings. These equations are taken from Eq. (B.23b) on page 143 of my thesis.";
DefConstantSymbol[kLambda,PrintAs->"\[Lambda]"];
DefConstantSymbol[kR1,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptR]\), \(1\)]\)"];
DefConstantSymbol[kR2,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptR]\), \(2\)]\)"];
DefConstantSymbol[kR3,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptR]\), \(3\)]\)"];
DefConstantSymbol[kR4,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptR]\), \(4\)]\)"];
DefConstantSymbol[kR5,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptR]\), \(5\)]\)"];
DefConstantSymbol[kR6,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptR]\), \(6\)]\)"];
DefConstantSymbol[kT1,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptT]\), \(1\)]\)"];
DefConstantSymbol[kT2,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptT]\), \(2\)]\)"];
DefConstantSymbol[kT3,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptT]\), \(3\)]\)"];
MikeKarananasEquations={Alp0==2*kLambda,Bet1==kT1/3+kT2/12+kLambda/4,Bet2==kT1/3-kT2/6+kLambda/2,Bet3==-kT1/3+2*kT3/3-kLambda,Alp1==kR6,Alp2==kR4+kR5,Alp3==kR4-kR5,Alp4==kR1/3+kR2/6,Alp5==2*kR1/3-2*kR2/3,Alp6==kR1/3+kR2/6-kR3};
DisplayExpression[(MikeKarananasEquations~MyRaggedBlock~5),EqnLabel->"MikeKarananas"];

Comment@"Now we introduce the cosmological couplings, and give them in terms of Mike's couplings. These equations are taken from Eq. (B.24a) on page 144 of my thesis.";
DefConstantSymbol[Sigma1,PrintAs->"\!\(\*SubscriptBox[\(\[Sigma]\), \(1\)]\)"];
DefConstantSymbol[Sigma2,PrintAs->"\!\(\*SubscriptBox[\(\[Sigma]\), \(2\)]\)"];
DefConstantSymbol[Sigma3,PrintAs->"\!\(\*SubscriptBox[\(\[Sigma]\), \(3\)]\)"];
DefConstantSymbol[Upsilon1,PrintAs->"\!\(\*SubscriptBox[\(\[Upsilon]\), \(1\)]\)"];
DefConstantSymbol[Upsilon2,PrintAs->"\!\(\*SubscriptBox[\(\[Upsilon]\), \(2\)]\)"];
DefConstantSymbol[MPl,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptCapitalM]\), \(Pl\)]\)"];
DefConstantSymbol[CConstant,PrintAs->"\[CapitalLambda]"];
MikeCosmologicalEquations={
	Sigma1==3*Alp1/2+Alp2/4+Alp3/4+Alp5/4-Alp6/2,
	Sigma2==3*Alp1/2+Alp2/2+Alp3/2+3*Alp4/2-Alp5/4+Alp6/2,
	Sigma3==3*Alp1/2+Alp2/2+Alp3/2+Alp4/2+Alp5/4+Alp6/2,
	Upsilon1==-2*Bet1+2*Bet2,
	Upsilon2==2*Bet1+Bet2+3*Bet3
};
DisplayExpression[(MikeCosmologicalEquations~MyRaggedBlock~5),EqnLabel->"MikeCosmological"];

Comment@"Now we define Case 2 from arXiv:1910.14197 in terms of Karananas' couplings. These equations are taken from Table 1 on page 5/6 of arXiv:1910.14197.";
Case2EquationsKarananas={kR1==0,kR3/2-kR4==0,kT1==0,kLambda==0,kR6==0};
DisplayExpression[(Case2EquationsKarananas~MyRaggedBlock~5),EqnLabel->"Case2"];

Comment@{"Now we use",Cref@{"MikeKarananas","Case2"}," define Case 2 from arXiv:1910.14197 in terms of Mike's couplings."};
KarananasToMike=First@Solve[MikeKarananasEquations,{kLambda,kR1,kR2,kR3,kR4,kR5,kR6,kT1,kT2,kT3}];
Case2EquationsMike=Case2EquationsKarananas/.KarananasToMike//FullSimplify;
DisplayExpression[(Case2EquationsMike~MyRaggedBlock~5),EqnLabel->"Case2Mike"];

Comment@{"Now we use",Cref@"Case2Mike"," define Case 2 from arXiv:1910.14197 as a set of rules."};
Case2RulesMike=First@Quiet@Solve[Case2EquationsMike,{Alp0,Alp1,Alp2,Alp3,Alp4,Alp5,Alp6,Bet1,Bet2,Bet3}];
DisplayExpression[(Case2RulesMike~MyRaggedBlock~5),EqnLabel->"Case2RulesMike"];

Comment@{"Now we use",Cref@{"MikeCosmological","Case2Mike"}," to consider the effect of Case 2 from arXiv:1910.14197 on the cosmological couplings."};
MikeCosmologicalEquations=MikeCosmologicalEquations/.Case2RulesMike//FullSimplify;
DisplayExpression@(MikeCosmologicalEquations~MyRaggedBlock~5);
Comment@"We see that only one of the \"special\" cosmological conditions actually follows from Case 2 from arXiv:1910.14197, i.e. the vanishing of the third sigma parameter. Next, we will explore the combined implications of Case 2 from arXiv:1910.14197 with the preferred cosmological interpretation.";

Comment@"Here are the preferred cosmological conditions, which are collected from comments in Chapters 2 and 3 of my thesis."; 
PreferredCosmologicalEquations={Sigma1==Sigma2,Sigma3==0,Upsilon1==Sigma1*CConstant/MPl^2,Upsilon2==-4/3*MPl^2};
DisplayExpression[(PreferredCosmologicalEquations~MyRaggedBlock~5),EqnLabel->"PreferredCosmological"];

Comment@{"Here are the combined implications of Case 2 from arXiv:1910.14197 with the preferred cosmological interpretation expressed as a set of rules. We obtain these by combining",Cref@{"PreferredCosmological","Case2Mike"},"."};
PreferredCosmologicalEquations=PreferredCosmologicalEquations/.Case2RulesMike//FullSimplify;
MikeCosmologicalToMike=First@Quiet@Solve[MikeCosmologicalEquations,{Sigma1,Sigma2,Sigma3,Upsilon1,Upsilon2}];
UnifiedRules=First@Quiet@Solve[Join[PreferredCosmologicalEquations/.MikeCosmologicalToMike,Case2EquationsMike],{Alp0,Alp1,Alp2,Alp3,Alp4,Alp5,Alp6,Bet1,Bet2,Bet3}];
DisplayExpression[(UnifiedRules~MyRaggedBlock~5),EqnLabel->"UnifiedRules"];

Supercomment@{"For Anthony. You will presumably wish to apply",Cref@"UnifiedRules"," to your setup, or perhaps the less restrictive",Cref@"Case2RulesMike","."};

GRRules={Alp0->-MPl^2/2,Alp1->0,Alp2->0,Alp3->0,Alp4->0,Alp5->0,Alp6->0,Bet1->0,Bet2->0,Bet3->0,q0->0};

PureGRRules={Alp0->-MPl^2/2,Alp1->0,Alp2->0,Alp3->0,Alp4->0,Alp5->0,Alp6->0,Bet1->0,Bet2->0,Bet3->0,q0->0,Zet->0,Thet->0,ChiCon->0,XiCon->0};

FlatCTEG=(UnifiedRules/.{CConstant->0})~Join~{CConstant->0};
