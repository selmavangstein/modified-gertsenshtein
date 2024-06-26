(* ::Package:: *)

BeginPackage["Testpackage`"]


(* ::Text:: *)
(*We can have text here as well I see. That might look worse*)


MainFunction::usage = 
	"MainFunction[x] does simple math."


Begin["Private`"]
MainFunction[x_]:=
Module[{y},
y=x^2;
y+x
]
End[]
EndPackage[]



