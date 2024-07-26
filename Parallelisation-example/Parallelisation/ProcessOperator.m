(* ::Package:: *)

(*===================*)
(*  ProcessOperator  *)
(*===================*)

ProcessOperator::identify="Studying the operator `1` now.";

ProcessOperator[InputFileName_]:=Module[{Expr},
		
	Expr=TheOperator;
	ProcessOperator::identify~Message~TheOperator;
	Expr=VarD[H[a,b],CD][Expr];
Expr];
