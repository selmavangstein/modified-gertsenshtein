(* ::Package:: *)

(*===================*)
(*  ProcessOperator  *)
(*===================*)

ProcessOperator::identify="Studying the operator `1` now.";

ProcessOperator[]:=Module[{Expr},
		
	Expr=TheOperator;
	ProcessOperator::identify~Message~TheOperator;
	Expr=VarD[H[q,r],CD][Expr];
Expr];
