(* ::Package:: *)

(*===================*)
(*  ProcessOperator  *)
(*===================*)

ProcessOperator::identify="Studying the operator `1` now.";

ProcessOperator[commands_List]:=Module[{Expr},
		
	Expr=TheOperator;
	ProcessOperator::identify~Message~TheOperator;
	Expr=Fold[#2[#1]&,Expr,commands];
Expr];
