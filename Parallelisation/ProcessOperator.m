(*===================*)
(*  ProcessOperator  *)
(*===================*)

ProcessOperator::identify="Studying the operator `1` now.";

ProcessOperator[InputFileName_]:=Module[{Expr},
		
	Expr=TheOperator;
	ProcessOperator::identify~Message~TheOperator;
	Expr=Expr<>"1";
Expr];
