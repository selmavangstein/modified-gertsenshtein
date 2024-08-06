(*=======================*)
(*  QuietParallelSubmit  *)
(*=======================*)

QuietParallelSubmit~SetAttributes~HoldAll;
QuietParallelSubmit[Expr_]:=ParallelSubmit@Block[{Print=Null&,PrintTemporary=Null&},Expr];
