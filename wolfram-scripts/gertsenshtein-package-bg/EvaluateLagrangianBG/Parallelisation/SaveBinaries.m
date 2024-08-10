(* ::Package:: *)

(*================*)
(*  SaveBinaries  *)
(*================*)

DumpSave[FileNameJoin[
	{$ThisDirectory,#<>".mx"}],#]&/@{"xAct`xTensor`",
	"xAct`xTensor`Private`","xAct`xCoba`",
	"xAct`xCoba`Private`",
	"TangentM`",
	"Global`"};
