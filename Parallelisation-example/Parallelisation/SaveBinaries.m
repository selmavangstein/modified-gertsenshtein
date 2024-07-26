(*================*)
(*  SaveBinaries  *)
(*================*)

DumpSave[FileNameJoin[
	{$ThisDirectory,#<>".mx"}],#]&/@{"xAct`xTensor`",
	"xAct`xTensor`Private`",
	"TangentM4`",
	"Global`"};
