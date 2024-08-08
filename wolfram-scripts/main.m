(*========*)
(*  main  *)
(*========*)

$ThisDirectory=If[NotebookDirectory[]==$Failed,Directory[],NotebookDirectory[],NotebookDirectory[]];
Get@FileNameJoin@{$ThisDirectory,"gertsenshtein-package-bg.m"};

resultsFileName="FullModel";

Section@"Analysis of GR";
Supercomment@"Analysis of GR without any constant-torsion background.";
StudySystem[GRRules];
Section@"Analysis of Case 2";
Supercomment@"Analysis of Case 2 as defined by Yun-Cherng.";
StudySystem[Case2RulesMike];
Section@"Analysis of CTEG";
Supercomment@"Analysis of constant-torsion emergent gravity."; 
StudySystem[UnifiedRules];
Section@"Analysis of flat CTEG";
Supercomment@"Analysis of constant-torsion emergent gravity without the cosmological constant."; 
StudySystem[FlatCTEG];
Supercomment@"This is the end of the script.";

Quit[];
