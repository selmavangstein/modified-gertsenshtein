(*=============*)
(*  DefFields  *)
(*=============*)

DefTensor[A[-a], M];
DefTensor[F[-a,-b],M,Antisymmetric[{-a,-b}]];
DefTensor[H[-a,-b],M,Symmetric[{-a,-b}],PrintAs->"\[ScriptH]"];
DefTensor[pertA[-a],M,PrintAs->"\[ScriptCapitalA]"];
DefTensor[pertF[-a,-b],M, Antisymmetric[{-a,-b}], PrintAs->"\[ScriptCapitalF]"];
DefTensor[pertT[a,-b,-c],M, Antisymmetric[{-b,-c}], PrintAs->"\[ScriptCapitalT]"];
DefTensor[Q[a],M];
DefTensor[pertQ[-a],M,PrintAs->"\[ScriptCapitalQ]"];
DefTensor[pertU[-a],M, PrintAs->"\[ScriptCapitalU]"];
