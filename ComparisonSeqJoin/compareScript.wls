#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Comparison script*)


(* ::Text:: *)
(*Script used to generate comparison (included in the associated article)*)


(* ::Title:: *)
(*Script code*)


(* ::Subchapter:: *)
(*Project root directory & other initialization*)


SetDirectory["~/GitHub/Joint_Baseline_PeakDeconv"];
On[Assert];


(* ::Subchapter:: *)
(*Some misc functions*)


(* ::Subsection:: *)
(*Generate ground truth*)


GenerateGroundTruthPeakList[]:=Block[{errorCode,file},
	file=CreateFile[];
	errorCode=Run["./generateSynthetic --what 2 > "<>file];
	If[errorCode!=0,Return[$Failed]];
	Return[Import[file,"csv"][[2;;-1,{1,2}]]];
];


GenerateGroundTruthSpectrum[c_?IntegerQ]:=Block[{errorCode,file},
	file=CreateFile[];
	errorCode=Run["./generateSynthetic --what 1 -c "<>ToString[c]<>" -n 0 > "<>file];
	If[errorCode!=0,Return[$Failed]];
	Return[AssociationThread[{"x","y","peakConv","baseline"}->Transpose[Import[file,"csv"][[2;;-1,1;;-2]]]]];
];


(* ::Subsection:: *)
(*Avoid some mistake... Be sure of the current ground truth*)


(* ::Text:: *)
(*CAVEAT: must be run once for initialization*)


Unprotect["TruePeakList"];
Unprotect["TrueSpectrum"];
Unprotect["TrueC"];
TrueC=1;
TruePeakList=GenerateGroundTruthPeakList[];
TrueSpectrum=GenerateGroundTruthSpectrum[TrueC];
Protect["TruePeakList"];
Protect["TrueSpectrum"];
Protect["TrueC"];


(* ::Subsection:: *)
(*Generate synthetic data*)


GenerateData[c_?IntegerQ,noise_,seed_?IntegerQ]:=Block[{errorCode,file},
	Assert[Abs[c]<=1]; (* c={-1,0,1} *)
	Assert[noise>=0];
	Assert[seed>=0];
	file=CreateFile[];
	errorCode=Run["./generateSynthetic -c "<>ToString[c]<>" -n "<>ToString[noise]<>" -s "<>ToString[seed]<>" > "<>file];
	If[errorCode!=0,Return[$Failed]];
	Return[file];
];


(* ::Text:: *)
(*Call me to generate new data*)


GenerateData[-1,0.2,0]


(* ::Subsection:: *)
(*Run a computation and return result stored into an association*)


ReadResult[filename_?StringQ]:=Block[{tmp},
	If[!FileExistsQ[filename],
		Print[filename<>" does not exist!"];
		Return[$Failed]
	];
	tmp=Transpose[Import[filename,"csv"][[2;;-1]]];
	Assert[Length[tmp]==6];
	Return[AssociationThread[{"x","y","peakDirac","baseline","peakConv","peakBaseline"}->tmp]];
];


ResultToPeakList[<|"x"->x_,"y"->y_,"peakDirac"->peakDirac_,"baseline"->baseline_,"peakConv"->peakConv_,"peakBaseline"->peakBaseline_|>]:=Transpose[Select[Transpose[{x,peakDirac}],#[[2]]>0&]]


PlotResult[<|"x"->x_,"y"->y_,"peakDirac"->peakDirac_,"baseline"->baseline_,"peakConv"->peakConv_,"peakBaseline"->peakBaseline_|>]:=
	ListLinePlot[{Transpose@{x,y},Transpose@{x,peakDirac},Transpose@{x,peakBaseline},Transpose@{x,baseline}},PlotRange->All];
PlotResult[filename_?StringQ]:=PlotResult[ReadResult[filename]];


SequentialDeconv[filename_?StringQ,nSNIP_?IntegerQ,\[Lambda]1_]:=Block[{errorCode,commandLine,fileOut,result},
	Assert[(nSNIP>2)&&(nSNIP<60)];
	Assert[(\[Lambda]1>=0)];
	Assert[FileExistsQ[filename]];
	fileOut=CreateFile[];
	commandLine="./sequencialDeconvolution -i "<>filename<>" -o "<>fileOut<>" --snip "<>ToString[nSNIP]<>" --lambda_1 "<>ToString[\[Lambda]1];
	errorCode=Run[commandLine];
	If[errorCode!=0,
		Print["The command line:\n",commandLine,"\nFAILED, error code:=",errorCode];
		Return[$Failed]
	];
	result=ReadResult[fileOut];
	DeleteFile[fileOut];
	Return[result];
];


JointDeconv[filename_?StringQ,\[Mu]_,\[Lambda]1_]:=Block[{errorCode,commandLine,fileOut,result},
	Assert[(\[Lambda]1>=0)];
	Assert[(\[Mu]>=0)];
	Assert[FileExistsQ[filename]];
	fileOut=CreateFile[];
	commandLine="./jointDeconvolution -i "<>filename<>" -o "<>fileOut<>" --mu "<>ToString[\[Mu]]<>" --lambda_1 "<>ToString[\[Lambda]1];
	errorCode=Run[commandLine];
		If[errorCode!=0,
		Print["The command line:\n",commandLine,"\nFAILED, error code:=",errorCode];
		Return[$Failed]
	];
	result=ReadResult[fileOut];
	DeleteFile[fileOut];
	Return[result];
];


(* ::Subsection:: *)
(*Compute error/quality factors*)


ComputeError[<|"x"->xComputed_,"y"->yComputed_,"peakDirac"->peakDiracComputed_,"baseline"->baselineComputed_,"peakConv"->peakConvComputed_,"peakBaseline"->peakBaselineComputed_|>]:=
	Norm[peakConvComputed-TrueSpectrum["peakConv"],2]/Norm[TrueSpectrum["peakConv"],2];



filename=GenerateData[-1,2,0];
res=SequentialDeconv[filename,40,0];
PlotResult[res]
ResultToPeakList[res]
ComputeError[res]


res2=JointDeconv[filename,1000,0.1];
PlotResult[res2]
ResultToPeakList[res2]
ComputeError[res2]


(* ::Subsection:: *)
(*Generate hyperparameters grid*)


hyperParameterGridSNIP=Range[20,40,2]
hyperParameterGridMU=Range[0,4500,250]
hyperParameterGridLAMBDA1=Range[0,4,.2]


hyperParameterGridSequential=Partition[Flatten[Outer[##&,hyperParameterGridSNIP,hyperParameterGridLAMBDA1]],2];
Length @ hyperParameterGridSequential
hyperParameterGridJoint=Partition[Flatten[Outer[##&,hyperParameterGridMU,hyperParameterGridLAMBDA1]],2];
Length @ hyperParameterGridJoint


(* ::Subsubsection:: *)
(*Return the best hyperparameter*)


BestJoinDeconv[filename_?StringQ]:=SortBy[ParallelMap[Function[{#[[1]],#[[2]],ComputeError[JointDeconv[filename,#[[1]],# [[2]]]]}],hyperParameterGridJoint],#[[3]]&][[1]]
BestSequentialDeconv[filename_?StringQ]:=SortBy[ParallelMap[Function[{#[[1]],#[[2]],ComputeError[SequentialDeconv[filename,#[[1]],# [[2]]]]}],hyperParameterGridSequential],#[[3]]&][[1]]


(* ::Text:: *)
(*Extend to newly generate data*)


BestJoinDeconv[c_,noise_,seed_]:=Block[{filename},
	filename=GenerateData[c,noise,seed];
	Return[BestJoinDeconv[filename]];
];
BestSequentialDeconv[c_,noise_,seed_]:=Block[{filename},
	filename=GenerateData[c,noise,seed];
	Return[BestSequentialDeconv[filename]];
];	


(* ::Subsection:: *)
(*Real computation (Take a lot of times!)*)


noiseGrid=Range[0,2,0.2]


RealComputationJoinDeconv[c_?IntegerQ,noiseGrid_?ListQ]:=Return[Map[{c,#,Table[BestJoinDeconv[c,#,seed],{seed,1,10}]}&,noiseGrid]];
RealComputationSequentialDeconv[c_?IntegerQ,noiseGrid_?ListQ]:=Return[Map[{c,#,Table[BestSequentialDeconv[c,#,seed],{seed,1,10}]}&,noiseGrid]];


(* ::Chapter:: *)
(*Result: caveat, takes a long time (>1 hour)!*)


(*
resultJoincPos=RealComputationJoinDeconv[1,noiseGrid];
resultSeqcPos=RealComputationSequentialDeconv[1,noiseGrid];

resultJoincNeg=RealComputationJoinDeconv[-1,noiseGrid];
resultSeqcNeg=RealComputationSequentialDeconv[-1,noiseGrid];

comparisonResult=<|"resultJoincPos"\[Rule]resultJoincPos,"resultSeqcPos"\[Rule]resultSeqcPos,"resultJoincNeg"\[Rule]resultJoincNeg,"resultSeqcNeg"\[Rule]resultSeqcNeg|>;
Save["comparisonResult.wl",comparisonResult]
*)


(* ::Subchapter:: *)
(*Here only the saved result!*)


(* ::Text:: *)
(*We reload the saved association <|"resultJoincPos"->resultJoincPos,"resultSeqcPos"->resultSeqcPos,"resultJoincNeg"->resultJoincNeg,"resultSeqcNeg"->resultSeqcNeg|>*)


Get["comparisonResult.wl"];
Keys[comparisonResult]


(* ::Subchapter:: *)
(*Here some functions to visualize results*)


ConstantArrayQ=Equal@@#&;

(* result format: a list of elements {c,noise,{{mu_or_snip,lambda_1,merit}, .... }} *)
FormatGlobalResult[resSeq_,resJoin_]:=Block[{i,n,assocError,plotError,assoc\[Lambda]1,plot\[Lambda]1,assocSNIP,plotSNIP,assocMu,plotMu,toMap},
(* assume c=constant *)
Assert[ConstantArrayQ @ Join[resSeq[[All,1]],resJoin[[All,1]]]];
(* same noise sample *)
Assert[resSeq[[All,2]]==resJoin[[All,2]]];

n=Length[resSeq];
Assert[n>0];

assocError=<||>;
assoc\[Lambda]1=<||>;
assocSNIP=<||>;
assocMu=<||>;
For[i=1,i<=n,i++,
assocError[ToString[resJoin[[i,2]]]]=<|"Sequential processing"->resSeq[[i,3]][[All,3]],"Proposed method"->resJoin[[i,3]][[All,3]]|>;
assoc\[Lambda]1[ToString[resJoin[[i,2]]]]=<|"Sequential processing"->resSeq[[i,3]][[All,2]],"Proposed method"->resJoin[[i,3]][[All,2]]|>;
assocSNIP[ToString[resJoin[[i,2]]]]=<|"Sequential processing"->resSeq[[i,3]][[All,1]]|>;
assocMu[ToString[resJoin[[i,2]]]]=<|"Proposed method"->resJoin[[i,3]][[All,1]]|>;
];
plotError=BoxWhiskerChart[assocError,ChartLegends->Placed[Automatic,{0.25,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(e\)]\)"},ImageSize->400];
plot\[Lambda]1=BoxWhiskerChart[assoc\[Lambda]1,ChartLegends->Placed[Automatic,{0.25,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\)"},ImageSize->400];
plotSNIP=BoxWhiskerChart[assocSNIP,ChartLegends->Placed[Automatic,{0.25,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\!\(\*SubscriptBox[\(m\), \(SNIP\)]\)"},ImageSize->400];
plotMu=BoxWhiskerChart[assocMu,ChartLegends->Placed[Automatic,{0.25,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\[Mu]"},ImageSize->400];

Return[{plotError,plot\[Lambda]1,plotSNIP,plotMu}];
];


FormatTeXGlobalResult[resSeq_,resJoin_]:=Block[{i,n,assoc,data,g1,g2,f,tabError,ToLaTex},
(* assume c=constant *)
Assert[ConstantArrayQ @ Join[resSeq[[All,1]],resJoin[[All,1]]]];
(* same noise sample *)
Assert[resSeq[[All,2]]==resJoin[[All,2]]];

n=Length[resSeq];
Assert[n>0];

g1[l_]:={N[Mean[l]],N[StandardDeviation[l]]};
g2[{m_,std_}]:=ToString[TeXForm[EngineeringForm[m,3]]]<>"\\pm "<>ToString[TeXForm[EngineeringForm[std,3]]];
f[l_]:=g2[g1[l]];

assoc=<|"c"->1,"\[Sigma]noise"->2, "\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(e\)]\)Seq"->3, "\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(e\)]\)Join"->4, "mSNIP"->5, "\[Mu]"->6, "\[Lambda]1Seq"->7, "\[Lambda]1Join"->8|>;
data={};
For[i=1,i<=n,i++,
data=Append[data,
{
ToString[resJoin[[i,1]]],
ToString[resJoin[[i,2]]],
f[resSeq[[i,3]][[All,3]]],
f[resJoin[[i,3]][[All,3]]],
f[resSeq[[i,3]][[All,1]]],
f[resJoin[[i,3]][[All,1]]],
f[resSeq[[i,3]][[All,2]]],
f[resJoin[[i,3]][[All,2]]]
}];
];
data=Map[Function["$"<>#<>"$"],data,{2}];
(*data=StringReplace[Map[StringJoin,r],"$$"->"$ & $"];*)
data=Transpose[data];

tabError=Transpose[{data[[2]],data[[3]],data[[4]]}];

ToLaTex[tab_]:=Block[{jsize,t1,t2,t3,t4},
jsize=Dimensions[tab][[2]];
t1=StringReplace[Map[StringJoin,tab],"$$"->"$ & $"]; (* & separator *)
t2=StringJoin[Map[Function["\n"<>#<>" \\\\ \\hline"],t1]]; (* ends line with \\ *)
(*t3=StringTake[t2,1;;-4]; remove final \\ *)
t3=t2;
(* 2 columns Error Join / Parallet *)
t4="\\begin{tabular}{"<>Table["|c",jsize]<>"|}"<>t3<>"\n\\end{tabular}\n";
Return[t4];
];

Return[ToLaTex[tabError]];
];


(* ::Subsection:: *)
(*Usage example*)


(* ::Text:: *)
(*Usage example: BoxWhiskerChart[FormatGlobalResult[resSeq,resJoin],ChartLegends->Automatic,ChartLabels->{Automatic,None}]*)


FormatGlobalResult[comparisonResult["resultSeqcPos"],comparisonResult["resultJoincPos"]]


FormatGlobalResult[comparisonResult["resultSeqcNeg"],comparisonResult["resultJoincNeg"]]


FormatTeXGlobalResult[comparisonResult["resultSeqcPos"],comparisonResult["resultJoincPos"]]


FormatTeXGlobalResult[comparisonResult["resultSeqcNeg"],comparisonResult["resultJoincNeg"]]


 


(* ::Chapter:: *)
(*"Manual" usage examples*)


(* ::Text:: *)
(*This were our first tests, we only keep them for illustratation purpose.*)


(* ::Section:: *)
(*Generate data*)


data=Import["!./generateSynthetic -s 2 | tee data3.txt","csv"][[2;;-1]];
ListLinePlot[data]


(* ::Section:: *)
(*Deconvolution using sequential approach*)


errorCode=Run["./sequencialDeconvolution data3.txt"];
Assert[errorCode==0]


deconvolvedData=Import["data3.txt.out","csv"][[2;;-1]];


ListLinePlot[deconvolvedData[[All,5]]]


deconvolvedData[[1]]


Run["ls qsf"]


(* ::Section:: *)
(*Generate synthetic data *)


Print[Directory[]]
