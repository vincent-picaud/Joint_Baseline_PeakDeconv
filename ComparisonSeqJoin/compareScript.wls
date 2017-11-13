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
hyperParameterGridMU=Range[100,3000,200]
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
c=-1;


RealComputationJoinDeconv[]:=Return[Map[{c,#,Table[BestJoinDeconv[c,#,seed],{seed,1,10}]}&,noiseGrid]];
RealComputationSequentialDeconv[]:=Return[Map[{c,#,Table[BestSequentialDeconv[c,#,seed],{seed,1,10}]}&,noiseGrid]];


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
