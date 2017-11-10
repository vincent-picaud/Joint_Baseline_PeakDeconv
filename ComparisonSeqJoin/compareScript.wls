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


GenerateData[c_?IntegerQ,noise_,seed_?IntegerQ]:=Block[{errorCode},
	Assert[Abs[c]<=1]; (* c={-1,0,1} *)
	Assert[noise>=0];
	Assert[seed>=0];
	errorCode=Run["./generateSynthetic -c "<>ToString[c]<>" -n "<>ToString[noise]<>" -s "<>ToString[seed]<>" > data.txt"];
	If[errorCode!=0,Return[$Failed]];
];


(* ::Text:: *)
(*Call me to generate new data*)


GenerateData[-1,0.2,0]


(* ::Subsection:: *)
(*Run a computation and return result stored into an association*)


ReadResult[filename_?StringQ]:=Block[{tmp},
	If[!FileExistsQ[filename],Return[$Failed]];
	tmp=Transpose[Import[filename,"csv"][[2;;-1]]];
	Assert[Length[tmp]==6];
	Return[AssociationThread[{"x","y","peakDirac","baseline","peakConv","peakBaseline"}->tmp]];
];


ResultToPeakList[<|"x"->x_,"y"->y_,"peakDirac"->peakDirac_,"baseline"->baseline_,"peakConv"->peakConv_,"peakBaseline"->peakBaseline_|>]:=Transpose[Select[Transpose[{x,peakDirac}],#[[2]]>0&]]


PlotResult[<|"x"->x_,"y"->y_,"peakDirac"->peakDirac_,"baseline"->baseline_,"peakConv"->peakConv_,"peakBaseline"->peakBaseline_|>]:=
	ListLinePlot[{Transpose@{x,y},Transpose@{x,peakDirac},Transpose@{x,peakBaseline},Transpose@{x,baseline}},PlotRange->All]


SequentialDeconv[nSNIP_?IntegerQ,\[Lambda]1_]:=Block[{errorCode},
	Assert[(nSNIP>2)&&(nSNIP<60)];
	Assert[(\[Lambda]1>=0)];
	errorCode=Run["./sequencialDeconvolution -i data.txt --snip "<>ToString[nSNIP]<>" --lambda_1 "<>ToString[\[Lambda]1]];
	If[errorCode!=0,Return[$Failed]];
	Return[ReadResult["data.txt.out"]]
];


JointDeconv[\[Mu]_,\[Lambda]1_]:=Block[{errorCode},
	Assert[(\[Lambda]1>=0)];
	Assert[(\[Mu]>=0)];
	errorCode=Run["./jointDeconvolution -i data.txt --mu "<>ToString[\[Mu]]<>" --lambda_1 "<>ToString[\[Lambda]1]];
	If[errorCode!=0,Return[$Failed]];
	Return[ReadResult["data.txt.out"]]
];


(* ::Subsection:: *)
(*Compute error/quality factors*)


ComputeError[<|"x"->xComputed_,"y"->yComputed_,"peakDirac"->peakDiracComputed_,"baseline"->baselineComputed_,"peakConv"->peakConvComputed_,"peakBaseline"->peakBaselineComputed_|>]:=
	Norm[peakConvComputed-TrueSpectrum["peakConv"],2]/Norm[TrueSpectrum["peakConv"],2];



res=SequentialDeconv[40,0];
PlotResult[res]
ResultToPeakList[res]
ComputeError[res]


ComputeError[res]


res2=JointDeconv[200,0];
PlotResult[res2]
ResultToPeakList[res2]
ComputeError[res2]


(* ::Subsection:: *)
(*Generate hyperparameters grid*)


hyperParameterGridSNIP=Range[20,40,2]
hyperParameterGridMU=Range[200,2400,200]
hyperParameterGridLAMBDA1=Range[0,4,2.2]


hyperParameterGridSequential=Partition[Flatten[Outer[##&,hyperParameterGridSNIP,hyperParameterGridLAMBDA1]],2];
Length @ hyperParameterGridSequential
hyperParameterGridJoint=Partition[Flatten[Outer[##&,hyperParameterGridMU,hyperParameterGridLAMBDA1]],2];
Length @ hyperParameterGridJoint


(* ::Subsubsection:: *)
(*Return the best hyperparameter*)


BestJoinDeconv[]:=SortBy[Map[Function[{#[[1]],#[[2]],ComputeError[JointDeconv[#[[1]],# [[2]]]]}],hyperParameterGridJoint],#[[3]]&][[1]]
BestSequentialDeconv[]:=SortBy[Map[Function[{#[[1]],#[[2]],ComputeError[SequentialDeconv[#[[1]],# [[2]]]]}],hyperParameterGridSequential],#[[3]]&][[1]]


(* ::Subsection:: *)
(*Complete with noise and seed*)


BestJoinDeconv[c_,noise_]:=TODO


(* ::Chapter:: *)
(*"Manual" usage examples*)


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
