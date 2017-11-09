#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Comparison script*)


(* ::Text:: *)
(*Script used to generate comparison (included in the associated article)*)


(* ::Subchapter:: *)
(*Project root directory & other initialization*)


SetDirectory["~/GitHub/Joint_Baseline_PeakDeconv"];
On[Assert];


(* ::Subchapter:: *)
(*Some misc function*)


(* ::Subsection:: *)
(*Generate synthetic data*)


GenerateData[c_?IntegerQ,noise_,seed_?IntegerQ]:=Block[{errorCode},
	Assert[Abs[c]<=1]; (* c={-1,0,1} *)
	Assert[noise>=0];
	Assert[seed>=0];
	errorCode=Run["./generateSynthetic -c "<>ToString[c]<>" -n "<>ToString[noise]<>" -s "<>ToString[seed]<>" > data.txt"];
	If[errorCode!=0,Return[$Failed]];
];


GenerateData[-1,0.2,40]


(* ::Subsection:: *)
(*Run a computation and return result stored into an association*)


ReadResult[filename_?StringQ]:=Block[{tmp},
	If[!FileExistsQ[filename],Return[$Failed]];
	tmp=Transpose[Import[filename,"csv"][[2;;-1]]];
	Assert[Length[tmp]==6];
	Return[AssociationThread[{"x","y","peakDirac","baseline","peakConv","peakBaseline"}->tmp]];
];


ResultToPeakList[<|"x"->x_,"y"->y_,"peakDirac"->peakDirac_,"baseline"->baseline_,"peakConv"->peakConv_,"peakBaseline"->peakBaseline_|>]:=Transpose[Select[Transpose[{x,peakDirac}],#[[2]]>0&]]


PlotResult[<|"x"->x_,"y"->y_,"peakDirac"->peakDirac_,"baseline"->baseline_,"peakConv"->peakConv_,"peakBaseline"->peakBaseline_|>]:=ListLinePlot[{Transpose@{x,y},Transpose@{x,peakDirac},Transpose@{x,peakBaseline},Transpose@{x,baseline}}]


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


res=SequentialDeconv[40,2];
PlotResult[res]


res=SequentialDeconv[40,2];
PlotResult[res]


res2=JointDeconv[2000,0.8];
PlotResult[res2]


res2=JointDeconv[2000,0.8];
PlotResult[res2]


 





Dimensions @ Transpose[Import["data3.txt.out","csv"][[2;;-1]]]


ReadResult["data3.txt.out"]


(* ::Subchapter:: *)
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
