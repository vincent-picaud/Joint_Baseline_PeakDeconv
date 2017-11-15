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
hyperParameterGridMU=Range[0,3200,200]
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
resJoincPos=RealComputationJoinDeconv[1,noiseGrid];
resSeqcPos=RealComputationSequentialDeconv[1,noiseGrid];

resJoincNeg=RealComputationJoinDeconv[-1,noiseGrid];
resSeqcNeg=RealComputationSequentialDeconv[-1,noiseGrid];
*)


(* ::Subchapter:: *)
(*Here only the saved result!*)


resJoin={{-1,0.`,{{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`},{800,0.25`,0.027626726056173532`}}},{-1,0.2`,{{600,0.25`,0.054785440597964724`},{1000,0.25`,0.03384211145430906`},{600,0.25`,0.05100527010341116`},{400,0.25`,0.04955803247260004`},{600,0.25`,0.03608135426464135`},{800,0.25`,0.0459452240690661`},{400,0.25`,0.040160347749914206`},{600,0.25`,0.060506055203091685`},{1000,0.25`,0.04608367348230634`},{1200,0.25`,0.033725619849677965`}}},{-1,0.4`,{{800,0.25`,0.09707313001695114`},{1800,0.25`,0.062271816737100985`},{1000,0.25`,0.08541113640149152`},{400,0.25`,0.08927608557079532`},{800,0.25`,0.07840725830744691`},{400,0.25`,0.07852547647001681`},{600,0.25`,0.0711226053758534`},{1000,0.5`,0.10171251723930821`},{1600,0.25`,0.0757331443155274`},{1400,0.25`,0.07619443088779847`}}},{-1,0.6000000000000001`,{{800,0.25`,0.13803461545945248`},{1600,0.25`,0.0980419466933216`},{800,0.25`,0.129022752853032`},{1800,0.5`,0.1359171786403417`},{1000,0.25`,0.10346706668736244`},{2000,0.5`,0.11097050499090855`},{400,0.25`,0.09916556151693266`},{2400,0.75`,0.14478074760200738`},{1000,0.25`,0.09982907905126671`},{1600,0.25`,0.12206948107219015`}}},{-1,0.8`,{{1000,0.25`,0.18470577313002126`},{2200,0.5`,0.11780429656459289`},{2000,0.75`,0.16155897396541344`},{2400,1.`,0.17708834498745402`},{1800,0.5`,0.14340886952986776`},{600,0.25`,0.15122858073459886`},{1000,0.25`,0.14296745951137213`},{2400,0.75`,0.17406092588960115`},{1000,0.25`,0.14704411058248915`},{2400,0.5`,0.14367305885710374`}}},{-1,1.`,{{2400,0.5`,0.2181985660486051`},{1400,0.25`,0.15008336084200094`},{1600,0.75`,0.17538990512708652`},{1600,0.5`,0.21228604311827853`},{2000,0.5`,0.18524138904328544`},{1600,0.25`,0.18457865097657283`},{600,0.25`,0.17880189803501056`},{2400,0.75`,0.19837917854154688`},{1000,0.25`,0.1763796898960712`},{1600,0.5`,0.15982714910283177`}}},{-1,1.2000000000000002`,{{2000,0.75`,0.2526598871364748`},{1200,0.25`,0.17910216873791568`},{600,0.5`,0.1981631532554938`},{1800,0.5`,0.253409425806479`},{2400,0.5`,0.20964719965668505`},{1800,0.25`,0.22461383162490342`},{800,0.25`,0.2133391161794948`},{2400,0.75`,0.22715082068787673`},{1600,0.25`,0.20543848215723234`},{2400,0.5`,0.18305744186576053`}}},{-1,1.4000000000000001`,{{2200,0.75`,0.2826123097295538`},{1000,0.25`,0.2096546590038972`},{1600,0.75`,0.21505273194242944`},{2400,0.5`,0.29034403627938066`},{2400,0.5`,0.24296501963905998`},{1200,0.25`,0.26684446998854155`},{1000,0.5`,0.24821314925969812`},{2000,1.`,0.25402714150923567`},{2400,0.5`,0.2358673864205374`},{800,0.5`,0.2062983720808907`}}},{-1,1.6`,{{2400,1.25`,0.3139636075172902`},{2400,0.5`,0.2392912154216494`},{1800,0.75`,0.23368279076825033`},{2200,1.`,0.31708565485505513`},{2400,0.5`,0.29990954463091984`},{1400,0.25`,0.30705877439360324`},{1000,0.5`,0.2605750380220697`},{1800,1.`,0.27920869252467023`},{2400,0.5`,0.25742253631937195`},{1200,0.5`,0.24089188025011082`}}},{-1,1.8`,{{2400,1.`,0.32885303867684856`},{2400,0.5`,0.260184518976902`},{800,0.75`,0.25196222570639243`},{2400,0.5`,0.36546248867103714`},{2400,0.5`,0.3172874287579224`},{1200,0.25`,0.34534802706476025`},{2200,0.75`,0.2884152961194753`},{2000,1.`,0.3029483670887953`},{2400,0.5`,0.27801752239552524`},{1000,0.5`,0.265933883642268`}}},{-1,2.`,{{2400,1.25`,0.35418360379267283`},{2400,0.5`,0.2845505181153679`},{800,0.75`,0.26013935270430355`},{2200,1.25`,0.37960195181641476`},{2200,0.5`,0.34742192821567347`},{2000,0.5`,0.3921053477060484`},{1000,0.5`,0.2991192056381303`},{2200,1.`,0.3240426263624634`},{2400,0.75`,0.2984083654744468`},{2400,0.75`,0.28439162873156093`}}}};


resSeq={{-1,0.`,{{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`},{36,0.5`,0.015000487918700582`}}},{-1,0.2`,{{38,1.`,0.06265547840312533`},{38,0.75`,0.03831739070262936`},{32,0.5`,0.03568844928091362`},{34,0.75`,0.053731700792370884`},{36,1.25`,0.0463391016026236`},{38,0.75`,0.04646925947449819`},{34,0.5`,0.04836851763985121`},{36,0.75`,0.0560846870023022`},{36,0.75`,0.05281656776587868`},{40,0.5`,0.03779613214015582`}}},{-1,0.4`,{{36,1.5`,0.10941542651611813`},{38,1.`,0.07205216910860672`},{30,0.75`,0.059631128234133246`},{32,1.25`,0.10280382344499858`},{36,1.25`,0.09049327521450401`},{40,0.75`,0.08457116879411655`},{32,0.75`,0.07977339746926074`},{32,1.5`,0.10281174546051774`},{36,1.25`,0.09988817932256921`},{40,1.`,0.0843245045948035`}}},{-1,0.6000000000000001`,{{32,1.`,0.15917125138339439`},{38,1.25`,0.10804257725157995`},{28,1.25`,0.10000367776725184`},{30,1.5`,0.15507767056984448`},{34,1.5`,0.12518458493723839`},{34,1.25`,0.11010987633092244`},{32,1.`,0.12657819415924174`},{30,1.5`,0.1555164575118535`},{34,1.5`,0.1428255022669766`},{38,0.5`,0.12685279829120197`}}},{-1,0.8`,{{32,1.25`,0.21060502591492672`},{36,1.5`,0.15435699960287635`},{26,1.`,0.14587397541154548`},{30,1.75`,0.21625006919776316`},{32,1.75`,0.1787426184652192`},{32,1.25`,0.14301746240815222`},{30,1.`,0.16888690868531891`},{28,1.`,0.19331777017881507`},{30,1.75`,0.18710635593525665`},{30,1.`,0.1580584716781213`}}},{-1,1.`,{{28,1.75`,0.2537966058518411`},{38,1.75`,0.1943572720921377`},{24,1.25`,0.18208708211682118`},{26,1.75`,0.2740587761980877`},{32,2.`,0.20969688579691328`},{32,1.25`,0.18030393974961204`},{32,1.25`,0.21554463018960346`},{26,1.25`,0.23785341587320932`},{36,2.`,0.24152327255047285`},{28,1.25`,0.17911318035804571`}}},{-1,1.2000000000000002`,{{26,1.5`,0.2938104045081801`},{34,2.25`,0.23440944303724476`},{24,1.5`,0.2057959405749864`},{28,2.`,0.33526046590635933`},{32,2.`,0.28614233952707596`},{30,1.25`,0.2211716823480597`},{24,1.25`,0.27505247824448914`},{26,1.25`,0.28291749709632635`},{30,2.25`,0.2944136329903671`},{28,1.5`,0.19426519408493162`}}},{-1,1.4000000000000001`,{{30,3.`,0.3203017303107421`},{32,2.`,0.2673155050835572`},{22,2.25`,0.226745079332672`},{24,2.`,0.3769184584829791`},{32,2.25`,0.3168171018089292`},{36,1.5`,0.2771669304454587`},{24,1.5`,0.30895280061902336`},{26,3.`,0.3103032098755064`},{30,2.25`,0.3244545473857885`},{28,1.75`,0.21733011482178258`}}},{-1,1.6`,{{34,4.25`,0.33126350735101867`},{32,2.25`,0.30029791609488105`},{22,2.25`,0.2439184839105416`},{26,2.25`,0.43147201649305605`},{26,2.25`,0.3392001953264343`},{26,1.`,0.33253984761379624`},{24,1.5`,0.33790813886141075`},{26,2.75`,0.35711888326749436`},{30,2.25`,0.3590819324885512`},{26,1.75`,0.263365459681749`}}},{-1,1.8`,{{26,2.75`,0.3658157962389096`},{30,2.`,0.3324801410330364`},{20,1.75`,0.26424708813567266`},{26,2.25`,0.48200011655386277`},{24,2.25`,0.3679375107919133`},{24,0.5`,0.37068656352338103`},{22,2.`,0.3655333026108554`},{26,2.75`,0.38541804492020415`},{32,2.5`,0.3859665155624473`},{24,2.`,0.2835867344861921`}}},{-1,2.`,{{24,2.25`,0.38537172211079834`},{30,2.`,0.36813373720238285`},{20,1.75`,0.2767327937546347`},{22,2.25`,0.5341054257934323`},{28,2.5`,0.39715492627713445`},{26,1.5`,0.40037799556573384`},{24,2.`,0.3970528350634238`},{26,2.75`,0.42076706101128275`},{34,2.5`,0.4247563312119269`},{22,2.`,0.2996623669083737`}}}};


(* result format: a list of elements {c,noise,{{mu_or_snip,lambda_1,merit}, .... }} *)
FormatGlobalResult[resSeq_,resJoin_]:=Block[{i,n,assocError,plotError,assoc\[Lambda]1,plot\[Lambda]1,assocSNIP,plotSNIP,assocMu,plotMu,f,g1,g2},
(* assume c=constant *)
Assert[ConstantArrayQ @ Join[resSeq[[All,1]],resJoin[[All,1]]]];
(* same noise sample *)
Assert[resSeq[[All,2]]==resJoin[[All,2]]];

n=Length[resSeq];
Assert[n>0];

f[g_]:=Block[{},
assocError=<||>;
assoc\[Lambda]1=<||>;
assocSNIP=<||>;
assocMu=<||>;
For[i=1,i<=n,i++,
assocError[ToString[resJoin[[i,2]]]]=<|"Sequential"->g[resSeq[[i,3]][[All,3]]],"Joint"->g[resJoin[[i,3]][[All,3]]]|>;
assoc\[Lambda]1[ToString[resJoin[[i,2]]]]=<|"Sequential"->g[resSeq[[i,3]][[All,2]]],"Joint"->g[resJoin[[i,3]][[All,2]]]|>;
assocSNIP[ToString[resJoin[[i,2]]]]=<|"Sequential"->g[resSeq[[i,3]][[All,1]]]|>;
assocMu[ToString[resJoin[[i,2]]]]=<|"Joint"->g[resJoin[[i,3]][[All,1]]]|>;
];
];

g1[l_]:=l;
f[g1];

plotError=BoxWhiskerChart[assocError,ChartLegends->Placed[Automatic,{0.15,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(e\)]\)"}];
plot\[Lambda]1=BoxWhiskerChart[assoc\[Lambda]1,ChartLegends->Placed[Automatic,{0.15,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\)"}];
plotSNIP=BoxWhiskerChart[assocSNIP,ChartLegends->Placed[Automatic,{0.15,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\!\(\*SubscriptBox[\(m\), \(SNIP\)]\)"}];
plotMu=BoxWhiskerChart[assocMu,ChartLegends->Placed[Automatic,{0.15,0.85}],ChartLabels->{Automatic,None},FrameLabel->{"\!\(\*SubscriptBox[\(\[Sigma]\), \(noise\)]\)","\[Mu]"}];

g2[l_]:={N[Mean[l]],N[StandardDeviation[l]]};
f[g2];

Return[{assocError,assoc\[Lambda]1,plotError,plot\[Lambda]1,assocSNIP,plotSNIP,assocMu,plotMu}];
];


(* ::Text:: *)
(*Usage example: BoxWhiskerChart[FormatGlobalResult[resSeq,resJoin],ChartLegends->Automatic,ChartLabels->{Automatic,None}]*)


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
