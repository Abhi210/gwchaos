(* ::Package:: *)

(* ::Title:: *)
(*KerrSpinOrbit subpackage of KerrSpinningSecondary*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrSpinOrbit`",
	{"ConstantsOfMotion`",
	"KerrGeoOrbit`"}];

KerrSpinOrbit::usage = "KerrSpinOrbit[a,p,e,x] returns a KerrSpinOrbitFunction[..] which stores the orbital trajectory and parameters.";
KerrSpinOrbitFunction::usage = "KerrSpinOrbitFunction[a,p,e,x,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Circular (Mino time) (Linearized in spin)*)


KerrSpinOrbitMino[a_,p_,0,1,{\[Sigma]para_,\[Sigma]perp_,\[Phi]s_},initPhases:{_,_,_,_}:{0,0,0,0}]:=Module[{M=1,consts,En,L,Q,assoc,v,A,B,\[CapitalUpsilon]s,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[Alpha]1,qt0,qr0,qz0,q\[Phi]0,t,r,\[Theta],\[Phi]},
En=KerrSpinEnergy[a,p,0,1,\[Sigma]para] ;
L=KerrSpinEnergy[a,p,0,1,\[Sigma]para] ;
Q=-2 a \[Sigma]para;
consts ={En,L,Q};
	
v=Sqrt[1/p];
\[CapitalUpsilon]\[Phi] =1/v Sqrt[1/(1-3v^2+2a v^3)]-((3 \[Sigma]para)/2) ((v^2) (1-a v)(1-2 v^2+a v^3) )/(1-3 v^2+2a v^3)^(3/2);
\[CapitalUpsilon]t = (1+a v^3)/(v^4 Sqrt[1-3 v^4+2 a v^3]) -(3 \[Sigma]para)/2 (v(1-a v)(1-2 a v^3+a^2 v^4))/(1-3 v^2+2 a v^3)^(3/2);
\[CapitalUpsilon]s = 1/v;
\[CapitalUpsilon]\[Theta] = (1/v) Sqrt[(1-4a v^3+3a^2 v^4)/(1-3 v^2 +2a v^3)];
\[Alpha]1=-(v (1-a v)Sqrt[(1-2v^2+a^2 v^4)])/Sqrt[1-3 v^2+2a v^3] ;
A[\[Lambda]_]=-((3\[Sigma]perp \[Alpha]1)/(2\[CapitalUpsilon]\[Theta]))(Sin[\[Phi]s+(\[CapitalUpsilon]s-\[CapitalUpsilon]\[Theta])\[Lambda]]/(\[CapitalUpsilon]s-\[CapitalUpsilon]\[Theta])+Sin[\[Phi]s+(\[CapitalUpsilon]s+\[CapitalUpsilon]\[Theta])\[Lambda]]/(\[CapitalUpsilon]s+\[CapitalUpsilon]\[Theta]));
B[\[Lambda]_]=(3\[Sigma]perp \[Alpha]1)/(2\[CapitalUpsilon]\[Theta]) (Cos[\[Phi]s+(\[CapitalUpsilon]s-\[CapitalUpsilon]\[Theta])\[Lambda]]/(\[CapitalUpsilon]s-\[CapitalUpsilon]\[Theta])+Cos[\[Phi]s+(\[CapitalUpsilon]s+\[CapitalUpsilon]\[Theta])\[Lambda]]/(\[CapitalUpsilon]s+\[CapitalUpsilon]\[Theta]));


t=Function[{Global`\[Lambda]}, Evaluate[\[CapitalUpsilon]t Global`\[Lambda] ], Listable];
r=Function[{Global`\[Lambda]}, Evaluate[p ], Listable];
\[Theta]=Function[{Global`\[Lambda]}, Evaluate[\[Pi]/2 + A[Global`\[Lambda]] Sin[\[CapitalUpsilon]\[Theta] Global`\[Lambda]]+B[Global`\[Lambda]] Cos[\[CapitalUpsilon]\[Theta] Global`\[Lambda]] ], Listable];
\[Phi]=Function[{Global`\[Lambda]}, Evaluate[ \[CapitalUpsilon]\[Phi] Global`\[Lambda]  ], Listable];


	assoc = Association[
	"a" -> a,
	"p" -> p,
	"e" -> 0,
	"Inclination" -> 1,
	"Parametrization"->"Mino", 
	"Energy" -> En, 
	"AngularMomentum" -> L, 
	"CarterConstant" -> Q, 
	"ConstantsOfMotion" -> consts,
	"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];

	KerrSpinOrbitFunction[a,p,0,1,{\[Sigma]para,\[Sigma]perp,\[Phi]s},assoc]

]


KerrSpinOrbitPhases[a_,p_,0,1,{\[Sigma]para_,\[Sigma]perp_,\[Phi]s_},initPhases:{_,_,_,_}:{0,0,0,0}]:=Module[{M=1,consts,En,L,Q,assoc,v,A,B,\[CapitalUpsilon]s,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[Alpha]1,qt0,qr0,qz0,q\[Phi]0,t,r,\[Theta],\[Phi]},
En=KerrSpinEnergy[a,p,0,1,\[Sigma]para] ;
L=KerrSpinEnergy[a,p,0,1,\[Sigma]para] ;
Q=-2 a \[Sigma]para;
consts ={En,L,Q};
	
v=Sqrt[1/p];
\[CapitalUpsilon]\[Phi] =1/v Sqrt[1/(1-3v^2+2a v^3)]-((3 \[Sigma]para)/2) ((v^2) (1-a v)(1-2 v^2+a v^3) )/(1-3 v^2+2a v^3)^(3/2);
\[CapitalUpsilon]t = (1+a v^3)/(v^4 Sqrt[1-3 v^4+2 a v^3]) -(3 \[Sigma]para)/2 (v(1-a v)(1-2 a v^3+a^2 v^4))/(1-3 v^2+2 a v^3)^(3/2);
\[CapitalUpsilon]s = 1/v;
\[CapitalUpsilon]\[Theta] = (1/v) Sqrt[(1-4a v^3+3a^2 v^4)/(1-3 v^2 +2a v^3)];
\[Alpha]1=-(v (1-a v)Sqrt[(1-2v^2+a^2 v^4)])/Sqrt[1-3 v^2+2a v^3] ;
A[q\[Theta]_,qs_]:=-((3\[Sigma]perp \[Alpha]1)/(2\[CapitalUpsilon]\[Theta]))(Sin[\[Phi]s+(qs-q\[Theta])]/(\[CapitalUpsilon]s-\[CapitalUpsilon]\[Theta])+Sin[\[Phi]s+(qs-q\[Theta])]/(\[CapitalUpsilon]s+\[CapitalUpsilon]\[Theta]));
B[q\[Theta]_,qs_]:=(3\[Sigma]perp \[Alpha]1)/(2\[CapitalUpsilon]\[Theta]) (Cos[\[Phi]s+(qs-q\[Theta])]/(\[CapitalUpsilon]s-\[CapitalUpsilon]\[Theta])+Cos[\[Phi]s+(qs-q\[Theta])]/(\[CapitalUpsilon]s+\[CapitalUpsilon]\[Theta]));


t=Function[{Global`qt}, Evaluate[ Global`qt], Listable];
r=Function[{Global`q\[Theta]}, Evaluate[p ], Listable];
\[Theta]=Function[{Global`q\[Theta],Global`qs}, Evaluate[\[Pi]/2 + A[Global`q\[Theta],Global`qs] Sin[Global`q\[Theta]]+B[Global`q\[Theta],Global`qs] Cos[Global`q\[Theta]] ], Listable];
\[Phi]=Function[{Global`q\[Phi]}, Evaluate[ Global`q\[Phi]  ], Listable];


	assoc = Association[
	"a" -> a,
	"p" -> p,
	"e" -> 0,
	"Inclination" -> 1,
	"Parametrization"->"Phases", 
	"Energy" -> En, 
	"AngularMomentum" -> L, 
	"CarterConstant" -> Q, 
	"ConstantsOfMotion" -> consts,
	"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];

	KerrSpinOrbitFunction[a,p,0,1,{\[Sigma]para,\[Sigma]perp,\[Phi]s},assoc]

]




(* ::Subsection:: *)
(*Circular (Mino time) (Not linearized in spin)*)


(* ::Subsection::Closed:: *)
(*Equatorial eccentric aligned (Darwin)*)


(* ::Text:: *)
(*Eccentric equatorial orbit with aligned spin calculated using DCT with code from KerrGeodesics package (Hopper, Forseth, Osburn, and Evans, PRD 92 (2015))*)


KerrSpinOrbitDarwinNonlinearized[a_, p_, e_, x_/;x^2==1, \[Sigma]_, initPhases:{_,_,_,_}:{0,0,0,0}] := 
Module[{M=1,En,Jz,\[Xi],\[CapitalSigma]s,\[CapitalDelta],Ps,jp,je,assoc,var,t0, \[Chi]0, \[Phi]0,r0,\[Theta]0,t,r,\[Theta],\[Phi],\[Chi],growthRateT,growthRatePh,\[CapitalOmega]r,\[CapitalOmega]\[Phi],
		\[Chi]r,NrMax,pg,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[Phi]C,tC,Pr,r0Sample,PrSample,dtd\[Chi],d\[Phi]d\[Chi],TVr,PVr},
	En = KerrSpinEnergy[a,p,e,x,\[Sigma],"Linear"->False];
	Jz = KerrSpinAngularMomentum[a,p,e,x,\[Sigma],"Linear"->False];
    \[Xi] = Jz-(a+\[Sigma])*En;
	
	(*Precision of sampling depends on precision of arguments*)
	pg=Min[{Precision[{a,p,e,x,\[Sigma]}],Precision[initPhases]}];
	
	(*Parameterize r in terms of Darwin parameter \[Chi]*)
	r0[chi_]:=p M/(1+e Cos[chi]);
	\[Theta]0[chi_?NumericQ]:=N[Pi/2,pg];
	\[Theta]0[chi_List]:=\[Theta]0[#]&/@chi;
	
	(* Expressions for dt/d\[Lambda] = TVr and d\[Phi]/d\[Lambda] = PVr *)
	
    \[CapitalSigma]s[rp_]:=rp^2*(1-M*\[Sigma]^2/rp^3);
    \[CapitalDelta][rp_]:=rp^2-2rp+a^2;
    Ps[rp_]:=((rp^2+a^2)+a*\[Sigma]/rp*(rp+M))*En-(a+M*\[Sigma]/rp)*Jz;
  
	TVr[rp_]:=a*(1+3M*\[Sigma]^2/rp/\[CapitalSigma]s[rp])*\[Xi]+(rp^2+a^2)/\[CapitalDelta][rp]*Ps[rp];
	PVr[rp_]:=(1+3M*\[Sigma]^2/rp/\[CapitalSigma]s[rp])*\[Xi]+a/\[CapitalDelta][rp]*Ps[rp];
	
	(* Sampling of radial position using evenly spaced values of Darwin parameter \[Chi]*)
	If[pg==$MachinePrecision,
		\[Chi]r[Nr_]:=N[Table[i Pi/(Nr-1),{i,0,Nr-1}]],
		\[Chi]r[Nr_]:=N[Table[i Pi/(Nr-1),{i,0,Nr-1}],1.5pg]
	];
	r0Sample[Nr_]:=r0[\[Chi]r[Nr]];
	
	(* Expression for d\[Lambda]/d\[Chi]*)
	jp = { 1 - En^2, -2, a^2 + 2*a*En*\[Xi] + \[Xi]^2, -2*( ( 1 - En^2 )*\[Sigma]^2 - En*\[Sigma]*\[Xi] + \[Xi]^2 ), 4*\[Sigma]^2, -2*a*\[Sigma]*( a*\[Sigma] + \[Xi]*( En*\[Sigma] + \[Xi] ) ), -\[Sigma]^2*( -( 1 - En )*\[Sigma] + \[Xi] )*( ( 1 + En )*\[Sigma] + \[Xi] ) }; (* coefficients for J(\[Psi]) *)
	je = { 1, 2, e^2 + 3, 4*( e^2 + 1), e^4 + 10*e^2 + 5, 2*( e^2 + 3 )*( 3*e^2 + 1 ), e^6 + 21*e^4 + 35*e^2 + 7 };(* coefficients for J(\[Psi]) *)
	Pr[\[Chi]_]:=Sqrt[1-e^2]/p/Sqrt[Sum[ Sum[ jp[[l+1]]*je[[k-l+1]]/(1 - e^2)^(k - l)/p^l,{l,0,k}]*(1 + e*Cos[\[Chi]])^k, {k,0,6}]];
	PrSample[Nr_]:=Pr[\[Chi]r[Nr]];
	
	(* Sampling of expressions for dt/d\[Chi] and d\[Phi]/d\[Chi] *)
	dtd\[Chi][Nr_]:=TVr[r0Sample[Nr]]PrSample[Nr];
	d\[Phi]d\[Chi][Nr_]:=PVr[r0Sample[Nr]]PrSample[Nr];
	
	(*Spectral integration of t and \[Phi] as functions of \[Chi]*)
	{growthRateT,\[CapitalDelta]tr}=KerrGeodesics`KerrGeoOrbit`Private`DarwinFastSpecIntegrateAndConvergenceCheck[dtd\[Chi]];
	{growthRatePh,\[CapitalDelta]\[Phi]r}=KerrGeodesics`KerrGeoOrbit`Private`DarwinFastSpecIntegrateAndConvergenceCheck[d\[Phi]d\[Chi]];
	
	(*Collect initial phases*)
	{t0, \[Chi]0, \[Phi]0} = {initPhases[[1]], initPhases[[2]], initPhases[[4]]};
	(* Find integration constants for t and \[Phi], so that t(\[Chi]=0)=t0 and \[Phi](\[Chi]=0)=\[Phi]0 *)
	\[Phi]C=\[CapitalDelta]\[Phi]r[\[Chi]0];
	tC=\[CapitalDelta]tr[\[Chi]0];
	
	t[\[Chi]_]:=\[CapitalDelta]tr[\[Chi]+\[Chi]0]+growthRateT \[Chi]+t0-tC;
	r[\[Chi]_]:=r0[\[Chi]+\[Chi]0];
	\[Theta][\[Chi]_]:=\[Theta]0[\[Chi]];
	\[Phi][\[Chi]_]:=\[CapitalDelta]\[Phi]r[\[Chi]+\[Chi]0]+growthRatePh \[Chi]+\[Phi]0-\[Phi]C;
	
	\[CapitalOmega]r = 1/growthRateT;
	\[CapitalOmega]\[Phi] = growthRatePh/growthRateT;
	
	assoc = Association[
		"Parametrization"->"Darwin",
		"Energy" -> En, 
		"AngularMomentum" -> Jz, 
		"ConstantsOfMotion" -> {"\[ScriptCapitalE]"->En,"\!\(\*SubscriptBox[\(\[ScriptCapitalJ]\), \(z\)]\)"->Jz},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"RadialFrequency" -> \[CapitalOmega]r,
		"AzimuthalFrequency" -> \[CapitalOmega]\[Phi],
		"a" -> a,
		"p" -> p,
		"e" -> e,
		"Inclination" -> x,
		"\[Sigma]" -> \[Sigma]
		];
		
	KerrSpinOrbitFunction[a,p,e,1,\[Sigma],assoc]
	
]


(* ::Section:: *)
(*KerrSpinOrbit and KerrSpinOrbitFunction*)


Options[KerrSpinOrbit] = {"Parametrization" -> "Mino","Linear"->True};
SyntaxInformation[KerrSpinOrbit] = {"ArgumentsPattern"->{_,OptionsPattern[]}};


KerrSpinOrbit[a_,p_,e_,x_,{\[Sigma]para_,\[Sigma]perp_,\[Phi]s_}, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param,lin},

param = OptionValue["Parametrization"];
lin = OptionValue["Linear"];

If[lin,

If[param=="Mino",Return[KerrSpinOrbitMino[a, p, e, x,{\[Sigma]para,\[Sigma]perp,\[Phi]s}, initPhases]]];
If[param=="Phases",  Return[KerrSpinOrbitPhases[a, p, e, x,{\[Sigma]para,\[Sigma]perp,\[Phi]s}, initPhases]]];

]

];

KerrSpinOrbit[a_,p_,e_,x_,\[Sigma]_, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param},

	param = OptionValue["Parametrization"];
	If[!OptionValue["Linear"],
		If[param=="Darwin",Return[KerrSpinOrbitDarwinNonlinearized[a, p, e, x, \[Sigma], initPhases]]];
	];

];


KerrSpinOrbitFunction /:
 MakeBoxes[kgof:KerrSpinOrbitFunction[a_, p_, e_, x_,{\[Sigma]para_,\[Sigma]perp_,\[Phi]s_}, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"p: ", p}], "  ",
                  BoxForm`SummaryItem[{"e: ", e}], "  ",
                  BoxForm`SummaryItem[{"x: ", x}]}],
             BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}]};
  extended = {BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}],
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["AngularMomentum"]}],
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["CarterConstant"]}]};
  BoxForm`ArrangeSummaryBox[
    KerrSpinOrbitFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];

KerrSpinOrbitFunction /:
 MakeBoxes[kgof:KerrSpinOrbitFunction[a_, p_, e_, x_, \[Sigma]_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"p: ", p}], "  ",
                  BoxForm`SummaryItem[{"e: ", e}], "  ",
                  BoxForm`SummaryItem[{"x: ", x}]}],
             BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}]};
  extended = {BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}],
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["AngularMomentum"]}](*,
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["CarterConstant"]}]*)};
  BoxForm`ArrangeSummaryBox[
    KerrSpinOrbitFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];


KerrSpinOrbitFunction[a_, p_, e_, x_, {\[Sigma]para_,\[Sigma]perp_,\[Phi]s_},assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrSpinOrbitFunction[a_, p_, e_, x_,{\[Sigma]para_,\[Sigma]perp_,\[Phi]s_}, assoc_][y_?StringQ] := assoc[y]

KerrSpinOrbitFunction[a_, p_, e_, x_, \[Sigma]_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrSpinOrbitFunction[a_, p_, e_, x_, \[Sigma]_, assoc_][y_?StringQ] := assoc[y]

Keys[g_KerrSpinOrbitFunction]^:=Keys[g[[5]]]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
