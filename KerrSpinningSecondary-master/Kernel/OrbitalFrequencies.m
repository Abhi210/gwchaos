(* ::Package:: *)

(* ::Title:: *)
(*OrbitalFrequencies subpackage of KerrSpinningSecondary*)


(* ::Section::Closed:: *)
(*Define usage for public functions*)


BeginPackage["KerrSpinningSecondary`OrbitalFrequencies`",
	{"KerrGeodesics`ConstantsOfMotion`",
	"KerrSpinningSecondary`ConstantsOfMotion`",
	"KerrGeodesics`OrbitalFrequencies`"}];

KerrSpinFrequencyCorrections::usage = "KerrSpinFrequencyCorrections[a, p, e, x] returns the spin corrections to geodesic orbital frequencies."
KerrSpinFrequencies::usage = "KerrSpinFrequencies[a, p, e, x, \[Sigma]] returns the orbital frequencies."

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Orbital Frequencies*)


(* ::Subsection:: *)
(*Schwarzschild*)


(* ::Subsection::Closed:: *)
(*Kerr*)


(* ::Subsubsection::Closed:: *)
(*Circular*)


KerrSpinMinoFrequenciesCorrectionTulczyjew[a_,p_,0,1,\[Sigma]_]:=Module[
	{M=1,v,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

	v=Sqrt[1/p];
    \[CapitalUpsilon]\[Phi] =1/v Sqrt[1/(1-3v^2+2a v^3)]-((3 \[Sigma])/2) ((v^2) (1-a v)(1-2 v^2+a v^3) )/(1-3 v^2+2a v^3)^(3/2);
    \[CapitalGamma] = (1+a v^3)/(v^4 Sqrt[1-3 v^4+2 a v^3]) -(3 \[Sigma])/2 (v(1-a v)(1-2 a v^3+a^2 v^4))/(1-3 v^2+2 a v^3)^(3/2);


	<| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
		"\[CapitalGamma]" -> \[CapitalGamma] 
	|>

]

KerrSpinMinoFrequenciesGeoPlusLinearTulczyjew[a_,p_,0,1,\[Sigma]_]:=KerrGeoAngularMomentum[a,p,0,1]+\[Sigma]*KerrSpinMinoFrequenciesCorrectionTulczyjew[a,p,0,1,\[Sigma]]


(* ::Subsubsection::Closed:: *)
(*Eccentric equatorial aligned*)


KerrSpinMinoFrequencyCorrectionsTulczyjew[a_?NumericQ,p_?NumericQ,e_?NumericQ,1]:=Module[{M=1,\[CapitalSigma],En,EnG,L,LG,RSKerr,\[CapitalUpsilon]rG,\[CapitalLambda]rG,r,drd\[Chi]r,OneOverRs,d\[Lambda]d\[Chi]r,d\[Lambda]d\[Chi]rexpand,d\[Lambda]d\[Chi]r\[Delta]S,\[Delta]\[CapitalLambda]r,\[Delta]\[CapitalUpsilon]r,\[CapitalGamma],\[Sigma]t,udt,EKillingTerm,LKillingTerm,ud\[Phi],d\[Phi]d\[Tau],d\[Phi]d\[Lambda],Integrand\[Phi],Integrand\[Phi]S,\[Delta]\[CapitalUpsilon]\[Phi], dtd\[Tau],dtd\[Lambda],Integrandt,IntegrandtS,\[Delta]\[CapitalGamma]},
    En=KerrSpinEnergy[a,p,e,1,\[Sigma]t,Linear->True] ;
    L=KerrSpinAngularMomentum[a,p,e,1,\[Sigma]t,Linear->True] ;
    EnG=KerrGeoEnergy[a,p,e,1] ;
    LG=KerrGeoAngularMomentum[a,p,e,1] ;
    RSKerr[r_]=  1/r (2 a^3 \[Sigma]t En^2+a^2 (2 r^2 En^2+r^3 (-1+En^2)-4 \[Sigma]t En L)+r^2 (2 r^2+r^3 (-1+En^2)-r L (-2 \[Sigma]t En+L)+2 L (-3 \[Sigma]t En+L))-2 a (-\[Sigma]t L^2+r^2 En (-3 \[Sigma]t En+2 L)));
    \[CapitalUpsilon]rG=KerrGeoFrequencies[a,p,e,1,Time->"Mino"][[1]];
    \[CapitalLambda]rG=(2 Pi)/\[CapitalUpsilon]rG;
    r[\[Chi]r_]:= p/(1 + e Cos[\[Chi]r]); 
    drd\[Chi]r[\[Chi]r_]:=Simplify[ D[r[\[Chi]r],\[Chi]r]];
    OneOverRs[\[Chi]r_]:=1/Sqrt[RSKerr[r[\[Chi]r]]];
    d\[Lambda]d\[Chi]r[\[Chi]r_]:= drd\[Chi]r[\[Chi]r]OneOverRs[\[Chi]r];
    d\[Lambda]d\[Chi]rexpand[\[Chi]r_]:=Normal[Series[d\[Lambda]d\[Chi]r[\[Chi]r],{\[Sigma]t,0,1}]];
    d\[Lambda]d\[Chi]r\[Delta]S[\[Chi]r_]=Coefficient[d\[Lambda]d\[Chi]rexpand[\[Chi]r],\[Sigma]t];
    \[Delta]\[CapitalLambda]r=2 NIntegrate[d\[Lambda]d\[Chi]r\[Delta]S[\[Chi]r],{\[Chi]r,0,\[Pi]}];
    \[Delta]\[CapitalUpsilon]r=-((2 \[Pi] \[Delta]\[CapitalLambda]r)/\[CapitalLambda]rG^2);
    
    udt[r_]:=-En+EKillingTerm[r]; ud\[Phi][r_]:=L+LKillingTerm[r];
    EKillingTerm[r_]:=-(((a EnG-LG)\[Sigma]t)/r^3);
    LKillingTerm[r_]:=-(((-a^2 EnG +a LG+EnG r^3) \[Sigma]t)/r^3);
    d\[Phi]d\[Tau][r_]:=(-((2 a udt[r])/(r (a^2+(-2+r) r)))+((-2+r) ud\[Phi][r])/(r (a^2+(-2+r) r))); 
    d\[Phi]d\[Lambda][r_]=\[CapitalSigma][r]d\[Phi]d\[Tau][r];
    \[CapitalSigma][r_]:=r^2;
    Integrand\[Phi][r_]=Normal[Series[(d\[Phi]d\[Lambda][r]drd\[Chi]r[\[Chi]r])/(Sqrt[RSKerr[r]](\[CapitalLambda]rG+ \[Sigma]t \[Delta]\[CapitalLambda]r)),{\[Sigma]t,0,1}]];
    Integrand\[Phi]S[r_]=Coefficient[Integrand\[Phi][r],\[Sigma]t];
    \[Delta]\[CapitalUpsilon]\[Phi]=2NIntegrate[Integrand\[Phi]S[r[\[Chi]r]] ,{\[Chi]r,0,\[Pi]}];
    
    dtd\[Tau][r_]=(-(((-a^2 (a^2+(-2+r) r)+(a^2+r^2)^2) udt[r])/(r^2 (a^2+(-2+r) r)))-(2 a ud\[Phi][r])/(r (a^2+(-2+r) r)));
    dtd\[Lambda][r_]=\[CapitalSigma][r]dtd\[Tau][r];
    Integrandt[r_]=Normal[Series[dtd\[Lambda][r]/(Sqrt[RSKerr[r]](\[CapitalLambda]rG+\[Sigma]t \[Delta]\[CapitalLambda]r)),{\[Sigma]t,0,1}]];
    IntegrandtS[r_]=Coefficient[Integrandt[r],\[Sigma]t];
    \[Delta]\[CapitalGamma]=2NIntegrate[IntegrandtS[r[\[Chi]r]] drd\[Chi]r[\[Chi]r],{\[Chi]r,0,\[Pi]}];

		<| "\!\(\*SubscriptBox[\(\[Delta]\[CapitalUpsilon]\), \(r\)]\)" -> \[Delta]\[CapitalUpsilon]r,
		"\!\(\*SubscriptBox[\(\[Delta]\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[Delta]\[CapitalUpsilon]\[Phi],
		"\[Delta]\[CapitalGamma]" -> \[Delta]\[CapitalGamma] 
	|>

]


KerrSpinBoyerLindquistFrequencyCorrectionsTulczyjew[a_?NumericQ,p_?NumericQ,e_?NumericQ,1]:=Module[{M=1,J,\[CapitalDelta],d\[CapitalDelta],Ps,Vt,V\[Phi],rp,En,Jz,Q,dEnds,dJzds,x,dxds,
    JDerivatives,dItds,dI\[Phi]ds,\[CapitalOmega]r,\[CapitalOmega]\[Phi],d\[CapitalOmega]rds,d\[CapitalOmega]\[Phi]ds,It,I\[Phi]},
  J[\[Chi]_]:=1-En^2+2*((1-En^2)/(1-e^2)-1/p)*(1+e*Cos[\[Chi]]);
  \[CapitalDelta][r_]:=r^2-2r+a^2;
  d\[CapitalDelta][r_]:=2r-2;
  Ps[r_]:=r^2*En-a*x;
  Vt[r_]:=a*x+(r^2+a^2)/\[CapitalDelta][r]*Ps[r];
  V\[Phi][r_]:=x+a/\[CapitalDelta][r]*Ps[r];
  rp[\[Chi]_]:=p/(1+e*Cos[\[Chi]]);
  JDerivatives[chi_]:=Module[{jp,je,djpds,djpdEn,djpdx,dJds,dJdEn,dJdx},
    jp={1-En^2,-2,a^2+2 a En x+x^2,-2 x^2,0,0,0};
    je={1,2,e^2+3,4(e^2+1),e^4+10e^2+5,2(e^2+3)(3e^2+1),e^6+21e^4+35e^2+7};
    djpds={0,0,0,2 En x,0,-2 a x^2,0};
    djpdEn={-2 En,0,2 a x,0,0,0,0};
    djpdx={0,0,2 a En+2 x,-4 x,0,0,0};
    dJds=Sum[Sum[djpds[[l+1]]*je[[k-l+1]]/(1-e^2)^(k-l)/p^l,{l,3,k}]*(1+e*Cos[chi])^k,{k,3,6}];
    dJdEn=Sum[Sum[djpdEn[[l+1]]*je[[k-l+1]]/(1-e^2)^(k-l)/p^l,{l,0,k}]*(1+e*Cos[chi])^k,{k,0,6}];
    dJdx=Sum[Sum[djpdx[[l+1]]*je[[k-l+1]]/(1-e^2)^(k-l)/p^l,{l,2,k}]*(1+e*Cos[chi])^k,{k,0,6}];
    <|
      "dJds"->dJds,
      "dJdEn"->dJdEn,
      "dJdx"->dJdx
    |>
  ];
  {En, Jz, Q} = Values[KerrGeoConstantsOfMotion[a,p,e,1]];
  dEnds = KerrSpinEnergyCorrection[a,p,e,1,"SSC"->"Tulczyjew"];
  dJzds = KerrSpinAngularMomentumCorrection[a,p,e,1,"SSC"->"Tulczyjew"];
  x=Jz-a*En;
  dxds=dJzds-a*dEnds-En;
  dItds = 2*Sqrt[1-e^2]/p*NIntegrate[(-((x (a^2+rp[\[Chi]]^2))/(rp[\[Chi]] \[CapitalDelta][rp[\[Chi]]]))+(rp[\[Chi]]^2 (a^2+rp[\[Chi]]^2))/\[CapitalDelta][rp[\[Chi]]]*dEnds+(a-(a (a^2+rp[\[Chi]]^2))/\[CapitalDelta][rp[\[Chi]]])*dxds)/Sqrt[J[\[Chi]]]-Vt[rp[\[Chi]]]/2/(J[\[Chi]])^(3/2)*(Values[JDerivatives[\[Chi]]] . {1,dEnds,dxds}),{\[Chi],0,Pi}];
  dI\[Phi]ds = 2*Sqrt[1-e^2]/p*NIntegrate[(-((a x)/(rp[\[Chi]] \[CapitalDelta][rp[\[Chi]]]))+(a rp[\[Chi]]^2)/\[CapitalDelta][rp[\[Chi]]]*dEnds+(1-a^2/\[CapitalDelta][rp[\[Chi]]])*dxds)/Sqrt[J[\[Chi]]]-V\[Phi][rp[\[Chi]]]/2/(J[\[Chi]])^(3/2)*(Values[JDerivatives[\[Chi]]] . {1,dEnds,dxds}),{\[Chi],0,Pi}];
  {\[CapitalOmega]r,\[CapitalOmega]\[Phi]} = Map[KerrGeoFrequencies[a,p,e,1],{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"}];
  It = 2*Pi/\[CapitalOmega]r;
  I\[Phi] = 2*Pi*\[CapitalOmega]\[Phi]/\[CapitalOmega]r;
  d\[CapitalOmega]rds = -\[CapitalOmega]r/It*dItds;
  d\[CapitalOmega]\[Phi]ds = \[CapitalOmega]\[Phi]*(1/I\[Phi]*dI\[Phi]ds-1/It*dItds);
  <|
    "\!\(\*SubscriptBox[\(\[Delta]\[CapitalOmega]\), \(r\)]\)"->d\[CapitalOmega]rds,
    "\!\(\*SubscriptBox[\(\[Delta]\[CapitalOmega]\), \(\[Phi]\)]\)"->d\[CapitalOmega]\[Phi]ds
  |>
]


KerrSpinMinoFrequenciesGeoPlusLinearTulczyjew[a_?NumericQ,p_?NumericQ,e_?NumericQ,1,\[Sigma]_?NumericQ]:=Module[{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},
\[CapitalUpsilon]r=KerrGeoFrequencies[a,p,e,1,Time->"Mino"][[1]]+\[Sigma]*KerrSpinMinoFrequencyCorrectionsTulczyjew[a,p,e,1][[1]];
\[CapitalUpsilon]\[Phi]=KerrGeoFrequencies[a,p,e,1,Time->"Mino"][[3]]+\[Sigma]*KerrSpinMinoFrequencyCorrectionsTulczyjew[a,p,e,1][[2]];
\[CapitalGamma]=KerrGeoFrequencies[a,p,e,1,Time->"Mino"][[4]]+\[Sigma]*KerrSpinMinoFrequencyCorrectionsTulczyjew[a,p,e,1][[3]];
		<| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
		"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
		"\[CapitalGamma]" -> \[CapitalGamma] 
	|>
]


KerrSpinMinoFrequenciesNonlinearTulczyjew[a_?NumericQ,p_?NumericQ,e_?NumericQ,1,\[Sigma]_?NumericQ]:=Module[
	{M=1,En,Jz,x,jp,je,J,\[CapitalSigma]s,\[CapitalDelta],Ps,Vt,V\[Phi],\[Chi],\[CapitalLambda]r,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},
	
	En=KerrSpinEnergy[a,p,e,1,\[Sigma],"Linear"->False];
	Jz=KerrSpinAngularMomentum[a,p,e,1,\[Sigma],"Linear"->False];
	x=Jz-(a+\[Sigma])*En;

	jp={1-En^2,-2,a^2+2a*En*x+x^2,-2*((1-En^2)*\[Sigma]^2-En*\[Sigma]*x+x^2),4\[Sigma]^2,-2a*\[Sigma]*(a*\[Sigma]+x*(En*\[Sigma]+x)),-\[Sigma]^2*(-(1-En)*\[Sigma]+x)*((1+En)*\[Sigma]+x)}; (* coefficients for J(\[Psi]) *)
	je={1,2,e^2+3,4(e^2+1),e^4+10e^2+5,2(e^2+3)(3e^2+1),e^6+21e^4+35e^2+7};(* coefficients for J(\[Psi]) *)

	J[\[Chi]_]:=Sum[Sum[jp[[l+1]]*je[[k-l+1]]/(1-e^2)^(k-l)/p^l,{l,0,k}]*(1+e*Cos[\[Chi]])^k,{k,0,6}];
	\[CapitalSigma]s[r_]:=r^2*(1-M*\[Sigma]^2/r^3);
	\[CapitalDelta][r_]:=r^2-2r+a^2;
	Ps[r_]:=((r^2+a^2)+a*\[Sigma]/r*(r+M))*En-(a+M*\[Sigma]/r)*Jz;
	Vt[r_]:=a*(1+3M*\[Sigma]^2/r/\[CapitalSigma]s[r])*(Jz-(a+\[Sigma])*En)+(r^2+a^2)/\[CapitalDelta][r]*Ps[r];
	V\[Phi][r_]:=(1+3M*\[Sigma]^2/r/\[CapitalSigma]s[r])*(Jz-(a+\[Sigma])*En)+a/\[CapitalDelta][r]*Ps[r];

	\[CapitalLambda]r = 2*Sqrt[1-e^2]/p*Quiet[NIntegrate[1/Sqrt[J[\[Chi]]],{\[Chi],0,Pi},WorkingPrecision->If[Precision[{a,p,e,\[Sigma]}]==MachinePrecision,MachinePrecision,Precision[{a,p,e,\[Sigma]}]-8]],{NIntegrate::precw}];

	\[CapitalUpsilon]r = 2*Pi/\[CapitalLambda]r;
	
	\[CapitalUpsilon]\[Phi] = 2*Sqrt[1-e^2]/p/\[CapitalLambda]r*Quiet[NIntegrate[V\[Phi][p*M/(1+e*Cos[\[Chi]])]/Sqrt[J[\[Chi]]],{\[Chi],0,Pi},WorkingPrecision->If[Precision[{a,p,e,\[Sigma]}]==MachinePrecision,MachinePrecision,Precision[{a,p,e,\[Sigma]}]-8]],{NIntegrate::precw}];
	
	\[CapitalGamma]  = 2*Sqrt[1-e^2]/p/\[CapitalLambda]r*Quiet[NIntegrate[Vt[p*M/(1+e*Cos[\[Chi]])]/Sqrt[J[\[Chi]]],{\[Chi],0,Pi},WorkingPrecision->If[Precision[{a,p,e,\[Sigma]}]==MachinePrecision,MachinePrecision,Precision[{a,p,e,\[Sigma]}]-8]],{NIntegrate::precw}];

	<| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
		"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
		"\[CapitalGamma]" -> \[CapitalGamma] 
	|>

]


KerrSpinBoyerLindquistFrequenciesNonlinearTulczyjew[a_?NumericQ,p_?NumericQ,e_?NumericQ,1,\[Sigma]_?NumericQ]:=Module[
	{M=1,En,Jz,x,jp,je,J,\[CapitalSigma]s,\[CapitalDelta],Ps,Vt,V\[Phi],\[Chi],\[CapitalDelta]\[Phi],T},

	En=KerrSpinEnergy[a,p,e,1,\[Sigma],"Linear"->False];
	Jz=KerrSpinAngularMomentum[a,p,e,1,\[Sigma],"Linear"->False];
	x=Jz-(a+\[Sigma])*En;

	jp={1-En^2,-2,a^2+2a*En*x+x^2,-2*((1-En^2)*\[Sigma]^2-En*\[Sigma]*x+x^2),4\[Sigma]^2,-2a*\[Sigma]*(a*\[Sigma]+x*(En*\[Sigma]+x)),-\[Sigma]^2*(-(1-En)*\[Sigma]+x)*((1+En)*\[Sigma]+x)}; (* coefficients for J(\[Psi]) *)
	je={1,2,e^2+3,4(e^2+1),e^4+10e^2+5,2(e^2+3)(3e^2+1),e^6+21e^4+35e^2+7};(* coefficients for J(\[Psi]) *)

	J[\[Chi]_]:=Sum[Sum[jp[[l+1]]*je[[k-l+1]]/(1-e^2)^(k-l)/p^l,{l,0,k}]*(1+e*Cos[\[Chi]])^k,{k,0,6}];
	\[CapitalSigma]s[r_]:=r^2*(1-M*\[Sigma]^2/r^3);
	\[CapitalDelta][r_]:=r^2-2r+a^2;
	Ps[r_]:=((r^2+a^2)+a*\[Sigma]/r*(r+M))*En-(a+M*\[Sigma]/r)*Jz;
	Vt[r_]:=a*(1+3M*\[Sigma]^2/r/\[CapitalSigma]s[r])*(Jz-(a+\[Sigma])*En)+(r^2+a^2)/\[CapitalDelta][r]*Ps[r];
	V\[Phi][r_]:=(1+3M*\[Sigma]^2/r/\[CapitalSigma]s[r])*(Jz-(a+\[Sigma])*En)+a/\[CapitalDelta][r]*Ps[r];

	\[CapitalDelta]\[Phi] = 2*Sqrt[1-e^2]/p*Quiet[NIntegrate[V\[Phi][p*M/(1+e*Cos[\[Chi]])]/Sqrt[J[\[Chi]]],{\[Chi],0,Pi},WorkingPrecision->If[Precision[{a,p,e,\[Sigma]}]==MachinePrecision,MachinePrecision,Precision[{a,p,e,\[Sigma]}]-8]],{NIntegrate::precw}];
	
	T  = 2*Sqrt[1-e^2]/p*Quiet[NIntegrate[Vt[p*M/(1+e*Cos[\[Chi]])]/Sqrt[J[\[Chi]]],{\[Chi],0,Pi},WorkingPrecision->If[Precision[{a,p,e,\[Sigma]}]==MachinePrecision,MachinePrecision,Precision[{a,p,e,\[Sigma]}]-8]],{NIntegrate::precw}];

	<| "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> 2*Pi/T,
		"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> \[CapitalDelta]\[Phi]/T
	|>
]


(* ::Subsection::Closed:: *)
(*Generic function for choosing between frequencies w.r.t different time coordinates*)


(* ::Subsubsection::Closed:: *)
(*\[Sigma] = 0*)


Options[KerrSpinFrequencies] = {"Time" -> "BoyerLindquist","SSC"->"Tulczyjew","Linear"->False};
SyntaxInformation[KerrSpinFrequencies] = {"ArgumentsPattern"->{_,_,_,_,_,OptionsPattern[]}};
KerrSpinFrequencies[a_,p_,e_,x_,\[Sigma]_/;PossibleZeroQ[\[Sigma]],OptionsPattern[]] := Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},


	If[OptionValue["Time"]=="Mino",
		Return[KerrGeoFrequencies[a,p,e,x,"Time"->"Mino"]]
	];

	If[OptionValue["Time"]=="BoyerLindquist" || OptionValue["Time"]=="BL", 
		Return[KerrGeoFrequencies[a,p,e,x,"Time"->"BoyerLindquist"]]
	];

	If[OptionValue["Time"]=="Proper",
		Return[KerrGeoFrequencies[a,p,e,x,"Time"->"Proper"]]
	];

]


(* ::Subsubsection::Closed:: *)
(*\[Sigma] != 0*)


KerrSpinFrequencies[a_,p_,e_,x_,\[Sigma]_,OptionsPattern[]] := Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

	If[OptionValue["SSC"]=="Tulczyjew",
	
		If[!OptionValue["Linear"],

			If[OptionValue["Time"]=="Mino",
				Return[KerrSpinMinoFrequenciesNonlinearTulczyjew[a,p,e,x,\[Sigma]]];
			];

			If[OptionValue["Time"]=="BoyerLindquist" || OptionValue["Time"]=="BL", 
				Return[KerrSpinBoyerLindquistFrequenciesNonlinearTulczyjew[a,p,e,x,\[Sigma]]]
			];

		];
		
		If[OptionValue["Linear"]=="True",
		
			If[OptionValue["Time"]=="Mino",
			Return[KerrSpinMinoFrequenciesGeoPlusLinearTulczyjew[a,p,e,x,\[Sigma]]];
		    ];
		];
	]

]


Options[KerrSpinFrequencyCorrections] = {"Time" -> "BoyerLindquist","SSC"->"Tulczyjew"}
SyntaxInformation[KerrSpinFrequencyCorrections] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};
KerrSpinFrequencyCorrections[a_,p_,e_,x_,OptionsPattern[]] := Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

	If[OptionValue["SSC"]=="Tulczyjew",
	
			If[OptionValue["Time"]=="Mino",
				Return[KerrSpinMinoFrequencyCorrectionsTulczyjew[a,p,e,x]];
			];

			If[OptionValue["Time"]=="BoyerLindquist" || OptionValue["Time"]=="BL", 
				Return[KerrSpinBoyerLindquistFrequencyCorrectionsTulczyjew[a,p,e,x]];
			];

	]

]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
