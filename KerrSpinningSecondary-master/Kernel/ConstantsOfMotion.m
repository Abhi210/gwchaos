(* ::Package:: *)

(* ::Title:: *)
(*ConstantOfMotion subpackage of KerrSpinningSecondary*)


BeginPackage["KerrSpinningSecondary`ConstantsOfMotion`",
	{"KerrGeodesics`ConstantsOfMotion`"}
	];

KerrSpinEnergyCorrection::usage = "KerrSpinEnergyCorrection[a,p,e,x,\[Sigma]] returns the linear part of the energy."
KerrSpinAngularMomentumCorrection::usage = "KerrSpinAngularMomentumCorrection[a,p,e,x,\[Sigma]] returns the linear part of the angular momentum."
KerrSpinEnergy::usage = "KerrSpinEnergy[a,p,e,x,\[Sigma]] returns the energy."
KerrSpinAngularMomentum::usage = "KerrSpinAngularMomentum[a,p,e,x,\[Sigma]] returns the angular momentum."

Begin["`Private`"];


(* ::Section::Closed:: *)
(*\[Sigma] = 0 limit*)


KerrSpinEnergy[a_,p_,e_,x_,\[Sigma]_/;PossibleZeroQ[\[Sigma]]]:=KerrGeoEnergy[a,p,e,x]
KerrSpinAngularMomentum[a_,p_,e_,x_,\[Sigma]_/;PossibleZeroQ[\[Sigma]]]:=KerrGeoAngularMomentum[a,p,e,x]


(* ::Section::Closed:: *)
(*Schwarzschild*)


(* ::Subsection::Closed:: *)
(*Circular equatorial aligned orbits*)


(* ::Text:: *)
(*E and L calculation taken from https://arxiv.org/pdf/1801.09616.pdf - Eq (2.15)*)


KerrSpinEnergyCorrectionTulczyjew[0,p_,0,1]:=1/Sqrt[1-3/p] (-p^(-5/2)/(2(1-3/p)))


KerrSpinAngularMomentumCorrectionTulczyjew[0,p_,0,1]:=1/Sqrt[1/p(1-3/p)] ( Sqrt[1/p] (1-2/p)/(1-3/p) (1-9/(2p)))


KerrSpinEnergyGeoPlusLinearTulczyjew[0,p_,0,1,\[Sigma]_]:=1/Sqrt[1-3/p] (1-2/p-\[Sigma] p^(-5/2)/(2(1-3/p)))


KerrSpinAngularMomentumGeoPlusLinearTulczyjew[0,p_,0,1,\[Sigma]_]:=1/Sqrt[1/p(1-3/p)] (1+\[Sigma] Sqrt[1/p] (1-2/p)/(1-3/p) (1-9/(2p)))


(* ::Subsection::Closed:: *)
(*Eccentric equatorial aligned*)


KerrSpinEnergyCorrectionTulczyjew[0,p_,e_,1]:=-(((-1+e^2)^2) /(2 (p (-3-e^2+p)^(3/2))))


KerrSpinAngularMomentumCorrectionTulczyjew[0,p_,e_,1]:=-(((9+3 e^2-2 p) Sqrt[(-4 e^2+(-2+p)^2) p])/(2 (p (-3-e^2+p)^(3/2))))


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Eccentric equatorial aligned*)


KerrSpinEnergyCorrectionTulczyjew[a_,p_,e_,1]:=Module[{sqrt,dsqrtds,num1,den,num2,dnum1ds,dnum2ds,ddends,En,dEnds},
	sqrt=Sqrt[16 p^9 (a^4 (-1+e^2)^2+(-4 e^2+(-2+p)^2) p^2+2 a^2 p (-2+p+e^2 (2+p)))];(* square root from Eq. (B1) and (B2) \[Sqrt](\[Epsilon]^2+\[Kappa]*\[Zeta]) for \[Sigma]=0 *)
	dsqrtds=1/2/sqrt*(16 a (3+e^2) p^7 (a^4 (-1+e^2)^2+(-4 e^2+(-2+p)^2) p^2+2 a^2 p (-2+p+e^2 (2+p))));(* derivative of the square root w.r.t. \[Sigma] \!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Sigma]\)]\(\[Sqrt]\((\[Epsilon]^2 + \[Kappa]*\[Zeta])\)\)\) for \[Sigma]=0 *)
	num1=-16 p^8 (-4 e^4 p-(-3+p) (-2+p)^2 p+a^2 (-1+e^2)^2 (-5+e^2+3 p)+e^2 p (-8+p^2))-8 a (-1+e^2)^2 p^3 sqrt;(* numerator in (B1) \[Kappa]*\[Rho]+2*\[Epsilon]*\[Sigma]-2*\[Sigma]*\[Sqrt](\[Epsilon]^2+\[Kappa]*\[Zeta]) *)
	num2=-16 a p^8 (a^2 (-1+e^2)^2 (-5+e^2+p)+p (12+e^4 (-4+p)-11 p+3 p^2+e^2 (-8-6 p+p^2)))+4 p^3 (-2 a^2 (-1+e^2)^2+p^2 (-3-e^2+p)) sqrt; (* numerator in (B2) \[Epsilon]*\[Rho]-2\[Kappa]*\[Eta]-\[Rho]*\[Sqrt](\[Epsilon]^2+\[Kappa]*\[Zeta]) *)
	den=16 p^9 (-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p); (* denominator in (B1) \[Rho]^2+4*\[Eta]*\[Sigma] *)
	dnum1ds=-4 p ((-1+e^2)^2 (2 a^2 (3+e^2)+p^2) sqrt+2 a p^2 (2 p^3 (a^2 (-1+e^2)^2 (-19+e^4+9 p+e^2 (-14+3 p))+p (40+e^6 (-8+p)-41 p+11 p^2+e^4 (-72-11 p+3 p^2)+e^2 (40-77 p+18 p^2)))+(-1+e^2)^2 dsqrtds)); (* derivaives of the denomiators and mnumerator *)
	dnum2ds=-8 p^6 (2 a^4 (-1+e^2)^2 (-19+e^4+e^2 (-14+p)+3 p)+p^3 (36+e^4 (-12+p)-43 p+17 p^2-2 p^3+3 e^2 (-8-2 p+p^2))+a^2 p (80-73 p+17 p^2+e^6 (-16+5 p)+3 e^4 (-48-17 p+3 p^2)+e^2 (80-137 p+6 p^2)))-8 a p (a^2 (-1+e^2)^2 (3+e^2)+2 (2+e^2+e^4) p^2) sqrt+4 p^3 (-2 a^2 (-1+e^2)^2+p^2 (-3-e^2+p)) dsqrtds;
	ddends=-64 a (3+e^2) p^7 (a^2 (-1+e^2)^2+p (-2+e^2 (-6+p)+p));
	En=Sqrt[num1/den];(* Formula for the energy (B1) *)
	1/2/En*(dnum1ds*den-num1*ddends)/den^2 (* derivatives of En from derivatives of (B1) *)
]

KerrSpinAngularMomentumCorrectionTulczyjew[a_,p_,e_,1]:=Module[{sqrt,dsqrtds,num1,den,num2,dnum1ds,dnum2ds,ddends,En,dEnds},
	sqrt=Sqrt[16 p^9 (a^4 (-1+e^2)^2+(-4 e^2+(-2+p)^2) p^2+2 a^2 p (-2+p+e^2 (2+p)))];(* square root from Eq. (B1) and (B2) \[Sqrt](\[Epsilon]^2+\[Kappa]*\[Zeta]) for \[Sigma]=0 *)
	dsqrtds=1/2/sqrt*(16 a (3+e^2) p^7 (a^4 (-1+e^2)^2+(-4 e^2+(-2+p)^2) p^2+2 a^2 p (-2+p+e^2 (2+p))));(* derivative of the square root w.r.t. \[Sigma] \!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Sigma]\)]\(\[Sqrt]\((\[Epsilon]^2 + \[Kappa]*\[Zeta])\)\)\) for \[Sigma]=0 *)
	num1=-16 p^8 (-4 e^4 p-(-3+p) (-2+p)^2 p+a^2 (-1+e^2)^2 (-5+e^2+3 p)+e^2 p (-8+p^2))-8 a (-1+e^2)^2 p^3 sqrt;(* numerator in (B1) \[Kappa]*\[Rho]+2*\[Epsilon]*\[Sigma]-2*\[Sigma]*\[Sqrt](\[Epsilon]^2+\[Kappa]*\[Zeta]) *)
	num2=-16 a p^8 (a^2 (-1+e^2)^2 (-5+e^2+p)+p (12+e^4 (-4+p)-11 p+3 p^2+e^2 (-8-6 p+p^2)))+4 p^3 (-2 a^2 (-1+e^2)^2+p^2 (-3-e^2+p)) sqrt; (* numerator in (B2) \[Epsilon]*\[Rho]-2\[Kappa]*\[Eta]-\[Rho]*\[Sqrt](\[Epsilon]^2+\[Kappa]*\[Zeta]) *)
	den=16 p^9 (-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p); (* denominator in (B1) \[Rho]^2+4*\[Eta]*\[Sigma] *)
	dnum1ds=-4 p ((-1+e^2)^2 (2 a^2 (3+e^2)+p^2) sqrt+2 a p^2 (2 p^3 (a^2 (-1+e^2)^2 (-19+e^4+9 p+e^2 (-14+3 p))+p (40+e^6 (-8+p)-41 p+11 p^2+e^4 (-72-11 p+3 p^2)+e^2 (40-77 p+18 p^2)))+(-1+e^2)^2 dsqrtds)); (* derivaives of the denomiators and mnumerator *)
	dnum2ds=-8 p^6 (2 a^4 (-1+e^2)^2 (-19+e^4+e^2 (-14+p)+3 p)+p^3 (36+e^4 (-12+p)-43 p+17 p^2-2 p^3+3 e^2 (-8-2 p+p^2))+a^2 p (80-73 p+17 p^2+e^6 (-16+5 p)+3 e^4 (-48-17 p+3 p^2)+e^2 (80-137 p+6 p^2)))-8 a p (a^2 (-1+e^2)^2 (3+e^2)+2 (2+e^2+e^4) p^2) sqrt+4 p^3 (-2 a^2 (-1+e^2)^2+p^2 (-3-e^2+p)) dsqrtds;
	ddends=-64 a (3+e^2) p^7 (a^2 (-1+e^2)^2+p (-2+e^2 (-6+p)+p));
	En=Sqrt[num1/den];(* Formula for the energy (B1) *)
	dEnds=1/2/En*(dnum1ds*den-num1*ddends)/den^2; (* derivatives of En from derivatives of (B1) *)
	-((En*num2*ddends+den*num2*dEnds-den*En*dnum2ds)/(den^2*En^2))(* derivatives of Jz from derivatives of (B2) *)
]


KerrSpinEnergyNonlinearTulczyjew[a_,p_,e_,x_/;x^2==1,\[Sigma]_]:=Module[{\[Kappa],\[Epsilon],\[Rho],\[Eta],\[Sigma]2,\[Zeta]},
  \[Kappa]=1/ (p^6) 4 (-p^10 (-2-2 e+p) (-2+2 e+p)+p^7 (-9+5 e^6-(-7+p) p-e^4 (3+(-7+p) p)+e^2 (7+2 p (1+p))) \[Sigma]^2+(-1+e^2)^2 (-4+p) p^4 (-3+e^4+2 p+2 e^2 (1+p)) \[Sigma]^4+(-1+e^2)^6 p \[Sigma]^6+a^3 (-1+e^2)^2 (3+e^2) p \[Sigma] (p^6+(-1+e^2)^3 \[Sigma]^4)+a p^2 \[Sigma] (p^6 (-8+5 p+10 e^2 p+e^4 (8+p))+4 (-1+e^2)^3 p^3 (-1+e^2+p) \[Sigma]^2+(-1+e^2)^5 (-4+p) \[Sigma]^4)+a^2 (-1+e^2)^2 (p^9-p^6 (3 (-2+p)+e^2 (2+p)) \[Sigma]^2+(3+e^2) p^3 (-2+2 e^4+p+3 e^2 p) \[Sigma]^4-(-1+e^2)^4 \[Sigma]^6));
  \[Epsilon]=-(1/ (p^6))2 (p^10 (e^2 (12+p)-(-4+p) (-3+2 p)) \[Sigma] - 2 (-1+e^2)^2 p^7 (6-6 e^2+(-4+p) p) \[Sigma]^3+(-1+e^2)^2 p^4 (24+e^4 (-24+p)+2 e^2 p (-7+2 p)+p (-19+4 p)) \[Sigma]^5+2 a^4 (3+e^2) p \[Sigma] ((-1+e^2)^2 p^6+(-1+e^2)^5 \[Sigma]^4)+a^2 p^2 \[Sigma] (p^6 (-16+13 p+14 e^2 p+e^4 (16+5 p))-2 (-1+e^2)^2 p^3 (-4-4 e^4+p (-2+3 p)+e^2 (8+p (2+p))) \[Sigma]^2+(-1+e^2)^2 (8+e^6 (-8+5 p)+p (-17+6 p)+3 e^4 (8+p (7+2 p))+e^2 (-24+p (-9+20 p))) \[Sigma]^4)+a^3 (-1+e^2)^2 (2 p^9+p^6 (e^2 (-4+p)+3 (4+p)) \[Sigma]^2+4 (-1+e^2) (3+4 e^2+e^4) p^3 \[Sigma]^4+(-1+e^2)^3 (2+e^2 (-2+p)+3 p) \[Sigma]^6)+a p (2 p^9 (-4+3 p+e^2 (4+p))+p^6 (-18+10 e^6+p (-2+5 p)+2 e^2 (7+5 p (2+p))+e^4 (-6+p (14+p))) \[Sigma]^2+4 (-1+e^2)^3 p^3 (-6+2 e^2 (-1+p)+p (2+p)) \[Sigma]^4+(-1+e^2)^5 (-2+2 e^2+(-4+p) p) \[Sigma]^6));
  \[Rho]=1/ (p^3) 4 ((3+e^2-p) p^8+2 a^3 (-1+e^2)^2 (3+e^2) p^4 \[Sigma]+4 a (2+e^2+e^4) p^6 \[Sigma]+(3+e^2) (1+3 e^2) p^5 \[Sigma]^2-a (-1+e^2)^2 p^2 (-9+3 e^4+p (-4+3 p)+e^2 (6+p (4+p))) \[Sigma]^3-(-1+e^2)^2 p^2 (-3+e^4+2 p+2 e^2 (1+p)) \[Sigma]^4+a^2 (-1+e^2)^2 (2 p^6+p^3 (e^2 (-4+p)+3 (4+p)) \[Sigma]^2-(-1+e^2)^3 \[Sigma]^4));
  \[Eta]=-(1/p^3)2 (2 a^4 (-1+e^2)^2 (3+e^2) p^4 \[Sigma]+(9+3 e^2-2 p) p^8 \[Sigma]-(-1+e^2)^2 p^6 \[Sigma]^3+a^2 p^2 \[Sigma] ((15+10 e^2+7 e^4) p^4-(-1+e^2)^2 (3 e^4-3 (3+4 p)+e^2 (6+4 p)) \[Sigma]^2)+a^3 (-1+e^2)^2 (2 p^6+p^3 (12+9 p+e^2 (-4+3 p)) \[Sigma]^2-(-1+e^2)^3 \[Sigma]^4)+a p^2 (2 (3+e^2) p^6+p^3 (6+9 p+2 e^2 (10+p)+e^4 (6+5 p)) \[Sigma]^2-(-1+e^2)^2 (-6+2 e^4+3 p^2+e^2 (4+p (8+p))) \[Sigma]^4));
  \[Sigma]2=(-1+e^2)^2/p^3 (4 a p^6+4 a^2 (3+e^2) p^4 \[Sigma]+2 p^6 \[Sigma]-2 a p^3 (3 (-4+p)+e^2 (4+p)) \[Sigma]^2-2 p^2 (-9+3 e^4+4 p+e^2 (6+4 p)) \[Sigma]^3-2 a (-1+e^2)^3 \[Sigma]^4);
  \[Zeta]=-(1/ (p^6))4 (a^5 (3+e^2) p \[Sigma] ((-1+e^2)^2 p^6+(-1+e^2)^5 \[Sigma]^4)+3 a p^4 \[Sigma] (p^6 (-4+3 p+e^2 (4+p))+4 (-1+e^2)^3 p^3 \[Sigma]^2+(-1+e^2)^3 (-8+e^2 (-8+p)+3 p) \[Sigma]^4)+4 a^3 p^2 \[Sigma] (p^6 (2 (-1+p)+e^2 p+e^4 (2+p))+(-1+e^2)^3 (-1+e^2-2 p) p^3 \[Sigma]^2+(-1+e^2)^3 (-1+e^4 (-1+p)+4 p+e^2 (2+7 p)) \[Sigma]^4)+p^4 (p^9+(3+e^2-p) (-4+p) p^6 \[Sigma]^2+p^3 (3+5 e^6+p-p^2+e^2 (47+2 (-9+p) p)+e^4 (9+p-p^2)) \[Sigma]^4+2 (-1+e^2)^2 (1+e^2) (-4 e^2+(-2+p)^2) \[Sigma]^6)+a^4 (-1+e^2)^2 (p^9+2 p^6 (e^2 (-1+p)+3 (1+p)) \[Sigma]^2+(3+e^2) p^3 (-2+2 e^4-p-3 e^2 p) \[Sigma]^4+(-1+e^2)^3 (1+e^2 (-1+p)+3 p) \[Sigma]^6)+a^2 p (2 p^9 (-2+p+e^2 (2+p))+p^6 (-9+5 e^6+3 p (-3+4 p)+e^4 (-3+p (7+4 p))+e^2 (7+2 p (9+8 p))) \[Sigma]^2-p^3 (e^8 (4-7 p)+3 p^2 (2+p)-3 (4+p)+e^6 p (2+p) (4+p)+e^2 (32+(26-5 p) p^2)+e^4 (-24+p (2+p (26+p)))) \[Sigma]^4+(-1+e^2)^2 ((-1+e^2)^4-4 (-1+e^2)^3 p+2 (-3-e^2+3 e^4+e^6) p^2+(3+e^2) (1+3 e^2) p^3) \[Sigma]^6));
  Sqrt[(\[Kappa]*\[Rho]+2*\[Sigma]2(\[Epsilon]-x*Sqrt[\[Epsilon]^2+\[Kappa]*\[Zeta]]))/(\[Rho]^2+4\[Eta]*\[Sigma]2)](* Energy *)
]

KerrSpinEnergyGeoPlusLinearTulczyjew[a_,p_,e_,1,\[Sigma]_]:=KerrGeoEnergy[a,p,e,1]+\[Sigma]*KerrSpinEnergyCorrectionTulczyjew[a,p,e,1]


KerrSpinAngularMomentumNonlinearTulczyjew[a_,p_,e_,x_/;x^2==1,\[Sigma]_]:=Module[{\[Kappa],\[Epsilon],\[Rho],\[Eta],\[Sigma]2,\[Zeta],energy},
  \[Kappa]=1/ (p^6) 4 (-p^10 (-2-2 e+p) (-2+2 e+p)+p^7 (-9+5 e^6-(-7+p) p-e^4 (3+(-7+p) p)+e^2 (7+2 p (1+p))) \[Sigma]^2+(-1+e^2)^2 (-4+p) p^4 (-3+e^4+2 p+2 e^2 (1+p)) \[Sigma]^4+(-1+e^2)^6 p \[Sigma]^6+a^3 (-1+e^2)^2 (3+e^2) p \[Sigma] (p^6+(-1+e^2)^3 \[Sigma]^4)+a p^2 \[Sigma] (p^6 (-8+5 p+10 e^2 p+e^4 (8+p))+4 (-1+e^2)^3 p^3 (-1+e^2+p) \[Sigma]^2+(-1+e^2)^5 (-4+p) \[Sigma]^4)+a^2 (-1+e^2)^2 (p^9-p^6 (3 (-2+p)+e^2 (2+p)) \[Sigma]^2+(3+e^2) p^3 (-2+2 e^4+p+3 e^2 p) \[Sigma]^4-(-1+e^2)^4 \[Sigma]^6));
  \[Epsilon]=-(1/ (p^6))2 (p^10 (e^2 (12+p)-(-4+p) (-3+2 p)) \[Sigma] - 2 (-1+e^2)^2 p^7 (6-6 e^2+(-4+p) p) \[Sigma]^3+(-1+e^2)^2 p^4 (24+e^4 (-24+p)+2 e^2 p (-7+2 p)+p (-19+4 p)) \[Sigma]^5+2 a^4 (3+e^2) p \[Sigma] ((-1+e^2)^2 p^6+(-1+e^2)^5 \[Sigma]^4)+a^2 p^2 \[Sigma] (p^6 (-16+13 p+14 e^2 p+e^4 (16+5 p))-2 (-1+e^2)^2 p^3 (-4-4 e^4+p (-2+3 p)+e^2 (8+p (2+p))) \[Sigma]^2+(-1+e^2)^2 (8+e^6 (-8+5 p)+p (-17+6 p)+3 e^4 (8+p (7+2 p))+e^2 (-24+p (-9+20 p))) \[Sigma]^4)+a^3 (-1+e^2)^2 (2 p^9+p^6 (e^2 (-4+p)+3 (4+p)) \[Sigma]^2+4 (-1+e^2) (3+4 e^2+e^4) p^3 \[Sigma]^4+(-1+e^2)^3 (2+e^2 (-2+p)+3 p) \[Sigma]^6)+a p (2 p^9 (-4+3 p+e^2 (4+p))+p^6 (-18+10 e^6+p (-2+5 p)+2 e^2 (7+5 p (2+p))+e^4 (-6+p (14+p))) \[Sigma]^2+4 (-1+e^2)^3 p^3 (-6+2 e^2 (-1+p)+p (2+p)) \[Sigma]^4+(-1+e^2)^5 (-2+2 e^2+(-4+p) p) \[Sigma]^6));
  \[Rho]=1/ (p^3) 4 ((3+e^2-p) p^8+2 a^3 (-1+e^2)^2 (3+e^2) p^4 \[Sigma]+4 a (2+e^2+e^4) p^6 \[Sigma]+(3+e^2) (1+3 e^2) p^5 \[Sigma]^2-a (-1+e^2)^2 p^2 (-9+3 e^4+p (-4+3 p)+e^2 (6+p (4+p))) \[Sigma]^3-(-1+e^2)^2 p^2 (-3+e^4+2 p+2 e^2 (1+p)) \[Sigma]^4+a^2 (-1+e^2)^2 (2 p^6+p^3 (e^2 (-4+p)+3 (4+p)) \[Sigma]^2-(-1+e^2)^3 \[Sigma]^4));
  \[Eta]=-(1/p^3)2 (2 a^4 (-1+e^2)^2 (3+e^2) p^4 \[Sigma]+(9+3 e^2-2 p) p^8 \[Sigma]-(-1+e^2)^2 p^6 \[Sigma]^3+a^2 p^2 \[Sigma] ((15+10 e^2+7 e^4) p^4-(-1+e^2)^2 (3 e^4-3 (3+4 p)+e^2 (6+4 p)) \[Sigma]^2)+a^3 (-1+e^2)^2 (2 p^6+p^3 (12+9 p+e^2 (-4+3 p)) \[Sigma]^2-(-1+e^2)^3 \[Sigma]^4)+a p^2 (2 (3+e^2) p^6+p^3 (6+9 p+2 e^2 (10+p)+e^4 (6+5 p)) \[Sigma]^2-(-1+e^2)^2 (-6+2 e^4+3 p^2+e^2 (4+p (8+p))) \[Sigma]^4));
  \[Sigma]2=(-1+e^2)^2/p^3 (4 a p^6+4 a^2 (3+e^2) p^4 \[Sigma]+2 p^6 \[Sigma]-2 a p^3 (3 (-4+p)+e^2 (4+p)) \[Sigma]^2-2 p^2 (-9+3 e^4+4 p+e^2 (6+4 p)) \[Sigma]^3-2 a (-1+e^2)^3 \[Sigma]^4);
  \[Zeta]=-(1/ (p^6))4 (a^5 (3+e^2) p \[Sigma] ((-1+e^2)^2 p^6+(-1+e^2)^5 \[Sigma]^4)+3 a p^4 \[Sigma] (p^6 (-4+3 p+e^2 (4+p))+4 (-1+e^2)^3 p^3 \[Sigma]^2+(-1+e^2)^3 (-8+e^2 (-8+p)+3 p) \[Sigma]^4)+4 a^3 p^2 \[Sigma] (p^6 (2 (-1+p)+e^2 p+e^4 (2+p))+(-1+e^2)^3 (-1+e^2-2 p) p^3 \[Sigma]^2+(-1+e^2)^3 (-1+e^4 (-1+p)+4 p+e^2 (2+7 p)) \[Sigma]^4)+p^4 (p^9+(3+e^2-p) (-4+p) p^6 \[Sigma]^2+p^3 (3+5 e^6+p-p^2+e^2 (47+2 (-9+p) p)+e^4 (9+p-p^2)) \[Sigma]^4+2 (-1+e^2)^2 (1+e^2) (-4 e^2+(-2+p)^2) \[Sigma]^6)+a^4 (-1+e^2)^2 (p^9+2 p^6 (e^2 (-1+p)+3 (1+p)) \[Sigma]^2+(3+e^2) p^3 (-2+2 e^4-p-3 e^2 p) \[Sigma]^4+(-1+e^2)^3 (1+e^2 (-1+p)+3 p) \[Sigma]^6)+a^2 p (2 p^9 (-2+p+e^2 (2+p))+p^6 (-9+5 e^6+3 p (-3+4 p)+e^4 (-3+p (7+4 p))+e^2 (7+2 p (9+8 p))) \[Sigma]^2-p^3 (e^8 (4-7 p)+3 p^2 (2+p)-3 (4+p)+e^6 p (2+p) (4+p)+e^2 (32+(26-5 p) p^2)+e^4 (-24+p (2+p (26+p)))) \[Sigma]^4+(-1+e^2)^2 ((-1+e^2)^4-4 (-1+e^2)^3 p+2 (-3-e^2+3 e^4+e^6) p^2+(3+e^2) (1+3 e^2) p^3) \[Sigma]^6));
  energy=KerrSpinEnergyNonlinearTulczyjew[a,p,e,x,\[Sigma]];(* Energy *)
  (-2\[Kappa]*\[Eta]+\[Rho](\[Epsilon]-x*Sqrt[\[Epsilon]^2+\[Kappa]*\[Zeta]]))/((\[Rho]^2+4\[Eta]*\[Sigma]2)*energy) (* angular momentum *)
]

KerrSpinAngularMomentumGeoPlusLinearTulczyjew[a_,p_,e_,1,\[Sigma]_]:=KerrGeoAngularMomentum[a,p,e,1]+\[Sigma]*KerrSpinAngularMomentumCorrectionTulczyjew[a,p,e,1]


(* ::Text:: *)
(*Some Kerr results can be found in, e.g., https://arxiv.org/pdf/1503.07060.pdf*)


Options[KerrSpinEnergyCorrection] = {"SSC" -> "Tulczyjew"}
KerrSpinEnergyCorrection[a_,p_,e_,x_,OptionsPattern[]] := Module[{ssc},

  ssc = OptionValue["SSC"];
  If[ssc == "Tulczyjew",
    Return[KerrSpinEnergyCorrectionTulczyjew[a,p,e,x]];
  ];
  Print["Unrecognized spin supplementary condition: " <> ssc];

]


Options[KerrSpinAngularMomentumCorrection] = {"SSC" -> "Tulczyjew"}
KerrSpinAngularMomentumCorrection[a_,p_,e_,x_,OptionsPattern[]] := Module[{ssc},

  ssc = OptionValue["SSC"];
  If[ssc == "Tulczyjew",
    Return[KerrSpinAngularMomentumCorrectionTulczyjew[a,p,e,x]];
  ];
  Print["Unrecognized spin supplementary condition: " <> ssc];

]


Options[KerrSpinEnergy] = {"SSC" -> "Tulczyjew","Linear"->True}
KerrSpinEnergy[a_,p_,e_,x_,\[Sigma]_,OptionsPattern[]] := Module[{ssc,lin},

  ssc = OptionValue["SSC"];
  lin = OptionValue["Linear"];
  
  If[ssc == "Tulczyjew",
    If[lin,
    Return[KerrSpinEnergyGeoPlusLinearTulczyjew[a,p,e,x,\[Sigma]]];,
    Return[KerrSpinEnergyNonlinearTulczyjew[a,p,e,x,\[Sigma]]];
    ]
  ];
  Print["Unrecognized spin supplementary condition: " <> ssc];

]

Options[KerrSpinAngularMomentum] = {"SSC" -> "Tulczyjew","Linear"->True}
KerrSpinAngularMomentum[a_,p_,e_,x_,\[Sigma]_,OptionsPattern[]] := Module[{ssc,lin},

  ssc = OptionValue["SSC"];
  lin = OptionValue["Linear"];
  
  If[ssc == "Tulczyjew",
    If[lin,
    Return[KerrSpinAngularMomentumGeoPlusLinearTulczyjew[a,p,e,x,\[Sigma]]];,
    Return[KerrSpinAngularMomentumNonlinearTulczyjew[a,p,e,x,\[Sigma]]];
    ]
  ];
  Print["Unrecognized spin supplementary condition: " <> ssc];

]


(* ::Section::Closed:: *)
(*Close package*)


End[];

EndPackage[];
