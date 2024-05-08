(****************************************************************)
(* KerrSpinAngularMomentum linear                               *)
(****************************************************************)
VerificationTest[
    KerrSpinAngularMomentum[0.9`20, 10, 0.1`20, 1, 0.1`20]
    ,
    3.5404817393316448786
    ,
    TestID->"KerrSpinAngularMomentum linear"
]

(****************************************************************)
(* KerrSpinAngularMomentum nonlinear                            *)
(****************************************************************)
VerificationTest[
    KerrSpinAngularMomentum[0.9`20, 10, 0.1`20, 1, 0.1`20, 
 "Linear" -> False]
    ,
    3.5404442223265530503
    ,
    TestID->"KerrSpinAngularMomentum nonlinear"
]

(****************************************************************)
(* KerrSpinAngularMomentum correction                           *)
(****************************************************************)
VerificationTest[
    KerrSpinAngularMomentumCorrection[0.9`20, 10, 0.1`20, 1]
    ,
    0.81907243297734421179
    ,
    TestID->"KerrSpinAngularMomentum correction"
]

