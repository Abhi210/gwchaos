(****************************************************************)
(* KerrSpinEnergy linear                                        *)
(****************************************************************)
VerificationTest[
    KerrSpinEnergy[0.9`20, 10, 0.1`20, 1, 0.1`20]
    ,
    0.9525434401597894423405
    ,
    TestID->"KerrSpinEnergy linear"
]

(****************************************************************)
(* KerrSpinEnergy nonlinear                                     *)
(****************************************************************)
VerificationTest[
    KerrSpinEnergy[0.9`20, 10, 0.1`20, 1, 0.1`20, 
 "Linear" -> False]
    ,
    0.952542996410938615799
    ,
    TestID->"KerrSpinEnergy nonlinear"
]

(****************************************************************)
(* KerrSpinEnergy correction                                    *)
(****************************************************************)
VerificationTest[
    KerrSpinEnergyCorrection[0.9`20, 10, 0.1`20, 1]
    ,
    -0.001534290476164244237
    ,
    TestID->"KerrSpinEnergy correction"
]

