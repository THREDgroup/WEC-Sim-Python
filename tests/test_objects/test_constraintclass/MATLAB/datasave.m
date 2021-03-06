body1_linearHydroRestCoef = body(1).hydroForce.linearHydroRestCoef 
body1_visDrag = body(1).hydroForce.visDrag
body1_linearDamping = body(1).hydroForce.linearDamping
body1_userDefinedFe = body(1).hydroForce.userDefinedFe
body1_re = body(1).hydroForce.fExt.re
body1_im = body(1).hydroForce.fExt.im
body1_md = body(1).hydroForce.fExt.md
body1_fAddedMass = body(1).hydroForce.fAddedMass
body1_fDamping = body(1).hydroForce.fDamping
body1_totDOF = body(1).hydroForce.totDOF
% irkb = body(1).hydroForce.irkb

body2_linearHydroRestCoef = body(2).hydroForce.linearHydroRestCoef 
body2_visDrag = body(2).hydroForce.visDrag
body2_linearDamping = body(2).hydroForce.linearDamping
body2_userDefinedFe = body(2).hydroForce.userDefinedFe
body2_re = body(2).hydroForce.fExt.re
body2_im = body(2).hydroForce.fExt.im
body2_md = body(2).hydroForce.fExt.md
body2_fAddedMass = body(2).hydroForce.fAddedMass
body2_fDamping = body(2).hydroForce.fDamping
body2_totDOF = body(2).hydroForce.totDOF

CTTime = simu.CTTime
w = waves.w
waveAmpTime = waves.waveAmpTime

save body_4_test/body1_linearHydroRestCoef.mat  body1_linearHydroRestCoef
save body_4_test/body1_visDrag.mat  body1_visDrag 
save body_4_test/body1_linearDamping.mat body1_linearDamping
save body_4_test/body1_userDefinedFe.mat  body1_userDefinedFe 
save body_4_test/body1_re.mat  body1_re
save body_4_test/body1_im.mat  body1_im 
save body_4_test/body1_md.mat  body1_md
save body_4_test/body1_fAddedMass.mat body1_fAddedMass
save body_4_test/body1_fDamping.mat  body1_fDamping
save body_4_test/body1_totDOF.mat  body1_totDOF

save body_4_test/body2_linearHydroRestCoef.mat  body2_linearHydroRestCoef
save body_4_test/body2_visDrag.mat  body2_visDrag 
save body_4_test/body2_linearDamping.mat body2_linearDamping
save body_4_test/body2_userDefinedFe.mat  body2_userDefinedFe 
save body_4_test/body2_re.mat  body2_re
save body_4_test/body2_im.mat  body2_im 
save body_4_test/body2_md.mat  body2_md 
save body_4_test/body2_fAddedMass.mat body2_fAddedMass
save body_4_test/body2_fDamping.mat  body2_fDamping
save body_4_test/body2_totDOF.mat  body2_totDOF

save body_4_test/CTTime.mat CTTime
save body_4_test/w.mat w
save body_4_test/waveAmpTime.mat waveAmpTime
