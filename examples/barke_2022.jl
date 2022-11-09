using AFMsimulations
# Trying to reproduce Ingo Sweep, powerpoint page 20 (LennardJones.pptx - Unibox)

d = 40.e-9
Q = 400. 
k = 0.7
f_0 = 50.e3
γ = 2.5e-14 # damping term for x^{-3} dependence 
A_free = 102.e-9
force = k * A_free/Q




#some parameters that dont matter in LJ-sweep anyway (need for construction of objects)
R = 10.e-9 
E_t = 130.e9
ν_t = 0.3

E_s = 1.e-9
ν_s = 0.3



