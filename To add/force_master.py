import src.force_calculations_v2 as fc
import src.moment_calculations as mc

def reaction_forces():
    FY = fc.Y_force(n_steps=10000, plot=False, info=False, ret=True)
    FZ = fc.Z_force(FY, n_steps=10000, plot=False, info=False, ret=True)
    return(FY,FZ)
    
def torque_calc():
    FY,FZ = reaction_forces()
    mc.moment_calc(FY,FZ)
    
print(reaction_forces())
# print(torque_calc())


"""
OUTPUT FORM:
    
[F1Y, -, F2Y, -, F3Y],[F1Z, R, F2Z, P, F3Z]

POSITIVE Y IS UP AND POSITIVE Z IS TOWARDS TRAILING EDGE, INCLUDING P AND R FOR CONSISTENCY!
"""


    

