import numpy as np
from advisoryParams import *
from safeable import *
from reluplexpy import Reluplex
    
# Constants for network normalization
R_h   = 16000.0  # ft
R_v  = 200.0     # ft/s
R_tau = 40.0     # s
mu_tau= 20.0     # s
R_out = 3.1023
mu_out = -0.432599379632

"""
Check safeable region given a previous RA and RA
Can also specify ownship climb rate ranges, decision time, pilot delay, and 
number of piecewise-linear approximations per curve
"""
def checkRegion(pra,ra,vmin=-100.,vmax=100.,deltaV=2.,eps=1,pd=0,overApprox=True, dt = 0.25,dti=0.0625,):
    v1 = vmin
    print("Checking pRA network: " + advName(pra) + ", Advisory: " + advName(ra))
    while v1<vmax:
        v2 = min(v1+deltaV,vmax)
        print("Speed Range: [%d, %d]"%(v1,v2))
        checkRegionSpeed(pra,ra,[v1,v2],eps,pd,overApprox,dt,dti)
        v1+=deltaV

"""
Check safeable region for a specific speed range.
Called by checkRegion
"""
def checkRegionSpeed(pra,ra,v,eps,pd,overApprox,dt,dti):
    # Get the bounds of the region we need to check
    bndsMin, bndsMax = getRegionToCheck(pra,ra,v,pd=pd,eps=eps,approx=True,overApprox=overApprox,dt=dt,dti=dti)
    
    # Split the bounds into pairs for small ranges of tau
    # Tau<6 is forced to be Tau=6
    bndsSplit = splitBnds(bndsMin + bndsMax)
    
    # Loop through bound pairs and check each slice
    for i,pair in enumerate(bndsSplit):
        checkSlice(pra,ra,v,pd,eps,pair,i)
   
"""
Check a small slice of the non-safeable region we are verifying. This function calls out to Reluplex 
and runs the query.
"""
def checkSlice(pra,ra,v,pd,eps,bnds,count=0):
    # Files paths for network, log file, and output file
    nnet_file_name = "./networks/VertCAS_pra%02d_v4_45HU_200.nnet" % (pra+1)
    fn = "reluplexResults_v1/pra%s_ra%s_vMin%d_vMax%d_pd%d_eps%d_%d.txt"%(advName(pra),advName(ra),v[0],v[1],pd,eps,count)
    marLog = "reluplexLogs_v1/pra%s_ra%s_vMin%d_vMax%d_pd%d_eps%d_%d.log"%(advName(pra),advName(ra),v[0],v[1],pd,eps,count)
    raInd = ra #np.where(np.array(ra_transitions[pra])==ra+1)[0][0]
    
    # Write header for output file
    file = open(fn,"w") 
    file.write("Pra: %s, Advisory: %s, Vel Range: [%d, %d]\nPilot Delay: %d, Decision Time: %d, Bound: %d\n" % (advName(pra),advName(ra),v[0],v[1],pd,eps,count))
    file.close()
    
    # Read network
    net = Reluplex.read_nnet(nnet_file_name)
    inputVars = net.inputVars[0]
    
    # Ownship climb rate bounds based on given range
    net.setLowerBound(inputVars[VOWN],v[0]/R_v)
    net.setUpperBound(inputVars[VOWN],v[1]/R_v)
    
    # Level flight for intruder
    net.setLowerBound(inputVars[VINT], 0.0)
    net.setUpperBound(inputVars[VINT], 0.0)
    
    # Tau bounds based on given range
    net.setLowerBound(inputVars[TAU],(bnds[0].getMinTau()-mu_tau)/R_tau)
    net.setUpperBound(inputVars[TAU],(bnds[0].getMaxTau()-mu_tau)/R_tau)
    
    # Apply each bound in set
    for bnd in bnds:
        
        # If the bound involves just H, set an uppr or lower bound
        if bnd.coeffs[1]==0 and bnd.coeffs[2]==0:
            if not bnd.isLower:
                net.setUpperBound(inputVars[H],bnd.getH_minTau()/R_h)
            else:
                net.setLowerBound(inputVars[H],bnd.getH_minTau()/R_h)
                
        # Otherwise, add the inequality
        elif not bnd.isLower:
            Reluplex.addInequality(net,[inputVars[H],inputVars[TAU]],[R_h,-bnd.coeffs[1]*R_tau],bnd.coeffs[0]-bnd.coeffs[1]*bnd.getMinTau() + mu_tau*bnd.coeffs[1])
            
        # Should not see any other kinds of bounds, not yet supported
        else:
            Reluplex.addInequality(net,[inputVars[H],inputVars[TAU]],[-R_h,bnd.coeffs[1]*R_tau],-bnd.coeffs[0]+bnd.coeffs[1]*bnd.getMinTau() - mu_tau*bnd.coeffs[1])
            
    # Output constraints
    outputVars = net.outputVars[0]

    # Property to be UNSAT: advisory is maximum cost
    # Inequalities take the form vars'weights <= scalar
    # The function inputs are the network, the variable vector, the weight vector, and the scalar
    for i in range(len(outputVars)):
        if i != raInd:
            Reluplex.addInequality(net, [outputVars[raInd], outputVars[i]], [-1, 1],0)

    # Solve
    #vals, stats = net.solve(marLog,verbose=False)
    stats,vals = net.solve(marLog,verbose=False)

    # Print out results and save to file
    file = open(fn,"a")
    if len(vals)==0:
        print("\nUNSAT: %d"%count)
        print("Pra: %s, RA: %s"%(advName(pra),advName(ra)))
        file.write("UNSAT\n")
    else:
        print("\nSAT: %d"%count)
        print("Pra: %s, RA: %s"%(advName(pra),advName(ra)))
        # Unnormalize input values
        print("Inputs:")
        print("h   = {}".format(vals[inputVars[0]]*R_h))
        print("dh0 = {}".format(vals[inputVars[1]]*R_v))
        print("dh1 = {}".format(vals[inputVars[2]]*R_v))
        print("tau = {}".format(vals[inputVars[3]]*R_tau + mu_tau))

        print("Outputs:")
        for i in range(len(outputVars)):
            print("{} = {}".format(advName(i), vals[outputVars[i]]*R_out + mu_out))

        file.write("SAT\n")
        # Unnormalize input values
        file.write("\nInputs:\n")
        file.write("h   = {}\n".format(vals[inputVars[0]]*R_h))
        file.write("dh0 = {}\n".format(vals[inputVars[1]]*R_v))
        file.write("dh1 = {}\n".format(vals[inputVars[2]]*R_v))
        file.write("tau = {}\n".format(vals[inputVars[3]]*R_tau + mu_tau))

        file.write("\nOutputs:\n")
        for i in range(len(outputVars)):
            file.write("{} = {}\n".format(advName(i), vals[outputVars[i]]*R_out + mu_out))
    file.close()