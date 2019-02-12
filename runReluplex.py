import sys
from reluplexFunctions import checkRegionConfig, checkRegion
from advisoryParams import *
import multiprocessing
import os
import time

'''
Check an advisory from config file
'''
def runConfig():
    cfg_file = "./reluplex.cfg"
    if len(sys.argv)>1:
        cfg_file = sys.argv[1]
    checkRegionConfig(cfg_file)

'''
Run convergence experiments on different threads
'''
def runExpConv(pra):
        
    ras         = getPossibleAdvisories(advInd(pra))
    overApprox  = False
    dts         = [0.125, 0.25, 0.5, 1, 2]
    threads     = []
    vmin        = float(-100.)
    vmax        = float(100.)
    deltaV      = float(2.)
    eps         = 1
    dti         = float(0.0625)
    pd          = 0
    pra         = advInd(pra)
    verbose     = False

    # Create Threads to run for the convergence problem
    for dt in dts:
        for ra in ras:

            dt = float(dt)
            # Build Folder Name
            if overApprox:
                s = "overApprox"
            else:
                s = "underApprox"

            # Make the folder to hold results
            folder = "res_%s_dt%.3f_dti0.0625_v4"%(s, dt)
            if not os.path.exists(folder):
                os.mkdir(folder)

            # Creating separate threads
            threads += [multiprocessing.Process(target=checkRegion,
                                args=(pra,ra,vmin,vmax,deltaV,eps,pd,overApprox,dt,dti,folder,verbose))]

    print("Running %d threads ...."%(len(threads)))

    # Start all threads
    for t in threads:
        t.start()

    # Wait until threads are complete
    for t in threads:
        t.join()

if __name__== "__main__":
    pra = "COC"
    if advInd(pra) == -1:
        print("ERROR: %s is a invalid advisory!"%pra)
        exit(1)
    start_time = time.time() # Start timer
    runExpConv(pra)
    print("------ Finished running in %s seconds ------"%(time.time() - start_time))




