
"""
Script to process through the results of Reluplex.
Test
"""

import os

"""
Generate dictionary of dictionaries for each pra containing
number of queries run for each ra
"""
def genDict(directory):
    adv         = {}

    for filename in os.listdir(directory):
        # Decompose filename
        split = filename.split('_',3)
        pra = split[0][3:]
        ra  = split[1][2:]
        
        if pra not in adv:
            adv[pra] = {}

        # Count number of results
        if ra in adv[pra]:
            adv[pra][ra] += 1
        else:
            adv[pra][ra] = 1
    
    return adv

"""
For each pra, return total number of queries run
"""
def praNum(directory,pra="COC"):
    adv = genDict(directory)

    if pra not in adv:
        print("ERROR: No queries run for Previous Advisory: %s"%pra)
        exit(1)

    # Iterate through pra and count total queries
    total = 0
    for ra in adv[pra]:
        total += adv[pra][ra]
    
    return total

"""
Total number of queries run (all pra's)
"""
def totalNum(directory):
    adv = genDict(directory)
    
    # Iterate through all pra's and count total queries
    total = 0
    for pra in adv:
        total += praNum(directory,pra)

    return(total)

"""
Read in inputs of SAT point from filename
"""
def readInputs(line):

    # Read in Inputs
    h   = line[0][line[0].find("=")+1:-1]
    tau = line[3][line[3].find("=")+1:-1]
    inputs = [float(h),float(tau)]
    return inputs

"""
Generate dict of number of SAT points
"""
def genDictSat(directory):
    adv = {}
    inputs = []

    counter=0
    for filename in os.listdir(directory):
        lines  = []
        # Decompose filename
        split = filename.split('_',3)
        pra = split[0][3:]
        ra  = split[1][2:]

        filename = os.path.join(directory,filename)
        with open(filename) as f:
            # Ignore first 2 lines
            next(f)
            next(f)
            
            # Read in SAT/UNSAT
            line = f.readline()

            # The query hasn't been run, skip
            if not len(line):
                print("Query not run for Pra: %s, ra: %s"%(pra,ra))
                continue
            
            # SAT point
            if line[0]=="S":
                # Ignore un-needed text
                next(f)
                next(f)

                # Read in input lines
                for i in range(4):
                    dat = f.readline()
                    lines += [dat]

                if pra not in adv:
                    adv[pra] = {}

                # Count SAT points
                if ra in adv[pra]:
                    adv[pra][ra] += 1
                else:
                    adv[pra][ra] = 1
                
                # Read in Inputs of SAT  
                dat = readInputs(lines)
                if dat[1] == 6.0:
                    counter += 1
                inputs += [dat]

    return adv, inputs, counter

"""
For each pra, return total number of SAT points
"""
def countSat(directory, pra="COC"):
    adv,inputs,offset = genDictSat(directory)

    # Iterate through all pra's and count total queries
    total = 0 
    for ra in adv[pra]:
        total += adv[pra][ra]

    return(total)

"""
Return total number of SAT points (all pra's)
"""
def totalSat(directory):
    adv,inputs,offset = genDictSat(directory)

    # Iterate through all pra's and count total queries
    total = 0
    for pra in adv:
        total += countSat(directory,pra)

    return(total)

def main(directory):
    print("Total queries run = %d"%totalNum(directory))
    print("Total number of SAT points = %d"%totalSat(directory))

if __name__ == "__main__":
    directory="res_underApprox_dt2.000_dti0.0625_v2"
    main(directory)







