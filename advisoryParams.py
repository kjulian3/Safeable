# CONSTANTS
H=0     # Variable index for h
VOWN = 1 # Variable for ownship climbrate
VINT = 2 # Variable for intruder climbrate
TAU=3   # Variable index for tau
HP=100  # Height of NMAC puck
G=32.2  # Graviational acceleration

# Advisory indices
COC=0
DNC=1
DND=2
DES1500 = 3
CL1500 = 4
SDES1500=5
SCL1500=6
SDES2500=7
SCL2500=8


def advisoryParams(actInd,eps=1,pd=0):
    params={}
    
    # Set w, vlo, and alo based on advisory
    params['alo'] = G/4
    params['w']   = 1
    params['vlo'] = 0
    if actInd>=SDES1500:
        params['alo']=G/3
    if actInd in [DNC,DES1500,SDES1500,SDES2500]:
        params['w'] =-1
    if actInd in [DES1500,CL1500,SDES1500,SCL1500]:
        params['vlo']=25*params['w']
    elif actInd in [SDES2500,SCL2500]:
        params['vlo']=2500./60.*params['w']
        
    params['vlor']=-2500./60.*params['w']
    params['vlos']= 2500./60.*params['w']
    params['ahi'] = G/2
    params['wr']  = - params['w']
    params['ws']  = params['w']
    params['eps'] = eps
    params['pd']  = pd
    
    # For COC, model advisory using an additional pilot delay and no decision time
    if actInd==COC:
        params['eps'] -= eps
        params['pd'] += eps
    
    params['name'] = advName(actInd)
    return params

def advName(adv):
    if adv==COC:
        return "COC"
    if adv==DNC:
        return 'DNC'
    elif adv==DND:
        return 'DND'
    elif adv==DES1500:
        return 'DES1500'
    elif adv==CL1500:
        return 'CL1500'
    elif adv==SDES1500:
        return 'SDES1500'
    elif adv==SCL1500:
        return 'SCL1500'
    elif adv==SDES2500:
        return 'SDES2500'
    elif adv==SCL2500:
        return 'SCL2500'
    else: 
        return "UNKNOWN"

def getPossibleAdvisories(pra):
    if pra<=DND:
        return [COC,DNC,DND,DES1500,CL1500]
    elif pra<=CL1500:
        return [COC,DNC,DND,DES1500,CL1500,SDES1500,SCL1500]
    else:
        return [COC,DNC,DND,DES1500,CL1500,SDES1500,SCL1500,SDES2500,SCL2500]