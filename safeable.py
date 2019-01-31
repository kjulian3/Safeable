import numpy as np
import matplotlib
#matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from advisoryParams import *

# CONSTANTS
H=0      # Variable index for h
VOWN = 1 # Variable for ownship climbrate
VINT = 2 # Variable for intruder climbrate
TAU=3    # Variable index for tau
HP=100   # Height of NMAC puck
G=32.2   # Graviational acceleration

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

'''
h = coeffs[0] + coeffs[1]*(t-minTau) + 0.5coeffs[2]*(t-minTau)^2
'''
class Bound(object):
    def __init__(self,minTau,maxTau,coeffs,isLower):
        self.minTau = minTau
        self.maxTau = maxTau
        self.coeffs = coeffs
        self.isLower = isLower
        
    '''
    Print representation
    '''
    def __repr__(self):
        rep = "Tau Range: [%.3f, %.3f]\n" % (self.minTau,self.maxTau)
        rep+= "Start/Final H: [%.2f, %.2f]\n" % (self.getH_minTau(),self.getH_maxTau())
        rep+= "Start/Final V: [%.2f, %.2f]\n" % (self.getV_minTau(), self.getV_maxTau())
        rep+= "Coefficients:  [%.2f, %.2f, %.2f]\n" % (self.coeffs[0],self.coeffs[1],self.coeffs[2])
        rep+= "Is Lower: " + str(self.isLower) + "\n\n"
        
        return rep
    
    def getLine(self, dt=0.01):
        n = int((self.maxTau-self.minTau)/dt)
        t = np.linspace(self.minTau,self.maxTau,n)
        h = np.array([self.getH(ti) for ti in t])
        return t,h
       
    def getMinTau(self):
        return self.minTau
    
    def getMaxTau(self):
        return self.maxTau
    
    def getH(self,t):
        return self.coeffs[0] + self.coeffs[1]*(t-self.minTau) + 0.5*self.coeffs[2]*(t-self.minTau)**2
   
    def getH_minTau(self):
        return self.getH(self.minTau)
    
    def getH_maxTau(self):
        return self.getH(self.maxTau)
    
    def getV(self,t):
        return self.coeffs[1] + self.coeffs[2]*(t-self.minTau)
    
    def getV_minTau(self):
        return self.getV(self.minTau)
    
    def getV_maxTau(self):
        return self.getV(self.maxTau)
    
    def setMinTau(self,newMinTau):
        diff = newMinTau - self.minTau
        self.coeffs[0] += self.coeffs[1]*diff + 0.5*self.coeffs[2]*diff**2
        self.coeffs[1] += self.coeffs[2]*diff
        self.minTau = newMinTau
        
    def setMaxTau(self,newMaxTau):
        self.maxTau = newMaxTau
        
    def reverse(self):
        self.isLower = not self.isLower
        
    def inRange(self, t):
        return self.minTau<= t and self.maxTau>t
    
    def getStandardCoeffs(self):
        a0, a1, a2 = self.coeffs
        ma = self.minTau
        return [a0-a1*ma+0.5*a2*ma**2, a1-a2*ma, 0.5*a2]
    
    def overlaps(self,bnd):
        c1 = self.getStandardCoeffs()
        c2 =  bnd.getStandardCoeffs()
        return np.all([c1[i]==c2[i] for i in range(len(c1))])
    
    def isLinear(self):
        return self.coeffs[2]==0.0
    
    def copy(self):
        return Bound(self.getMinTau(),self.getMaxTau(),[self.coeffs[0],self.coeffs[1],self.coeffs[2]],self.isLower)
    
    def getTimeAtH(self,h):
        a0, a1, a2 = self.coeffs
        ma = self.minTau
        c = a0 - a1*ma+0.5*a2*ma**2 - h
        b = a1 - a2*ma
        a = 0.5*a2
        disc = b**2-4*a*c
        if disc<0:
            return []
        elif disc==0:
            if self.inRange(-b/2/a): return [-b/2/a]
            return []
        times = []
        t1 = (-b-np.sqrt(disc))/2/a
        t2 = (-b+np.sqrt(disc))/2/a
        if self.inRange(t1):
            times+=[t1]
        if self.inRange(t2):
            times+=[t2]
        return times
    
    def satisfies(self,h,t):
        if t<self.minTau or t>self.maxTau: return 0
        hBound = self.getH(t)
        if self.isLower and hBound<=h: return 1
        if not self.isLower and h<=hBound: return 1
        return 0
        
def satisfiesBounds(bnds,h,t):
    return np.sum([bnd.satisfies(h,t) for bnd in bnds])==2

def intersectTimes(bound1, bound2):
    a0, a1, a2 = bound1.coeffs
    ma = bound1.minTau
    b0, b1, b2 = bound2.coeffs
    mb = bound2.minTau
    c = (a0-b0) - (a1*ma - b1*mb) + 0.5*(a2*ma**2 - b2*mb**2)
    b = (a1-b1) - (a2*ma - b2*mb)
    a = 0.5*(a2-b2)
    if a==0:
        if b==0:
            return []
        t = -c/b
        if bound1.inRange(t) and bound2.inRange(t):
            return [t]
        return []
    disc = b**2-4*a*c
    if disc<0:
        return []
    times = []
    t1 = (-b-np.sqrt(disc))/2/a
    t2 = (-b+np.sqrt(disc))/2/a
    if bound1.inRange(t1) and bound2.inRange(t1):
        times+=[t1]
    if bound1.inRange(t2) and bound2.inRange(t2):
        times+=[t2]
    return times
        
        
def truncateBounds(boundsMin,boundsMax):
    minInd = 0
    maxInd = 0
    while minInd < len(boundsMin) and maxInd < len(boundsMax):
        tInt = intersectTimes(boundsMin[minInd],boundsMax[maxInd])
        if len(tInt)>0:
            tInt = tInt[0]
            boundsMin[minInd].setMaxTau(tInt)
            boundsMax[maxInd].setMaxTau(tInt)
            return boundsMin[:minInd+1], boundsMax[:maxInd+1], tInt
        
        if boundsMin[minInd].getMaxTau() < boundsMax[maxInd].getMaxTau():
            minInd+=1
        else:
            maxInd+=1
    return boundsMin, boundsMax, -1

def truncBounds(boundsMin,boundsMax,t):
    minInd = 0
    while minInd < len(boundsMin) and not boundsMin[minInd].inRange(t):
        minInd +=1
    maxInd = 0
    while maxInd < len(boundsMax) and not boundsMax[maxInd].inRange(t):
        maxInd +=1
    boundsMin[minInd].setMaxTau(t)
    boundsMax[maxInd].setMaxTau(t)
    return boundsMin[:minInd+1], boundsMax[:maxInd+1]
        
    
def findExtrema(bounds,increase):
    ext = []
    for bnd in bounds:
        if bnd.coeffs[2]!=0.0:
            a = 0.5*bnd.coeffs[2]
            if np.sign(a)==increase:
                b = bnd.coeffs[1]-bnd.coeffs[2]*bnd.getMinTau()
                t = -b/2/a
                if bnd.inRange(t):
                    ext += [(t,bnd.getH(t))]
    return ext

def flattenBound(bnd,ext,increase):
    ext_t,ext_h = ext
    if bnd.isLinear():
        if bnd.coeffs[1]==0.0:
            if bnd.coeffs[0]*increase <= ext_h*increase:
                return [bnd]
            else:
                return [Bound(bnd.minTau,bnd.maxTau,[ext_h,0,0],bnd.isLower)]
        t = (ext_h-bnd.coeffs[0]+bnd.coeffs[1]*bnd.getMinTau())/bnd.coeffs[1]
        if t>= bnd.getMaxTau():
            if bnd.getH_maxTau()*increase <= increase*ext_h:
                return [bnd]
            else:
                return [Bound(bnd.minTau,bnd.maxTau,[ext_h,0,0],bnd.isLower)]
                
        elif t<=bnd.getMinTau(): return [Bound(bnd.minTau,bnd.maxTau,[ext_h,0,0],bnd.isLower)]

        keepFirst = np.sign(coeffs[1]==increase)
        if keepFirst:
            bnd.setMaxTau(t)
            return [bnd,Bound(t,bnd.maxTau,[ext_h,0,0],bnd.isLower)]
        bnd.setMinTau(t)
        return [Bound(bnd.minTau,t,[ext_h,0,0],bnd.isLower),bnd]
    else:
        if increase != np.sign(bnd.coeffs[2]):
            times = bnd.getTimeAtH(ext_h)
            if len(times)==0:
                return [Bound(bnd.minTau,bnd.maxTau,[ext_h,0,0],bnd.isLower)]
            t = times[0]
            t2 = bnd.getMaxTau()
            bnd.setMaxTau(t)
            return [bnd,Bound(t,t2,[ext_h,0,0],bnd.isLower)]
        
        if bnd.inRange(ext_t):
            mt = bnd.minTau
            bnd.setMinTau(ext_t)
            return [Bound(mt,ext_t,[ext_h,0,0],bnd.isLower),bnd]
        else:
            return [Bound(bnd.minTau,bnd.maxTau,[ext_h,0,0],bnd.isLower)]
        

def monotonic(bounds,increase): 
    ext = findExtrema(bounds,increase)
    if len(ext)==0: 
        if np.sign(bounds[0].getV_minTau())==increase:
            return bounds
        else:
            h = bounds[-1].getH_maxTau()
            mt = bounds[-1].getMaxTau()
            isLower = bounds[-1].isLower
            return [Bound(0,mt,[h,0,0],isLower)]
    
    finalBounds = []
    extInd = 0
    bndInd = 0
    while bndInd <len(bounds) and extInd<len(ext):
        finalBounds+=flattenBound(bounds[bndInd],ext[extInd],increase)
        if bounds[bndInd].inRange(ext[extInd][0]):
            extInd+=1
        bndInd+=1
    finalBounds += bounds[bndInd:]
    finalBounds = combineBounds(finalBounds)
    return finalBounds

def tryCombine(bnd1,bnd2):
    h1 = bnd1.getH_maxTau()
    h2 = bnd2.getH_minTau()
    v1 = bnd1.getV_maxTau()
    v2 = bnd2.getV_minTau()
    a1 = bnd1.coeffs[2]
    a2 = bnd2.coeffs[2]
    t1 = bnd1.getMaxTau()
    t2 = bnd2.getMinTau()
    l1 = bnd1.isLower
    l2 = bnd2.isLower
    
    if h1==h2 and v1==v2 and a1==a2 and t1==t2 and l1==l2:
        bnd1.setMaxTau(bnd2.getMaxTau())
        return True, [bnd1]
    return False, [bnd2]
    
def combineBounds(bnds):
    finalBounds = []
    bInd = 1
    finalBounds += [bnds[0]]
    while bInd<len(bnds):
        success, nb = tryCombine(finalBounds[-1],bnds[bInd])
        if success:
            finalBounds[-1] = nb[0]
        else:
            finalBounds += nb
        bInd+=1
    return finalBounds

'''
Use the "Inner" method for bounding a portion of the trajectory.
Returns a list of Bounds to bound the trajectory
'''
def innerApproximate(bnd,dt,dti):
    bounds = []
    minTau = bnd.getMinTau()
    maxTau = bnd.getMaxTau()
    
    t0 = minTau
    hSeg_0 = bnd.getH(t0)
    hSeg_1 = bnd.getH(t0+dt)
    
    i=0
    num = int(dt/dti)
    while t0+dti< maxTau:
        if i%num==0:
            hSeg_0 = bnd.getH(t0)
            hSeg_1 = bnd.getH(t0+dt)
            if t0+dt>maxTau:
                hSeg_1 = hSeg_0 + dt*(bnd.getH_maxTau()-hSeg_0)/(maxTau-t0)
            i=0
            
        t1 = t0 + dti
        h0 = hSeg_0*(num-i)/num + hSeg_1*(i/num)
        i+=1
        h1 = hSeg_0*(num-i)/num + hSeg_1*(i/num)
        
        bounds += [Bound(t0,t1,[h0,(h1-h0)/dti,0],bnd.isLower)]
        t0=t1
    if len(bounds)==0:
        h0 = bnd.getH_minTau()
        h1 = bnd.getH_maxTau()
        bounds+= [Bound(minTau,maxTau,[h0,(h1-h0)/(maxTau-minTau),0],bnd.isLower)]
    elif t0<maxTau:
        t1 = maxTau
        h0 = bounds[-1].getH_maxTau()
        h1 = bnd.getH_maxTau()
        bounds+= [Bound(t0,t1,[h0,(h1-h0)/(t1-t0),0],bnd.isLower)]
    return bounds

def outerApproximate(bnd,dt,dti):
    bounds = []
    minTau = bnd.getMinTau()
    maxTau = bnd.getMaxTau()
    
    t0 = minTau
    dt/=2
    num = int(dt/dti)
    i=0
    awayFromCurve=True
    while t0+dti< maxTau:
        if i%num==0:
            if awayFromCurve:
                v = bnd.getV(t0)
                h = bnd.getH(t0)
                tMin = t0
            else:
                v = bnd.getV(t0+dt)
                h = bnd.getH(t0+dt)
                tMin = t0+dt
            i=0
            awayFromCurve = not awayFromCurve
            
        t1 = t0 + dti
        bounds += [Bound(tMin,t1,[h,v,0],bnd.isLower)]
        bounds[-1].setMinTau(t0)
        t0=t1
        i+=1
    if len(bounds)==0:
        h = bnd.getH_minTau()
        v = bnd.getV_minTau()
        bounds+= [Bound(minTau,maxTau,[h,v,0],bnd.isLower)]
    elif t0<maxTau:
        t1 = maxTau
        h0 = bounds[-1].getH_maxTau()
        h1 = bnd.getH_maxTau()
        bounds+= [Bound(t0,t1,[h0,(h1-h0)/(t1-t0),0],bnd.isLower)]
    return bounds

def overApproximate(bnds,dt,dti):
    finalBounds = []
    for bnd in bnds:
        if bnd.isLower:
            if bnd.coeffs[2]<0: finalBounds += innerApproximate(bnd,dt,dti)
            else: finalBounds += outerApproximate(bnd,dt,dti)
        else:
            if bnd.coeffs[2]<0: finalBounds += outerApproximate(bnd,dt,dti)
            else: finalBounds += innerApproximate(bnd,dt,dti)
    return finalBounds
        
def underApproximate(bnds,dt,dti):
    finalBounds = []
    for bnd in bnds:
        if bnd.isLower:
            if bnd.coeffs[2]<0: finalBounds += outerApproximate(bnd,dt,dti)
            else: finalBounds += innerApproximate(bnd,dt,dti)
        else:
            if bnd.coeffs[2]<0: finalBounds += innerApproximate(bnd,dt,dti)
            else: finalBounds += outerApproximate(bnd,dt,dti)
    return finalBounds
        
    
    
'''
Compute nominal trajectory. Also returns the ownship climbrate at the end of the trajectory
'''
def getNominalTraj(v,vlo,a,tauMin,tauMax,hInit,isLower=True):
    if tauMax<=tauMin:
        return []
    
    # Compute time to end of trajectory
    tauFinal = (vlo-v)/a
    
    # If tauFinal is negative, then the ownship velocity is already compliant
    if tauFinal<=0:
        return [Bound(tauMin,tauMax,[hInit,v,0],isLower)]
    
    # Compute trajectory as ownship follows advisory
    tauMax_adv = min(tauFinal+tauMin,tauMax)
    traj = [Bound(tauMin,tauMax_adv,[hInit,v,a],isLower)]
    
    
    # If there is more time remaining after we reach compliance with advisory,
    # Assume that we follow along at vlo for remaining time
    if tauMax> tauMax_adv:
        tauMin= tauMax_adv
        hInit = traj[0].getH_maxTau()
        v     = traj[0].getV_maxTau()
        a     = 0
        traj+= [Bound(tauMin,tauMax,[hInit,v,a],isLower)]
        
    return traj

'''
Returns bounds that define the safeable region.
The direction of the bounds defines the unsafe region for an advisory
Input is a params dictionary and velocity range. Optional number of 
piece-wise linear approximations per curve can be given as n
'''
def getSafeable(params,v,worstCase=True):
    hInitLower = 0
    hInitUpper= 0
    boundListLo = []
    boundListHi = []
    trajTau = np.array([])
    trajH_upper = np.array([])
    trajH_lower = np.array([])
    
    # Get variables from params dict:
    w    = params['w']
    vlo  = params['vlo']
    alo  = params['alo']
    ahi  = params['ahi']
    wr   = params['wr']
    vlor = params['vlor']
    ws   = params['ws']
    vlos = params['vlos']
    eps  = params['eps']
    pd   = params['pd']
    
    # How long to generate trajectories
    tf = 15 + eps+pd
    vmin,vmax=v
    
    # Simulate trajectory and get bounds for an initial pilot delay
    trajMin = []
    trajMax = []
    boundMin = []
    boundMax = []
    
    trajMin  += getNominalTraj(vmin,vmin,w*alo,0,pd,hInitLower,isLower=True)
    trajMax  += getNominalTraj(vmax,vmax,w*alo,0,pd,hInitUpper,isLower=False)
    boundMin += getNominalTraj(vmin,vmin,w*alo,0,pd,hInitLower-HP,isLower=True)
    boundMax += getNominalTraj(vmax,vmax,w*alo,0,pd,hInitUpper+HP,isLower=False)
        
    # Update initial altitudes of trajectories
    if len(trajMin)>0:
        hInitLower = trajMin[-1].getH_maxTau()
        hInitUpper = trajMax[-1].getH_maxTau()
        
    # Compute bounds for epsilon period and strengthen/reversal period
    # Assume acceleration at G/3
    tauMin = pd
    tauMax = eps+pd
    if w>0:
        trajMin  += getNominalTraj(vmin,vlo,w*alo,tauMin,tauMax,hInitLower,isLower=True)
        trajMax  += getNominalTraj(vmax,vlo,w*ahi,tauMin,tauMax,hInitUpper,isLower=False)
        boundMin += getNominalTraj(vmin,vlo,w*alo,tauMin,tauMax,hInitLower-HP,isLower=True)
        boundMax += getNominalTraj(vmax,vlo,w*ahi,tauMin,tauMax,hInitUpper+HP,isLower=False) 
    else:
        trajMax  += getNominalTraj(vmax,vlo,w*alo,tauMin,tauMax,hInitUpper,isLower=False)
        trajMin  += getNominalTraj(vmin,vlo,w*ahi,tauMin,tauMax,hInitLower,isLower=True)
        boundMax += getNominalTraj(vmax,vlo,w*alo,tauMin,tauMax,hInitUpper+HP,isLower=False)
        boundMin += getNominalTraj(vmin,vlo,w*ahi,tauMin,tauMax,hInitLower-HP,isLower=True)
        
        
    tauMin = eps+pd
    tauMax = tf
    hInitLower = trajMin[-1].getH_maxTau()
    hInitUpper = trajMax[-1].getH_maxTau()
    vLower = trajMin[-1].getV_maxTau()
    vUpper = trajMax[-1].getV_maxTau()
    if w>0:
        trajMin  += getNominalTraj(vLower,vlos,ws*G/3,tauMin,tauMax,hInitLower,isLower=True)
        trajMax  += getNominalTraj(vUpper,vlor,wr*G/3,tauMin,tauMax,hInitUpper,isLower=False)
        boundMin += getNominalTraj(vLower,vlos,ws*G/3,tauMin,tauMax,hInitLower-HP,isLower=True)
        boundMax += getNominalTraj(vUpper,vlor,wr*G/3,tauMin,tauMax,hInitUpper+HP,isLower=False)
    else:
        trajMax  += getNominalTraj(vUpper,vlos,ws*G/3,tauMin,tauMax,hInitUpper,isLower=False)
        trajMin  += getNominalTraj(vLower,vlor,wr*G/3,tauMin,tauMax,hInitLower,isLower=True)
        boundMax += getNominalTraj(vUpper,vlos,ws*G/3,tauMin,tauMax,hInitUpper+HP,isLower=False)
        boundMin += getNominalTraj(vLower,vlor,wr*G/3,tauMin,tauMax,hInitLower-HP,isLower=True)
        
    
    ## Clean up the bounds
    boundMin, boundMax, tTrunc = truncateBounds(boundMin,boundMax)
    trajMin, trajMax = truncBounds(trajMin, trajMax, tTrunc)
    
    ## Assume worst case when tau=0
    if worstCase:
        boundMin = monotonic(boundMin,increase=1)
        boundMax = monotonic(boundMax,increase=-1)
    
    ## Return list of bounds and two trajectories
    return trajMin,trajMax,boundMin,boundMax


'''
Remove bounds from bnds1 and bnds2 where the boundaries overlap
'''
def removeOverlap(bnds1,bnds2):
    finalBnds1 = []; finalBnds2 = []
    ind1 = 0; ind2 = 0;
    while ind1<len(bnds1) and ind2<len(bnds2):
        if bnds1[ind1].overlaps(bnds2[ind2]):
            if bnds1[ind1].getMaxTau()==bnds2[ind2].getMaxTau():
                ind1+=1; ind2+=1
            elif bnds1[ind1].getMaxTau()<bnds2[ind2].getMaxTau():
                bnds2[ind2].setMinTau(bnds1[ind1].getMaxTau())
                ind1+=1
            else:
                bnds1[ind1].setMinTau(bnds2[ind2].getMaxTau())
                ind2+=1
        else:
            finalBnds1+= [bnds1[ind1]]
            finalBnds2+= [bnds2[ind2]]
            ind1+=1; ind2+=1
    finalBnds1+=bnds1[ind1:]
    finalBnds2+=bnds2[ind2:]
    return finalBnds1, finalBnds2
            

def getRegionToCheck(pra,ra,v,pd=0,eps=1,approx=False,overApprox=True,dt=1.0,dti=0.125):
    advisory = advisoryParams(ra,pd=pd,eps=eps)
    _,_,boundMin,boundMax = getSafeable(advisory,v)
    uBoundMax, uBoundMin  = getUnsafeRegion(pra,v,pd,eps)
    for bnd in uBoundMin: bnd.reverse()
    for bnd in uBoundMax: bnd.reverse()
        
    boundMin, uBoundMax = removeOverlap(boundMin,uBoundMax)
    boundMax, uBoundMin = removeOverlap(boundMax,uBoundMin)
    
    if approx:
        return linearizeRegion(boundMin+uBoundMin, boundMax+uBoundMax, overApprox,dt,dti)
    return boundMin+uBoundMin, boundMax+uBoundMax

def linearizeRegion(boundsMin,boundsMax,overApprox=True,dt=1.0,dti=0.125):
    if overApprox: return overApproximate(boundsMin,dt,dti), overApproximate(boundsMax,dt,dti)
    return underApproximate(boundsMin,dt,dti), underApproximate(boundsMax,dt,dti)

'''
Get the region where all advisories are unsafeable
Given a previous RA (pra) and velocity range, this function computes the safeable region for all
possible RAs and reduces the bounds to the minimum set of bounds where all advisories are unsafe
'''
def getUnsafeRegion(pra,v,pd=0,eps=1):
    # Loop through each possible advisory
    possAdv = getPossibleAdvisories(pra)
    _,_,boundMin,_ = getSafeable(advisoryParams(possAdv[-1],pd=pd,eps=eps),v)
    _,_,_,boundMax = getSafeable(advisoryParams(possAdv[-2],pd=pd,eps=eps),v)          
    boundMin, boundMax,_ = truncateBounds(boundMin,boundMax)
    return boundMin,boundMax

def plotAll(pra,v,pd=0,eps=1,name="",withTraj=False,axes=None,plotUnsafe=False):
    if axes==None:
        fig = plt.figure(figsize=(7,5),dpi=160,)
        axes = fig.add_subplot(111)
        
    cs = ['brown','r','b','g','c','m','orange','purple','lime']
    for ra in getPossibleAdvisories(pra):
        adv = advisoryParams(ra,pd=pd,eps=eps)
        trajMin,trajMax,boundMin,boundMax = getSafeable(adv,v)
        axes = plotBounds(boundMin,boundMax,trajMin,trajMax,withTraj=False,axes=axes,c=cs[ra],label=advName(ra))
    if plotUnsafe:
        boundMin,boundMax = getUnsafeRegion(pra,v,pd,eps)
        axes = plotBounds(boundMin,boundMax,withTraj=False,axes=axes,c='k--',label='Unsafe')
    
    
    axes.legend()
    return axes
        
def plotGetRegionToCheck(pra,ra,v,pd=0,eps=1,name="",label="",c='k',axes=None, approx=False,overApprox=True,dt=1.0,dti=0.125):
    if axes==None:
        fig = plt.figure(figsize=(7,5),dpi=160,)
        axes = fig.add_subplot(111)
        
    boundMin,boundMax = getRegionToCheck(pra,ra,v,pd,eps,approx=approx,overApprox=overApprox,dt=dt,dti=dti)
    axes = plotBounds(boundMin,boundMax,axes=axes,c=c,name=name,label=label)
    return axes    

'''
Plot the given set of bounds and trajectories
Can pass in a set of axes to add the plot to, otherwise a new set of axes will be made
'''
def plotBounds(boundsMin,boundsMax,trajMin=[],trajMax=[],name="",v=[],withTraj=True,c='k',axes=None,label=None):
    
    # Generate axes
    if axes==None:
        fig = plt.figure(figsize=(7,5),dpi=160,)
        axes = fig.add_subplot(111)
        
    # Plot bounds
    for b in boundsMin:
        t,h = b.getLine()
        axes.plot(t,h,c,label=label)
        label=None # We don't need to label every segment in bound set, just the first one
        
    for b in boundsMax:
        t,h = b.getLine()
        axes.plot(t,h,c,label=label)
        label=None # We don't need to label every segment in bound set, just the first one
        
    # Plot trajectories
    if withTraj:
        for b in trajMin:
            t,h = b.getLine()
            axes.plot(t,h,'b--')

        for b in trajMax:
            t,h = b.getLine()
            axes.plot(t,h,'b--')
    
    # Format axes labels and title
    plt.ylabel(r'$h$ (ft)')
    plt.xlabel(r'$\tau$ (sec)')
    title=""
    if name:
        title=name
    if len(v)==2:
        if len(title)>0:
            title+=", "
        title += "Climbrates: [%.1f ft/s, %.1f ft/s]"%(v[0],v[1])
    plt.title(title)
    return axes

'''
Plot a color plot of the safeable advisory (green), the unsafeable region (red), and the region where the state is safeable but not with this advisory (white). The white region is where we wish to check if the advisory is ever given
'''
def plotSafeableAdvisory(pra,ra,v,pd=0,eps=1,nbin=250,withTraj=True):
    
    cmap = ListedColormap(['green', 'white', 'red'])
    cmap_bnds = [0,1,2]
    
    # Get safeable bounds
    advisory = advisoryParams(ra,pd=pd,eps=eps)
    trajMin,trajMax,boundMin,boundMax = getSafeable(advisory,v)
        
    # Get bounds for unsafeable region
    uBoundMin,uBoundMax = getUnsafeRegion(pra,v,pd=pd,eps=eps)

    # Select bounds for plot
    minH = boundMin[0].getH_minTau()
    maxH = boundMax[0].getH_minTau()
    maxTau = boundMax[-1].getMaxTau()
    extent = 0, maxTau+1, minH-50, maxH+50

    # Generate a heat map using the bounds for each set
    heatMap = np.zeros((nbin,nbin))
    for i,h in enumerate(np.linspace(extent[2],extent[3],nbin)):
        for j,tau in enumerate(np.linspace(extent[0],extent[1],nbin)):

            # Due to matplotlib's axes layout, we have to reverse the x axis
            if satisfiesBounds(uBoundMin+uBoundMax,h,tau):
                heatMap[nbin-1-i,j]=1
            elif satisfiesBounds(boundMin+boundMax,h,tau):
                heatMap[nbin-1-i,j]=2

    # Plot the heat map and bounds
    fig = plt.figure(figsize=(7,5),dpi=160,)
    axes = fig.add_subplot(111)

    axes.imshow(heatMap, cmap=cmap,interpolation='nearest',extent=extent,aspect='auto')
    plotBounds(boundMin,boundMax,trajMin=trajMin,trajMax=trajMax,v=v,c='k',withTraj=withTraj,axes=axes)
    plotBounds(uBoundMin,uBoundMax,name="pRA: "+advName(pra)+", RA: " + advName(ra),v=v,c='k',withTraj=False,axes=axes)
    
    
    
    
    
"""
This process splits the bound list into pairs of bounds for small region of tau, where each pair has an upper and lower bound.

This function enforces tau<=6 to be tau=6, and redundant regions are merged when possible
"""
def splitBnds(bnds):
    bListAll = []
    
    # Get all unqiue taus
    taus = np.unique([np.ceil(10e6*b.getMinTau())/10e6 for b in bnds]+[np.ceil(10e6*b.getMinTau())/10e6 for b in bnds]+[6.0])
    # Loop through taus and add upper and lower bound pairs
    for tauMin,tauMax in zip(taus[:-1],taus[1:]):
        pairs=[]
        relevantBnds = [b for b in bnds if b.inRange(tauMin)]
        
        # If four bounds share a tau, there are actually two regions to check! This sorts the bounds and 
        # pairs them up to define the two regions
        if len(relevantBnds)==4:
            relevantBnds = sorted(relevantBnds, key=lambda k: -k.getH(tauMin*0.5+tauMax*0.5))
            pairs += [(relevantBnds[:2])]
            pairs += [(relevantBnds[2:])]
            
        # If there are two, then we have the pair, assuming one is an upper bound and one is lower
        elif len(relevantBnds)==2:
            relevantBnds = sorted(relevantBnds, key=lambda k: -k.getH(tauMin*0.5+tauMax*0.5))
            pairs += [(relevantBnds)]
            
        # Add bound pairs to list and enforce tau=6 to be the minimum allowed tau
        for pair in pairs:
            bList = []
            for b in pair:
                if tauMax<=6.0:
                    if b.isLower:
                        bList+=[Bound(6.0,6.0,[min(b.getH(tauMin),b.getH(tauMax)),0,0],b.isLower)]
                    else:
                        bList+=[Bound(6.0,6.0,[max(b.getH(tauMin),b.getH(tauMax)),0,0],b.isLower)]
                else:
                    b2 = b.copy()
                    b2.setMinTau(tauMin)
                    b2.setMaxTau(tauMax)
                    bList += [b2]
            bListAll += [(bList)]
            
    # Merge regions when possible for tau=6
    mergePairs(bListAll)
    return bListAll

"""
Loops through bounds and merges regions when possible to remove redundancy
"""
def mergePairs(bndPairs):
    originalNum = len(bndPairs)+1
    newNum = len(bndPairs)
    
    # Continue trying to combine bounds until there is no change
    while newNum<originalNum:
        i=1
        
        # Loop through pairs and try to merge them
        # If a bound can be merged, replace the first bound 
        # with the new bound and remove the redundant second bound
        while i<len(bndPairs):
            pair = bndPairs[i]
            for j,b in enumerate(bndPairs[:i]):
                newPair = tryMerge(b,pair)
                if newPair is not None:
                    bndPairs.remove(pair)
                    bndPairs[j]=newPair
                    i-=1
                    break;
            i+=1
        originalNum = newNum
        newNum = len(bndPairs)
        
"""
Attempts to merge pairs of bounds when tau=6 for min and max

A new bound is returned if the merge is possible, otherwise returns None
"""
def tryMerge(pair1,pair2):
    u1 = pair1[0]
    u2 = pair2[0]
    l1 = pair1[1]
    l2 = pair2[1]
    
    if u1.getMaxTau()==6.0 and u2.getMaxTau()==6.0:
        
        # If the regions overlap, then the regions can be merged
        if u1.getH_maxTau()>=u2.getH_maxTau() and l1.getH_maxTau()<=u2.getH_maxTau():
            if l1.getH_maxTau()>l2.getH_maxTau():
                l1 = l2.copy()
            return [u1,l1]
        elif u1.getH_maxTau()<u2.getH_maxTau() and u1.getH_maxTau()>l2.getH_maxTau():
            u1 = u2.copy()
            if l1.getH_maxTau()>l2.getH_maxTau():
                l1 = l2.copy()
            return [u1,l1]
    return None