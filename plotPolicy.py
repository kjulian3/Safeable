import numpy as np
import h5py
from safeable import *
from advisoryParams import *
from NNet import *

HS     = np.concatenate([np.linspace(-8000,-4000,5),np.linspace(-3000,-1250,8),np.linspace(-1000,-800,3),np.linspace(-700,-150,12),np.linspace(-100,100,9),np.linspace(150,700,12),np.linspace(800,1000,3),np.linspace(1250,3000,8),np.linspace(4000,8000,5)])

VOWNS  = np.concatenate([np.linspace(-100,-60,5),np.linspace(-50,-35,4),np.linspace(-30,30,21),np.linspace(35,50,4),np.linspace(60,100,5)])
VINTS  = VOWNS
TAUS   = np.linspace(0,40,41)
RESPS  = np.array([0.0,1.0])


'''
Plot safeable region over the neural network policy
The average of the ownship velocity range is used to compute policy heat map
SAT points are highlighted in RED
SAT highlighting can be turned of by setting highlightSAT to False
'''
def plotSafeablePolicy(pra,ra,v,pd=0,eps=1,nbin=250,highlightSAT=True,axes=None,pltBnds=True,useTable=False,drawLegend=True,dt=0.125,dti=0.0625):
    # Get safeable bounds
    bndsMin,bndsMax = getRegionToCheck(pra,ra,v,pd=pd,eps=eps,dt=dt,dti=dti)
    
    # Select bounds for plot
    minH = bndsMin[0].getH_minTau()
    maxH = bndsMax[0].getH_minTau()
    maxTau = max([b.getMaxTau() for b in bndsMax])
    extent = 0, maxTau+1, minH-50, maxH+150

    # Plot the heat map and bounds
    if axes is None:
        fig = plt.figure(figsize=(7,5),dpi=160,)
        axes = fig.add_subplot(111)

    # Pass in bounds to policy plotter if we want to highlight SAT points
    if highlightSAT:
        plotPolicy(pra,extent,np.mean(v),nbin=nbin,drawLegend=drawLegend,axes=axes,ra=ra,bounds=bndsMin+bndsMax,useTable=useTable)
    else:
        plotPolicy(pra,extent,np.mean(v),nbin=nbin,drawLegend=drawLegend,axes=axes,useTable=useTable)

    # Plot bounds
    if pltBnds:
        plotBounds(bndsMin,bndsMax,name="pRA: "+advName(pra)+", RA: " + advName(ra),v=v,c='k', axes=axes)

'''
Compute heatmap of policy
'''
def getPolicy(pra,extent,vown,vint=0.,nbin=250,ra=-1,bounds=[],useTable=False):
    if useTable:
        return getPolicy_Table(pra,extent,vown,vint,nbin,ra,bounds)
    return getPolicy_NNet(pra,extent,vown,vint,nbin,ra,bounds)

def getPolicy_NNet(pra,extent,vown,vint=0.,nbin=250,ra=-1,bounds=[]):
    
    # Load network of given pRA into Python
    nnet_file_name = "networks/VertCAS_pra%02d_v4_45HU_200.nnet" % (pra+1)
    net = NNet(nnet_file_name)
    
    tauMin,tauMax,hMin,hMax = extent
    
    # Get safeable bounds
    possAdv = getPossibleAdvisories(pra)
    
    # Generate a heat map using the bounds for each set
    # Use meshgrid to define array of all inputs to network
    hVec = np.linspace(hMin,hMax,nbin)
    tauVec = np.linspace(tauMin,tauMax,nbin)
    hMesh,tauMesh = np.meshgrid(hVec,tauVec)
    hMesh = hMesh.reshape(nbin**2,1); tauMesh = tauMesh.reshape(nbin**2,1)
    vownMesh = np.ones((nbin**2,1))*vown 
    vintMesh = np.ones((nbin**2,1))*vint 
    netIn = np.concatenate((hMesh,vownMesh,vintMesh,tauMesh),axis=1)
    
    # Evaluate network on all inputs
    netOut = net.evaluate_network_multiple(netIn)
    
    # Convert outputs to best advisories
    bestAdv = np.argmax(netOut,axis=1)
    heatMap = np.array([possAdv[i] for i in bestAdv])
    
    # Highlight SAT points if bounds given
    if len(bounds)>0:
        for i in range(nbin**2):
            if heatMap[i]==ra and satisfiesBounds(bounds,netIn[i,0],netIn[i,3]):
                heatMap[i]=10
                
    # Reshape the map and flip around the axes for plotting purposes
    heatMap = heatMap.reshape((nbin,nbin)).T[::-1]
    return heatMap


def getIndex(h,vown,vint,tau):
    h_ind = np.argmin(abs(h-HS))
    vown_ind = np.argmin(abs(vown-VOWNS))
    vint_ind = np.argmin(abs(vint-VINTS))
    tau_ind  = np.argmin(abs(tau-TAUS))
    return h_ind +len(HS)*vown_ind + len(HS)*len(VOWNS)*vint_ind + len(HS)*len(VOWNS)*len(VINTS)*tau_ind

def getPolicy_Table(pra,extent,vown,vint=0.,nbin=250,ra=-1,bounds=[]):
    tauMin,tauMax,hMin,hMax = extent
    possAdv = getPossibleAdvisories(pra)
    table_folder = "/raid/kyle/kyle/VertCAS"
    with h5py.File(table_folder+"/VertCAS_TrainingData_v2_%02d.h5"%(pra+1), "r") as f:
        table = f['y'].value
    
    policy = np.zeros((nbin,nbin),int)
    for i,tau in enumerate(np.linspace(tauMin,tauMax,nbin)):
        for j,h in enumerate(np.linspace(hMin,hMax,nbin)):
            policy[i,j] = possAdv[np.argmax(table[getIndex(h,vown,vint,tau)])]
            
            # Highlight SAT points if bounds given
            if len(bounds)>0 and policy[i,j]==ra and satisfiesBounds(bounds,{H:h,TAU:tau}):
                policy[i,j]=10
    
    return policy.T[::-1]


'''
Plot policy heat map
'''
def plotPolicy(pra,extent,vown,vint=0.,nbin=250,drawLegend=True,axes=None,ra=-1,bounds=[],useTable=False):
    cmap = ListedColormap(['white',       # COC
                           'cyan',        # DNC
                           'lightgreen',  # DND
                           'dodgerblue',  # DES1500
                           'lime',        # CL1500
                           'blue',        # SDES1500
                           'forestgreen', # SCL1500
                           'navy',        # SDES2500
                           'darkgreen',   # SCL2500
                           'red',         # Non-compliant (SAT) point color
                          ])

    # Get heatmap and plot on given axes or new axes if none given
    heatMap = getPolicy(pra,extent,vown,vint,nbin,ra,bounds,useTable)
    if axes is None:
        fig = plt.figure(figsize=(7,5),dpi=160,)
        axes = fig.add_subplot(111)
        
    # Show the heatmap and draw legend
    axes.imshow(heatMap, cmap=cmap,interpolation='nearest',extent=extent,aspect='auto',vmin=0, vmax=len(cmap.colors))
    if drawLegend: writeColorfulLegend(axes,cmap.colors)

'''
Draw colorful legend on given axes
'''
def writeColorfulLegend(axes,colors):
    p = matplotlib.patches.Rectangle(
    (1.01, 0.14), 0.18, 0.8,
    fill=True,facecolor="brown",alpha=0.5,edgecolor='k',linewidth=2, transform=axes.transAxes, clip_on=False
    )
    axes.add_patch(p)
    axes.text(1.03, 0.92, "LEGEND", transform=axes.transAxes, fontsize=12,
        verticalalignment='top', color='k')
    for i, c in enumerate(colors):
        if i==9:
            adv="SAT"
        else:
            adv = advName(i)
        axes.text(1.03, 0.85-i*0.06, adv, transform=axes.transAxes, fontsize=10,
        verticalalignment='top', color=c)