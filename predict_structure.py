import numpy as np
from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C ,WhiteKernel as Wht,Matern as matk
from sklearn.gaussian_process.kernels import RationalQuadratic as expker
from sklearn.metrics import mean_squared_error as MSError
from scipy.stats import norm
from base64 import b16encode
import warnings
warnings.filterwarnings("ignore")

inputmap=dict()
ninputmap=dict()
totfea_atom=2              #total number of atoms per layer
Natomsinpayer=6           # number of atoms in 3 layer
natom_layer=Natomsinpayer*totfea_atom   #total number of features
Nbandpoint=30             # number of points at which GP model is build

#input paramters
inputfile_name="3-layer-band_structure.txt"    #file name of the input data
train_test_split=0.60                    #split between training and test set
Nruns = 1

#create input feature vector of the given n-layer heterostructure
def createinputmap(inputmap,ninputmap,totfea_atom):
    #define the eletronegetivity and ionization potential of each atoms
    inputmap['Mo'] = [2.16,684.3,190.0]
    inputmap['W'] = [2.36,770.0,193.0]
    inputmap['S'] = [2.58,999.6,88.8]
    inputmap['Se'] = [2.55,941.0,103.0]
    inputmap['Te'] = [2.10,869.3,123.0]

    #normalize the input features by (tt-xmax)/(xmax-xmin)
    Xmax = np.empty(totfea_atom,dtype=float)
    Xmin = np.empty(totfea_atom, dtype=float)
    Xmean= np.empty(totfea_atom,dtype=float)
    Xstd = np.empty(totfea_atom,dtype=float)
    Xmax.fill(0.0)
    Xmin.fill(10000.0)
    Xmean.fill(0.0)
    Xstd.fill(0.0)
    nfeatures=0
    for keys in inputmap:
        nfeatures+=1
        for ii in range(0,totfea_atom):
            if Xmax[ii] < inputmap[keys][ii]: Xmax[ii]=inputmap[keys][ii]
            if Xmin[ii] > inputmap[keys][ii]: Xmin[ii]=inputmap[keys][ii]
            Xmean[ii]+=inputmap[keys][ii]
    for ii in range(0,totfea_atom):
        Xmean[ii]=Xmean[ii]/float(nfeatures)
    for keys in inputmap:
        for ii in range(0, totfea_atom):
            Xstd[ii]+=(inputmap[keys][ii]- Xmean[ii])*(inputmap[keys][ii]- Xmean[ii])
    for ii in range(0, totfea_atom):
        Xstd[ii]=np.sqrt(Xstd[ii]/float(nfeatures))
    print("Xmax and  Xmin: ",Xmax,Xmin)
    print("Xmean  and Xstd: ",Xmean,Xstd)
    for keys in inputmap:
        ninputmap[keys]=list()
        for ii in range(0, totfea_atom):
            ninputmap[keys].append((inputmap[keys][ii]-Xmin[ii])/(Xmax[ii]-Xmin[ii]))   # normalized by by (tt-xmax)/(xmax-xmin)

#read input data
def readinput(filename,natom_layer,Natomsinpayer,Nbandpoint):
    inputfile=open(filename,'r')
    itag=0
    count=-1
    ndata=0
    for lines in inputfile:
        if itag==0:
            ndata=int(lines)
            Xdata = np.ndarray(shape=(ndata, natom_layer), dtype=float)
            Xinfo =  np.chararray(ndata, itemsize=20)
            Ytopdata = np.ndarray(shape=(ndata,Nbandpoint),dtype=float)
            Ybotdata = np.ndarray(shape=(ndata, Nbandpoint), dtype=float)
            itag=1
        else :
            lines = lines.replace("\n", "").split()
            structname=str()
            count+=1
            for ii in range(0,Natomsinpayer):
                jj=lines[ii]
                if (ii >0): structname = structname + '-' + jj
                if(ii==0): structname=jj
                Xdata[count][2 * ii] = ninputmap[jj][0]
                Xdata[count][2 * ii + 1] = ninputmap[jj][1]
                Xinfo[count]=structname
            for itop in range(Natomsinpayer,Natomsinpayer+Nbandpoint):
                Ytopdata[count][itop-Natomsinpayer]=float(lines[itop])
            for ibot in range(Natomsinpayer+Nbandpoint,Natomsinpayer+Nbandpoint+Nbandpoint):
                Ybotdata[count][ibot - Natomsinpayer-Nbandpoint] = float(lines[ibot])
#            print("structname: ",structname,lines)
#            print("datatop : ",lines[Natomsinpayer:Natomsinpayer+Nbandpoint])
#            print("databottom: ",lines[Natomsinpayer+Nbandpoint:lines.__len__()])
    return Xdata,Ytopdata,Ybotdata,Xinfo,ndata

#read x-axis value of the 30 points at which gp model is build
def readxaxisval(filename,Nbandpoint):
    inputfile = open(filename, 'r')
    Xdata = np.empty(Nbandpoint, dtype=float)
    count=-1
#    Xdata=range(0,Nbandpoint)
    for val in inputfile:
        count+=1
        Xdata[count]=float(val)
    return Xdata

def plotbandstructure(XX,structure,YY1,YY2):
    fig = plt.figure(figsize=(7, 7))
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('font', weight='bold')
    plt.plot(XX, YY1, 'b-', linewidth=3.5, label=u'HOMO')
    plt.plot(XX, YY2, 'r-', linewidth=3.5, label=u'LOMO')
    plt.title(structure, fontsize=20, fontweight='bold')
    plt.show()

#make polt of the CBM/VBM for the test set inside a folder Bandstructure
def plotbandtest(XX,structure,YY1true,YY2true,YY1predict,YY2predict,YY1sigma,YY2sigma,myid):
    fig = plt.figure(figsize=(10,10))
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('font', weight='bold')
    plt.plot(XX, YY1true, 'b-', linewidth=3.5, label=u'CBM-Ground-Truth')
    plt.plot(XX, YY1predict, 'c--', linewidth=3.5)
    #plt.fill(np.concatenate([XX, XX[::-1]]),np.concatenate([YY1predict - 1.9600 * YY1sigma, (YY1predict + 1.9600 * YY1sigma)[::-1]]), alpha=.3, fc='y', ec='None',label='95% confidence interval')
    plt.plot(XX, YY2true, 'r-', linewidth=3.5, label=u'VBM-Ground-Truth')
    plt.plot(XX, YY2predict, 'm--', linewidth=3.5)
    #plt.fill(np.concatenate([XX, XX[::-1]]),np.concatenate([YY2predict - 1.9600 * YY2sigma, (YY2predict + 1.9600 * YY2sigma)[::-1]]), alpha=.3, fc='g',ec='None', label='95% confidence interval')
    plt.title(structure, fontsize=20, fontweight='bold')
    plt.legend(loc='upper right', bbox_to_anchor=(0.28, 1.16), ncol=1, fancybox=True, shadow=True, prop={'size': 14})
    plt.xlabel('Wave Vector',fontsize=20, fontweight='bold')
    plt.ylabel('Energy(eV)',fontsize=20, fontweight='bold')
    imagefile = "Bandstructure/Strucuture" + str(myid+1)
    plt.savefig(imagefile)
#    plt.show()

#Build GP regression model
def gpregression(Xtrain,Ytrain,Nfeature):
    cmean=[1.0]*Nfeature
    cbound=[[1e-3, 10000]]*Nfeature
#    kernel = C(1.0, [1e-3, 1e3]) * RBF(cmean, cbound)
    kernel = C(1.0, (1e-3, 1e3)) * matk(cmean, cbound, 2.5)+ Wht(1.0, (1e-3, 1e3))

    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=50, normalize_y=False)
    gp.fit(Xtrain, Ytrain)
    print("initial parameters:", kernel)
    print("optimal parameters A:", gp.kernel_, "likelihood:", gp.log_marginal_likelihood(gp.kernel_.theta))
#    print("likelihood:", gp.log_marginal_likelihood(gp.kernel_.theta))
    return gp

def gprediction(gpnetwork,xtest):
    y_pred, sigma = gpnetwork.predict(xtest, return_std=True)
    return y_pred, sigma

#------- Program Starts from here -------------
createinputmap(inputmap,ninputmap,totfea_atom)
Xdata,Ytopdata,Ybotdata,Xinfo,ndata=readinput(inputfile_name,natom_layer,Natomsinpayer,Nbandpoint)
xaxisval=readxaxisval("xaxisvalue.txt",Nbandpoint)

for ii in range(0,ndata):
    materialname = Xinfo[ii].decode('utf-8')
    print("Xinfo: ",materialname)
 #   plotbandstructure(xaxisval,materialname,Ytopdata[ii,:].ravel(),Ybotdata[ii,:].ravel())

# Make  a regression model for each of the 60 points
ntrain=int(train_test_split*ndata)
ntest=ndata-ntrain
print("Total training and Test Data: ",ntrain,ntest)
for ii in range(0,Nruns):
    dataset = np.random.permutation(ndata)
    a1data=np.empty(ntrain, dtype=int)
    a2data=np.empty(ntest, dtype=int)
    a1data[:]=dataset[0:ntrain]
    a2data[:]=dataset[ntrain:ndata]
    #create the training set
    Xtrain = np.ndarray(shape=(ntrain, natom_layer), dtype=float)
    Xtraininfo = np.chararray(ntrain, itemsize=20)
    Ytoptrain = np.ndarray(shape=(ntrain, Nbandpoint), dtype=float)
    Ybottrain = np.ndarray(shape=(ntrain, Nbandpoint), dtype=float)
    Xtrain[0:ntrain, :] = Xdata[a1data, :]
    Xtraininfo[0:ntrain] = Xinfo[a1data]
    Ytoptrain[0:ntrain,:] = Ytopdata[a1data,:]
    Ybottrain[0:ntrain, :] = Ybotdata[a1data, :]
    #create the test set
    Xtest = np.ndarray(shape=(ntest, natom_layer), dtype=float)
    Xtestinfo = np.chararray(ntest, itemsize=20)
    Ytoptest = np.ndarray(shape=(ntest, Nbandpoint), dtype=float)
    Ytoppredict = np.ndarray(shape=(ntest, Nbandpoint), dtype=float)
    Ytopsigma = np.ndarray(shape=(ntest, Nbandpoint), dtype=float)
    Ybottest = np.ndarray(shape=(ntest, Nbandpoint), dtype=float)
    Ybotpredict = np.ndarray(shape=(ntest, Nbandpoint), dtype=float)
    Ybotsigma = np.ndarray(shape=(ntest, Nbandpoint), dtype=float)
    Xtest[0:ntest, :] = Xdata[a2data, :]
    Xtestinfo[0:ntest] = Xinfo[a2data]
    Ytoptest[0:ntest,:] = Ytopdata[a2data,:]
    Ybottest[0:ntest, :] = Ybotdata[a2data, :]
    gptopmodel=list()
    gpbotmodel=list()
    for jj in range(0,Nbandpoint):
        Ytemp1=Ytoptrain[:,jj]
        Ytemp2 = Ybottrain[:, jj]
        print("Point number: ",jj)
        gptemp1 = gpregression(Xtrain,Ytemp1.ravel(),natom_layer)
        gptemp2 = gpregression(Xtrain, Ytemp2.ravel(), natom_layer)
        gptopmodel.append(gptemp1)
        gpbotmodel.append(gptemp2)
        del gptemp1,gptemp2
    #Test the Model
    for jj in range(0,Nbandpoint):
        Ytoppredict[:,jj], Ytopsigma[:,jj] = gprediction(gptopmodel[jj], Xtest)
        Ybotpredict[:,jj], Ybotsigma[:,jj] = gprediction(gpbotmodel[jj], Xtest)
    #Make plots of the all the test cases
    for kk in range(0,ntest):
        materialname = Xtestinfo[kk].decode('utf-8')
        plotbandtest(xaxisval, materialname, Ytoptest[kk, :].ravel(), Ybottest[kk, :].ravel(),Ytoppredict[kk, :].ravel(), Ybotpredict[kk, :].ravel(),Ytopsigma[kk, :].ravel(), Ybotsigma[kk, :].ravel(),kk)
