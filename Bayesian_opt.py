import sys
import numpy as np
from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C ,WhiteKernel as Wht,Matern as matk
from sklearn.gaussian_process.kernels import RationalQuadratic as expker
from sklearn.metrics import mean_squared_error as MSError
from scipy.stats import norm

inputmap=dict()
ninputmap=dict()
totfea_atom=2              #total number of atoms per layer
n_3layer_atoms=6           # number of atoms in 3 layer
natom_layer=n_3layer_atoms*totfea_atom   #total number of features
Niteration = 30            # number of iteration in a given Bayesian  Optimization
#input parameters
train_test_split=0.10                    # initial sampled data in a given Bayesian  Optimization run
Nruns = 1                                # total number of Bayesian  Optimization runs


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
#            ninputmap[keys].append((inputmap[keys][ii]-Xmean[ii])/Xstd[ii])
    #print the final keys:
    for keys in inputmap:
        print("key :", keys,inputmap[keys])
    for keys in ninputmap:
        print("nkey :", keys, ninputmap[keys])

#read input data
def readinput(filename,natom_layer):
    inputfile=open(filename,'r')
    dataset=list()
    itag=0
    count=-1
    ndata=0
    for lines in inputfile:
        if itag==0:
            ndata=int(lines)
            Xdata = np.ndarray(shape=(ndata, natom_layer), dtype=float)
            Xinfo =  np.chararray(ndata, itemsize=20)
            Ydata = np.empty(ndata,dtype=float)
            itag=1
        else :
            lines = lines.replace("\n", "").split()
            structname=str()
            count+=1
            for ii in range(0,lines.__len__()-1):
                jj=lines[ii]
                if (ii > 0) : structname = structname + '-' + jj
                else: structname=jj
 #               print(3*ii,3*ii+1,3*ii+2,jj,inputmap[jj][0],inputmap[jj][1],inputmap[jj][2])
                Xdata[count][2 * ii] = ninputmap[jj][0]
                Xdata[count][2 * ii + 1] = ninputmap[jj][1]
                Xinfo[count]=structname
            Ydata[count] = float(lines[lines.__len__() - 1])
            print("structname: ",structname,lines,lines.__len__())
    #print the entire dataset
#    for ii in range(0,ndata):
#        print("data: ",ii,Xdata[ii][:],Ydata[ii])
    return Xdata,Ydata,Xinfo,ndata

#building a gaussian process regression model
def gpregression(Xtrain,Ytrain,Nfeature):
    cmean=[1.0]*Nfeature
    cbound=[[1e-3, 1000]]*Nfeature
    kernel = C(1.0, (1e-3, 1e3)) * matk(cmean, cbound, 1.5)

    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=40, normalize_y=False)
    gp.fit(Xtrain, Ytrain)
    return gp

#predict result using GP regression model
def gprediction(gpnetwork,xtest):
    y_pred, sigma = gpnetwork.predict(xtest, return_std=True)
    return y_pred, sigma

#compute expected improvement 
def expectedimprovement(xdata,gpnetwork,ybest,itag,epsilon):
    ye_pred, esigma = gprediction(gpnetwork, xdata)
    expI = np.empty(ye_pred.size, dtype=float)
    for ii in range(0,ye_pred.size):
        if esigma[ii] > 0:
            zzval=itag*(ye_pred[ii]-ybest)/float(esigma[ii])
            expI[ii]=itag*(ye_pred[ii]-ybest-epsilon)*norm.cdf(zzval)+esigma[ii]*norm.pdf(zzval)
        else:
            expI[ii]=0.0
    return expI

#Bayesian optimization run
def numberofopt(Xdata,Ydata,Xinfo,ndata,natom_layer,totfea_atom):
    itag = 1
    epsilon = 0.1
    ntrain = int(train_test_split * ndata)
    nremain = ndata - ntrain
    dataset = np.random.permutation(ndata)
    a1data = np.empty(ntrain, dtype=int)
    a2data = np.empty(nremain, dtype=int)
    a1data[:] = dataset[0:ntrain]
    a2data[:] = dataset[ntrain:ndata]
    # info for the initial training set
    Xtrain = np.ndarray(shape=(ntrain, natom_layer), dtype=float)
    Xtraininfo = np.chararray(ntrain, itemsize=20)
    Ytrain = np.empty(ntrain, dtype=float)
    Xtrain[0:ntrain, :] = Xdata[a1data, :]
    Xtraininfo[0:ntrain] = Xinfo[a1data]
    Ytrain[0:ntrain] = Ydata[a1data]
    yopttval = np.max(Ytrain)
    xoptval = Xtraininfo[np.argmax(Ytrain)]
    yoptstep=0
    yopinit = yopttval
    xoptint = xoptval
    # info for the remaining data set
    Xremain = np.ndarray(shape=(nremain, natom_layer), dtype=float)
    Xremaininfo = np.chararray(nremain, itemsize=20)
    Yremain = np.empty(nremain, dtype=float)
    Xremain[0:nremain, :] = Xdata[a2data, :]
    Xremaininfo[0:nremain] = Xinfo[a2data]
    Yremain[0:nremain] = Ydata[a2data]
    print("Xremain: ", Xremain.shape)
    print("Yremain:", Yremain.shape)
    print("Xremaininfo: ", Xremaininfo.shape)
    print("Initial max value 0th run : ", xoptval, yopttval)
    print("Total number of inital training points: ", ntrain)
    # print("Xtrain: ",Xtrain)
    for ii in range(0, Niteration):
        if ii > int(0.5*Niteration):
            epsilon=0.01
            print("updated epsilon: ",epsilon)
        gpnetwork = gpregression(Xtrain, Ytrain, natom_layer)
        yt_pred, tsigma = gprediction(gpnetwork, Xtrain)
        ybest = np.max(yt_pred)
        ybestloc = np.argmax(yt_pred)
        print("current Best in iteration ii", ii + 1, " is ", ybest, "for the structure: ", Xtraininfo[ybestloc])
        if yopttval < ybest:
            yopttval = ybest
            xoptval = Xtraininfo[ybestloc]
        print("Best Strucutre so far", yopttval, "for the structure: ", xoptval)
        expI = expectedimprovement(Xremain, gpnetwork, ybest, itag, epsilon)
        expImax = np.max(expI)
        expimaxloc = np.argmax(expI)
        print("Next Structure to evaluate has expI ", expImax, "for the structure: ", Xremaininfo[expimaxloc],
              "has Y: ", Yremain[expimaxloc])
        xnew = np.append(Xtrain, Xremain[expimaxloc]).reshape(-1, natom_layer)
        xnewinfo = np.append(Xtraininfo, Xremaininfo[expimaxloc])
        ynew = np.append(Ytrain, Yremain[expimaxloc])
        xrnew = np.delete(Xremain, expimaxloc, 0)
        xrnewinfo = np.delete(Xremaininfo, expimaxloc)
        yrnew = np.delete(Yremain, expimaxloc)
        if ii==0:
            Xexplored=Xremaininfo[expimaxloc]
            Yexplored=Yremain[expimaxloc]
        else:
            Xexploredtemp=np.append(Xexplored, Xremaininfo[expimaxloc])
            Yexploredtemp=np.append(Yexplored, Yremain[expimaxloc])
            del Xexplored,Yexplored
            Xexplored=Xexploredtemp
            Yexplored=Yexploredtemp
        #    print("Xremain info: ",xrnew.shape,yrnew.shape,xrnewinfo.shape)
        del Xtrain, Ytrain, Xremaininfo, gpnetwork
        Xtrain = xnew
        Xtraininfo = xnewinfo
        Ytrain = ynew
        Xremain = xrnew
        Xremaininfo = xrnewinfo
        Yremain = yrnew
        del xnew, xnewinfo, ynew, xrnew, xrnewinfo, yrnew

    if not yopinit==yopttval:
         yoptstep=np.argmax(Yexplored)+1
    else:
        yoptstep=0
    dataorder = np.argsort(Yexplored)
    Yexploredtemp=Yexplored[dataorder]
    Xexploredtemp = Xexplored[dataorder]
    print(Yexplored)
    Xbest=Xexploredtemp[Niteration-3:Niteration]
    Ybest=Yexploredtemp[Niteration - 3:Niteration]
    print("\n")
    print("Initial Best Strucuture: ", xoptint, "has value: ", yopinit)
    print("Final Optimal Strucuture: ", xoptval, "has value: ", yopttval,"in step: ",yoptstep)
    print("Final Best Structure 1st: ",Xbest[2],"has value: ",  Ybest[2])
    print("Final Best Structure 2st: ", Xbest[1],"has value: ", Ybest[1])
    print("Final Best Structure 2st: ", Xbest[0],"has value: ", Ybest[0])
    return xoptint,yopinit,xoptval,yopttval


#------- Program Starts from here -------------
createinputmap(inputmap,ninputmap,totfea_atom)
in_file=sys.argv[1]
Xdata,Ydata,Xinfo,ndata=readinput(in_file,natom_layer)
print("Original Training X and Y :",np.shape(Xdata),np.shape(Xdata))

Xinitguess = np.chararray(Nruns, itemsize=20)
Yinitguess = np.empty(Nruns, dtype=float)
Xoptimal = np.chararray(Nruns, itemsize=20)
Yoptimal = np.empty(Nruns, dtype=float)

for ii in range(0,Nruns):
    Xinitguess[ii], Yinitguess[ii], Xoptimal[ii], Yoptimal[ii] =numberofopt(Xdata, Ydata, Xinfo, ndata, natom_layer, totfea_atom)

print("\n-----Final Result------\n")
for ii in range(0,Nruns):
    print("Initial Best Strucuture:   ", Xinitguess[ii], "  has value:  ", Yinitguess[ii],"  Final Optimal Strucuture:   ", Xoptimal[ii], "  has value: ", Yoptimal[ii])
