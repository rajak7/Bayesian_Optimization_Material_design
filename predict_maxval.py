import numpy as np
from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C ,WhiteKernel as Wht,Matern as matk
from sklearn.gaussian_process.kernels import RationalQuadratic as expker
from sklearn.metrics import mean_squared_error as MSError


inputmap=dict()
ninputmap=dict()
totfea_atom=2
natom_layer=6*totfea_atom

def createinputmap(inputmap,ninputmap,totfea_atom):
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
#            ninputmap[keys].append((inputmap[keys][ii]-Xmin[ii])/(Xmax[ii]-Xmin[ii]))   # normalized by by (tt-xmax)/(xmax-xmin)
            ninputmap[keys].append((inputmap[keys][ii]-Xmean[ii])/Xstd[ii])
    #print the final keys:
#    for keys in inputmap:
#        print("key :", keys,inputmap[keys])
#    for keys in ninputmap:
#        print("nkey :", keys, ninputmap[keys])



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
#            Ydata = np.ndarray(shape=(ndata, 1), dtype=float)
            Ydata = np.empty(ndata,dtype=float)
            itag=1
        else :
            lines = lines.replace("\n", "").split()
#            print(lines,lines.__len__())
            count+=1
            for ii in range(0,lines.__len__()-1):
                jj=lines[ii]
 #               print(3*ii,3*ii+1,3*ii+2,jj,inputmap[jj][0],inputmap[jj][1],inputmap[jj][2])
                Xdata[count][2 * ii] = ninputmap[jj][0]
                Xdata[count][2 * ii + 1] = ninputmap[jj][1]
#                Xdata[count][3 * ii + 2] = ninputmap[jj][2]
            Ydata[count] = float(lines[lines.__len__() - 1])

    #print the entire dataset
#    for ii in range(0,ndata):
#        print("data: ",ii,Xdata[ii][:],Ydata[ii])
    return Xdata,Ydata,ndata


def gpregression(Xtrain,Ytrain,Xtest,Ytest,ntrain,ntest):
    print("regression")
    cmean=[1.0]*12
    cbound=[[1e-3, 1000]]*12
#    kernel = C(1.0, [1e-3, 1e3]) * RBF(cmean, cbound)
    kernel = C(1.0, (1e-3, 1e3)) * matk(cmean, cbound, 1.5)+ Wht(1.0, (1e-3, 1e3))
#    kernel = C(1.0, (1e-3, 1e3)) * matk(1, (1e-05, 1000.0), 2.5) + Wht(1.0, (1e-3, 1e3))+ C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
#    kernel = C(1.0, (1e-3, 1e3)) * matk(1, (1e-05, 1000.0), 2.5)+C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))

    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=100, normalize_y=False)
    gp.fit(Xtrain, Ytrain)
    print("initial parameters:", kernel)
    print("optimal parameters:", gp.kernel_, "likelihood:", gp.log_marginal_likelihood(gp.kernel_.theta))
    y_pred, sigma = gp.predict(Xtest, return_std=True)
    dataorder=np.argsort(Ytest)
    tYest=Ytest[dataorder]
    ty_pred=y_pred[dataorder]
    tsigma=sigma[dataorder]
    del Ytest,y_pred,sigma
    Ytest=tYest
    y_pred=ty_pred
    sigma=tsigma
    toterr=0.0
    for val in range(0,ntest):
#        print("Prediction: ",Ytest[val],"  ",y_pred[val]," ",sigma[val])
        toterr+=np.abs(Ytest[val]-y_pred[val])
    print("toterr  prediction loss : ",toterr,toterr/float(ntest))
    fig = plt.figure(figsize=(14,10))
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('font', weight='bold')
    xxdummy=range(ntest)
    plt.plot(xxdummy, Ytest, 'r-', linewidth=3.5, label=u'True Value')
    plt.plot(xxdummy, Ytest, 'r.', markersize=20)
    plt.plot(xxdummy, y_pred, 'b--', linewidth=3.5, label=u'Prediction')
    plt.plot(xxdummy, y_pred, 'b.', markersize=20)
    plt.fill(np.concatenate([xxdummy, xxdummy[::-1]]),np.concatenate([y_pred - 1.9600 * sigma,(y_pred + 1.9600 * sigma)[::-1]]),alpha=.5, fc='y', ec='None', label='95% confidence interval')
    plt.xlabel('tri-layer structure',fontsize=40, fontweight='bold')
    plt.ylabel('Band GAP',fontsize=40, fontweight='bold')
    plt.legend(loc='upper left', ncol=1, fancybox=True, shadow=True, prop={'size': 20})
#    plt.legend(loc='upper left')
    plt.title("TEST DATA",fontsize=40,fontweight='bold')
    #-----training set-----
    yt_pred, tsigma = gp.predict(Xtrain, return_std=True)
#    for val in range(0,ntrain):
#        print("Training set: ",Ytrain[val],"  ",yt_pred[val]," ",tsigma[val])
    print("Total training errror: ",np.sqrt(MSError(Ytrain,yt_pred)))
    print("Total prediction errror: ", np.sqrt(MSError(Ytest,y_pred)))
#    xxtdummy=range(ntrain)
#    plt.plot(xxtdummy, Ytrain, 'r-', markersize=10, label=u'Observations')
#    plt.plot(xxtdummy, Ytrain, 'r.', markersize=10)
#    plt.plot(xxtdummy, yt_pred, 'b-', markersize=10, label=u'Prediction')
#    plt.plot(xxtdummy, yt_pred, 'b.', markersize=10)
#    plt.fill(np.concatenate([xxtdummy, xxtdummy[::-1]]),np.concatenate([yt_pred - 1.9600 * tsigma,(yt_pred + 1.9600 * tsigma)[::-1]]),alpha=.8, fc='b', ec='None', label='95% confidence interval')
#    plt.xlabel('$x$')
#    plt.ylabel('$f(x)$')
#    plt.legend(loc='upper left')
#    plt.title("Training data")
    plt.ylim(-0.6,1.6)
    plt.show()
#    plt.savefig('fig1a.png')
#    plt.close()
    return

#------- Program Starts from here -------------

createinputmap(inputmap,ninputmap,totfea_atom)
Xdata,Ydata,ndata=readinput("177data.txt",natom_layer)
print("Original Training and Y :",np.shape(Xdata),np.shape(Ydata))
print("Transpose Training and Y : ",np.shape(np.transpose(Xdata)),np.shape(np.transpose(Ydata)))
print("Original Training and Y :",np.shape(Xdata),np.shape(Ydata))

ntrain=int(0.40*ndata)
ntest=ndata-ntrain
print("Total training and Test Data: ",ntrain,ntest)
for ii in range(0,1):
    dataset=np.random.permutation(ndata)
    a1data=np.empty(ntrain, dtype=int)
    a2data=np.empty(ntest, dtype=int)
    a1data[:]=dataset[0:ntrain]
    a2data[:]=dataset[ntrain:ndata]
#    print("ii=",ii," : ",dataset)
#    print("a1: ", a1data,a1data.__len__())
#    print("a2: ", a2data,a2data.__len__())
    Xtrain=np.ndarray(shape=(ntrain, natom_layer), dtype=float)
    Ytrain = np.empty(ntrain, dtype=float)
    Xtest = np.ndarray(shape=(ntest, natom_layer), dtype=float)
    Ytest = np.empty(ntest, dtype=float)
    for itrain in range(0,ntrain):
        mm=a1data[itrain]
        Xtrain[itrain][:]=Xdata[mm][:]
        Ytrain[itrain]=Ydata[mm]
    for itest in range(0,ntest):
        mm = a2data[itest]
        Xtest[itest][:]=Xdata[mm][:]
        Ytest[itest]=Ydata[mm]
    gpregression(Xtrain,Ytrain,Xtest,Ytest,ntrain,ntest)
    del Xtrain,Ytrain
    del Xtest,Ytest



