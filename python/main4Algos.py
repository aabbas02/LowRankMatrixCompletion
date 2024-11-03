import numpy as np
import time
from dask.distributed import Client, LocalCluster
from altMin import altMinGD, altMinFullyDist, altMinFedCol
from functions4Algos import factGD, altGDMinFedDaskScttrSparseNew
from utils import plotErrAgnstTime, plotTimeAgnstWrkrs, getDataPtr, fedSVD
import random
import os
from scipy.sparse.linalg import svds
from dask.distributed import Client
from scipy.sparse import csr_matrix
import scipy
import pickle, coiled
if __name__ == '__main__':  
    r = 2
    n = 500
    q = 1000
    eta_c = 1.0
    p = 0.1
    ID  = np.random.randint(low=0, high=999)
    numWrkrs_ = [5,10]
    MC = 2
    T = 1001
    lim = 1e-13
    T_in = 10
    #
    #timeAltMinFedRow = -1*np.ones( (MC,len(numWrkrs_)) )
    #timeFedVAltMinRow = -1*np.ones((MC,len(numWrkrs_)))
    #timeFedUAltMinRow = -1*np.ones((MC,len(numWrkrs_)))
    #timeAltMin_ = -1*np.ones((len(numWrkrs_),T))
    #SDAltMin_ = -1*np.ones((len(numWrkrs_),T))
    #---
    timeAltMinFedCol = -1*np.ones( (MC,len(numWrkrs_)) )
    #timeFedVAltMinRow = -1*np.ones((MC,len(numWrkrs_)))
    #timeFedUAltMinRow = -1*np.ones((MC,len(numWrkrs_)))
    timeAltMin_ = -1*np.ones((len(numWrkrs_),T))
    SDAltMin_ = -1*np.ones((len(numWrkrs_),T))

    #---
    timeAltGDMinFedSparse = -1*np.ones( (MC,len(numWrkrs_)) )
    timeCntrAltGDMinSparse = -1*np.ones( (MC,len(numWrkrs_)))
    timeFedAltGDMinSparse =  -1*np.ones ( (MC,len(numWrkrs_)))
    timeAltGDMin_ = -1*np.ones((len(numWrkrs_),T))
    SDAltGDMin_ = -1*np.ones((len(numWrkrs_),T))    
    #---
    timeAltMinGD = -1*np.ones( (MC,len(numWrkrs_)) )
    timeFedVAltMinGD = -1*np.ones((MC,len(numWrkrs_)))
    timeFedUAltMinGD = -1*np.ones((MC,len(numWrkrs_)))
    timeAltMinGD_ = -1*np.ones((len(numWrkrs_),T))
    SDAltMinGD_ = -1*np.ones((len(numWrkrs_),T))
    #--- 
    timeFactGD = -1*np.ones( (MC,len(numWrkrs_)) )
    timeFedFactGD = -1*np.ones((MC,len(numWrkrs_)))
    timeCntrlFactGD= -1*np.ones((MC,len(numWrkrs_)))
    timeFactGD_ = -1*np.ones((len(numWrkrs_),T))
    SDFactGD_ = -1*np.ones((len(numWrkrs_),T))
    #--------------------------------------------
    Ustr = np.linalg.qr(np.random.randn(n, r))[0]
    for i,numWrkrs in enumerate(numWrkrs_):
        #cluster = coiled.Cluster(
        #n_workers=numWrkrs,
        #region="us-east-2",
        #backend_options = {"zone_name": "us-east-2b"},
        #spot_policy="on-demand",
        #scheduler_vm_types="c7gn.4xlarge",
        #worker_vm_types="c7gn.2xlarge")
        #client = cluster.get_client()
        #----------------------------------------------
        cluster = LocalCluster(name='lll',n_workers=numWrkrs)
        client = Client(cluster)
        #----------------------------------------------
        timeAltMinMC = -1*np.ones((MC,T))
        SDAltMinMC = -1*np.ones((MC,T))    
        timeAltGDMinMC = -1*np.ones((MC,T))
        SDAltGDMinMC = -1*np.ones((MC,T))    
        timeAltMinGDMC = -1*np.ones((MC,T))
        SDAltMinGDMC = -1*np.ones((MC,T))
        timeFactGDMC = -1*np.ones((MC,T))
        SDFactGDMC = -1*np.ones((MC,T))
        inputList = [0]*numWrkrs
        colsPerWrkr = [q//numWrkrs]*numWrkrs
        colsPerWrkr[0] = int(colsPerWrkr[0])  + int(q - (q//numWrkrs)*numWrkrs)
        inputList = [[n,colsPerWrkr[j],r,p,0,Ustr] for j in range(numWrkrs)]
        offset = 0
        for j in range(numWrkrs):
            inputList[j][-2] = offset
            offset = offset + colsPerWrkr[j]
        for mc in range(MC):
            a = time.time()
            dataPtrs = client.map(getDataPtr,inputList)
            if mc < MC:
                U0_init,S,Vptr_ = fedSVD(dataPtrs,n,q,r,numWrkrs,client,T = 15)
                S = S/p
            b = time.time()
            print(f"Time taken indexing and Init SVD = {b - a:.3f} seconds. MC = {mc}.")
            #  run code once to make sure that the code is distributed to the workers 
            #  distributing the code to ther workers incurs extra time, which makes the timing reported by th first MC run inaccruate 
            if mc == 0: 
                Tlcl = 2
                limLcl = 1
                [a,b,c,d] = altMinFedCol(r, q, Ustr, Tlcl, dataPtrs, numWrkrs,client, limLcl, U0_init)
                [a,b,c,d,e,f,g,h] = altGDMinFedDaskScttrSparseNew(r, eta_c, Ustr,  2*Tlcl, p, dataPtrs,numWrkrs, client, limLcl, U0_init,S)
                [a,b,c,d,e,f,g,h] = factGD(q,r,  Ustr, Tlcl, p, dataPtrs,numWrkrs, client, limLcl, U0_init,np.diag(S),Vptr_)
                [a,b,c,d,e,f,g,h] = altMinGD(r, eta_c, Ustr,  2*Tlcl, p,dataPtrs, numWrkrs, client, limLcl, U0_init, S, T_in)

            #[SD_,tme_, tmeC_, tmeF_,tme,tmeFU,tmeFV,idx] = altMinFullyDist(r,q, Ustr, T,dataPtrs
            #                                                                      , numWrkrs,client, lim, U0_init)

            #timeAltMinFedRow[mc,i] = tme
            #timeFedUAltMinRow[mc,i] = tmeFU
            #timeFedVAltMinRow[mc,i] = tmeFV
            #timeAltMinMC[mc,:idx] = tme_
            #SDAltMinMC[mc,:idx] = SD_
            #print(f"altMinFed-FullyDist. time = {tme:.2f}s, timeFed-U = {tmeFU:.2f}s, timeFed-V = {tmeFV:.2f}s,  SD = {SD_[-1]:.2E},   Workers  =  {numWrkrs}, n = {n}, q = {q}, r = {r}, p = {p}", flush = True)

            #client.restart()
            [SD_,tme_,tme,idx] = altMinFedCol(r, q, Ustr, T, dataPtrs, numWrkrs,client, lim, U0_init)
            #print(SD_)
            timeAltMinFedCol[mc,i] = tme 
            timeAltMinMC[mc,:idx] = tme_
            SDAltMinMC[mc,:idx] = SD_
            print(f"altMinFedCol. time = {tme:.2f}s,   SD = {SD_[-1]:.2E},   Workers  =  {numWrkrs}, n = {n}, q = {q}, r = {r}, p = {p}", flush = True)
            print(f"------")
            #
            [SD_,tme_,tmeC_,tmeF_,tme,tmeC,tmeF,idx] = altGDMinFedDaskScttrSparseNew(r, eta_c, Ustr,  2*T, p,
                                                                                 dataPtrs,numWrkrs, client, lim, U0_init,S)
            timeAltGDMinFedSparse[mc,i] = tme
            timeCntrAltGDMinSparse[mc,i] = tmeC
            timeFedAltGDMinSparse[mc,i] = tmeF
            timeAltGDMinMC[mc,:idx] = tme_
            SDAltGDMinMC[mc,:idx] = SD_
            print(f"altGDMinSparse. time = {tme:.2f}s,  timeCenter = {tmeC:.2f}s,  timeFed = {tmeF:.2f}s,  SD = {SD_[-1]:.2E}, Workers  =  {numWrkrs}, n = {n}, q = {q}, r = {r}, p = {p}", flush = True)
            print(f"-----")
            #
            [SD_,tme_,tmeC_,tmeF_,tme,tmeC,tmeF,idx] = factGD(q,r,  Ustr, T, p, dataPtrs,
                                                            numWrkrs, client, lim, U0_init,np.diag(S),Vptr_)
            timeFactGD[mc,i] = tme
            timeCntrlFactGD[mc,i] = tmeC
            timeFedFactGD[mc,i] = tmeF
            timeFactGDMC[mc,:idx] = tme_
            SDFactGDMC[mc,:idx] = SD_
            print(f"factGD. time = {tme:.2f}s,  timeCenter = {tmeC:.2f}s, timeFed = {tmeF:.2f}s,  SD = {SD_[-1]:.2E}, Workers  =  {numWrkrs}, n = {n}, q = {q}, r = {r}, p = {p}", flush = True)
            print(f"-----")
            #
            #[SD_,tme_,tmeFU_,tmeFV_,tme,tmeFU,tmeFV,idx] = altMinGD(r, eta_c, Ustr,  2*T, p,
            #                                                         dataPtrs, numWrkrs, client, lim, U0_init, S, T_in)
            #timeAltMinGD[mc,i] = tme
            #timeFedUAltMinGD[mc,i] = tmeFU
            #timeFedVAltMinGD[mc,i] = tmeFV
            #timeAltMinGDMC[mc,:idx] = tme_
            #SDAltMinGDMC[mc,:idx] = SD_
            #print(f"altMinGD. time = {tme:.2f}s,  timeFedU = {tmeFU:.2f}s, timeFedV = {tmeFV:.2f}s, SD = {SD_[-1]:.2E}, Workers  =  {numWrkrs},n = {n}, q = {q}, r = {r}, p = {p}", flush = True)
            #print(f"-----")
            #
            #
            client.restart()
            #Save variables
            #varDictSmall = {"n":n,"q":q,"r":r,"p":p,"MC":mc,"ID":ID,
            #"timeFactGD":timeFactGD[:mc,i],
            #"timeAltGDMinFedSparse":timeAltGDMinFedSparse[:mc,i],
            #"timeAltMinFedRow":timeAltMinFedCol[:mc,i],
            #"timeAltMinGD":timeAltMinGD[:mc,i],
            #"numWrkrs":numWrkrs}
            #scipy.io.savemat(f"data/n_{n}_q_{q}_r_{r}_p_{p}_numWrks_{numWrkrs}_ID_{ID}.mat", varDictSmall)
            # Store data (serialize)
            #with open(f"data/n_{n}_q_{q}_r_{r}_p_{p}_numWrks_{numWrkrs}_ID_{ID}.pickle", 'wb') as handle:
                #pickle.dump(varDictSmall, handle, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"--------------------------------------------------------")
        cluster.close()
        client.close()
        client.shutdown()
        # Averaging
        timeAltMin_[i,:] =  np.sum(timeAltMinMC,axis=0)/MC
        SDAltMin_ [i,:]=  np.sum(SDAltMinMC,axis=0)/MC
        #---
        timeAltGDMin_[i,:] = np.sum(timeAltGDMinMC,axis=0)/MC
        SDAltGDMin_[i,:] = np.sum(SDAltGDMinMC,axis=0)/MC   
        #----
        timeAltMinGD_[i,:] = np.sum(timeAltMinGDMC,axis=0)/MC
        SDAltMinGD_[i,:]= np.sum(SDAltMinGDMC,axis=0)/MC
        #---
        timeFactGD_[i,:] = np.sum(timeFactGDMC,axis=0)/MC
        SDFactGD_[i,:]= np.sum(SDFactGDMC,axis=0)/MC
        # Plotting
        LAltMin = np.where(SDAltMin_[i,:] <= 1e-10)[0][0] + 1
        for l in range(LAltMin):
            timeAltMin_[i,l+1] = timeAltMin_[i,l+1] + timeAltMin_[i,l]
        LAltGDMin = np.where(SDAltGDMin_[i,:] <= 1e-10)[0][0] + 1
        for l in range(LAltGDMin):
            timeAltGDMin_[i,l+1] = timeAltGDMin_[i,l+1] + timeAltGDMin_[i,l]    
        LFactGD = np.where(SDFactGD_[i,:] <= 1e-10)[0][0] + 1
        for l in range(LFactGD):
            timeFactGD_[i,l+1] = timeFactGD_[i,l+1] + timeFactGD_[i,l]
        LAltMinGD = np.where(SDAltMinGD_[i,:] <= 1e-10)[0][0] + 1
        for l in range(LAltMinGD):
            timeAltMinGD_[i,l+1] = timeAltMinGD_[i,l+1] + timeAltMinGD_[i,l]
        plotErrAgnstTime(numWrkrs,n,q,r,p,ID,lim,MC,
                        timeAltMin_[i,:LAltMin+1] ,SDAltMin_[i,:LAltMin+1],
                        timeAltGDMin_[i,:LAltGDMin+1], SDAltGDMin_[i,:LAltGDMin+1], 
                        timeFactGD_[i,:LFactGD+1], SDFactGD_[i,:LFactGD+1], 
                        timeAltMinGD_[i,:LAltMinGD+1], SDAltMinGD_[i,:LAltMinGD+1],T_in)
        #----------------------------------------------------------
        print(f"timeFactGD ={np.sum(timeFactGD,axis=0)/MC}", flush = True)
        #---
        print(f"timeAltGDMinFedSparse ={np.sum(timeAltGDMinFedSparse,axis=0)/MC}", flush = True)
        #---
        print(f"timeAltMinFedCol = {np.sum(timeAltMinFedCol, axis=0)/MC}", flush = True)
        #---
        print(f"timeAltMinGD ={np.sum(timeAltMinGD,axis=0)/MC}", flush = True)
        #---
    #---------------------------------------------------------------
    # Average (Mean) time to convergence
    timeAltMinFedColAvg = np.sum(timeAltMinFedCol, axis=0)/MC
    timeAltGDMinFedSparseAvg = np.sum(timeAltGDMinFedSparse, axis=0)/MC
    timeAltMinGDAvg = np.sum(timeAltMinGD, axis=0)/MC
    timeFactGDAvg = np.sum(timeFactGD, axis=0)/MC
    #-------------------------------------------------------------
    plotTimeAgnstWrkrs(np.array(numWrkrs_),  n,q,r,p,ID,lim,MC,
                    timeAltMinFed = timeAltMinFedColAvg, 
                    timeAltGDMinFedSparse = timeAltGDMinFedSparseAvg,
                    timeAltMinGD = timeAltMinGDAvg,
                    timeFactGD = timeFactGDAvg, T_in  = T_in, median = 0)
    # Median time to convergence
    timeAltMinFedColMed = np.median(timeAltMinFedCol, axis = 0)
    timeAltGDMinFedSparseMed = np.median(timeAltGDMinFedSparse,axis = 0)
    timeAltMinGDMed = np.median(timeAltMinGD, axis = 0)
    timeFactGDMed = np.median(timeFactGD, axis = 0)
    # -------------------------------------------------------
    plotTimeAgnstWrkrs(np.array(numWrkrs_),  n,q,r,p,ID,lim,MC,
                timeAltMinFed = timeAltMinFedColMed, 
                timeAltGDMinFedSparse = timeAltGDMinFedSparseMed,
                timeAltMinGD = timeAltMinGDMed,
                timeFactGD = timeFactGDMed, T_in  = T_in, median = 1)
    #------------------------------------------------------------
    varDict = {"n":n,"q":q,"r":r,"p":p,"MC":MC,"ID":ID,
                "timeAltMin_":timeAltMin_,"SDAltMin_":SDAltMin_,
                "timeAltGDMin_":timeAltGDMin_, "SDAltGDMin_":SDAltGDMin_,
                "timeAltMinGD_":timeAltMinGD_,"SDAltMinGD_":SDAltMinGD_,
                "timeFactGD_":timeFactGD_,"SDFactGD_":SDFactGD_,
                "timeFactGD":timeFactGD,
                "timeAltGDMinFedSparse":timeAltGDMinFedSparse,
                "timeAltMinFedRow":timeAltMinFedCol,
                "timeAltMinGD":timeAltMinGD,
                "numWrkrs_":numWrkrs_}
