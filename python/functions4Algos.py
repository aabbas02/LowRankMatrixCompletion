import numpy as np
import time
from dask.distributed import Client, LocalCluster
import random
import os
from scipy.sparse.linalg import svds
from dask.distributed import Client
from scipy.sparse import csr_matrix
from utils import plotErrAgnstTime, plotTimeAgnstWrkrs
import copy

def altGDMinFedDaskScttrSparseNew(r,eta_c, Ustr, T, p, wrkrDataPtrs, numWrkrs, client, lim, U,S):
    n = U.shape[0]
    eta = eta_c / (S[-1] ** 2 * p)
    # this estimates mu as max j such that ||u^j||_2 = mu sqrt(r/n)  - implies no truncation at all for U_init
    mu = np.max(np.linalg.norm(U, axis=1)) * np.sqrt(n / r) 
    U = U * np.minimum(np.ones((n, 1)),
                        mu*np.sqrt(r/n)/ np.sqrt(np.sum(U**2, axis=1)).reshape(-1, 1))
    U,_ = np.linalg.qr(U)
    SDVals = np.zeros(T + 1)
    SDVals[0] = np.linalg.norm(Ustr - U @ (U.T @ Ustr), 'fro')
    timeArr = np.zeros(T + 1)    
    timeArrC = np.zeros(T + 1)
    timeArrF = np.zeros(T + 1)
    gradient = np.zeros( (n, r) )
    for i in range(T):
        tStart = time.time()
        gradient = 0*gradient
        newList = [[U,wrkrDataPtrs[j]] for j in range(numWrkrs)]
        grdnts = client.map(wrkrJGradScttrSparseNew, newList)
        grdnts = client.gather(grdnts)
        gradient = np.sum(grdnts,axis=0)
        tStartC = time.time()
        U = U - eta * gradient
        U, _ = np.linalg.qr(U, mode='reduced')
        tEnd = time.time()
        tFed = tStartC - tStart
        tCntr = tEnd - tStartC
        timeArr[i + 1] = tEnd - tStart
        timeArrC[i + 1] = tCntr
        timeArrF[i + 1] = tFed
        SD = np.linalg.norm(Ustr - U @ (U.T@ Ustr), 'fro')
        #print(f"AltGDMin SD = {SD:.5E}")
        SDVals[i + 1] = SD
        if SD < lim:
            break  
    idx = np.where(SDVals < lim)[0][0]+1
    tme = np.sum(timeArr)
    tmeC = np.sum(timeArrC)
    tmeF = np.sum(timeArrF)
    return SDVals[:idx], timeArr[:idx], timeArrC[:idx], timeArrF[:idx],tme,tmeC,tmeF,idx
#-----------------------------------------------------
def wrkrJGradScttrSparseNew(inputList):
    U = inputList[0]
    rowIdx_jk = inputList[1][2]
    XcolJ = inputList[1][1]
    rowsFlat = inputList[1][4]
    colsFlat = inputList[1][5]
    Xlin = inputList[1][6]
    n,r = U.shape
    temp = np.asarray( [ np.linalg.lstsq(U[rowIdx_jk[k], :], XcolJ[k], rcond=None)[0] for k in range(len(rowIdx_jk))] ) #r x q/L
    diff = np.sum(U[rowsFlat,:]*temp[colsFlat,:],axis=1) - Xlin
    M = csr_matrix((diff, (rowsFlat,colsFlat)),[n,len(XcolJ)])
    return M@temp

def factGD(q, r, Ustr, T, p, wrkrDataPtrs, numWrkrs, client, lim, U,S,Vptr_):
    n = Ustr.shape[0]
    SDVals = np.zeros(T + 1)
    SDVals[0] = np.linalg.norm(Ustr - U @ (U.T @ Ustr), 'fro')
    U = U[:, :r] @ np.sqrt(S[:r, :r]) #replace by just S if S is r x r matrix
    # Step size selection
    step_const = 0.75
    steplength = step_const / S[-1, -1]
    # Incoherence parameter mu  
    norm_U = np.max(np.linalg.norm(U, axis=1)) * np.sqrt(n / r) 
    #norm_B = np.max(np.linalg.norm(V, axis=1)) * np.sqrt(q / r)  # 2 
    #mu = max(norm_U, norm_B)
    mu = norm_U 
    # Initialization of Uzero and Vzero after Projection 
    const1 = np.sqrt(4 * mu * r / n) * S[-1, -1]
    const2 = np.sqrt(4 * mu * r / q) * S[-1, -1]
    U = U * np.minimum(np.ones((n, 1)), const1 / np.sqrt(np.sum(U**2, axis=1)).reshape(-1, 1))
    timeArr = np.zeros(T + 1)
    timeArrC = np.zeros(T + 1)
    timeArrF = np.zeros(T + 1)
    # Pacakge worker data
    colsPerWrkr = [q//numWrkrs]*numWrkrs
    colsPerWrkr[0] = int(colsPerWrkr[0])  + int(q - (q//numWrkrs)*numWrkrs)
    List = [[p,steplength,const2] for j in range(numWrkrs)]
    listFuture = client.scatter(List,broadcast = True)
    #-----------
    for i in range(T):
        tStart = time.time()
        newList = [[U,Vptr_[j],wrkrDataPtrs[j]] for j in range(numWrkrs)]
        # First Exchange
        diffM  = client.map(getDiff,newList)
        newList = [[U,diffM[j],Vptr_[j],wrkrDataPtrs[j]] for j in range(numWrkrs)]
        outputs_ = client.map(factGD1,newList)
        outputs_ = client.gather(outputs_)
        T1 = np.zeros((n,r))
        T3 = np.zeros((r,r))
        for j in range(numWrkrs):
            T1 = T1 + outputs_[j][0]
            T3 = T3 + outputs_[j][1]
        # Center 
        tStartCenter = time.time()
        T4 = U.T@U
        diffNrmBlnc = T4 - T3
        Unew = U - steplength*T1/p - steplength/16*U@(diffNrmBlnc)  
        Unew = Unew * np.minimum(np.ones((n, 1)), const1 / np.sqrt(np.sum(Unew**2, axis=1)).reshape(-1, 1))
        tEndCenter = time.time() 
        # Second Exchange - update VJ locally
        newList = [[U,diffM[j],diffNrmBlnc,Vptr_[j],listFuture[j]] for j in range(numWrkrs)]
        Vptr_ = client.map(factGD2,newList)
        tEnd = time.time()
        timeArr[i + 1]  = tEnd - tStart
        timeArrC[i + 1] = tEndCenter - tStartCenter
        timeArrF[i + 1] = tEnd - tStart - (tEndCenter - tStartCenter)
        Uproj, _ = np.linalg.qr(Unew, mode='reduced')
        SD = np.linalg.norm(Ustr - Uproj @ (Uproj.T@ Ustr), 'fro')
        #print(f" Fact GD  SD = {SD:.3E}")
        U = Unew
        SDVals[i + 1] = SD
        if SD < lim:
            break
    idx = np.where(SDVals < lim)[0][0]+1
    tme = np.sum(timeArr)
    tmeC = np.sum(timeArrC)
    tmeF = np.sum(timeArrF)
    return SDVals[:idx], timeArr[:idx], timeArrC[:idx],timeArrF[:idx],tme,tmeC,tmeF, idx

def getDiff(inputList):
    U = inputList[0]
    VJ = inputList[1] # VJ is q/numWrkrs x r
    XcolJ = inputList[2][1]
    rowsFlat = inputList[2][4]
    colsFlat = inputList[2][5]
    Xlin = inputList[2][6]
    n = U.shape[0]
    diff = np.sum(U[rowsFlat,:]*VJ[colsFlat,:],axis=1) - Xlin
    M = csr_matrix((diff, (rowsFlat,colsFlat)),[n,len(XcolJ)])
    return M

def factGD1(inputList):
    diffM = inputList[1]
    VJ = inputList[2]
    gradU =  diffM@VJ
    nrmBlnc = VJ.T@VJ
    return gradU, nrmBlnc

def factGD2(inputList):
    U = inputList[0]
    diffM = inputList[1]
    diffNrmBlnc = -inputList[2] # the minus makes it V'V - U'U
    VJ = inputList[3]
    p = inputList[4][-3]
    stepLength = inputList[4][-2]
    const2 = inputList[4][-1]
    T2 = diffM.T@U
    VJ = VJ - stepLength*T2/p - (stepLength/16)*VJ@diffNrmBlnc
    VJ = VJ * np.minimum(np.ones((VJ.shape[0], 1)), const2 / np.sqrt(np.sum(VJ**2, axis=1)).reshape(-1, 1))
    return VJ    
