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

def altMinGD(r, eta_c, Ustr, T, p, wrkrDataPtrs, numWrkrs, client, lim, U,S, T_in):
    n = U.shape[0]
    # this estimates mu as max j such that ||u^j||_2 = mu sqrt(r/n)  - implies no truncation at all for U_init
    mu = np.max(np.linalg.norm(U, axis=1)) * np.sqrt(n / r) 
    U = U * np.minimum(np.ones((n, 1)),
                        mu*np.sqrt(r/n)/ np.sqrt(np.sum(U**2, axis=1)).reshape(-1, 1))
    eta = eta_c / (S[-1] ** 2 * p)
    SDVals = np.zeros(T + 1)
    SDVals[0] = np.linalg.norm(Ustr - U @ (U.T @ Ustr), 'fro')
    timeArr = np.zeros(T + 1)    
    timeArrFU = np.zeros(T + 1)
    timeArrFV = np.zeros(T + 1)
    gradient = np.zeros( (n, r))
    for i in range(T):
        tStart = time.time()
        newList = [[U,wrkrDataPtrs[j]] for j in range(numWrkrs)]
        temp = client.map(wrkrJLstSqrsScttrCol,newList)
        tStartU = time.time()
        for _ in range(T_in):
            newList = [[U,temp[j],wrkrDataPtrs[j]] for j in range(numWrkrs)]
            grdnts = client.map(wrkrJAltMinGD, newList)
            grdnts = client.gather(grdnts)
            gradient = np.sum(grdnts,axis=0)
            U = U - eta * gradient
        tEnd = time.time()
        #print(f"Iter = {i}. Gradient Norm = {np.linalg.norm(gradient,"fro"):.9E}")
        tFedU =  tEnd - tStartU
        timeArr[i + 1] =   tEnd - tStart
        timeArrFU[i + 1] = tFedU
        timeArrFV[i + 1] = tFedU
        Uproj, _ = np.linalg.qr(U, mode='reduced')
        SD = np.linalg.norm(Ustr - Uproj @ (Uproj.T@ Ustr), 'fro')
        SDVals[i + 1] = SD
        #print(f"AltMinGD. SD = {SD}")
        if SD < lim:
            break  
    
    idx = np.where(SDVals < lim)[0][0]+1
    tme = np.sum(timeArr)
    tmeFU = np.sum(timeArrFU)
    tmeFV = np.sum(timeArrFV)        
    return SDVals[:idx], timeArr[:idx], timeArrFU[:idx], timeArrFV[:idx],tme,tmeFU,tmeFV,idx
#-----------------------------------------------------
def wrkrJLstSqrsScttrCol(inputList):
    U = inputList[0]
    rowIdx_jk = inputList[1][2]
    XcolJ = inputList[1][1]
    _,r = U.shape
    V = np.zeros((r,len(XcolJ)))
    V = np.asarray([np.linalg.lstsq(U[rowIdx_jk[j], :], XcolJ[j], rcond=None)[0] for j in range(len(XcolJ))]).T
    return V
#----
def wrkrJAltMinGD(inputList):
    U = inputList[0]
    temp = inputList[1]
    XcolJ = inputList[2][1]
    rowsFlat = inputList[2][4]
    colsFlat = inputList[2][5]
    Xlin = inputList[2][6]
    n,r = U.shape
    diff = np.sum(U[rowsFlat,:]*temp[:,colsFlat].T,axis=1) - Xlin
    M = csr_matrix((diff, (rowsFlat,colsFlat)),[n,len(XcolJ)])
    return M@temp.T

def altMinFullyDist(r, q,Ustr, T, wrkrDataPtrs, numWrkrs,client, lim, Uoriginal):
    U = copy.deepcopy(Uoriginal)
    n = U.shape[0]
    # this estimates mu as max j such that ||u^j||_2 = mu sqrt(r/n)  - implies no truncation at all for U_init
    mu = np.max(np.linalg.norm(U, axis=1)) * np.sqrt(n / r) 
    U = U * np.minimum(np.ones((n, 1)),
                        mu*np.sqrt(r/n)/ np.sqrt(np.sum(U**2, axis=1)).reshape(-1, 1))
    #q = len(Xcol)
    SDVals = np.zeros(T+1)
    SDVals[0] = np.linalg.norm(Ustr - U @ (U.T @ Ustr), 'fro')
    timeArr = np.zeros(T+1)
    timeArrFV = np.zeros(T + 1)
    timeArrFU = np.zeros(T + 1)
    V = np.zeros((r, q))
    colsPerWrkr = [q//numWrkrs]*numWrkrs
    colsPerWrkr[0] = int(colsPerWrkr[0])  + int(q - (q//numWrkrs)*numWrkrs)
    #--------------
    a = time.time()
    rowDataPtrs = client.map(getRowData,wrkrDataPtrs)
    rowData = client.gather(rowDataPtrs)
    b = time.time()
    # The following for loop stiches row indices and values from all workers; 
    # each worker sends row data for n rows in the corresponding qL columns
    for j in range(numWrkrs): 
        if j > 0:
            colIdx = [np.append(colIdx[j_],rowData[j][0][j_] ) for j_ in range(n)]
            Xrow = [np.append(Xrow[j_], rowData[j][1][j_]) for j_ in range(n)]
        else:
            colIdx = rowData[j][0] # list of length n - each item is a list of indices - indices are global, i.e., range from 0 to n
            Xrow = rowData[j][1] # list of length n 
    rowsPerWrkr = [n//numWrkrs]*numWrkrs
    rowsPerWrkr[0] = int(rowsPerWrkr[0])  + int(n - (n//numWrkrs)*numWrkrs)
    strtR = 0
    List = [0]*numWrkrs
    for j in range(numWrkrs): 
        colIdx_jk = colIdx[strtR:strtR + rowsPerWrkr[j]]
        #print(colIdx_jk)
        XrowJ = Xrow[strtR:strtR + rowsPerWrkr[j]]
        List[j] = [colIdx_jk, XrowJ]
        strtR = strtR + rowsPerWrkr[j]
    c = time.time()
    listFutureRow = client.scatter(List, broadcast = False)
    d = time.time()
    offset = b - a + d - c
    print(f"Time to gather and scatter AltMin Fed. Fully (Dist.) = {offset:.3f}s. Gather Time = {b - a:.3f}. Scatter Time = {d - c:.3f}s")
    timeArr[0] = offset
    for i in range(T):
        tStart = time.time()
        newList = [[U,wrkrDataPtrs[j]] for j in range(numWrkrs)] 
        result = client.map(wrkrJLstSqrsScttrCol,newList)
        result = client.gather(result)
        strt = 0
        for j in range(numWrkrs):
            V[:,strt : strt + colsPerWrkr[j]] = result[j]
            strt = strt + colsPerWrkr[j]
        tStart2 = time.time()
        #newList = [[V,listFuture[j]] for j in range(numWrkrs)]
        newList = [[V, listFutureRow[j]] for j in range(numWrkrs)]
        result = client.map(wrkrJLstSqrsScttrRow,newList)
        result = client.gather(result)
        strtR = 0
        for j in range(numWrkrs):
            U[strtR:strtR + rowsPerWrkr[j], :] = result[j]
            strtR = strtR + rowsPerWrkr[j]
        tEnd = time.time()
        tFedU = tStart2 - tStart
        tFedV = tEnd - tStart2
        timeArr[i + 1] =  tEnd - tStart
        timeArrFU[i + 1] =  tFedU
        timeArrFV[i + 1] = tFedV
        Uproj, _ = np.linalg.qr(U, mode='reduced')
        SD = np.linalg.norm(Ustr - Uproj @ (Uproj.T@ Ustr), 'fro')
        SDVals[i + 1] = SD
        if SD < lim:
            break
    idx = np.where(SDVals < lim)[0][0]+1
    tme = np.sum(timeArr)
    tmeFU = np.sum(timeArrFU)
    tmeFV = np.sum(timeArrFV)        
    return SDVals[:idx], timeArr[:idx], timeArrFU[:idx],  timeArrFV[:idx],tme,tmeFU,tmeFV,idx

def getRowData(inputList):
    return inputList[3], inputList[0]

def wrkrJLstSqrsScttrRow(inputList):
    V = inputList[0]
    colIdx_jk = inputList[1][0]
    XrowJ = inputList[1][1]
    r,_ = V.shape
    U = np.zeros((r,len(XrowJ)))
    U = np.asarray( [np.linalg.lstsq(V[:, colIdx_jk[j]].T, XrowJ[j], rcond=None)[0].T.flatten() for j in range(len(XrowJ))] )
    return U
#-------------------------------------------------------------------------------
def altMinFedCol(r, q, Ustr, T, wrkrDataPtrs, numWrkrs,client, lim, Uoriginal):
    U = copy.deepcopy(Uoriginal)
    n = U.shape[0]
    # this estimates mu as max j such that ||u^j||_2 = mu sqrt(r/n)  - implies no truncation at all for U_init
    mu = np.max(np.linalg.norm(U, axis=1)) * np.sqrt(n / r) 
    U = U * np.minimum(np.ones((n, 1)),
                        mu*np.sqrt(r/n)/ np.sqrt(np.sum(U**2, axis=1)).reshape(-1, 1))
    #q = len(Xcol)
    SDVals = np.zeros(T+1)
    SDVals[0] = np.linalg.norm(Ustr - U @ (U.T @ Ustr), 'fro')
    timeArr = np.zeros(T+1)
    #timeArrFV = np.zeros(T + 1)
    #timeArrFU = np.zeros(T + 1)
    V = np.zeros((r, q))
    colsPerWrkr = [q//numWrkrs]*numWrkrs
    colsPerWrkr[0] = int(colsPerWrkr[0])  + int(q - (q//numWrkrs)*numWrkrs)
    #--------------
    a = time.time()
    rowDataPtrs = client.map(getRowData,wrkrDataPtrs)
    rowData = client.gather(rowDataPtrs)
    b = time.time()
    # The following for loop stiches row indices and values from all workers; 
    # each worker sends row data for n rows in the corresponding qL columns
    for j in range(numWrkrs): 
        if j > 0:
            colIdx = [np.append(colIdx[j_],rowData[j][0][j_] ) for j_ in range(n)]
            Xrow = [np.append(Xrow[j_], rowData[j][1][j_]) for j_ in range(n)]
        else:
            colIdx = rowData[j][0] # list of length n - each item is a list of indices - indices are global, i.e., range from 0 to n
            Xrow = rowData[j][1] # list of length n 
    rowsPerWrkr = [n//numWrkrs]*numWrkrs
    rowsPerWrkr[0] = int(rowsPerWrkr[0])  + int(n - (n//numWrkrs)*numWrkrs)
    strtR = 0
    List = [0]*numWrkrs
    for j in range(numWrkrs): 
        colIdx_jk = colIdx[strtR:strtR + rowsPerWrkr[j]]
        #print(colIdx_jk)
        XrowJ = Xrow[strtR:strtR + rowsPerWrkr[j]]
        List[j] = [colIdx_jk, XrowJ]
        strtR = strtR + rowsPerWrkr[j]
    c = time.time()
    #listFutureRow = client.scatter(List, broadcast = False)
    d = time.time()
    offset = b - a + d - c
    print(f"Time to gather and scatter AltMin Fed. Fully (Dist.) = {offset:.3f}s. Gather Time = {b - a:.3f}. Scatter Time = {d - c:.3f}s")
    #timeArr[0] = offset
    timeArr[0] = 0
    for i in range(T):
        tStart = time.time()
        newList = [[U,wrkrDataPtrs[j]] for j in range(numWrkrs)] 
        result = client.map(wrkrJLstSqrsScttrCol,newList)
        result = client.gather(result)
        strt = 0
        for j in range(numWrkrs):
            V[:,strt : strt + colsPerWrkr[j]] = result[j]
            strt = strt + colsPerWrkr[j]
        #tStart2 = time.time()
        U = np.asarray([np.linalg.lstsq(V[:, colIdx[j]].T, Xrow[j], rcond=None)[0].T.flatten() for j in range(n)])
        tEnd = time.time()
        timeArr[i + 1] =  tEnd - tStart
        Uproj, _ = np.linalg.qr(U, mode='reduced')
        SD = np.linalg.norm((np.eye(n) - Uproj @ Uproj.T) @ Ustr, 'fro')
        SDVals[i + 1] = SD
        if SD < lim:
            break
    idx = np.where(SDVals < lim)[0][0]+1
    tme = np.sum(timeArr)
    return SDVals[:idx],timeArr[:idx], tme,idx
#--------------------------------
#def wrkrJLstSqrsScttr(inputList):
#    U = inputList[0]
#    rowIdx_jk = inputList[1][0]
#    XcolJ = inputList[1][1]
#    n,r = U.shape
#    V = np.zeros((r,len(XcolJ)))
#    V = np.asarray([np.linalg.lstsq(U[rowIdx_jk[j], :], XcolJ[j], rcond=None)[0] for j in range(len(XcolJ))]).T
#    return V