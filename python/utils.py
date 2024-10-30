import matplotlib.pyplot as plt
import pickle
import scipy.io
import numpy as np
from scipy.sparse import csr_matrix
import time

def saveVariables(varDict,n,q,r,p,ID,MC):
    scipy.io.savemat(f"n_{n}_q_{q}_r_{r}_p_{p}_MC_{MC}_ID_{ID}.mat", varDict)

def plotTimeAgnstWrkrs(numWrkrs_,  n,q,r,p,ID,lim,MC, timeFactGD = -1*np.ones(1000),
                timeAltGDMinFedSparse = -1*np.ones(1000), timeAltGDMinFedDaskScttr = -1*np.ones(1000),
                timeAltMinCntrl = -1*np.ones(1000), timeAltMinFed = -1*np.ones(1000), timeAltMinGD = -1*np.ones(1000), T_in = 10
               ):
    fig1, _ = plt.subplots(nrows=1,ncols=1)
    # Plotting the data
    L = numWrkrs_.shape[0]
    plt.plot(numWrkrs_[timeFactGD[:L] > 0], timeFactGD[timeFactGD > 0], label='FactGD', marker="o")
    plt.plot(numWrkrs_[timeAltGDMinFedSparse[:L] > 0], timeAltGDMinFedSparse[timeAltGDMinFedSparse > 0], label = 'AltGDMin', marker="^")
    plt.plot(numWrkrs_[timeAltMinFed[:L] > 0], timeAltMinFed[timeAltMinFed > 0], label='AltMin(FullyDist.)', marker="s")
    plt.plot(numWrkrs_[timeAltMinGD[:L] > 0], timeAltMinGD[timeAltMinGD > 0], label='AltMin(GD. Prvt.)', marker="*")
    plt.xlabel('Number of workers')
    plt.ylabel(r"Time taken for $SD_F(\hat \mathbf{U},\mathbf{U}^*) \leq$" f"{lim:.2E} / seconds")
    plt.title(f"n = {n}, q = {q}, r = {r}, p = {p}.")
    plt.legend()
    plt.grid()
    plt.xticks(numWrkrs_)
    #plt.show()
    # Save the figure
    fig1.savefig(f"figures/tmeAgnstWrkrs_n_{n}_q_{q}_r_{r}_p_{p}_T_in_{T_in}_MC_{MC}_ID_{ID}_.pdf", bbox_inches='tight')


def plotErrAgnstTime(numWrkrs,n,q,r,p,ID,lim,MC,
                    timeAltMinFullyDist,  SDAltMinFullyDist,
                    timeAltGDMin, SDAltGDMin, 
                    timeFactGD, SDFactGD, 
                    timeAltMinGDPrvt, SDAltMinGDPrvt,T_in
                    ):
    fig1, _ = plt.subplots()
    plt.semilogy(timeFactGD, SDFactGD, label='FactGD', marker="o")
    plt.semilogy(timeAltGDMin, SDAltGDMin, label='AltGDMin', marker="^")
    plt.semilogy(timeAltMinFullyDist, SDAltMinFullyDist, label='AltMin(FullyDist.)', marker="s")
    plt.semilogy(timeAltMinGDPrvt, SDAltMinGDPrvt, label='AltMin(GD.Prvt.)', marker="*")
    # Setting labels and title
    plt.xlabel(f"Time/s")
    plt.ylabel(r"$SD_F(\hat \mathbf{U},\mathbf{U}^*)$")
    plt.title(f"n = {n}, q = {q}, r = {r}, p = {p}, Workers = {numWrkrs}.")
    plt.legend()
    plt.grid()
    #plt.show()
    # Save the figure
    fig1.savefig(f"figures/ErrAgnstTime_n_{n}_q_{q}_r_{r}_p_{p}__Wrkrs_{numWrkrs}_MC_{MC}_T_in_{T_in}_ID_{ID}_.pdf", bbox_inches='tight')


def plotandSaveErrAgnstIter(timeSDAltGDMinFed,  SDAltGDMinFed,
                            timeSDAltGDMinFedDask, SDAltGDMinFedDask, 
                            timeSDAltGDMinFedDaskScttr, SDAltGDMinFedDaskScttr, 
                            timeSDAltMinCntrl, SDAltMinCntrl,
                            timeSDAltMinFed, SDAltMinFed,
                            timeSDAltMinFedDaskScttr, SDAltMinFedDaskScttr,n,q,r,p,ID):
    fig1, _ = plt.subplots()
    # Plotting the data
    plt.semilogy(SDAltGDMinFed[SDAltGDMinFed>0], label='AltGDMin (Fed.)', marker="o")
    plt.semilogy(SDAltGDMinFedDask[SDAltGDMinFedDask>0], label='AltGDMin (Dask.)', marker="x")
    plt.semilogy(SDAltGDMinFedDaskScttr[SDAltGDMinFedDaskScttr > 0], label='AltGDMin (Dask-Scttr)', marker="v")
    #plt.plot(timeSDAltMinCntrl, SDAltMinCntrl, label='AltMin (Cntrl.)', marker="v")
    plt.semilogy(SDAltMinFed[SDAltMinFed > 0], label='AltMin (Fed.)', marker="^")
    plt.semilogy(SDAltMinFedDaskScttr [SDAltMinFedDaskScttr > 0], label='AltMin (Dask-Scttr)', marker="h")
    # Setting labels and title
    plt.xlabel(f"Iter")
    plt.ylabel(r"$SD_F(\hat \mathbf{U},\mathbf{U}^*)$")
    plt.title(f"n = {n}, q = {q}, r = {r}, p = {p}.")
    # Adding a legend
    plt.legend()
    # Display the plot
    plt.grid()
    #plt.xticks(numWrkrs_)
    #plt.show()
    # Save the figure
    fig1.savefig(f"figures/ErrAgnstIter_n_{n}_q_{q}_r_{r}_p_{p}_ID_{ID}_.pdf", bbox_inches='tight')

def plotandSaveAllTime(numWrkrs_, timeAltGDMinFed,timeCntrAltGDMinFed, timeFedAltGDMinFed,
                timeAltGDMinFedDask,timeAltGDMinFedDaskScttr,
                timeAltMinCntrl, timeAltMinFed, timeAltMinFedDaskScttr,
                n,q,r,p,ID,MC):
    fig1, _ = plt.subplots()
    # Plotting the data
    plt.plot(numWrkrs_[timeAltGDMinFed > 0], timeAltGDMinFed[timeAltGDMinFed > 0], label='AltGDMin (Fed.)', marker="o")
    plt.plot(numWrkrs_[timeCntrAltGDMinFed > 0], timeCntrAltGDMinFed[timeCntrAltGDMinFed > 0], label='AltGDMin (Fed.)-Cntr', marker = "o", linestyle = '--')
    plt.plot(numWrkrs_[timeFedAltGDMinFed > 0], timeFedAltGDMinFed[timeFedAltGDMinFed > 0], label='AltGDMin (Fed.)-Nodes', marker="o", linestyle = "--")
    #
    plt.plot(numWrkrs_[timeAltGDMinFed > 0], timeAltGDMinFed[timeAltGDMinFed > 0], label='AltGDMin (Fed.)', marker="o")
    plt.plot(numWrkrs_[timeAltGDMinFedDask > 0], timeAltGDMinFedDask[timeAltGDMinFedDask > 0], label='AltGDMin (Dask.)', marker="x")
    plt.plot(numWrkrs_[timeAltGDMinFedDaskScttr > 0], timeAltGDMinFedDaskScttr[timeAltGDMinFedDaskScttr > 0], label='AltGDMin (Dask-Scttr)', marker="o")
    plt.plot(numWrkrs_[timeAltMinCntrl > 0], timeAltMinCntrl[timeAltMinCntrl > 0], label='AltMin (Cntrl.)', marker="v")
    plt.plot(numWrkrs_[timeAltMinFed > 0], timeAltMinFed[timeAltMinFed > 0], label='AltMin (Fed.)', marker="^")
    plt.plot(numWrkrs_[timeAltMinFedDaskScttr > 0], timeAltMinFedDaskScttr[timeAltMinFedDaskScttr > 0], label='AltMin (Dask-Scttr)', marker="h")
    plt.xlabel('Number of workers')
    plt.ylabel(r"Time taken for $SD_F(\hat \mathbf{U},\mathbf{U}^*) \leq 10^{-10}$")
    plt.title(f"n = {n}, q = {q}, r = {r}, p = {p}.")
    plt.legend()
    plt.grid()
    plt.xticks(numWrkrs_)
    fig1.savefig(f"figures/All_Time_NEW_n_{n}_q_{q}_r_{r}_p_{p}_MC_{MC}_ID_{ID}_.pdf", bbox_inches='tight')

def getDataPtr(inputs_):
    n = inputs_[0]
    qL = inputs_[1]
    r = inputs_[2]
    p = inputs_[3]
    offset = inputs_[4]
    U = inputs_[5]
    Bl = np.random.randn(r,qL)
    Xl = U@Bl
    idx = np.random.choice(n * (qL), size=int(p * n * (qL) ), replace=False)
    row, col = np.unravel_index(idx, (n, qL))
    rowIdx = [row[col == j] for j in range(qL)]
    Xcol = [Xl[rowIdx[j], j] for j in range(qL)]
    colIdx = [col[row == i] for i in range(n)]
    Xrow = [Xl[i, col_idx].reshape(-1, 1) for i, col_idx in enumerate(colIdx, start=0)]
    colIdx = [col_idx + offset for col_idx in colIdx]
    rowsFlat = []
    colsFlat = []
    Xlin = []
    for k, rowIdx_k in enumerate(rowIdx):
        numRows = len(rowIdx_k)
        rowsFlat.extend(rowIdx_k)
        colsFlat.extend([k] * numRows)
        Xlin.extend(Xcol[k])
    # Xrow is a list of length n. Each item in the list is a numpy array of the observed row entries.
    # Xcol is a list of length qL. Each item in the list is a numpy array of the observed column entries.
    # rowIdx is a list of lenth  qL. The items are the observed row indices in qL columns. Indices are values ranging from 0 to n-1.
    # colIdx is a list of length n.  The items are the observed column indices in the n rows. Indices are values ranging from 0 + offset to qL - 1 + offset.
    # rowFlat is a vectorization of rowIdx. 
    # colsFlat are local indices in that they take values from 0 to qL - 1. That is, no offset. 
    return Xrow, Xcol, rowIdx, colIdx, rowsFlat, colsFlat, Xlin


def fedSVD(dataPtrs,n,q,r,numWrkrs,client,T = 15):
    U = np.random.randn(n,r)
    for i in range(T):
        inputList = [[U,dataPtrs[j]] for j in range(numWrkrs)]
        U_ = client.map(fedSVDgetU, inputList)
        U_ = client.gather(U_)
        U = np.sum(U_,axis=0)
        if i == T - 1:
             S = np.sqrt(np.linalg.norm(U, axis=0))
        U = np.linalg.qr(U,mode='reduced')[0]
    idx = np.argsort(S) # ascending sort, that is, idx[-1] is the largest element.
    U = U[:,idx]
    S = S[idx]
    mu = np.max(np.linalg.norm(U, axis=1)) * np.sqrt(n / r) 
    const2 = np.sqrt(4 * mu * r / q) * S[-1]
    inputList = [[U,const2,dataPtrs[j]] for j in range(numWrkrs)]
    Vptr_ = client.map(fedSVDgetV,inputList) #Vptr_ is of size numWrkrs x qL x r
    return U,S,Vptr_

def fedSVDgetU(inputList):
    U = inputList[0]
    Xlin = inputList[1][-1]
    rowsFlat = inputList[1][4]
    colsFlat = inputList[1][5]
    qL = len(inputList[1][1])
    n = U.shape[0]
    M = csr_matrix((Xlin, (rowsFlat,colsFlat)),[n,qL])
    tmp = M.T@U
    return M@tmp

def fedSVDgetV(inputList):
    U = inputList[0]
    const2 = inputList[1]
    Xlin = inputList[2][-1]
    rowsFlat = inputList[2][4]
    colsFlat = inputList[2][5]
    qL = len(inputList[2][1])
    n = U.shape[0]
    M = csr_matrix((Xlin, (rowsFlat,colsFlat)),[n,qL])
    VJ = (U.T@M).T
    VJ = VJ * np.minimum(np.ones((VJ.shape[0], 1)), const2 / np.sqrt(np.sum(VJ**2, axis=1)).reshape(-1, 1))
    return VJ