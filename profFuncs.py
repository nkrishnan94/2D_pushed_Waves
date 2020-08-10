import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
from scipy import fft, ifft


def graphTot(profArr):
    fig,ax=plt.subplots()
    
    profArr_=np.copy(profArr)
    my_cmap = cm.get_cmap('viridis')
    my_cmap.set_bad('grey')
    normalize = cm.colors.LogNorm(vmin=0, vmax=np.max(np.sum(profArr_,axis=2)))
    #smap = cm.ScalarMappable(norm=normalize, cmap=my_cmap,no)
    #smap = cm.get_cmap('Purples',norm=normalize)
    


    im = ax.imshow(np.sum(profArr_,axis=2).T,cmap=my_cmap
                   ,norm=normalize
                   ,vmin=1)
    ax.set_xticks([])
    ax.set_yticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    cbar_ticks= [1,10,100,1000,10000]
    cbar = fig.colorbar(im,ticks=cbar_ticks,orientation='horizontal',cax=cax)
    cbar.ax.set_xticklabels(cbar_ticks,fontsize=12) 
    cbar.ax.set_xlabel('Population Total',fontsize=15)
    plt.show()
   
def graphFreq(freqArr):
    fig,ax=plt.subplots()
    
    freqArr_=np.copy(freqArr)
    my_cmap = cm.get_cmap('viridis')
    my_cmap.set_bad('grey')
    #normalize = cm.colors.LogNorm(vmin=0, vmax=np.max(np.sum(prof_arr,axis=2)))
    #smap = cm.ScalarMappable(norm=normalize, cmap=my_cmap,no)
    #smap = cm.get_cmap('Purples',norm=normalize)
    


    im = ax.imshow(freqArr_.T,cmap=my_cmap
                   #,norm=normalize
                   ,vmin=0)
    ax.set_xticks([])
    ax.set_yticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    #cbar_ticks= [1,10,100,1000,10000]
    cbar = fig.colorbar(im
                        #t,icks=cbar_ticks
                        ,orientation='horizontal',cax=cax)
    #cbar.ax.set_xticklabels(cbar_ticks,fontsize=12) 
    cbar.ax.tick_params(labelsize=12) 
    cbar.ax.set_xlabel('$f_a$',fontsize=15)
    plt.show()
   

def profFileToArray(file,ndemeh,ndemew):
    arr = np.loadtxt(file,delimiter =  ', ')
    arr_=np.copy(arr)
    if len(file)!=0:
        
        return arr_[:,2:].reshape((ndemeh,ndemew,2))
    
    else:
        return False
    
def arrToFreq(arr):
    arr_=np.copy(arr)
    return arr_[:,:,0]/(arr_[:,:,1]+arr_[:,:,0])

def getLastRow(arrFreq):
    arrFreq_ = np.copy(arrFreq)
    max_all = np.sum(np.all(~np.isnan(arrFreq_),axis=1))-1
    f_x = arrFreq_[max_all,:]
    return f_x

def getFrontLength(arrFreq):
    arrFreq_ = np.copy(arrFreq)
    width = arrFreq_.shape[1]
    max_deme = np.zeros(width)
    front =np.argwhere(~np.isnan(arrFreq_))
    for i in range(width):
        col = front[:,1]==i
        max_deme[i] = np.max(front[col,:][:,0])
    return(width+np.sum(np.abs(np.diff(max_deme))))
    
def getFront(arrFreq):
    arrFreq_ = np.copy(arrFreq)
    width = arrFreq_.shape[1]
    max_deme = np.zeros(width)
    front = np.argwhere(~np.isnan(arrFreq_))
    for i in range(width):
        col = front[:,1]==i
        max_col = np.max(front[col,:][:,0])
        max_deme[i] = arrFreq_[max_col,i]
    return max_deme

def getSectorLengths(arrFreq):
    arrFreq_ = np.copy(arrFreq)
    width = arrFreq_.shape[1]
    max_deme = np.zeros(width)
    max_col = np.zeros(width)
    sec_size = np.zeros(width)
    front = np.argwhere(~np.isnan(arrFreq_))
    for i in range(width):
        col = front[:,1]==i
        max_col[i] = np.max(front[col,:][:,0])
        max_deme[i] = arrFreq_[int(max_col[i]),i]
        sec_size[i] = np.argmax(((arrFreq_[:int(max_col[i]),i] == max_deme[i])*1==0)[::-1])+1
    return max(sec_size)
    
def sectorCountMed(arr_row):
    arr_row_ = np.copy(arr_row)
    diffs = np.zeros(50)
    diff_ = np.abs(np.diff(arr_row_))
    for i,l in enumerate(np.linspace(0,1,50)):

        diffs[i] = np.sum(diff_>l)
    return np.median(diffs)

def sectorCountMean(arr_row):
    arr_row_ = np.copy(arr_row)
    diffs = np.zeros(50)
    diff_ = np.abs(np.diff(arr_row_))
    for i,l in enumerate(np.linspace(0,1,50)):

        diffs[i] = np.sum(diff_>l)
    return np.mean(diffs)


def signalFilter(arrNum,thresh):
    arrNum= np.copy(arr_row)
    
    return 
def rowTrim(arr_row,split_thresh):
    arr_row_=np.copy(arr_row)
    diff1 = np.diff(arr_row_)
    if len(np.argwhere(diff1!=0).flatten())!=0:
        size= np.concatenate((1+np.argwhere(diff1!=0)[0],
                    np.diff(np.argwhere(diff1!=0).flatten()),
                    len(diff1)-np.argwhere(diff1!=0)[-1]) )


        inds = np.cumsum(size)-1
        split_inds = np.argwhere(size<split_thresh).flatten()
        cnt = 0 
        size_d = 0
        for i,ind in enumerate(inds[split_inds]):


            ind+=size_d

            curr_size = size[split_inds][cnt]

            if ind-curr_size>0:
                arr_row_ = np.concatenate((arr_row_[:ind-curr_size],arr_row_[ind+1:]))
            else: 
                arr_row_= arr_row_[ind+1:]

            size_d-=curr_size
            cnt+=1
    return arr_row_
    

def RowBinary(f_x, thresh):
    f_x_=np.copy(f_x)
    f_x_[f_x>thresh]=1
    f_x_[f_x_<(1-thresh)]=0
    f_x_[((f_x_ != 0) & (f_x_ != 1))]=.5
    return f_x_

#def anyRowBinary(arrFreq,row_ind, thresh):

#    f_x = arrFreq[rowind,:]
#    f_x[f_x>thresh]=1
#    f_x[f_x<(1-thresh)]=0
#    f_x[((f_x != 0) & (f_x != 1))]=.5
#    return f_x


#def sectorCount(arr_row):
#    arr_row_=np.copy(arr_row)
#    return np.sum(np.diff(arr_row_)!=0)+1

def sectorCount(arr_row):
    arr_row_=np.copy(arr_row)
    inds = np.argwhere(np.diff(arr_row_)!=0).flatten()

    if len(inds)>0:
        counts = np.bincount((arr_row_[inds]*2).astype(int),minlength=3)
    else:
        counts = np.array([1,0,0])
    return counts


def sectorProp(arr_row):
    arr_row_=np.copy(arr_row)
    counts = np.bincount((arr_row_*2).astype(int),minlength=2)
    return counts/np.sum(counts)
    
def demePropMix(arr_row):
    arr_row_=np.copy(arr_row)
    return np.sum((arr_row_ != 0) & (arr_row!=1))/len(arr_row_)

    
def corrFunc(arr_row):
    arr_row_=np.copy(arr_row)
    width =len(arr_row_)
    c_d =np.zeros(int(width/2))
    for i in range(int(width/2)):
        prod = 0
        for j in range(width):
            #print(f_x[i]*f_x[(i+1)%40])
            prod +=(arr_row_[j]*arr_row_[(j+i)%width])
        c_d[i] = prod/width
        
    return c_d


def corrSpikeCount(arr_row):
    arr_row_=np.copy(arr_row)
    return np.sum(np.diff(np.diff(fft(corrFunc(arr_row_))))<0) +1

def corrSpikes(arr_row):
    arr_row_=np.copy(arr_row)
    co = np.diff(fft(corrFunc(arr_row_)))
    sums=0
    for i in range(len(co)-1):
        if co[i]*co[i+1]<0 and (co[i]>0):
            sums +=1
        
    return sums  
    