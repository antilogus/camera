#!/usr/bin/env python
######################################################################
## Filename:      file_list.py
## Version:       V3
## Author:        P.Antilogus
##
## HISTORY
##
## Version        when                who           Comment  
## -----------   -------------------  ------------  ----------------
##  V1                               P.Antilogus    Initial version 
##  V2            18 Dec 2017        P.Antilogus    Frozen version for some time
##  V3            13 Oct 2020        P.Antilogus    Implement the long standing correction
##  Done in V3 :
##  13-10-20    - fix bug in the exact fraction of the amplifier used (error in the usage of  xoff and yoff)
##  13-10-20    - the error(=sigma)  on the variance 'evar' should refer to the <mean> value , so we should divide by the number of entry 
######################################################################
import glob
import os
try:
    import pyfits
except:
    import astropy.io.fits as pyfits
from optparse import OptionParser
import numpy as np
import pickle
#
class File_List :
    # Handel ( select et all ) file list from LSST ccd test bench 
    def __init__(self,dirall=['/Users/antilog/scratch/20160901'],Pickle=False,fkey={},nskip=0,nkeep=-1):
        # dirall : a list of directory/file   to read in  : the file header will be used to select or not the file if fkey is set  or the file will be read from the content of the pickle file (if Pickle is True ) , fkey will also be used for selection. 
        # fkey : {'selection_name' : {'fits_header_name':{'key':value} , ... }  : a triple  dictionary of { 'selection_name' : {header : {key:value}} } , to select a file 
        self.nkept=0
        self.all_file=[]
        self.clap=[]
        self.selection=[]
        self.fkey=fkey
        self.directory=sorted(dirall)
        self.stat={}
        self.nseen=0
        # loop on header 
        if Pickle : 
            self.all_file_from_pickle(dirall=dirall,fkey=fkey,nskip=nskip,nkeep=nkeep)
        else : 
            self.all_file_from_dir(dirall=dirall,fkey=fkey,nskip=nskip,nkeep=nkeep)
        return
    #
    def all_file_from_pickle(self,dirall,fkey,nskip,nkeep) :
        # fill all_file from the header of the files in dirall . 
        for dirname in self.directory :
            # have we allready all the needed file ? 
            if nkeep> 0 and self.nkept==nkeep : 
                break  
            # build the list of file for this directory 
            if (len(os.path.splitext(dirname)[1])>0) :
                file_list=glob.glob(dirname)
            else :
                file_list=glob.glob("%s/*.pkl" % (dirname))
            file_list.sort()
            # loop on files to select them if needed  
            for pickle_file  in file_list :
                # open the pickle file
                input=open(pickle_file,'rb')
                file=pickle.load(input)
                #
                for i_cur in range(len(file)) : 
                    #filename=file[i_cur].filename  
                    dname=file[i_cur].dir  
                    fname=file[i_cur].filename
                    clap=file[i_cur].clap
                    keep=True
                    if len(fkey)==0 :  
                        selection='Main'
                    else :
                        #
                        # is there any extra selection based on Header , key , value  ? 
                        for selection, sel_id in fkey.items() :
                            local_keep=False
                            # select of file name if any
                            if 'first' in sel_id.keys() :
                                if  fname < sel_id['first'] : continue
                            if 'last' in sel_id.keys() :
                                 if  fname > sel_id['last'] : continue
                            local_keep=True
                            # test key (=)  key+ (=>) key- (<=)
                            if 'key' in sel_id.keys() : 
                                for header, key_cur in sel_id['key'].items() :
                                    if not ( local_keep ) : break  
                                    for key, value in key_cur.items() :
                                        if ( key in file[i_cur].header[header] ) :
                                            if (file[i_cur].header[header][key] is None) or ( file[i_cur].header[header][key]!=value ):
                                                local_keep=False
                                                break 
                                        else : 
                                            local_keep=False
                                            break
                            # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                            if not(local_keep) : continue 
                            #
                            if 'key+' in sel_id.keys() : 
                                for header, key_cur in sel_id['key+'].items() :
                                    if not ( local_keep ) : break  
                                    for key, value in key_cur.items() :
                                        if ( key in file[i_cur].header[header] ) :
                                            if (file[i_cur].header[header][key] is None) or ( file[i_cur].header[header][key]<value ):
                                                local_keep=False
                                                break 
                                        else : 
                                            local_keep=False
                                            break
                            # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                            if not(local_keep) : continue
                            #
                            if 'key-' in sel_id.keys() : 
                                for header, key_cur in sel_id['key-'].items() :
                                    if not ( local_keep ) : break  
                                    for key, value in key_cur.items() :
                                        if ( key in file[i_cur].header[header] ) :
                                            if (file[i_cur].header[header][key] is None) or ( file[i_cur].header[header][key]>value ):
                                                local_keep=False
                                                break 
                                        else : 
                                            local_keep=False
                                            break
                            # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                            if local_keep : break 
                        keep=local_keep
                        #
                    if (keep) :
                        self.nseen+=1
                        if self.nseen>nskip :
                            self.all_file.append(dname+'/'+fname)
                            self.clap.append(clap)
                            self.selection.append(selection)
                            print ('%d : Selected %s File %s ' % (self.nkept,selection,fname) )
                            self.nkept+=1
                            if self.nkept==nkeep and nkeep > 0  :
                                # we selected the number of files requested 
                                break
        return
    
    def all_file_from_dir(self,dirall,fkey,nskip,nkeep) :
        # fill all_file from the header of the files in dirall . 
        for dirname in self.directory :
            # have we allready all the needed file ? 
            if nkeep> 0 and self.nkept==nkeep : 
                break  
            # build the list of file for this directory 
            if (len(os.path.splitext(dirname)[1])>0) :
                file_list=glob.glob(dirname)
            else :
                file_list=glob.glob("%s/*.fz" % (dirname))
            file_list.sort()
            # loop on files to select them if needed  
            for filename  in file_list :
                #
                keep=True
                if len(fkey)==0 :  
                    selection='Main'
                else : 
                    fits_is_open=False
                    dname, fname=os.path.split((filename))
                    # is there any extra selection based on Header , key , value  ? 
                    for selection, sel_id in fkey.items() :
                        local_keep=False
                        # select of file name if any
                        if 'first' in sel_id.keys() :
                            if  fname < sel_id['first'] : continue
                        if 'last' in sel_id.keys() :
                            if  fname > sel_id['last'] : continue
                        local_keep=True
                        if 'key' in sel_id.keys() : 
                            if not(fits_is_open) :
                                fitsfile=pyfits.open(filename)
                                fits_is_open=True
                            for header, key_cur in sel_id['key'].items() :
                                if not ( local_keep ) : break  
                                for key, value in key_cur.items() :
                                    if ( key in fitsfile[header].header ) :
                                        if ( fitsfile[header].header[key]!=value ):
                                            local_keep=False
                                            break 
                                    else : 
                                        local_keep=False
                                        break
                        # this file is ok for the current selection ,  remark : a file is selected in the first selection that is compatible with 
                        if local_keep : break 
                    keep=local_keep
                    if fits_is_open :
                        if keep :
                            exptime=fitsfile['Primary'].header['EXPTIME']
                        fitsfile.close()
                        del fitsfile
                #
                if (keep) :
                    self.nseen+=1
                    if self.nseen>nskip :
                        self.all_file.append(filename)
                        self.selection.append(selection)
                        # to be updated with a call to clap
                        self.clap.append(exptime)
                        print ('%d : Selected %s File %s ' % (self.nkept,selection,filename) )
                        self.nkept+=1
                        if self.nkept==nkeep and nkeep > 0  :
                            # we selected the number of files requested 
                            break
        return
    #
    def diff_pair(self,selection,clip_sig=4.,ncov_x=2,ncov_y=2,nx=5,ny=20,ckeep=range(16),subover=True,xf_r=10,xl_r=524,yf_r=0,yl_r=2002):
        # Cliping in number of sigma for the variance and mean
        # default : clip_sig=4.
        # Size of the cov matrix 
        # default : ncov_x=2
        #           ncov_y=2
        # The cov matrix will be computed in nx*ny bin for each  amplifier , not just 1 cov matrix per amplifier ( amplifier is subdivided in nx*ny boxes)    
        # default : nx=5
        #            ny=20 
        # size of the ccd
        # x : column , y : line , first & last+1 
        # ckeep : list of channel to study 
        # subover = True   ==> subtract overscan before taking the difference 
        # xf=10    : first column to use
        # xl=524   : last column+1 to use 
        # yf=0     : first line to use 
        # yl=2002  : last line + 1 to use 
        # xoverscan in number of column since the end ( e2v 0:576   overscan 524:576 ==> 531:576 : -45 sound good / away from CTE et all     
        xover=-45
        # in y exclusion region , and box size 
        yoff=10
        #
        yf=yf_r+yoff
        yl=yl_r-yoff
        ystep=int((yl-yf)/ny)
        # list of File to process 
        file_to_consider=[i for i,j in enumerate(self.selection) if j==selection ] 
        nread=len(file_to_consider)
        # Temporary strorage before averaging over all files 
        temp_var=np.zeros((16,ny,nx,ncov_y,ncov_x,int(nread/2)))
        temp_evar=np.zeros((16,ny,nx,ncov_y,ncov_x,int(nread/2)))
        temp_nentry=np.zeros((16,ny,nx,ncov_y,ncov_x,int(nread/2)))
        temp_mean=np.zeros((16,ny,nx,int(nread/2)))
        temp_clap=np.zeros((int(nread/2)))
        #
        file_cur=[]
        pair=0
        npair=0
        for ifile in file_to_consider :
            fitsfile=pyfits.open(self.all_file[ifile])
            file_cur.append(fitsfile)
            over=fitsfile[13].data[100:120,xover:].mean()
            print ('%d : Selected File %s (flux = %f , overscan = %f ) ' % (npair,self.all_file[ifile],fitsfile[13].data[100:120,10:120].mean()-over,over))
            if pair==1 :
                # proceed with the processing of the difference
                # loop on channels 
                for k in ckeep :
                    # x : column
                    if ( k in  [0,7,8,15] ) :
                        # exclude edge , for edge amplifier   
                        xoff=25
                    else : 
                        xoff=10
                    xf=xf_r+xoff
                    xl=xl_r-xoff
                    xstep=int((xl-xf)/nx)
                    #
                    ix=0
                    iy=0
                    # loop on the different sub-boxes of a channel
                    for y in range(yf,yl,ystep) :
                        if subover :
                            over=[]
                            for i in range(2) : 
                                # for each pair of file compute the overscan for the considered lines  :
                                over.append(file_cur[i][k+1].data[y:y+ystep,xover:].mean())
                        for x in range(xf,xl,xstep) :
                            im=[]
                            imf=[]
                            for i in range(2) : 
                                # for each pair of file compute the mean and sigma in each sub-box  :
                                mr=file_cur[i][k+1].data[y:y+ystep,x:x+xstep].mean()
                                stdr=file_cur[i][k+1].data[y:y+ystep,x:x+xstep].std()
                                # fill each pixel value and a mask to reject events at more than "clip_sigma" sigma in this local box
                                if subover :
                                    im.append(np.ma.array(file_cur[i][k+1].data[y:y+ystep,x:x+xstep]-over[i],mask=abs(mr-file_cur[i][k+1].data[y:y+ystep,x:x+xstep])>(clip_sig*stdr)))
                                else : 
                                    im.append(np.ma.array(file_cur[i][k+1].data[y:y+ystep,x:x+xstep],mask=abs(mr-file_cur[i][k+1].data[y:y+ystep,x:x+xstep])>(clip_sig*stdr)))
                                # compute the masked mean for this box & file  
                                imf.append(im[i].mean())
                            # compute the flux mean of the 2 files 
                            temp_mean[k,iy,ix,npair]=(imf[0]+imf[1])/2.
                            # compute the difference of the 2 files , reweigthing at the flux mean 
                            wd=(im[0]/imf[0]-im[1]/imf[1])*temp_mean[k,iy,ix,npair]
                            # compute the various variance and corvariance 
                            #diff_std=wd.std()
                            #diff_var=diff_std**2
                            for i in range(ncov_y) : 
                                for j in range (ncov_x) : 
                                    v1=wd[i:ystep,j:xstep]*wd[0:ystep-i,0:xstep-j]
                                    if i != 0 and j != 0 :
                                        n1=v1.count()
                                        v2=wd[i:ystep,0:xstep-j]*wd[0:ystep-i,j:xstep]
                                        n2=v2.count()
                                        temp_nentry[k,iy,ix,i,j,npair]=n1+n2
                                        #weighted average of the mean in function of the statistic in each subf sample 
                                        temp_var[k,iy,ix,i,j,npair]=(v1.mean()*n1+v2.mean()*n2)/temp_nentry[k,iy,ix,i,j,npair]
                                        #also weight the sigma , but the sigma refer to the mean of var not the sample dispersion , so let's / sqrt(nb_entries)  
                                        temp_evar[k,iy,ix,i,j,npair]=(v1.std()*n1+v2.std()*n2)/temp_nentry[k,iy,ix,i,j,npair]/np.sqrt(temp_nentry[k,iy,ix,i,j,npair])
                                    else:
                                        temp_nentry[k,iy,ix,i,j,npair]=v1.count()
                                        temp_var[k,iy,ix,i,j,npair]=v1.mean()
                                        temp_evar[k,iy,ix,i,j,npair]=v1.std()/np.sqrt(temp_nentry[k,iy,ix,i,j,npair])
                            ix+=1
                        #
                        ix=0
                        iy+=1
                # file clap value
                temp_clap[npair]=(self.clap[ifile]+self.clap[ifile-1])/2.
                # close the files ,
                for i in range(2) : 
                    # remark : due to the del the index 0 is needed for all ... ;-) 
                    file_cur[0].close()
                    del file_cur[0]
                # reset/set assocaited  "pair" flag & counters
                file_cur=[]
                pair=0
                npair+=1
            else : 
                pair+=1
        # average over the pairs of files 
        self.stat[selection]={'mean':temp_mean[:,:,:,0:npair].mean(axis=3),
                              'var':temp_var[:,:,:,:,:,0:npair].mean(axis=5),
                              'evar':temp_evar[:,:,:,:,:,0:npair].mean(axis=5),
                              'nentry':temp_nentry[:,:,:,:,:,0:npair].sum(axis=5),
                              'npair':npair,
                              'clap':temp_clap.mean(),
                              'sclap':temp_clap.std()}
        return
    def no_diff(self,selection,clip_sig=4.,ncov_x=2,ncov_y=2,nx=5,ny=20,ckeep=range(16)) :
        # Cliping in number of sigma for the variance and mean
        # default : clip_sig=4.
        # Size of the cov matrix 
        # default : ncov_x=2
        #           ncov_y=2
        # The cov matrix will be computed in nx*ny bin for each  amplifier , not just 1 cov matrix per amplifier ( amplifier is subdivided in nx*ny boxes)    
        # default : nx=5
        #            ny=20 
        # size of the ccd
        # x : column , y : line , first & last+1 
        # same function as diff pair , without diff and without overscan subtraction : used for pulse image 
        xf=10
        xl=524
        yf=0
        yl=2002
        # xoverscan in number of column since the end ( e2v 0:576   overscan 524:576 ==> 531:576 : -45 sound good / away from CTE et all     
        xover=-45
        # in y exclusion region , and box size 
        yoff=10
        ystep=int((yl-yf-2*yoff)/ny)
        # list of File to process 
        file_to_consider=[i for i,j in enumerate(self.selection) if j==selection ] 
        nread=len(file_to_consider)
        # Temporary strorage before averaging over all files 
        temp_var=np.zeros((16,ny,nx,ncov_y,ncov_x,nread))
        temp_evar=np.zeros((16,ny,nx,ncov_y,ncov_x,nread))
        temp_nentry=np.zeros((16,ny,nx,ncov_y,ncov_x,nread))
        temp_mean=np.zeros((16,ny,nx,nread))
        temp_clap=np.zeros((nread))
        #
        file_cur=[]
        pair=0
        npair=0
        for ifile in file_to_consider :
            fitsfile=pyfits.open(self.all_file[ifile])
            over=fitsfile[13].data[100:120,xover:].mean()
            print ('%d : Selected File %s (signal (no over sub)= %f , overscan = %f ) ' % (npair,self.all_file[ifile],fitsfile[13].data[100:120,10:120].mean(),over))
            # loop on channels 
            for k in ckeep :
                # x : column
                if ( k in  [0,7,8,15] ) :
                    # exclude edge , for edge amplifier   
                    xoff=25
                else : 
                    xoff=10
                xstep=int((xl-xf-2*xoff)/nx)
                #
                ix=0
                iy=0
                # loop on the different sub-boxes of a channel
                for y in range(yf,yl-ystep,ystep) :
                    for x in range(xf,xl-xstep,xstep) :
                        # for each pair of file compute the mean and sigma in each sub-box  :
                        mr=fitsfile[k+1].data[y:y+ystep,x:x+xstep].mean()
                        stdr=fitsfile[k+1].data[y:y+ystep,x:x+xstep].std()
                        # fill each pixel value and a mask to reject events at more than "clip_sigma" sigma in this local box
                        im=np.ma.array(fitsfile[k+1].data[y:y+ystep,x:x+xstep],mask=abs(mr-fitsfile[k+1].data[y:y+ystep,x:x+xstep])>(clip_sig*stdr))
                        # compute the masked mean for this box & file  
                        imf=im.mean()
                        # compute the flux means 
                        temp_mean[k,iy,ix,npair]=imf
                        # substract the mean to each pixel 
                        wd=(im-imf)
                        # compute the various variance and corvariance 
                        #diff_std=wd.std()
                        #diff_var=diff_std**2
                        for i in range(ncov_y) : 
                            for j in range (ncov_x) : 
                                v1=wd[i:ystep,j:xstep]*wd[0:ystep-i,0:xstep-j]
                                if i != 0 and j != 0 :
                                    n1=v1.count()
                                    v2=wd[i:ystep,0:xstep-j]*wd[0:ystep-i,j:xstep]
                                    n2=v2.count()
                                    temp_nentry[k,iy,ix,i,j,npair]=n1+n2
                                    #weighted average of the mean in function of the statistic in each subf sample 
                                    temp_var[k,iy,ix,i,j,npair]=(v1.mean()*n1+v2.mean()*n2)/temp_nentry[k,iy,ix,i,j,npair]
                                    temp_evar[k,iy,ix,i,j,npair]=(v1.std()*n1+v2.std()*n2)/temp_nentry[k,iy,ix,i,j,npair]/np.sqrt(temp_nentry[k,iy,ix,i,j,npair])
                                else:
                                    temp_nentry[k,iy,ix,i,j,npair]=v1.count()
                                    temp_var[k,iy,ix,i,j,npair]=v1.mean()
                                    temp_evar[k,iy,ix,i,j,npair]=v1.std()/np.sqrt(temp_nentry[k,iy,ix,i,j,npair])
                        ix+=1
                    #
                    ix=0
                    iy+=1
            # file the equivalent of the clap value
            temp_clap[npair]=float(self.selection[ifile][6:])
            # close the files ,
            fitsfile.close()
            del fitsfile
            # reset/set assocaited  "pair" flag & counters
            npair+=1
        # average over the pairs of files 
        self.stat[selection]={'mean':temp_mean[:,:,:,0:npair].mean(axis=3),
                              'var':temp_var[:,:,:,:,:,0:npair].mean(axis=5),
                              'evar':temp_evar[:,:,:,:,:,0:npair].mean(axis=5),
                              'nentry':temp_nentry[:,:,:,:,:,0:npair].sum(axis=5),
                              'npair':npair,
                              'clap':temp_clap.mean(),
                              'sclap':temp_clap.std()}
        return
    #
    def diff_box(self,selection,clip_sig=4.,ncov_x=2,ncov_y=2,nx=5,ny=20,ckeep=range(16)) :
        # Cliping in number of sigma for the variance and mean
        # default : clip_sig=4.
        # Size of the cov matrix 
        # default : ncov_x=2
        #           ncov_y=2
        # The cov matrix will be computed in nx*ny bin for each  amplifier , not just 1 cov matrix per amplifier ( amplifier is subdivided in nx*ny boxes)    
        # default : nx=5
        #            ny=20 
        # size of the ccd
        # x : column , y : line , first & last+1 
        xf=10
        xl=524
        yf=0
        yl=2002
        # xoverscan in number of column since the end ( e2v 0:576   overscan 524:576 ==> 531:576 : -45 sound good / away from CTE et all     
        xover=-45
        # in y exclusion region , and box size 
        yoff=10
        ystep=int((yl-yf-2*yoff)/ny)
        # list of File to process 
        file_to_consider=[i for i,j in enumerate(self.selection) if j==selection ] 
        nread=len(file_to_consider)
        # Temporary strorage before averaging over all files 
        temp_vcor=np.zeros((16,2,ny,nx,1+2*(ncov_y-1),1+2*(ncov_x-1),int(nread/2)))
        temp_nvcor=np.zeros((16,2,ny,nx,1+2*(ncov_y-1),1+2*(ncov_x-1),int(nread/2)))

        #
        file_cur=[]
        pair=0
        npair=0
        for ifile in file_to_consider :
            fitsfile=pyfits.open(self.all_file[ifile])
            file_cur.append(fitsfile)
            over=fitsfile[13].data[100:120,xover:].mean()
            print ('%d : Selected File %s (flux = %f , overscan = %f ) ' % (npair,self.all_file[ifile],fitsfile[13].data[100:120,10:120].mean()-over,over))
            if pair==1 :
                # proceed with the processing of the difference
                # loop on channels 
                for k in ckeep :
                    # x : column
                    if ( k in  [0,7,8,15] ) :
                        # exclude edge , for edge amplifier   
                        xoff=25
                    else : 
                        xoff=10
                    xstep=int((xl-xf-2*xoff)/nx)
                    #
                    ix=0
                    iy=0
                    # loop on the different sub-boxes of a channel
                    for y in range(yf,yl-ystep,ystep) :
                        over=[]
                        for i in range(2) : 
                            # for each pair of file compute the overscan for the considered lines  :
                            over.append(file_cur[i][k+1].data[y:y+ystep,xover:].mean())
                        for x in range(xf,xl-xstep,xstep) :
                            im=[]
                            imf=[]
                            for i in range(2) : 
                                # for each pair of file compute the mean and sigma in each sub-box  :
                                mr=file_cur[i][k+1].data[y:y+ystep,x:x+xstep].mean()
                                stdr=file_cur[i][k+1].data[y:y+ystep,x:x+xstep].std()
                                # fill each pixel value and a mask to reject events at more than "clip_sigma" sigma in this local box
                                im.append(np.ma.array(file_cur[i][k+1].data[y:y+ystep,x:x+xstep]-over[i],mask=abs(mr-file_cur[i][k+1].data[y:y+ystep,x:x+xstep])>(clip_sig*stdr)))
                                # compute the masked mean for this box & file  
                                imf.append(im[i].mean())
                            # compute the difference of the 2 files , reweigthing at the flux mean 
                            wd=(im[0]/imf[0]-im[1]/imf[1])*(imf[0]+imf[1])/2.
                            stdr=wd.std()
                            # compute the various variance and corvariance 
                            # Compute the vcor matrix , for flucuation above 2 sigma
                            # Flag all the pixel arround a bad pixel to exclude them all as a potentiel "good fluctuation location" ) , needed to be sure to be able to measure the flux in all the pixel arround a selected one ( in cov_y+1 cov_y-1  , idem with x ) 
                            wd_extra=np.ma.copy(wd) 
                            in_mask=np.ma.getmask(wd_extra)
                            ll=[[ii,jj] for ii in range(ystep) for jj in range(xstep) if in_mask[ii,jj]]
                            for i,j in ll : 
                                for ii in range(max(0,i-ncov_y+1),min(ystep,i+ncov_y),1) :
                                    for jj in range(max(0,j-ncov_x+1),min(xstep,j+ncov_x),1) : 
                                         wd_extra[ii,jj]=np.ma.masked
                            # look for all the pixel above / below 2 sigma (and selected as good one : a test on a masked pixel return always false)  
                            ll=[[ii,jj] for ii in range(ncov_y-1,ystep-ncov_y+1) for jj in range(ncov_x-1,xstep-ncov_x+1) if abs(wd_extra[ii,jj])>2*stdr]
                            # so in ll only index of wd above/bellow 2 sigam , with valid pixel in the range of the cov matrix are listed 
                            for i,j in ll : 
                                ssign=0
                                if wd[i,j] > 0 : ssign=1 
                                for ii in range(-ncov_y+1,ncov_y,1) :
                                    for jj in range(-ncov_x+1,ncov_x,1) :
                                      temp_vcor[k,ssign,iy,ix,ncov_y-1+ii,ncov_x-1+jj,npair]+=wd[ii+i,jj+j]
                                      temp_nvcor[k,ssign,iy,ix,ncov_y-1+ii,ncov_x-1+jj,npair]+=1
                            ix+=1
                        #
                        ix=0
                        iy+=1
                # close the files ,
                for i in range(2) : 
                    # remark : due to the del the index 0 is needed for all ... ;-) 
                    file_cur[0].close()
                    del file_cur[0]
                # reset/set assocaited  "pair" flag & counters
                file_cur=[]
                pair=0
                npair+=1
            else : 
                pair+=1
        # average over the pairs of files 
        nvcor=temp_nvcor[:,:,:,:,:,:,0:npair].sum(axis=6)
        vcor=temp_vcor[:,:,:,:,:,:,0:npair].sum(axis=6)
        vcor_normed=np.where(nvcor == 0., vcor, vcor/nvcor)
        self.stat[selection]={'vcor':vcor_normed,
                              'nvcor':nvcor,
                              'npair':npair}
        return
                        
if __name__ == "__main__" :
    usage = 'usage: %prog [options] arg '
    parser = OptionParser(usage=usage)
    #year default is current year
    parser.add_option('-f','--fkey', type='string', help='select from main header for files with fits key Name and Value {"selection_name_id":{"fits_header_name":{"header_key":"Value","header_key":"Value1"..}}} ) [%default]', default='')
    parser.add_option('-n', '--nkeep', type='int', help='Number of files to keep/process (0 = no limit)  [%default]', default=0)
    parser.add_option('-s', '--nskip', type='int', help='Number of files to skip , before starting select them for processing [%default]', default=0)

    opts,directory = parser.parse_args()

    if opts.fkey == '' :
        fkey={}
    else :
        fkey=eval(opts.fkey)
    #
    data=File_List(dirall=directory,fkey=fkey,nskip=opts.nskip,nkeep=opts.nkeep)

