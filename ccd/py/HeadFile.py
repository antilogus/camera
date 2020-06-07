#!/usr/bin/env python
try:
    import pyfits
except : 
    import astropy.io.fits as pyfits
import numpy as np
#
class HeadFile:
    """ Dump in object fits header """ 
    def __init__(self,filename,dname,fname,date,header,fkey) :
        # dirname : can be a directory name or a file with * at the moment it ends by .fz                                                                        
        #  ex : /Users/antilog/scratch/REB_DATA/20160513/linearity                                                                                               
        #  or : /Users/antilog/scratch/REB_DATA/20160513/linearity/reb3*.fz                                                                                      
        # first_col = Start of the signal to analyse ( before = pre-scan )                                                                                       
        # sum = key word for on the fly sum of the exposure (#0 ) , if = 1  sum each exposure , 2 sum odd and even exposure , 3 ...                       
        # keep  = .ture. ==> Keep the individual image                                                       
        self.dir=dname
        self.filename=fname
        self.date=date
        self.run=date
        self.comment='no-comment'
        # read the file and fill the header information 
        fitsfile=pyfits.open(filename)
        for ihead in range(len(header)):
            if ihead == 0 : 
                self.header={header[ihead]:header_dict(fitsfile,header[ihead],fkey[ihead])}
            else : 
                self.header.update({header[ihead]:header_dict(fitsfile,header[ihead],fkey[ihead])})
        # estimate the clap flux
        self.clap=GetCLAP(fitsfile)
        #
        fitsfile.close()
        del fitsfile
        return

    def print_summary(self) :
        # Print a summary of the file content based on the header 
        # 
        to_print='{:22s} '.format(self.filename)
        to_print+='{:12s} '.format(self.header['Primary']['TESTTYPE'])
        to_print+='{:12s} '.format(self.header['Primary']['IMGTYPE'])
        to_print+='{:7.2f} '.format(self.header['Primary']['EXPTIME'])
        try :
            to_print+='{:11.0f} '.format(self.clap)
        except :
            to_print+='{:11s} '.format('None')        
        if self.header['CHAN_12']['AVLIGHT'] != None  : 
            to_print+='{:6d} '.format(int(self.header['CHAN_12']['AVLIGHT']))
            to_print+='{:7.1f} '.format(self.header['CHAN_12']['STDLIGHT'])
        elif self.header['CHAN_12']['AVERAGE'] != None : 
            to_print+='{:6d} '.format(int(self.header['CHAN_12']['AVERAGE']))
            to_print+='{:7s} '.format('None')
        else :
            to_print+='{:6s} '.format('None')
            to_print+='{:7s} '.format('None')
        if  self.header['BSS']['VOLTAGE']  != None  :   
            to_print+='{:5.1f} '.format(self.header['BSS']['VOLTAGE'])
        else :
            to_print+='{:5s} '.format('None')  
        if  self.header['CCD_COND']['PU']  != None  :   
            to_print+='{:5.2f} '.format(self.header['CCD_COND']['PU'])
        else :
            to_print+='{:5s} '.format('None')  
        if  self.header['LASER']['POW_CH2']  != None  :   
            to_print+='{:4.1f} '.format(self.header['LASER']['POW_CH2'])
        else :
            to_print+='{:4s} '.format('None')  
        if self.header['CORNERSTONE']['WVLGTH']  != None  :
            to_print+='{:6.1f} '.format(self.header['CORNERSTONE']['WVLGTH'])
        else :
             to_print+='{:6s} '.format('None')  
        if self.header['QTH']['POWER'] != None : 
            to_print+='{:5.1f} '.format(self.header['QTH']['POWER'])
        else :
             to_print+='{:5s} '.format('None') 
        if self.header['XEHG']['POWER'] != None : 
            to_print+='{:5.1f} '.format(self.header['XEHG']['POWER'])
        else :
             to_print+='{:5s} '.format('None') 
        if self.header['XYZ']['XPOS'] !=None :
            to_print+='{:6.2f} '.format(self.header['XYZ']['XPOS'])
            to_print+='{:6.2f} '.format(self.header['XYZ']['APOS'])
        else :
            to_print+='{:6s} '.format('None')
            to_print+='{:6s} '.format('None')
        if self.header['SEQUENCER']['DATASUM']  !=None  :    
            to_print+='{:11s} '.format(self.header['SEQUENCER']['DATASUM'])
        else :
            to_print+='{:11s} '.format('None')
        to_print+='{:s} \n'.format(self.dir)
        return to_print
    def print_head_summary(self) :
        # Print a summary of the file content based on the header 
        # 
        to_print='#{:22s}'.format('Filename')
        to_print+='{:12s} '.format('TestType')
        to_print+='{:12s} '.format('ImgType')
        to_print+='{:7s} '.format('ExTime')
        to_print+='{:11s} '.format('CLAP_Flux')
        to_print+='{:6s} '.format('Ch12Fx')
        to_print+='{:7s} '.format('Ch12Si')
        to_print+='{:5s} '.format('V-BSS')
        to_print+='{:5s} '.format('Vup//')
        to_print+='{:4s} '.format('LAS2')
        to_print+='{:6s} '.format('LambMo')
        to_print+='{:5s} '.format('QTHPo')
        to_print+='{:5s} '.format('ArcPo')
        to_print+='{:6s} '.format('Xpos')
        to_print+='{:6s} '.format('Apos')
        to_print+='{:11s} '.format('SeqDsum')
        to_print+='Directory for day {:s} \n'.format(self.date)
        return to_print 
 
def header_dict(file,header,key) :
    fkey={}
    for i in range(len(key)) :
        try :
            val=file[header].header[key[i]]
        except:
            val=None
        fkey.update({key[i]:val})
    return fkey    

def GetCLAP(fitsfile):
    try:
        data=(fitsfile['CLAP'].data).astype('float64')
        # measure the reference flux : 
        ref=np.ma.array(data,mask=(data<-10) + (data>30)  )
        ref_mean=ref.mean()
        ref_std=ref.std()
        ref=np.ma.array(data,mask=abs(data-ref_mean)>4*ref_std)
        ref_mean=ref.mean()
        ref_std=ref.std()
        # measure the signal flux : 
        val=np.ma.array(data,mask=data>ref_mean-3*ref_std)
        m1=val.mean()
        s1=val.std()
        val2=np.ma.array(data,mask=abs(data-m1)>3*s1)
        m2=val2.mean()
        s2=val2.std()
        val3=np.ma.array(data,mask=abs(data-m2)>4*s2)
        # identify the noise (= 1 pixle glitch ) and replace it by the surrounding mean 
        l=np.ma.clump_masked(val3)
        smax=len(data)-1
        for s in l :
            index=s.__getattribute__('start')
            if index==s.__getattribute__('stop')-1 : 
                # if you pass teh edge it mean that the glitch is at the edge , so you take the closest pixel only [1] or [smax-1] 
                data[index]=(data[max(index-1,1)]+data[min(index+1,smax-1)])/2.
            # iem in the ref aread
        l=np.ma.clump_masked(ref)
        for s in l :
            index=s.__getattribute__('start')
            if index==s.__getattribute__('stop')-1 : 
                data[index]=(data[max(index-1,1)]+data[min(index+1,smax-1)])/2.
        # Compute the integrated flux
        clap=(data-ref.mean()).sum()
        # compute the exposure time 
        # to be done 
        return float(clap)
    except:
        return None
    return
