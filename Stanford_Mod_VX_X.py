"""
To be filled later
"""
import os, sys, pylab, math, warnings, datetime, time, shutil, string, gc, inspect, warnings, glob, pickle
import numpy as np

import scipy as sp
from scipy.io import loadmat
from scipy.optimize import curve_fit as cfit
from scipy.odr import ODR, Model, Data, RealData
from scipy.misc import factorial as fac
from scipy.signal import periodogram
from scipy.stats import skew
from scipy.linalg import toeplitz
from scipy.interpolate import interp1d

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import patches as patches

from mpl_toolkits.mplot3d import Axes3D
#from scipy import ndimage as nd
#import fitFunctions
#import fitVoigt


def analyzeData(noiseStr, pulseStr, templateStr = None, dpath = None, spath = None, fullpath = False,
                voltage = 'N/A', pw = 'N/a', **kwords):
    """
    Fill Later
    """
    # Keywords to set the input parameters for processTemplate()
    if 'tIter' in kwords.keys():
        tIter = kwords['tIter']
    else:
        tIter = 4
    if 'tendBin' in kwords.keys():
        tendBin = kwords['tendBin']
    else:
        tendBin = 1000
    if 'ttrigBin' in kwords.keys():
        ttrigBin = kwords['ttrigBin']
    else:
        ttrigBin = 150
    if 'Full' in kwords.keys():
        Full = kwords['Full']
    else:
        Full = True


    if templateStr is None:
        templateStr = pulseStr
    if dpath is None:
        dpath = '/Applications/Data-Processed/*/trace_data/{0}/{0}_*.mat'
    if spath is None:
        spath = '/Applications/Data-Processed/Pickles/{0}'
    print('Noise Series: ',noiseStr)
    print('Template Series: ',templateStr)
    if isinstance(pulseStr, str):
        print('Pulse Series: ',pulseStr)
    else:
        print('Pulse Series Starting at: ', '_'.join(pulseStr[0].split('/')[-1].split('_')[:2]))

    # Data file names
    if not fullpath:
        if isinstance(pulseStr, str):
            pulsefiles      = glob.glob(dpath.format(pulseStr))
            pulsefiles.sort(key = lambda line: int(line.split('_')[-1].split('.')[0]))
        elif isinstance(pulseStr, list):
            pulsefiles      = list()
            for line in pulseStr:
                pulsefiles      = (pulsefiles
                                   + sorted(glob.glob(dpath.format(line)),
                                            key = lambda line: int(line.split('_')[-1].split('.')[0])))
            if templateStr == pulseStr:
                templateStr = pulseStr[0]
            pulseStr = pulseStr[0]
    else:
        pulsefiles = pulseStr
        pulseStr = '_'.join(pulseStr[0].split('/')[-1].split('_')[:2])
        print(pulseStr)
    #Think these two lines are obsolete...
    #noisefiles          = glob.glob(dpath.format(noiseStr))
    #templatefiles       = glob.glob(dpath.format(templateStr))

    # Data save names
    noiseSaveFile       = 'PSDs_'+noiseStr+'.pkl'
    templateSaveFile    = 'Templates_'+templateStr+'.pkl'
    OFSaveFile          = 'OFResults_'+pulseStr+'.pkl'
    rawDataSaveFile     = 'Pulses_'+pulseStr+'.pkl'
    
    # Save settings
    saveNoise           = True
    saveTemplate        = True

    scaleTrig=False
    psInt=int(pulseStr[0:6])
    if(psInt > 180208):
        scaleTrig = False #True

    print('Loading Data...')
    if pulsefiles == []:
        print('No files in list')
        print(dpath.format(pulseStr))
        raise
    data = loadDataset(pulsefiles, scaleTrig = scaleTrig)

    Fs=data['Info']['Fs']

    voltage=str(data['Info']['Voltage'])
    pw=str(data['Info']['LaserPulseWidth'])
    lp=str(data['Info']['LaserPower'])
    atten=str(data['Info']['LaserAttenuation'])

    #forcing template to be fixed if not triggered on the laser
    if(not data['Info']['LaserTrigger']):
        templateSaveFile='Templates_180119_2256.pkl'

    print('Template: '+templateSaveFile)

    if os.path.exists(os.path.join(spath.format('PSDs/'), noiseSaveFile)):
        print('Using existing PSDs from '+noiseSaveFile)
        saveNoise=False
        with open(os.path.join(spath.format('PSDs/'), noiseSaveFile), 'rb') as output:
            res=pickle.load(output)
            PSD=res['PSD']
    else:
        print('No existing PSDs from '+os.path.join(spath.format('PSDs/'), noiseSaveFile))
        #get noise PSDs
        PSD=dict()
    
        fInterp=np.fft.rfftfreq(len(data['Pulses']['A'][0]),d=1/Fs)
    
        #Processing PSDs
        print('Processing PSDs')
        for chan in ['A','B','Total']:
            #get PSD for channel from noise traces
            noiseResults = processNoise(data['Noise'][chan], makePlots = False, autoCut = True,
                                        fs = Fs, fileStr = pulseStr+'_'+chan, saveResults = False,
                                        slopeCut = None)
            os.system('mkdir -p plots')
            #os.system('mv Noise*.png plots/')
            freq=noiseResults['f']
        
            #re-sample PSD to match PSD binning of signal
            fPSD=interp1d(freq,noiseResults['psdMean'])
            PSD[chan]=fPSD(fInterp)

    if(saveNoise):
        print('...Saving Noise: ' + os.path.join(spath.format('PSDs/'), noiseSaveFile))
        noiseResults={'freq':fInterp,'PSD':PSD}
        with open(os.path.join(spath.format('PSDs/'), noiseSaveFile), 'wb') as output:
            pickle.dump(noiseResults, output)

    if(os.path.exists(os.path.join(spath.format('Templates/'), templateSaveFile))):
        print('Using existing templates from '+templateSaveFile)
        saveTemplate=False
        with open(os.path.join(spath.format('Templates/'), templateSaveFile), 'rb') as output:
            template=pickle.load(output)
    else:
        print('No existing template from '+templateSaveFile)
        print('Creating Templates')
        #get pulses and create templates
        template=dict()
        for chan in ['A','B','Total']:
            template[chan] = processTemplate(data['Pulses'][chan],
                                             trigBin = ttrigBin,
                                             endBin = tendBin,
                                             cutIters = tIter)
            #normalize template amplitude to 1
            template[chan]-=np.mean(template[chan][0:100])
            template[chan]/=np.max(np.abs(template[chan]))

    if(saveTemplate):
        print('...Saving Templates')
        with open(os.path.join(spath.format('Templates/'), templateSaveFile), 'wb') as output:
            pickle.dump(template, output)

    for traceType in ['Noise','Pulses']:
    
        print('Running Optimal Filter for '+traceType)
        #process OF with templates and noise generated above
        resultsOF=dict()
        
        resultsOF['Info']=dict()
        for k,v in data['Info'].items():
            resultsOF['Info'][k]=v
    
        if(traceType == 'Noise'):
            extraStr='_Randoms'
            OFSaveFile='OFResults_'+pulseStr+'_Randoms.pkl'
        else:
            extraStr=''
            OFSaveFile='OFResults_'+pulseStr+'.pkl'
        if Full:
            resultsOF['TrigAmp']=np.max(data[traceType]['T'],axis=1)-np.min(data[traceType]['T'],axis=1)
            resultsOF['TrigTime']=np.argmax(data[traceType]['T'],axis=1)*(1.0/Fs)
            resultsOF['Triggered']=(resultsOF['TrigAmp']>1e-5)

            resultsOF['time']=data['Info']['time']
            resultsOF['StartTime']=min(resultsOF['time'])
            resultsOF['time']-=resultsOF['StartTime']
    
        tLen=len(template['Total'])
        trigpt=np.argmax(np.abs(template['Total']))
        ptsAfterTrig=tLen-trigpt
        time=np.arange(0,tLen)/Fs
    
        for chan in ['A','B','Total']:
            print('...',chan)
        
            chanPSD=PSD[chan]*1e-12
            #estimate noise variance from PSD
            varPSD=(np.fft.irfft(chanPSD**2)*Fs)[0]
            if Full:
                resultsOF['NoiseVariance_'+chan]=varPSD
        
                resultsOF['Sum_I_'+chan]=list()
                resultsOF['Sum_I2_'+chan]=list()
                resultsOF['OFP_A1_'+chan]=list()
                resultsOF['OFP_t1_'+chan]=list()
                resultsOF['OFP_A2_'+chan]=list()
                resultsOF['OFP_t2_'+chan]=list()
                resultsOF['OFP_chisq_'+chan]=list()
                resultsOF['OFP_chisqTime_'+chan]=list()
            
                resultsOF['Slope_'+chan]=list()
                resultsOF['Mean_'+chan]=list()
                resultsOF['MeanBase_'+chan]=list()
                resultsOF['Range_'+chan]=list()
            
            resultsOF['OF0_'+chan]=list()
            resultsOF['OF0_chisq_'+chan]=list()
            resultsOF['OF0_chisqTime_'+chan]=list()
        
            resultsOF['OF_'+chan]=list()
            resultsOF['OF_time_'+chan]=list()
            resultsOF['OF_chisq_'+chan]=list()
            resultsOF['OF_chisqTime_'+chan]=list()

            chanTemplate=template[chan]
            trigpt=np.argmax(np.abs(chanTemplate))
        
            count=0
            for trace in data[traceType][chan]:
                if Full:
                    resultsOF['Slope_'+chan].append(slope(time,trace,removeMeans=True))
                    resultsOF['Mean_'+chan].append(np.mean(trace))
                    resultsOF['MeanBase_'+chan].append(np.mean(trace[0:(trigpt-20)]))
                    resultsOF['Range_'+chan].append(np.ptp(trace))
            
                    resultsOF['Sum_I_'+chan].append(np.sum((trace-resultsOF['MeanBase_'+chan][-1])))
                    resultsOF['Sum_I2_'+chan].append(np.sum((trace-resultsOF['MeanBase_'+chan][-1])**2.0))
            
                count=count+1
                #calculate optimum filter amplitude in microamps
                a,t0,x = OptimumFilterAmplitude(trace, chanTemplate, chanPSD, Fs,normalize=False,withDelay=False)
                resultsOF['OF0_'+chan].append(a*1e6) #store amplitude in micro-amps
                resultsOF['OF0_chisq_'+chan].append(x)
            
                #compute no-delay time domain chi-square
                xt = np.var(trace - a*chanTemplate)/varPSD
                resultsOF['OF0_chisqTime_'+chan].append(xt)
            
                a,t,x = OptimumFilterAmplitude(trace, chanTemplate, chanPSD, Fs,normalize=False,withDelay=True)
                resultsOF['OF_'+chan].append(a*1e6)
                resultsOF['OF_time_'+chan].append(t)
                resultsOF['OF_chisq_'+chan].append(x)
            
                offset = int(t*Fs)
                if(offset > ptsAfterTrig):
                    offset-=tLen
                fit=chanTemplate.copy()
                fit=np.roll(fit,offset)
                fit*=a
            
                xt=np.var(trace-fit)/varPSD
                resultsOF['OF_chisqTime_'+chan].append(xt)

                if Full:
                    #calculate pileup OF
                    a1,a2,t1,t2,x = OptimumFilterAmplitude_Pileup(trace, chanTemplate,
                                                              chanPSD, Fs, downSample=8)
            
                    resultsOF['OFP_A1_'+chan].append(a1*1e6)
                    resultsOF['OFP_A2_'+chan].append(a2*1e6)
                    resultsOF['OFP_t1_'+chan].append(t1)
                    resultsOF['OFP_t2_'+chan].append(t2)
                    resultsOF['OFP_chisq_'+chan].append(x)
            
                    offset1=int(t1*Fs)
                    offset2=int(t2*Fs)
            
                    if(offset1 > ptsAfterTrig):
                        offset1-=tLen
                    if(offset2 > ptsAfterTrig):
                        offset2-=tLen
            
                    fit1=np.roll(chanTemplate.copy(),offset1)
                    fit1*=a1
            
                    fit2=np.roll(chanTemplate.copy(),offset2)
                    fit2*=a2
            
                    fit=fit1+fit2
                    xt=np.var(trace-fit)/varPSD
                    resultsOF['OFP_chisqTime_'+chan].append(xt)

            if Full:
                resultsOF['Slope_'+chan]=np.array(resultsOF['Slope_'+chan])
                resultsOF['Mean_'+chan]=np.array(resultsOF['Mean_'+chan])
                resultsOF['MeanBase_'+chan]=np.array(resultsOF['MeanBase_'+chan])
                resultsOF['Range_'+chan]=np.array(resultsOF['Range_'+chan])
                resultsOF['OFP_A1_'+chan]=np.array(resultsOF['OFP_A1_'+chan])
                resultsOF['OFP_t1_'+chan]=np.array(resultsOF['OFP_t1_'+chan])
                resultsOF['OFP_A2_'+chan]=np.array(resultsOF['OFP_A2_'+chan])
                resultsOF['OFP_t2_'+chan]=np.array(resultsOF['OFP_t2_'+chan])
                resultsOF['OFP_'+chan]=resultsOF['OFP_A1_'+chan]+resultsOF['OFP_A2_'+chan]
                resultsOF['OFP_chisq_'+chan]=np.array(resultsOF['OFP_chisq_'+chan])
                resultsOF['OFP_chisqTime_'+chan]=np.array(resultsOF['OFP_chisqTime_'+chan])
            
            resultsOF['OF0_'+chan]=np.array(resultsOF['OF0_'+chan])
            resultsOF['OF0_chisq_'+chan]=np.array(resultsOF['OF0_chisq_'+chan])
            resultsOF['OF0_chisqTime_'+chan]=np.array(resultsOF['OF0_chisqTime_'+chan])
        
            resultsOF['OF_'+chan]=np.array(resultsOF['OF_'+chan])
            resultsOF['OF_time_'+chan]=np.array(resultsOF['OF_time_'+chan])
            resultsOF['OF_chisq_'+chan]=np.array(resultsOF['OF_chisq_'+chan])
            resultsOF['OF_chisqTime_'+chan]=np.array(resultsOF['OF_chisqTime_'+chan])
        


        resultsOF['OF0_Sum']=resultsOF['OF0_A']+resultsOF['OF0_B']
        resultsOF['OF_Sum']=resultsOF['OF_A']+resultsOF['OF_B']
        if Full:
            resultsOF['OFP_A1_Sum']=resultsOF['OFP_A1_A']+resultsOF['OFP_A1_B']
            resultsOF['OFP_A2_Sum']=resultsOF['OFP_A2_A']+resultsOF['OFP_A2_B']
            resultsOF['OFP_Sum']=resultsOF['OFP_A']+resultsOF['OFP_B']
    del data # Clear memory space so help with saving maybe(?)
    print('Saving Results')
    #save results in pkl files for later use
    print(os.path.join(spath.format('OF_Results/'), OFSaveFile))
    with open(os.path.join(spath.format('OF_Results/'), OFSaveFile), 'wb') as output:
        pickle.dump(resultsOF, output)
    return None

def applyCut(inds, array):
    """
        Fill Later
    """
    for a, line in enumerate(array):
        array[a] = line[inds]
    return array

def getChannels(filelist, verbose = False):
    """
        Fill Later
    """
    if filelist == []:
        print('There are no files in list')
        raise
    if(type(filelist) == str):
        return getChannelsSingleFile(filelist,verbose=verbose)
    else:
        res1 = getChannelsSingleFile(filelist[0],verbose=verbose)
        combined = dict()
        combined['A'] = [res1['A']]
        combined['B'] = [res1['B']]
        combined['Total'] = [res1['Total']]
        combined['T'] = [res1['T']]
        combined['dVdI'] = res1['dVdI']
        combined['Fs'] = res1['Fs']
        combined['prop'] = res1['prop']
        combined['time'] = [res1['time']]
        combined['ai0'] = [res1['ai0']]
        combined['ai1'] = [res1['ai1']]
        combined['ai2'] = [res1['ai2']]
        #combined['pulsenum'] = [res1['pulsenum']]
        #combined['filename'] = [res1['filename']]
        combined['file'] = dict(zip(list(range(len(filelist))), filelist))
        if 'ai3' in res1.keys():
            combined['ai3'] = [res1['ai3']]
        
        for i in range(1,len(filelist)):
            try:
                res = getChannelsSingleFile(filelist[i],verbose=verbose)
                combined['A'].append(res['A'])
                combined['B'].append(res['B'])
                combined['Total'].append(res['Total'])
                combined['T'].append(res['T'])
                combined['time'].append(res['time'])
                combined['ai0'].append(res['ai0'])
                combined['ai1'].append(res['ai1'])
                combined['ai2'].append(res['ai2'])
                #combined['pulsenum'].append(res['pulsenum'])
                #combined['filename'].append(res['filename'])
                if 'ai3' in res.keys():
                    combined['ai3'].append(res['ai3'])
            except:
                print('Skipping '+filelist[i])

        combined['A'] = np.concatenate(combined['A'])
        combined['B'] = np.concatenate(combined['B'])
        combined['Total'] = np.concatenate(combined['Total'])
        combined['T'] = np.concatenate(combined['T'])
        combined['time'] = np.concatenate(combined['time'])
        combined['ai0'] = np.concatenate(combined['ai0'])
        combined['ai1'] = np.concatenate(combined['ai1'])
        combined['ai2'] = np.concatenate(combined['ai2'])
        #combined['pulsenum'] = np.concatenate(combined['pulsenum'])
        #combined['filename'] = np.concatenate(combined['filename'])
        if 'ai3' in combined.keys():
            combined['ai3'] = np.concatenate(combined['ai3'])
        else:
            print('ai3 is unavailable', end = '\r')
        #print(combined['T'].shape)
        combined['filenum'] = len(filelist)
        return combined

def getChannelsSingleFile(filename, verbose = False):
    """
        Function for opening a .mat file from the Stanford DAQ and returns a dictionary
        that contains the data.
        
        Parameters
        ----------
        filename : str
        The filename that will be opened. Should be a Stanford DAQ .mat file.
        
        Returns
        -------
        res : dict
        A dictionary that has all of the needed data taken from a Stanford DAQ
        .mat file.
        
        """
    
    res = loadmat(filename, squeeze_me = False)
    prop = res['exp_prop']
    data = res['data_post']
    
    exp_prop = dict()
    for line in prop.dtype.names:
        try:
            val = prop[line][0][0][0]
        except IndexError:
            val = 'Nothing'
        if type(val) is str:
            exp_prop[line] = val
        elif val.size == 1:
            exp_prop[line] = val[0]
        else:
            exp_prop[line] = np.array(val, dtype = 'f')

    gains = np.array(prop['SRS'][0][0][0], dtype = 'f')
    rfbs = np.array(prop['Rfb'][0][0][0], dtype = 'f')
    turns = np.array(prop['turn_ratio'][0][0][0], dtype = 'f')
    fs = float(prop['sample_rate'][0][0][0])
    minnum = min(len(gains), len(rfbs), len(turns))

    if 'daqrange' in prop.dtype.names:
        #The factor of 2 is because the range is +/- (added on 190905)
        convert = 2*prop['daqrange'][0][0][0][0]/2**12
    else:
        convert = 1
    
    ch1 = data[:,:,0]*convert
    ch2 = data[:,:,1]*convert
    try:
        trig = data[:,:,2]*convert
    except IndexError:
        trig = np.array([])
    ai0 = ch1[:]
    ai1 = ch2[:]
    ai2 = trig[:]
    try:
        ai3 = data[:, :, 3]*convert
    except:
        pass

    try:
        ttable  = np.array([24*3600.0, 3600.0, 60.0, 1.0])
        reltime = res['t_rel_trig'].squeeze()
        abstime = res['t_abs_trig'].squeeze()
        timestamp = abstime[:,2:].dot(ttable)+reltime
    except:
        timestamp = np.arange(0,len(ch1))
    
    dvdi = turns[:minnum]*rfbs[:minnum]*gains[:minnum]
    didv = 1.0/dvdi
    
    res = dict()
    res['A'] = ch1*didv[0]
    res['B'] = ch2*didv[1]
    res['Total'] = res['A']+res['B']
    res['T'] = trig
    res['dVdI'] = dvdi
    res['Fs'] = fs
    res['prop'] = prop
    res['filenum'] = 1
    res['time'] = timestamp
    res['exp_prop'] = exp_prop
    res['ai0'] = ai0
    res['ai1'] = ai1
    res['ai2'] = ai2
    try:
        res['ai3'] = ai3
    except:
        pass
    return res

"""
def getChannelsSingleFile(filename, verbose = False):
    if(verbose):
        print('Loading',filename, end = '\r')
    res = loadmat(filename, squeeze_me = False)
    prop = res['exp_prop']
    data = res['data_post']

    exp_prop = dict()
    for line in prop.dtype.names:
        try:
            val     = prop[line][0][0][0]
        except IndexError:
            val     = 'Nothing'
        if type(val) is str:
            exp_prop[line] = val
        elif val.size == 1:
            exp_prop[line] = val[0]
        else:
            exp_prop[line] = np.array(val, dtype = 'f')
            
    Gains = np.array(prop['SRS'][0][0][0], dtype = 'f')
    Rfbs = np.array(prop['Rfb'][0][0][0], dtype = 'f')
    Turns = np.array(prop['turn_ratio'][0][0][0], dtype = 'f')
    Fs = float(prop['sample_rate'][0][0][0])
    minnum = min(len(Gains), len(Rfbs), len(Turns))

    if 'daqrange' in prop.dtype.names:
        convert = prop['daqrange'][0][0][0][0]/2**12
    else:
        convert = 1
        
    ch1 = data[:,:,0]*convert
    ch2 = data[:,:,1]*convert
    try:
        trig = data[:,:,2]*convert
    except IndexError:
        trig = np.array([])
    ai0 = ch1[:]
    ai1 = ch2[:]
    ai2 = trig[:]
    try:
        ai3 = data[:, :, 3]*convert
    except:
        pass

    try:
        ttable  = np.array([24*3600.0,3600.0,60.0,1.0])
        reltime = res['t_rel_trig'].squeeze()
        abstime = res['t_abs_trig'].squeeze()
        timestamp = abstime[:,2:].dot(ttable)+reltime
    except:
        timestamp = np.arange(0,len(ch1))
        if(verbose):
            print('No Timestamps Found')

    dVdI        = Turns[:minnum]*Rfbs[:minnum]*Gains[:minnum]
    dIdV        = 1.0/dVdI
    
    res = dict()
    res['A'] = ch1*dIdV[0]
    res['B'] = ch2*dIdV[1]
    res['Total'] = res['A']+res['B']
    res['T'] = trig
    res['dVdI'] = dVdI
    res['Fs'] = Fs
    res['prop'] = prop
    res['filenum'] = 1
    res['time'] = timestamp
    res['exp_prop'] = exp_prop
    res['ai0'] = ai0
    res['ai1'] = ai1
    res['ai2'] = ai2
    res['pulsenum'] = np.arange(len(ai0))
    res['filename'] = np.array([filename.split('/')[-1]]*len(ai0))
    try:
        res['ai3'] = ai3
    except:
        pass
    return res
"""
def loadDataset(filelist, scaleTrig = True, **kwords):
    
    # keywords for tailoring the division of a trace into 'noise' and 'pulse'
    
    
    traces = getChannels(filelist, verbose=True)
    tlen=len(traces['A'][0])
    
    dataset=dict()
    dataset['Pulses']=dict()
    dataset['Noise']=dict()
    dataset['Info']=dict()
    
    chans=['A','B','T','Total']
    exclude = ['ai0', 'ai1', 'ai2', 'ai3']
    for chan in chans:
        
        if(chan == 'T'):
            #need to remove offset which tracks MC temperature
            trigTrace=traces['T']
            dataset['Info']['R_MC']=-1.0*np.mean(trigTrace,axis=1)
            trigChan=np.abs(np.add(trigTrace.transpose(),dataset['Info']['R_MC'])).transpose()
            print(np.shape(trigChan))
            #if offset present, need to scale back up due to attenuation
            if(scaleTrig):
                trigChan=trigChan*6e2
            dataset['Noise'][chan]=trigChan[:,0:int(tlen/2)]
            dataset['Pulses'][chan]=trigChan[:,int(tlen/2):]
        else:
            dataset['Noise'][chan]=traces[chan][:,0:int(tlen/2)]
            #dataset['Pulses'][chan]=traces[chan][:,int(tlen/2) - 100:-100]
            dataset['Pulses'][chan]=traces[chan][:,int(tlen/2):]
    for k in traces.keys():
        if k not in chans and k not in exclude:
            dataset['Info'][k]=traces[k]


    #if('R_MC' in dataset['Info'].keys()):
    #    dataset['Info']['Temperature']=temperature(1.0*dataset['Info']['R_MC']*dataset['Info']['dVdI'][2]*-1e3)
    
    try:
        dataset['Info']['V_S1']=dataset['Info']['prop']['hvs1'][0][0][0][0]
        dataset['Info']['V_S2']=dataset['Info']['prop']['hvs2'][0][0][0][0]
        dataset['Info']['Voltage']=dataset['Info']['V_S2']+dataset['Info']['V_S1']
    except:
        print('No Voltage Fields')
        dataset['Info']['V_S1']='NA'
        dataset['Info']['V_S2']='NA'
        dataset['Info']['Voltage']='NA'
    
    try:
        dataset['Info']['LaserRate']=dataset['Info']['prop']['ltr'][0][0][0]
    except:
        print('No Laser Rate Field')
        dataset['Info']['LaserRate']='NA'
    
    try:
        dataset['Info']['LaserPulseWidth']=dataset['Info']['prop']['lwd'][0][0][0]
    except:
        print('No Laser PW Field')
        dataset['Info']['LaserPulseWidth']='NA'
    
    try:
        dataset['Info']['LaserPower']=dataset['Info']['prop']['lpk'][0][0][0]
    except:
        print('No Laser Power Field')
        dataset['Info']['LaserPower']='NA'
    
    try:
        dataset['Info']['LaserAttenuation']=dataset['Info']['prop']['laseratten'][0][0][0][0]
    except:
        print('No Attenuation Field')
        dataset['Info']['LaserAttenuation']='NA'
    
    try:
        dataset['Info']['TriggerThreshold']=dataset['Info']['prop']['threshtrig'][0][0][0][0]
    except:
        print('No Trigger Threshold Field')
        dataset['Info']['TriggerThreshold']='NA'
    
    try:
        dataset['Info']['TriggerChannel']=dataset['Info']['prop']['trigchan'][0][0][0][0]
    except:
        print('No Trigger Channel Field')
        dataset['Info']['TriggerChannel']=0
    
    if(dataset['Info']['TriggerChannel'] == 3):
        dataset['Info']['LaserTrigger']=True
    else:
        dataset['Info']['LaserTrigger']=False

    try:
        dataset['Info']['Temp']=dataset['Info']['prop']['T_MC'][0][0][0][0]
    except:
        print('No Temperature Field')
        dataset['Info']['Temp']='NA'

    return dataset

def massOFA(Signal, Template, NoisePSD, Fs, withDelay = True, normalize = False, coupling = 'AC'):
    """
    Takes a signal trace and calculates the amplitude assuming the profile of template and noise from NoisePSD
    """
    dt = 1.0/Fs
    Ns = float(len(Signal))
    T = Ns*dt
    dnu = 1.0/T

    if(normalize):
        tRange = max(Template) - min(Template)
        Template = Template*1.0/tRange
        #Template /= tRange

    #take one-sided fft of Signal and Template
    Sf = np.fft.rfft(Signal, axis = 1)
    Tf = np.fft.rfft(Template)

    #check for compatibility between PSD and fft
    if(len(NoisePSD) != len(Sf[0])):
        raise ValueError("PSD length incompatible with signal size")
    
    #take squared noise PSD
    J = NoisePSD**2.0*dnu
    
    #If AC coupled, the 0 component of the PSD is non-sensical
    #If DC coupled, ignoring the DC component will still give the correct amplitude
    if(coupling == 'AC'):
        J[0] = np.inf

    #find optimum filter and norm
    OF = Tf.conjugate()/J
    #Norm = np.real(OF.dot(Tf)) # As written this is not normalized...
    Norm = np.sqrt(np.real(OF.dot(OF.conjugate()))) # FP Replacement
    OFp = OF/Norm

    Sfilt = OFp*Sf

    #compute OF with delay
    if(withDelay):
        #have to correct for numpy rfft convention by multiplying by N/2
        At = np.fft.irfft(Sfilt, axis = 1)*Ns/2.0
        
        #signal pary of chi-square
        chi0 = np.real(np.sum(Sf*Sf.conjugate()/J, axis = 1))*dt/Ns
        
        #fitting part of chi-square
        chit = (At**2)*Norm*dt/Ns
        
        #sum parts of chi-square
        chi = (np.ones((len(chit[0]), len(chit)))*chi0).transpose() - chit

        #find time of best-fit
        bInd = np.argmin(chi, axis = 1)
        A = At[np.arange(len(bInd)),bInd]
        Xr = chi[np.arange(len(bInd)),bInd]
        t0 = bInd*dt - T*(bInd*dt - T == 0).astype(int)

    #compute OF amplitude no delay
    else:
        A = np.real(np.sum(Sfilt, axis = 1))
        t0 = 0.0
    
        #compute reduced chi-square, normal chi-square divided by number of points
        delta = Sf - Tf*(np.ones((len(Sf[0]), len(Sf)))*A).transpose()
        Xr = np.real(np.sum((delta.conjugate()*(delta/J)), axis = 1))*dt/Ns

    return A,t0,Xr

def OptimumFilterAmplitude(Signal, Template, NoisePSD, Fs, withDelay = True, normalize = False, coupling = 'AC'):
    """
    Takes a signal trace and calculates the amplitude assuming the profile of template and noise from NoisePSD
    """
    dt = 1.0/Fs
    Ns = float(len(Signal))
    T = Ns*dt
    dnu = 1.0/T

    if(normalize):
        tRange = max(Template) - min(Template)
        Template = Template*1.0/tRange
        #Template /= tRange

    #take one-sided fft of Signal and Template
    Sf = np.fft.rfft(Signal)
    Tf = np.fft.rfft(Template)

    #check for compatibility between PSD and fft
    if(len(NoisePSD) != len(Sf)):
        raise ValueError("PSD length incompatible with signal size")
    
    #take squared noise PSD
    J = NoisePSD**2.0*dnu
    
    #If AC coupled, the 0 component of the PSD is non-sensical
    #If DC coupled, ignoring the DC component will still give the correct amplitude
    if(coupling == 'AC'):
        J[0] = np.inf

    #find optimum filter and norm
    OF = Tf.conjugate()/J
    #Norm = np.real(OF.dot(Tf)) # As written this is not normalized...
    Norm = np.sqrt(np.real(OF.dot(OF.conjugate()))) # FP Replacement
    OFp = OF/Norm

    Sfilt = OFp*Sf
    #global trouble
    #trouble = [Sf, Tf, NoisePSD, Fs, J, OF, Norm, OFp]
    #raise ValueError

    #compute OF with delay
    if(withDelay):
        #have to correct for numpy rfft convention by multiplying by N/2
        At = np.fft.irfft(Sfilt)*Ns/2.0
        
        #signal pary of chi-square
        chi0 = np.real(np.dot(Sf.conjugate()/J, Sf))*dt/Ns
        
        #fitting part of chi-square
        chit = (At**2)*Norm*dt/Ns
        
        #sum parts of chi-square
        chi = chi0 - chit
        
        #find time of best-fit
        bInd = np.argmin(chi)
        A = At[bInd]
        Xr = chi[bInd]
        t0 = bInd*dt

        if(t0 == T):
            t0 -= T

    #compute OF amplitude no delay
    else:
        A = np.real(np.sum(Sfilt))
        t0 = 0.0
    
        #compute reduced chi-square, normal chi-square divided by number of points
        delta = Sf - Tf*A
        Xr = np.real(np.dot(delta/J,delta.conjugate()))*dt/Ns

    return A,t0,Xr

def OptimumFilterAmplitude_Pileup(Signal, Template, NoisePSD, Fs, downSample = 10):

    dt = 1.0/Fs
    Ns = float(len(Signal))
    T = Ns*dt
    dnu = 1.0/T

    #take one-sided fft of Signal and Template
    Sf = np.fft.rfft(Signal)
    Tf = np.fft.rfft(Template)

    #check for compatibility between PSD and fft
    if(len(NoisePSD) != len(Sf)):
        raise ValueError("PSD length incompatible with signal size")
    
    #take squared noise PSD
    J = NoisePSD**2.0*dnu
    
    #TEMPORARY: NEED TO SWITCH TO FULL FFT
    J[0] = np.inf

    #find optimum filter and norm
    OF = Tf.conjugate()/J
    Norm = np.real(OF.dot(Tf))
    OFp = OF/Norm

    #filter template and trace
    Sfilt = OFp*Sf
    Tfilt = OFp*Tf

    #compute OF with delay
    
    #have to correct for numpy rfft convention by multiplying by N/2
    At = np.fft.irfft(Sfilt)*Ns/2.0
    Gt = np.real(np.fft.irfft(Tfilt))*Ns/2.0
        
    #signal part of chi-square
    chi0 = np.real(np.dot(Sf.conjugate()/J,Sf))*dt/Ns

    #construct matrices for t0 in the row and delta t in the column
    ds = int(downSample)
    At0 = timeMatrix(At[::ds])
    Atd = timeShiftMatrix(At[::ds])
    GM = (timeMatrix(Gt[::ds])).transpose()
    
    #compute full solution
    A2 = (Atd - At0*GM)/(1.0 - GM**2 + 1e-10)
    A1 = At0-A2*GM
    
    #compute chi-square
    chit = (A1**2 + A2**2 + 2*A1*A2*GM)*Norm*dt/Ns
        
    #sum parts of chi-square
    chi = chi0 - chit
        
    #find time of best-fit
    dti, ti = np.unravel_index(np.argmin(chi), np.shape(chi))
    A1s = A1[dti,ti]
    A2s = A2[dti,ti]
    Xr = chi[dti,ti]
    dtS = float(downSample)*dt
    t1 = float(ti)*dtS
    t2 = t1 + float(dti)*dtS

    #keep times in domain
    if(t1 >= T):
        t1 -= T
    if(t2 >= T):
        t2 -= T

    #first time is first amplitude
    if(t2 < t1):
        t1t = t1
        t1 = t2
        t2 = t1t

        A1t = A1
        A1 = A2
        A2 = A1t

    return A1s, A2s, t1, t2, Xr

def processNoise(rawTraces, traceGain = 1.0, fs = 625e3, fit = False, autoCut = False,
                 cutIters = 2, fileStr = '', makePlots = False, saveResults = True,
                 slopeCut = None, meanCutHigh = None, meanCutLow = None):

    #converting sampling rate to time step                                                                                                                                       
    dt = (1/fs)*1e6 #time in microseconds                                                                                                                                          

    #get trace x values, scale to microseconds                                                                                                                                   
    xvals = np.arange(0, len(rawTraces[0]))
    xvals_scaled = xvals*dt
    
    #for storing results
    tt = 1e6*rawTraces/traceGain
    ranges = tt.max(1) - tt.min(1)
    means = np.mean(tt, axis = 1)
    slopes = slope(xvals, tt)
    freq, ps = periodogram(tt, fs = fs, axis = 1)
    ps = np.sqrt(ps)*1e6

    traces = tt - (np.ones((len(tt[0]), len(tt)))*slopes).transpose()*xvals
    skewnesses = np.abs(skew(traces, axis = 1))
    freq, psc = periodogram(traces, fs = fs, axis = 1)
    psc = np.sqrt(psc)*1e6

    if(slopeCut != None):
        inds    = np.abs(slopes) < slopeCut/1e6
        [means, ranges, slopes, skewnesses, traces, ps, psc] = applyCut(inds, [means, ranges, slopes, skewnesses, traces, ps, psc])
 
    if(meanCutHigh != None):
        inds    = means < meanCutHigh
        [means, ranges, slopes, skewnesses, traces, ps, psc] = applyCut(inds, [means, ranges, slopes, skewnesses, traces, ps, psc])
 
    if(meanCutLow != None):
        inds    = means > meanCutLow
        [means, ranges, slopes, skewnesses, traces, ps, psc] = applyCut(inds, [means, ranges, slopes, skewnesses, traces, ps, psc])
    
    if(autoCut):
        for i in range(0, cutIters):
            #make cuts
            inds    = (removeOutliers(means) & removeOutliers(ranges) &
                       removeOutliers(slopes) & removeOutliers(skewnesses))
            [means, ranges, slopes, skewnesses, traces, ps, psc] = applyCut(inds, [means, ranges, slopes, skewnesses, traces, ps, psc])                                                                       
        print('Acceptance: ', float(len(slopes))/float(len(rawTraces)))
        
    #store results                                                                                                                                                                                               
    psd_corr_mean       = np.mean(psc, axis = 0)
    psd_corr_err        = stdComplex(psc, axis = 0)
    tmean               = np.mean(traces, axis = 0)
    f2, psd_correlated  = periodogram(tmean, fs = fs)
    psd_correlated      = np.sqrt(psd_correlated)*1e6
    offset              = np.mean(means)
    saveData = {'f':freq, 'psdMean':psd_corr_mean,
                'psdErr':psd_corr_err, 'psdCorr':psd_correlated,
                't':xvals_scaled, 'traceMean':tmean, 'Is0':offset}
    
    if(saveResults):
        pkl_file = open('Noise' + fileStr + '.pkl','wb')
        pickle.dump(saveData, pkl_file)
        pkl_file.close()
    return saveData

def processTemplate(rawTraces, autocut = True, cutIters = 2,trigBin = 450, endBin = 2000):
    """
        Inputs:
        rawTraces   - np.array of traces used to calculate the template only upto 10K will be used
        autocut     - True, default,
    
    """
    x = np.arange(0,len(rawTraces[0]))
    nonTrace = np.logical_or(x < trigBin, x > endBin)
    xNT = x[nonTrace]
    xSlope= xNT - np.mean(xNT)
    
    #for storing results
    tt = 1e6*rawTraces[:10000,nonTrace] # Current in microAmps, use only up to 10000 traces
    ranges = tt.max(1) - tt.min(1)
    means = np.mean(tt, axis = 1)
    slopes = slope(xSlope, tt - (np.ones((len(tt[0]), len(tt)))*means).transpose(), removeMeans = False)
    skewnesses = np.abs(skew(tt, axis = 1))
    traces = rawTraces[:10000]

    if(autocut):
        for i in range(0,cutIters):
            #make cuts
            inds    = (removeOutliers(means) & removeOutliers(ranges) &
                       removeOutliers(slopes) & removeOutliers(skewnesses))
            [means, ranges, slopes, skewnesses] = applyCut(inds, [means, ranges, slopes, skewnesses])
            traces = traces[inds, :]
    tmean = np.mean(means)
    template = np.mean(traces, axis = 0)
    
    plt.figure(0)
    for trace in rawTraces:
        plt.plot(trace, alpha = 0.05)    
    plt.plot(template, color = 'black')
    plt.show(block = False)
    print('Acceptance',float(len(traces))/float(len(rawTraces)))
    return template

def slope(x, y, removeMeans = True):
    """Computes the maximum likelihood slope of a set of x and y points
    
    Arguments:
    - x: Array of real-valued independent variables
    - y: Array of real-valued dependent variables

    Optional Keyword Arguments:
    - removeMeans (default True) - whether the mean needs to be removed from x and y. Set to false if mean has already been removed from the 

    Return Value:
    - slope: Maximum likelihood slope estimate calculated as
                 sum((y-<y>)*(x-<x>))/sum((x-<x>)^2)

    """
    if(removeMeans):
        xTemp = x - np.mean(x)
        if len(np.array(y).shape) == 1:
            y = np.array([y])
        yTemp = y - (np.ones((len(y[0]), len(y)))*np.mean(y, axis = 1)).transpose()
    else:
        xTemp = x
        yTemp = y
    return xTemp.dot(yTemp.transpose())/xTemp.dot(xTemp)

def stdComplex(x, axis = 0):
    """
    Return complex standard deviation (individually computed for real and imaginary components)
    """
    rstd = np.std(x.real, axis = axis)
    istd = np.std(x.imag, axis = axis)
    return rstd + 1.0j*istd

def removeOutliers(x, maxiter = 20, skewTarget = 0.05):
    """Function to return indices of inlying points
    Tries to symmetrize cut distribution to remove asymmetric outliers
    Arguments: 
    - x: Array of real-valued variables from which to remove outliers
    
    Optional Keyword Arguments:
    - maxiter (default 20): Maximum number of iterations to continue to minimize skewness
    - skewTarget (default 0.05): Desired residual skewness of distribution
    
    Returns:
    - inds: boolean indicies indicating which values to select/reject, same length as x
    """
    i = 1
    inds = (x != np.inf)
    sk = skew(x[inds])
    while(sk > skewTarget):
        dmed = x - np.median(x[inds])
        dist = np.min([abs(min(dmed)), abs(max(dmed))])
        inds = inds & (abs(dmed) < dist)
        sk = skew(x[inds])
        if(i > maxiter):
            break
        i += 1
    return inds

def timeMatrix(trace):
    return np.full([len(trace), len(trace)], trace)

def timeShiftMatrix(trace):
    x = np.flip(trace, 0)
    xf = np.flip(x, -1)
    xf = np.roll(xf, 1)
    xShift = toeplitz(xf, x)
    xM = np.flip(xShift, 1)
    return xM


def fitfun(x, s, c0, c1, n, L, Lupp, L1):
    L0 = (1 - L1 - Lupp)
    dx = x[1] - x[0]
    y = L*(L0*gzero(x, c0, s, n) +
           Lupp*dx/2/(c1 - c0)*(erf((x - c0)/np.sqrt(2*s**2)) -
                                erf((x-c1)/np.sqrt(2*s**2))) +
           L1*gaussian(x, c1, s))
    return y

# These functions are intended to perform the trapping and impact ionization fits...
def Amijn(A, m, i, j, n):
    A1, Alow, Aupp = A
    numerator = (A1**i)*(Alow**j)*Aupp**(m - i - j)*(-1)**(m - i - n)*(m - i)*fac(m)
    denominator = fac(i)*fac(j)*fac(m - i - j)*fac(n)*fac(m - i - n)
    return numerator/denominator

def h(x, m, A):
    A1, Alow, Aupp = A
    theta = np.heaviside
    adx = np.sum((m - 1 <= x) & (x < m))
    out = np.zeros(x.size)
    #out[x == m] = A1**m
    out[(m - 1 <= x) & (x < m)] = m*Alow*A1**(m-1)
    out[(m <= x) & (x < m + 1)] = m*Aupp*A1**(m-1)
    for i in range(0, m - 2 + 1): # the plus one is intended to run m - 2
        for j in range(0, m - i + 1): # the plus one is intended to run m - i
            for n in range(1, m - i + 1):# the plus one is intended to run 1 to m - i
                out = out + Amijn(A, m, i, j, n)*(n - x)**(m - i - 1)*theta(n - x, 1)*theta(x - c[m-j], 0)
    return out*(1 - A1**m)/sum(out)

def hfitgood(x, m, A, c):
    A1, Alow, Aupp = A
    theta = np.heaviside
    out = np.zeros(x.size)
    dx = x[1] - x[0]
    #out[x == m] = A1**m
    out[(c[m-1] <= x) & (x < c[m])] = m*Alow*A1**(m-1)
    out[(c[m] <= x) & (x < c[m+1])] = m*Aupp*A1**(m-1)
    for i in range(m - 2, -1, -1): # the minus one is intended to run 0
        for j in range(m - i, -1, -1):
            for n in range(m - j + 1, 2*m - i - j + 1):# the plus one is intended to run m-j + 1 to 2*m - i - j
                temp = Amijn(A, m, i, j, n)*(c[n] - x)**(m - i - 1)*theta(c[n] - x, 1)*theta(x - c[m-j], 0)
                out = out + temp
    return out*(1 - A1**m)/sum(out)

def hfit(x, m, A, c):
    A1, Alow, Aupp = A
    theta = np.heaviside
    out = np.zeros(x.size)
    dx = x[1] - x[0]
    norm = A1**m + m*Alow*A1**(m-1)*(c[m] - c[m-1]) + m*Aupp*A1**(m-1)*(c[m+1] - c[m])
    out = m*Alow*A1**(m-1)*theta(c[m] - x, 0)*theta(x - c[m-1], 1)
    out = out + m*Aupp*A1**(m-1)*theta(c[m+1] - x, 0)*theta(x - c[m], 1)
    for i in range(0, m - 2 + 1): # the plus one is intended to run m - 2
        for j in range(0, m - i + 1): # the plus one is intended to run m - i
            for n in range(1, m - i + 1):# the plus one is intended to run 1 to m - i
                temp = Amijn(A, m, i, j, n)*(c[n+m-j] - x)**(m-i-1)*theta(c[n+m-j] - x, 1)*theta(x - c[m-j], 0.5)
                out = out + temp
                norm = norm + Amijn(A, m, i, j, n)/(m-i)*(c[n+m-j] - c[m-j])**(m-i)
    return out*dx/norm


def hnorm(m, A, c):
    A1, Alow, Aupp = A
    norm = A1**m + m*Alow*A1**(m-1)*(c[m] - c[m-1]) + m*Aupp*A1**(m-1)*(c[m+1] - c[m])
    for i in range(0, m - 2 + 1): # the plus one is intended to run m - 2
        for j in range(0, m - i + 1): # the plus one is intended to run m - i
            for n in range(1, m - i + 1):# the plus one is intended to run 1 to m - i
                norm = norm + Amijn(A, m, i, j, n)/(m-i)*(c[n+m-j] - c[m-j])**(m-i)
    return norm

def goodfunform(x, s, c0, c1, c2, c3, c4, c5, n, L, Alow, Aupp, lam):
    numpeaks = 10
    dx = x[1] - x[0]
    cent = [c0, c1, c2, c3, c4, c5] + list(np.arange(1, 2*numpeaks) + c5)
    A1 = 1 - Alow - Aupp
    norm = [hnorm(m = value, A = np.array([A1, Alow, Aupp]), c = cent) for value in range(1, numpeaks)]
    amp = np.array([np.exp(-lam)*lam**(value)/fac(value) for value in range(numpeaks)])
    yvalues = [hfit(x, value, np.array([A1, Alow, Aupp]), cent) for value in range(1, numpeaks)]
    peakform = np.sum([yvalues[value - 1]*amp[value] for value in range(1, numpeaks)], axis = 0)
    peakform = np.convolve(gaussian(np.arange(-2, 2, dx), 0, s), peakform, 'same')
    peakform = peakform + np.sum([gaussian(x, cent[value - 1], s)*amp[value]*A1**value/norm[value - 1]
                                  for value in range(1, numpeaks)], axis = 0)
    return L*(peakform + amp[0]*gzero(x, c0, s, n))

def funform(x, s, c0, c1, c2, c3, c4, c5, n, L, Lbulk, Lsurf, Alow, Aupp, lam):
    numpeaks = 5
    cent = [c0, c1, c2, c3, c4, c5] + list(np.arange(1, 2*numpeaks) + c5)
    dx = x[1] - x[0]
    A1 = 1 - Alow - Aupp
    L0 = 1 - Lbulk - Lsurf
    theta = np.heaviside
    adx = np.sum((cent[0] <= x) & (x < cent[1]))
    
    # Amplitudes of the Gaussian
    amp = [np.exp(-lam)*lam**(value)/fac(value) for value in range(numpeaks)]
    L = L/sum(amp) # removes artifact of missing counts by redistributing evenly over all other energies
    
    # Events that 'pile-up' with a 0 event:
    yvalues = [hfit(x, value, np.array([A1, Alow, Aupp]), cent) for value in range(1, numpeaks)]
    peakform = L0*np.sum([yvalues[value - 1]*amp[value] for value in range(1, numpeaks)], axis = 0)
    
    #one = np.sum(peakform, axis = 0)/L0
    #print('One', np.sum(peakform, axis = 0), one*L0)
    
    # Events that 'Pile-up' with a surface event:
    yvalues = [hfit(x - cent[0], value, np.array([A1, Alow, Aupp]), cent) for value in range(1, numpeaks)]
    peakform = peakform + Lsurf*np.sum([yvalues[value - 1]*amp[value] for value in range(1, numpeaks)], axis = 0)
    
    #print('Two', np.sum(peakform, axis = 0), one*(L0 + Lsurf))
    
    # Events that 'Pile-up' with a bulk event:
    yvalues = np.zeros(x.size)
    for value in range(1, numpeaks):
        temp = np.zeros(x.size)
        temp[(cent[value - 1] <= x) & (x < cent[value])] = 1
        yvalues = yvalues + amp[value]*temp/sum(temp)

    peakform = peakform + Lbulk*yvalues

    #print('Three', np.sum(peakform, axis = 0), one*(L0 + Lsurf + Lbulk))
    norm = [hnorm(m = value, A = np.array([A1, Alow, Aupp]), c = cent) for value in range(1, numpeaks)]
    peakform = np.convolve(gaussian(np.arange(-2, 2, dx), 0, s), peakform, 'same')
    peakform = (peakform + np.sum([gaussian(x, cent[value - 1], s)*amp[value]*A1**value/norm[value-1]*L0
                                   for value in range(1, numpeaks)], axis = 0)
                + np.sum([gaussian(x-cent[0],cent[value-1], s)*amp[value]*A1**value/norm[value-1]*Lsurf
                          for value in range(1, numpeaks)], axis = 0))
        
#print('Four', np.sum(peakform, axis = 0), (L0 + Lbulk + Lsurf)*sum(amp[1:]))
    return L*(peakform + fitfun(x, s, c0, c1, n, amp[0], Lbulk, Lsurf))
