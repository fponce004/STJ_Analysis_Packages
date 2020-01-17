"""
To be filled later
"""
import os, sys, pylab, math, warnings, datetime, time, shutil, string, gc, inspect
import warnings, fitFunctions, fitPulse, fitVoigt, glob

import numpy as np
import scipy as sp
from scipy.optimize import curve_fit as cfit
from scipy.odr import ODR, Model, Data, RealData
from scipy.misc import factorial as fac
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import scipy.ndimage as nd
import scipy.signal as sig
#from Stanford_Mod_VX_X import *

class plotFile(dict):
    def __init__(self, files = None, paths = None):
        self.files = {}
        self.paths = {}
        if (files != None and paths != None) and (
            not isinstance(paths, str) and len(paths) != len(files)):
            print('Input error, paths must be either singular' +
                  'or array/list of equal size to files')
            raise
        if files != None:
            if isinstance(files, (list, np.ndarray)):
                for L, line in enumerate(files):
                    self.files['File_{}'.format(L)] = line
            elif isinstance(files, dict):
                self.files = files
            else:
                print('Input Error for variable "files".')
                raise
            if paths != None and (isinstance(paths, str) or len(paths) == 1):
                for line in enumerate(files):
                    self.paths['File_{}'.format(L)] = path
            else:
                for L, line in enumerate(paths):
                    self.paths['File_{}'.format(L)] = line

class plotdata(dict):
    def __init__(self, file, path):
        self.file = file
        self.path = path
        self.abscissa = {} # This is the x-axis
        self.ordinate = {} # This is the y-axis
    def keys(self):
        return sorted(list(self.abscissa.keys()))
    def xkeys(self):
        temp = sorted(list(self.abscissa.items()), key = lambda line: line[0])
        return reindex(temp)[1]
    def ykeys(self):
        temp = sorted(list(self.ordinate.items()), key = lambda line: line[0])
        return reindex(temp)[1]
    def plot(self, xaxis, yaxis, file = None):
        pass
        
class dataFile(dict):
    """
    For later updates change 'object' above to 'dict' in order to have this class
    inherit the properties of the dict class.
    """
    ######################################################################
    # Sets up the object with either a user defined list of titles, self
    # generated titles or given dictionary.
    ######################################################################
    def __init__(self, N = 0, path = None, file = None):
        if isinstance(N, (list, np.ndarray)):
            dictionary = {}
            for line in N:
                dictionary[line] = None
        elif isinstance(N, int):
            dictionary = {}
            if N != 0:
                for n in range(N):
                    dictionary['entry' + str(n)] = None
        elif isinstance(N, dict):
            dictionary = N
        else:
            print('Input error for class object construction.\n')
            raise
        self.dictionary = dictionary
        self.path = path
        self.file = file

    ######################################################################
    # Adds two entries in the current dictionary together and creates a
    # new entry in the dictionary with their combined titles
    ######################################################################
    def add(self, keyA, keyB, Newkey = None):
        valA = keyA if isinstance(keyA, (int, float)) else np.asarray(self.get(keyA))
        valB = keyB if isinstance(keyB, (int, float)) else np.asarray(self.get(keyA))
        Newkey = Newkey if Newkey else '({A} + {B})'.format(A = keyA, B = keyB)
        self.dictionary[Newkey] = valA + valB

    ######################################################################
    # 
    ######################################################################
    def append(self, key, value):
        if not any(key == el for el in self.keys()):
            self.set(key, value)
        elif self.get(key) == None:
            self.set(key, value)
        else:
            if isinstance(self.get(key), list):
                temp = self.get(key)
            elif (isinstance(self.get(key), int) or isinstance(self.get(key), float)
                 or isinstance(self.get(key), str)):
                temp = [self.get(key)]
            elif not (isinstance(self.get(key), int) or isinstance(self.get(key), float)
                 or isinstance(self.get(key), str)):
                temp = self.get(key)
                temp = temp.tolist()
            temp.append(value)
            self.set(key, temp)
        
    ######################################################################
    # Peforms a fit to laser calibration data
    ######################################################################
    def calib(self, keyX, keyY, newKey = None, Auto = False, peaks = None,
                   LLD = None, ULD = None, sigma = None, Energy = None,
                   Plot = True):
        global npeaks
        if Auto:
            Plot = False
        if not newKey:
            newKey = 'Calibrated'
        if not LLD:
            LLD = 0
        if not ULD:
            ULD = len(self.get(keyX))
        if not Energy:
            Energy = 3.5
        if not sigma and Auto:
            sigma = 10
        elif not sigma and not Auto:
            while True:
                try:
                    sigma = float(input('Please enter a value for sigma.\n'))
                    break
                except(ValueError):
                    print('Please enter an integar or float.\n')
        if not peaks:
            peaks = findPeak(self.get(keyX), self.get(keyY))
            if not Auto:
                print('Located peaks\n', peaks)
                vertical(peaks, self.get(keyY))
                self.plot(keyX, keyY)
                plt.show()
                peaks = str2num(input('Enter peaks to use.\n'))
        npeaks = len(peaks)
        conv = np.polyfit(peaks, range(1, len(peaks) + 1), 1)
        A = [self.get(keyY)[peak] for peak in peaks]
        mu = np.polyval(conv, sum([A[p]*peak/sum(A) for p, peak in enumerate(peaks)]))
        try:
            pos, poserror = cfit(compFun04, np.asarray([np.polyval(conv, x) for x in peaks]),
                                 np.asarray(A), p0 = [mu, max(A)],
                                 sigma = [np.sqrt(el) for el in A])
            out, outerror = cfit(compFun01, self.get(keyX)[LLD:ULD], self.get(keyY)[LLD:ULD],
                                 p0 = [1/conv[0], -conv[1]/conv[0], 10, pos[1], pos[0]],
                                 maxfev = 10000)
            self.set(newKey, [np.polyval([Energy/out[0], -Energy*out[1]/out[0]], x)
                                for x in self.get(keyX)])            
            if Plot:
                plt.plot(self.get(keyX), self.get(keyY), linewidth = 5, label = 'Data')
                plt.plot(self.get(keyX),
                         compFun01(np.asarray(self.get(keyX)), *out),
                         linewidth = 3, label = 'Fit')
                plt.plot(self.get(keyX),
                         compFun04(np.polyval(conv, np.asarray(self.get(keyX))),
                                              out[4], out[3]),
                         linewidth = 3, label = 'Poisson')
                plt.legend()
                plt.show()
            if not Auto:
                print('M_cent (Ch/Unit), B_off (Ch), sigma (Ch), amp, mu (Ch)\n', out)
            return out.tolist(), [row[r] for r, row in enumerate(outerror)]
        except(RuntimeError):
            print('Could not fit data.\n', sys.exc_info())
        
        
    ######################################################################
    # Convert one element to another using a polynomial, p = [pn,...p0]
    # for pn*x**n + ... + p0
    ######################################################################
    def convert(self, key, p, newKey = None, Offset = False):
        if not Offset:
            p[-1] = 0
        if not newKey:
            newKey = key + ' (eV)'
        self.set(newKey, [sigFigN(el, 4) for el in np.polyval(p, self.get(key))])

    ######################################################################
    # Calculate the conversion from Ch to eV output is [m, b], [dm, db]
    ######################################################################    
    def conversion(self, keyX, keyY, keyW = None, Newkey = None):
        if any(keyW == el for el in self.keys()):
            weights = self.get(keyW)
        elif isinstance(keyW, list):
            weights = keyW
        else:
            weights = [ 1 for el in self.get(keyY)]
        line, lcov = cfit(linear, np.asarray(self.get(keyX)),
                          np.asarray(self.get(keyY)), sigma = weights)
        if len(self.get(keyY)) == 2:
            lcov = np.zeros((2, 2))
        self.set('p1 (eV/Ch)', [sigFigN(line[0])]*len(self.get(keyY)))
        self.set('p0 (eV)', [sigFigN(line[1])]*len(self.get(keyY)))
        self.set('dp1 (eV/Ch)', [sigFigN(lcov[0][0])]*len(self.get(keyY)))
        self.set('dp0 (eV)', [sigFigN(lcov[1][1])]*len(self.get(keyY)))
        return [sigFigN(el) for el in line], [sigFigN(lcov[0][0]), sigFigN(lcov[1][1])]

    ######################################################################
    # Divides two entries in the current dictionary together and creates a
    # new entry in the dictionary with their combined titles
    ######################################################################
    def div(self, keyA, keyB, Newkey = None):
        valA = keyA if isinstance(keyA, (int, float)) else np.asarray(self.get(keyA))
        valB = keyB if isinstance(keyB, (int, float)) else np.asarray(self.get(keyA))
        Newkey = Newkey if Newkey else '({A} / {B})'.format(A = keyA, B = keyB)
        self.dictionary[Newkey] = valA/valB

    ######################################################################
    # Display a range of rows in the data.
    ######################################################################
    def display(self, rows = None):
        if all(isinstance(self.get(key), list) for key in self.keys()):
            entire = reindex([[key] + self.get(key) for key in self.keys()])
            if rows == None and len(entire[0]) >= 2:
                rows = [0, 2]
            for l, line in enumerate(entire[rows[0]:rows[1]]):
                if l == 0:
                    print('', '\t', line)
                else:
                    print(l - 1, '\t', line)
        else:
            print('Not all entries in object are lists.')

    ######################################################################
    # Delete a row of data
    ######################################################################
    def delrow(self, row = None, Auto = False):
        if row == None:
            print('Please determine the row to permanently delete first.\n')
        else:
            entire = reindex([[key] + self.get(key) for key in self.keys()])
            if not Auto:
                print(entire[0])
                print(entire[row + 1])
                delete = input('Confirm that the above is the row to remove.\n')
            else:
                delete = 'y'
            if delete == 'y':
                print(len(entire))
                del entire[row + 1]
                print(len(entire))
                for column in reindex(entire):
                    self.set(column[0], column[1:])

    ######################################################################
    # Extract all the columns from the file as a list of lists
    ######################################################################    
    def extract(self, Columns = None):
        if Columns:
            try:
                out = [[key] + list(self.get(key)) for key in Columns]
            except(KeyError):
                out = [[key] + list(self.get(key)) for key in self.keys()]
        else:
            out = [[key] + list(self.get(key)) for key in self.keys()]
        return out
            
    ######################################################################
    # Retrieves information stored in key
    ######################################################################
    def get(self, key):
        return self.dictionary[key]

    ######################################################################
    # Retrieves keys (keywords) for this dictionary as a list
    ######################################################################
    def keys(self):
        return sorted(list(self.dictionary.keys()))

    ######################################################################
    # Merge all the loaded files into a single file and save it
    ######################################################################
    def merge(self, fileName = None, Path = None, counter = True,
              remove = False, rename = True):
        if type(self) != dataFile:
            print('The data type: {} \nis not currently coded to merge!'.format(type(self)))
            raise
        folder = self.keys()
        for f, file in enumerate(folder):
            titles = self.get(file).keys()
            for col in titles:
                if not counter and rename:
                    self.get(file).rename(col, '{}_{}'.format(file.replace('.txt', ''), col))
                elif rename:
                    self.get(file).rename(col, '{}_{}'.format(index_gen(f), col))
        tofile = self.get(folder[0]).extract()
        for file in folder[1:]:
            tofile = tofile + self.get(file).extract()
        fileName = fileName if fileName else self.keys()[0].replace('.', '_merge.')
        fileName = fileName if fileName.endswith('.txt') else fileName + '.txt'
        Path = Path if Path else self.path
        save(data2str(reindex(tofile, Save = True)), fileName, Path, Verbose = True)
        if remove:
            for file in folder:
                os.remove(os.path.join(Path, file))

    ######################################################################
    # Multiplies two entries in the current dictionary together and
    # creates a new entry in the dictionary with their combined titles
    ######################################################################
    def mul(self, keyA, keyB, Newkey = None):
        valA = keyA if isinstance(keyA, (int, float)) else np.asarray(self.get(keyA))
        valB = keyB if isinstance(keyB, (int, float)) else np.asarray(self.get(keyA))
        Newkey = Newkey if Newkey else '({A} * {B})'.format(A = keyA, B = keyB)
        self.dictionary[Newkey] = valA * valB

    ######################################################################
    # Calculates the centroid and sigma of the input data
    ######################################################################
    def peakInfo(self, keyX, keyY, keyY_p = None, keyY_pp = None, units = None, Auto = True):
        if units == None:
               units = ' (Ch)'
        if not any(keyY_p == key for key in self.keys()):
            keyY_p = keyY + '_p' 
            self.set(keyY_p, firstDer(self.get(keyX), self.get(keyY)))
        else:
            keyY_p = keyY + '_p' 
        if not any(keyY + '_pp' == key for key in self.keys()):
            keyY_pp = keyY + '_pp'
            self.set(keyY_pp, secondDer(self.get(keyX), self.get(keyY)))
        else:
            keyY_pp = keyY + '_pp'
        if not any('PeakInfo' == key for key in self.keys()):
            self.set('PeakInfo', dataFile())
        if not any('peaks' == key for key in self.get('PeakInfo').keys()) or Auto:
            self.get('PeakInfo').set('peaks', findPeak(self.get(keyY), self.get(keyY)))        
        if not Auto:
            fit = analyze(self.get(keyX), self.get(keyY),
                          self.get('PeakInfo').get('peaks'))
            columns = ['peaks', 'C (Ch)', 'dC (Ch)',
                       'S (Ch)', 'dS (Ch)', 'b (Ch)', 'A']
            for c, col in enumerate(columns):
                self.get('PeakInfo').set(col, fit[c])
        elif Auto:
            for p, peak in enumerate(self.get('PeakInfo').get('peaks')):    
                try:
                    while True:
                        temp, vari = dataCurveFit(self.get(keyX), self.get(keyY),
                                                  param = [peak, 1, 1, 1], Auto = Auto,
                                                  Plot = not Auto)
                except(UnboundLocalError, TypeError):
                    del self.get('PeakInfo').get('peaks')[p]
                    print('Could not fit ' + str(peak) + ' at .peakInfo()\n', sys.exc_info())
                if temp != [1, 1, 1, 1]:
                    self.get('PeakInfo').append('C' + units, sigFigN(temp[0]))
                    self.get('PeakInfo').append('dC' + units, sigFigN(vari[0][0]))
                    self.get('PeakInfo').append('S' + units, sigFigN(temp[1]))
                    self.get('PeakInfo').append('dS' + units, sigFigN(vari[1][1]))
                    self.get('PeakInfo').append('b', sigFigN(temp[2]))
                    self.get('PeakInfo').append('A', sigFigN(temp[3]))
                else:
                    del self.get('PeakInfo').get('peaks')[p]
                    print('Could not fit to ' + str(peak) + ' peak') 

    #######################################################################
    # Plot columns in this dictionary
    ######################################################################
    def plot(self, keyX, keyY, keyErrorxy = [None, None], Symbol = '-',
             plotkw = {}, Hold = True, Plot = True, Type = 'Linear'):
        if Type == 'Linear':
            plt.plot(self.dictionary[keyX], self.dictionary[keyY], Symbol,**plotkw)
        elif Type =='semilogy':
            plt.semilogy(self.dictionary[keyX], self.dictionary[keyY], Symbol, **plotkw)
        elif Type =='semilogx':
            plt.semilogx(self.dictionary[keyX], self.dictionary[keyY], Symbol, **plotkw)
        elif Type =='loglog':
            plt.loglog(self.dictionary[keyX], self.dictionary[keyY], Symbol, **plotkw)
        elif Type =='errorbar':
            if keyErrorxy[0]:
                xError = self.get[keyErrorxy[0]]
            if keyErrorxy[1]:
                yError = self.get[keyErrorxy[1]]
            plt.errorbar(self.dictionary[keyX], self.dictionary[keyY],
                         fmt = Symbol, yerr = yError, xerr = xError, **plotkw)
        plt.xlabel(keyX)
        plt.ylabel(keyY)
        if any('label' == key for key in plotkw):
            plt.legend()
        if Plot == True:
            plt.show(block = Hold)

    ######################################################################
    # Perform a polynomial fit on a subset of the data
    ######################################################################    
    def poly(self, keyX, keyY, Order = 1, Newkey = None):
        self.plot(keyX, keyY, Symbol = 'o')
        mask = [[0, len(self.get(keyY))]]
        while True:
            while True:
                change = mask
                while True:
                    print('Current mask:', mask)
                    action = input('Enter i[#,#] to change interval, m[#,#] to add a mask, c to clear all masks or e to exit.\n')
                    if action == 'e':
                        mask = change
                        break
                    elif action[0] == 'c':
                            mask = [mask[0]]  
                    elif action[0] == 'i':
                        try:
                            change[0] = str2num(action[1:])
                        except(ValueError):
                            print('Error at masking lines\n')
                            print("Unexpected error:", sys.exc_info())
                    elif action[0] == 'm':
                        try:
                            change.append(str2num(action[1:]))
                        except(ValueError):
                            print('Error at masking lines\n')
                            print("Unexpected error:", sys.exc_info())
                break
            try:
                for line in mask:
                    if line == mask[0]:
                        exclude = [element for e, element in enumerate(self.get(keyX))
                                   if element < line[0] or element > line[1]]
                        xval, rem = removepts(self.get(keyX)[:], ret = True, spe = exclude)
                        yval = removepts(self.get(keyY)[:], rem)
                    else:
                        special = [int(element) for element in
                                   np.linspace(line[0], line[1], line[1] - line[0] + 1)]
                        xval, rem = removepts(xval, ret = True, spe = special)
                        yval = removepts(yval, rem)  
                p, pcov = np.polyfit(xval, yval, Order, cov = True)
                self.plot(keyX, keyY, Symbol = 'o', Plot = False)
                plt.plot(xval, yval, 'ko')
                plt.plot(self.get(keyX), np.polyval(p, self.get(keyX)), '--b', linewidth = 3)
                plt.xlabel(keyX, size = 30)
                plt.ylabel(keyY, size = 30)
                plt.show()
                if 'e' == input('Enter e to exit if acceptable.\n'):
                    break     
            except(RuntimeError):
                print('RuntimeError', sys.exc_info())
            except(AttributeError):
                print('AttributeError', sys.exc_info())
            except(TypeError):
                print('TypeError', sys.exc_info())
        if Newkey == None:
            Newkey = 'Lfit '+ keyY + ' vs ' + keyX
        self.dictionary[Newkey] = np.polyval(p, self.dictionary[keyX]).tolist()

    ######################################################################
    # Rename a key
    ######################################################################
    def rename(self, oldkey, newkey):
        self.dictionary[newkey] = self.dictionary[oldkey]
        if newkey != oldkey:
            self.remove(oldkey)

    ######################################################################
    # Remove an existing key
    ######################################################################
    def remove(self, key):
        del self.dictionary[key]

    ##################################################
    # Save data in class to ascii format
    ##################################################
    def save(self, fileName = None, Path = None, Columns = None, Verbose = True):
        keys = Columns if Columns else self.keys()
        Path = Path if Path else self.path
        fileName = fileName if fileName else self.file
        data = self.extract(keys)
        save(data2str(reindex(data, Save = True)), fileName, Path, Verbose)

    ######################################################################
    # Set the entry for a key or create new key
    ######################################################################
    def set(self, key, value):
        self.dictionary[key] = value

    ##################################################
    # Subtracts two entries in the current dictionary
    # together and creates a new entry in the
    # dictionary with their combined titles
    ##################################################
    def sub(self, keyA, keyB, Newkey = None):
        valA = keyA if isinstance(keyA, (int, float)) else np.asarray(self.get(keyA))
        valB = keyB if isinstance(keyB, (int, float)) else np.asarray(self.get(keyA))
        Newkey = Newkey if Newkey else '({A} - {B})'.format(A = keyA, B = keyB)
        self.dictionary[Newkey] = valA - valB

    ######################################################################
    # Converts to np array
    ######################################################################
    def toarray(self, key):
        self.set(key, np.asarray(self.get(key)))
                
    ######################################################################
    # Converts to list
    ######################################################################
    def tolist(self, key):
        self.set(key, self.get(key).tolist())

################################################################################
################################################################################
################################################################################
class poly(object):
    r"""
    This class is intended to generate and evaluate polynomials identical to
    the np.polyval function:
    y = p0*x**n + p1*x**(n-1) + ... + pn
    """
    
    def __init__(self, dataX = None, dataY = None, dataW = None, p = None,
                 order = 1, error = False, fit = True):
        if p == None and (dataX == None or dataY == None):
            print('Please enter data to fit or coeff')
            raise
        if fit:
            dataX = array(dataX)
            dataY = array(dataY)
            dataW = array(dataW)
            p = np.polyfit(dataX, dataY, order, w = dataW, cov = error)
        if error:
            self.coeff = p[0]
            self.error = np.sqrt(np.diag(p[1]))
        else:
            self.val = array(p)
            self.err = np.zeros(len(p))

    def eval(self, dataX):
        return np.polyval(self.coeff, array(dataX))

class gaussian(object):
    r"""
    This class is intended to generate and evaluate the normalized gaussian function:
    y = area/sqrt(2*pi*sig**2)*exp(-(x-center)**2/(2*sig**2))
    with
    p = [center, sig, area]
    """
    
    def __init__(self, dataX = None, dataY = None, dataW = None, p = None,
                 order = 1, error = False, fit = True):
        if fit:
            dataX = array(dataX)
            dataY = array(dataY)
            dataW = array(dataW)
            center = sum(dataX*dataY)/sum(dataY)
            sigma = np.sqrt(sum(dataY*(dataX - center)**2)/sum(dataY))
            area = sum((dataY[:-1] + dataY[1:])*(dataX[1:] - dataX[:-1])/2)
            p, pcov = cfit(self.eval, dataX, dataY, p0 = [center, sigma, area],
                           sigma = dataW, absolute_sigma = dataW.size)
            self.val = p
            self.err = np.sqrt(np.diag(pcov))

    def eval(self, dataX, center = 0, sigma = 0, area = 0):
        p = np.asarray([center, sigma, area])
        if not all(p):
            p = self.val
        return p[2]/np.sqrt(2*np.pi*p[1]**2)*np.exp(-(x-p[0])**2/(2*p[1]**2))

################################################################################
# Begin none class definitions
################################################################################
        ########################################################################
######### Perform Gaussian fit to dataY with dataX
        ########################################################################
def alpha_ij(i_0, i_1, j_0, j_1):
    """
    This function is the matrix generator for the calibaration fit
    """
    return round((max(min(i_1, j_1) - max(i_0, j_0), 0))/(i_1 - i_0),5)

def analyze(dataX, dataY, peaks = np.asarray([]), **keys):
    r"""
    This program is intended as a 'robust' peak fitter.
    Parameters:
    -----------
    dataX : array or list with the x component of the data
    dataY : array or list with the y component of the data
    peaks : <optional> array or list with the location of the peaks to be fitted
            They will be calculated from the data if not given.
    Returns:
    --------
    output : list of np.ndarrays with the values of the best fit to the peaks with
             the input peak followed by the fitting parameters alternating with their
             total error added in quadrature then the standard error; for compFun03:
             peak, C, dCtot, dC, S, dStot, dS, b, dbtot, db, A, dAtot, dA
             for gaussian fits with constant offset note here that the A is area not
             amplitude.
    keywords:
    Auto : Boolean,  default = False changes this program to run with self calculated
           param for the given data, currently untested
    findPeak : dictionary {'baseNoise':###} for setting the basenoise in findPeak, default
               is None
    function : handle for the function the data will be fit to, the default is compFun03
               a gaussian with a constant offset
    """

    if keys:
        try:
            Auto = keys['Auto']
        except(KeyError):
            Auto = False
        try:
            bN = keys['findPeak']['base']
        except(KeyError):
            bN = 1
        try:
            per = keys['findPeak']['per']
        except(KeyError):
            per = 5
        try:
            function = keys['function']
        except(KeyError):
            function = compFun03
    else:
        Auto = False
        bN = 1
        per = 5
        function = compFun03

    fit_param = inspect.getargspec(function)[0][1:]
    if not isinstance(dataX, np.ndarray):
        dataX = np.asarray(dataX)
    if not isinstance(dataY, np.ndarray):
        dataY = np.asarray(dataY)
    if not isinstance(peaks, np.ndarray):
        peaks = np.asarray(peaks)
    
    #C, dC, S, dS, b, A = [],[],[],[],[],[]
    output = []
    ##################################################
    # Allows user to specify which peaks to use and
    # enter new ones if necessary.
    ##################################################
    Yp =firstDer(dataX, dataY)
    Ypp = secondDer(dataX, dataY)
    
    if (not Auto and not peaks.size):
        hold = findPeak(dataX, dataY)#, base = bN, per = per)
        while True:
            vertical(hold, dataY)
            plt.plot(dataX, dataY)
            print(hold)
            plt.show()
            action = input('Enter new list of peaks or e to exit.\n')
            if action == 'e':
                peaks = np.asarray(hold)
                del hold
                break
            else:
                try:
                    hold = str2num(action)
                except(ValueError):
                    print('Error, please enter list or e.\n', sys.exc_info())
    ##################################################
    # Fit to peaks
    ##################################################
    for peak in peaks:
        param = [peak, 1, 1, dataY[int(peak)]]
        mask = [findInterval(peak, firDer = Yp, secDer = Ypp)]
        hold = []
        while True:
            param, psig = dataCurveFit(dataX, dataY, param,
                                       mask, dfspe = 0, Plot = True, Fun = function)
            hold.append(param)
            print('Initial guess:\n' +
                  '[{}]'.format(('{}, '*len(fit_param))[:-2]).format(*fit_param))
            print(param)
            print('Current mask:\n[[interval], [mask 1], ...]')
            print(mask)
            while True:
                if function == compFun03 and param[1] < 0:
                    param[1] = -param[1]
                    hold.pop()
                    break
                else:
                    change = input('Enter i[#, #] for new interval, ' +
                                   'p[#, #] to change initial parameters, ' +
                                   'm[#, #] to mask data, ' +
                                   'mc[#] to clear specific mask, ' +
                                   'mc to clear all masks, ' +
                                   'rc[#] to specific systematic error value, ' +
                                   'rc to clear all systematic error values, ' +
                                   'dis to display all current parameters, ' + 
                                   'or e to exit.\n')
                if change == 'e':
                    break
                elif change == 'dis':
                    print('Initial guess:\n[{}]' +
                          ''.format(('{},'*len(fit_param))[:-2]).format(*fit_param))
                    print(param)
                    print('Current mask:\n[[interval], [mask 1], ...]')
                    print(mask)
                    print('First column of the values to determine systematic error')
                    dis(reindex(hold)[0])
                elif change.startswith('rc['):
                    try:
                        del hold[int(change[3:-1])]
                    except(ValueError, IndexError):
                        print('Error clearing systematic error entry')
                        print("Unexpected error:", sys.exc_info())
                elif change == 'rc':
                    hold = []
                elif change.startswith('m['):
                    try:
                        mask.append(str2num(change[1:]))
                        if len(mask[-1]) !=2:
                            mask.pop()
                        break
                    except(ValueError, IndexError):
                        print('Error at masking lines')
                        print("Unexpected error:", sys.exc_info())
                elif change.startswith('mc['):
                    try:
                        del mask[int(change[3:-1])]
                        break
                    except(ValueError, IndexError):
                        print('Error clearing mask entry')
                        print("Unexpected error:", sys.exc_info())                        
                elif change == 'mc':
                    mask = [mask[0]]
                    break
                elif change.startswith('i['):
                    try:
                        mask[0] = str2num(change[1:])
                        break
                    except(TypeError, ValueError):
                        print('Error entering new interval')
                        print('Unexpected error:', sys.exc_info())
                elif change.startswith('p['):
                    try:
                        param = str2num(change[1:])
                        break
                    except(TypeError, ValueError):
                        print('Error entering new parameters')
                        print('Unexpected error:', sys.exc_info())
            if change == 'e':
                #hold = np.asarray([np.std(line) for line in reindex(hold)])
                hold = 0 # This removes any consideration of systematic error on range
                ptot = np.sqrt(psig**2 + hold**2)
                #output.append([peak] +
                #              list(np.asarray(list(zip(param, ptot, psig))).flatten()))
                output.append([peak] +
                              list(np.asarray(list(zip(param, psig))).flatten()))
                break
    return reindex(output)

def array(data):
    """
    This program checks if the input data is in an array format and if not converts it
    Parameters:
    -----------
    data : List or array type

    Return:
    -------
    data : array type
    """
    if not isinstance(data, (list, np.ndarray)):
        print('Input type is {} this takes list or array types only.'.format(type(data)))
        raise
    return np.asarray(data)

        ########################################################################
######### Calculate 'base' noise
        ########################################################################
def basenoise(dataY, interval = None, time = None):
    if interval == None:
        interval = [0, len(dataY)]
    return sum(dataY[interval[0]: interval[1]])/(interval[1] - interval[0])

def bin_line(E, p0, p1):
    """
    This function takes an energy value E and converts it to the corresponding
    bin based on p1, p0.
    [p0] : [eV/ch]
    [p1] : [eV]
    """
    return E/p0 - p1/p0

def bin_quad(E, q0, q1, q2):
    """
    This function takes an energy value E and converts it to the corresponding
    bin based on q2, q1, q0.
    [q0] : [eV/ch**2]
    [q1] : [eV/ch]
    [q2] : [eV]
    """
    return (np.sqrt(q1**2 + 4*(E - q2)*q0) - q1)/(2*q0)

def boxfilter(paths):
    """
    This function will load the given data and process using a box filter
    currently this very basic.
    """
    # Load file names in folder
    folder = sorted(glob.glob(os.path.join(paths, '*')),
                    key = lambda line: int(line.split('_')[-1][:-4]))

    # Load data from folder and set basic run parameters
    host = getChannels(folder)
    nsamples = host['prop']['samples'][0][0][0][0]
    pretrig = host['prop']['pretrig'][0][0][0][0]
    
    # Process the 'voltage' unit data ai0 and ai1 with a box filter
    gap = 25    # Must be odd
    length = 55
    output = [np.array([pulse_filter(host['ai0'][val], gap, length)
                        for val in range(len(host['ai0']))]),
              np.array([pulse_filter(host['ai1'][val], gap, length)
                        for val in range(len(host['ai1']))])]

    # Calculate ranges to determine pulse properties
    out = np.array([sum(line)/len(line) for line in output])
    exloc = [[np.where(max(line) == line)[0][0], np.where(min(line) == line)[0][0]]
             for line in out]

    # Set left side
    ltag = [(0.2*line[exloc[l][0]] < line) &
            (line < 0.8*line[exloc[l][0]]) &
            (np.arange(nsamples) < exloc[l][0]) for l, line in enumerate(out)]
    lpol = [np.polyfit(np.arange(nsamples)[ltag[l]], line[ltag[l]], 1)
            for l, line in enumerate(out)]
    left = np.array([int(-line[1]/line[0]) for line in lpol]).astype(int)

    # Set midpoint
    mtag = [(0.8*line[exloc[l][1]] < line) &
            (line < 0.8*line[exloc[l][0]]) &
            (exloc[l][0] < np.arange(nsamples)) &
            (np.arange(nsamples) < exloc[l][1]) for l, line in enumerate(out)]

    mpol = [np.polyfit(np.arange(nsamples)[mtag[l]], line[mtag[l]], 1)
            for l, line in enumerate(out)]
    middle = np.array([int(-line[1]/line[0]) for line in mpol]).astype(int)

    # Set right side
    rtag = [(0.10*line[exloc[l][1]] < line) &
            (line < 0.05*line[exloc[l][1]]) &
            (exloc[l][1] < np.arange(nsamples)) for l, line in enumerate(out)]
    rpol = [np.polyfit(np.arange(nsamples)[rtag[l]], line[rtag[l]], 1)
            for l, line in enumerate(out)]
    right = np.array([int(-line[1]/line[0]) for line in rpol]).astype(int)

    # Calculate a response for each event
    forhistogram = np.array([np.array([sum(line[left[val]:middle[val]]) -
                                       sum(line[middle[val]:right[val]])
                                       for line in output[val]])
                             for val in range(len(output))])
    counts, end = np.histogram(sum(forhistogram)/np.sqrt(forhistogram.size),
                               bins = 2000, range = (-1, 1))
    center = (end[:-1] + end[1:])/2
    
    # Save results
    data = dataFile()
    data.file = '_'.join(folder[0].split('/')[-1].split('_')[:-1]) + '_Histograms.txt'
    data.path = '/'.join(folder[0].split('/')[:4]) + '/Analysis/'
    data.set('000_Bin_center', np.round(center, 5).astype(str))
    data.set('001_Joint_hist', counts)
    data.set('002_ai0_hist', np.histogram(forhistogram[0], bins = 1000, range = (-0.5, 0.5))[0])
    data.set('003_ai1_hist', np.histogram(forhistogram[1], bins = 1000, range = (-0.5, 0.5))[0])
    data.set('004_ai0_area', forhistogram[0])
    data.set('005_ai1_area', forhistogram[1])
    data.save()
    return data

def bruteforce(xval, yval, xerr = None, yerr = None, fun = np.polyval, p0 = [1, 1]):
    output = np.array([])
    outerr = np.array([])
    model = Model(fun)
    if xerr != None:
        for b in range(2*xval.size + 1):
            ytemp = yval
            corr = np.zeros(xval.size)
            if b < xval.size:
                corr[b] = 1
            elif b == 2*xval.size:
                pass
            elif b > xval.size:
                corr[b%xval.size] = -1
            data = RealData(xval, yval + corr*yerr, sx = xerr, sy = yerr)
            pro = ODR(data, model, beta0 = p0)
            pro.set_job(fit_type = 0)
            pro_out = pro.run()
            if b == 2*xval.size:
                output = np.insert(output, 0, pro_out.beta)
                outerr = np.insert(outerr, 0, pro_out.sd_beta)
            else:
                output = np.insert(output, output.size, pro_out.beta)
                outerr = np.insert(outerr, outerr.size, pro_out.sd_beta)
    if yerr != None:
        for b in range(2*yval.size + 1):
            ytemp = yval
            corr = np.zeros(yval.size)
            if b < yval.size:
                corr[b] = 1
            elif b == 2*yval.size:
                pass
            elif b > yval.size:
                corr[b%yval.size] = -1
            data = RealData(xval, yval + corr*yerr, sx = xerr, sy = yerr)
            pro = ODR(data, model, beta0 = [1, 1])
            pro.set_job(fit_type = 0)
            pro_out = pro.run()
            if xerr != None and b != 2*yval.size:
                output = np.insert(output, output.size, pro_out.beta)
                outerr = np.insert(outerr, outerr.size, pro_out.sd_beta)
            elif xerr == None:
                output = np.insert(output, 0, pro_out.beta)
                outerr = np.insert(outerr, 0, pro_out.sd_beta)
    return output, outerr

def calib_data(pu_file = None, laser_file = None, calib_file = None, time = None, rang = (0, 100), bins = 1000):
    if not isinstance(pu_file, dataFile):
        print('Please select the Pu file you wish to calibrate.')
        pu_file = impData()
        pu_file = pu_file.get(pu_file.keys()[0])
    if not isinstance(laser_file, dataFile):
        print('Please select the corresponding Laser data.')
        laser_file = impData()
        laser_file = laser_file.get(laser_file.keys()[0])
    if not isinstance(calib_file, dataFile):
        print('Please select the corresponding Calib data.')
        calib_file = impData()
        calib_file = calib_file.get(calib_file.keys()[0])
  
    column = calib_file.keys()
    column.sort(key = lambda line: int(line.split('_')[0]))
    order = round(len(column[2:-9])/3 - 1)
    num_files = len(calib_file.get('00_File'))

    pfiles = pu_file.keys()    
    lfiles = laser_file.keys()
    instance = 0
    fit = ['Line', 'Quad']
    while True:
        path = os.path.join(pu_file.path, 'Analysis_{}_{}'.format(fit[order-1], instance))
        if not os.path.exists(path):
            os.makedirs(path)
            break
        elif os.path.exists(path):
            instance += 1
        else:
            raise
    pu_file.path = path
    pu_file.file = 'Analysis_' + pu_file.file
    laser_file.path = path
    laser_file.file = 'Analysis_' + laser_file.file

    inverse = inv_line if order == 1 else inv_quad
    num2ev = np.asarray([np.ones(num_files)/3.5**(order - p) for p in range(order + 1)])
    pval = np.asarray([np.asarray(line[1:])
                       for line in calib_file.extract()[2:-9:3]])*num2ev
    pval = np.asarray(list(zip(*pval)))
    perr = np.sqrt(np.asarray([np.asarray(line[1:])
                               for line in calib_file.extract()[3:-9:3]])**2 +
                   np.asarray([np.asarray(line[1:])
                               for line in calib_file.extract()[4:-9:3]])**2)*num2ev
    perr = np.asarray(list(zip(*perr)))
    time = time if time else int(input('Enter the acquisition time '+
                                       'of a single file in minutes.\n'))
    while True:
        runs = []
        for f, file in enumerate(pfiles[:num_files]):
            puse = (pval[f] + pval[f+1])/2
            euse = np.sqrt(perr[f]**2 + perr[f+1]**2)/2
            file_en, file_err = inverse(np.arange(1, pu_file.get(file).size + 1),
                                        0, puse, euse)
            pu_file.set('Energy (eV)_{}'.format(f), file_en)
            pu_file.set('dE (eV)_{}'.format(f), file_err)
            laser_file.set('Energy (eV)_{}'.format(f), file_en)
            laser_file.set('dE (eV)_{}'.format(f), file_err)
            hist, end = np.histogram(file_en, bins, rang, weights = pu_file.get(file))
            runs.append(hist/time)
        energy = (end[1:] + end[:-1])/2
        plt.plot(energy, sum(runs), linewidth = 5)
        plt.show()
        action = input('What now?\n')
        if action.startswith('e'):
            break
        try:
            hold = str2num(action)
            rang = hold[:2]
            bins = hold[-1]
        except:
            print('Error: {}'.format(sys.exc_info()))
    out = dataFile()
    out.path = path
    out.file = pu_file.file.replace('.txt', '_{}Bins_{}-{}eV.txt'.format(bins, *rang))
    out.set('00_Energy (eV)', energy)
    out.set('00_Intensity_Sum', sum(runs))
    for f, file in enumerate(pfiles[:num_files]):
        out.set(file, runs[f])
    
    pcol = pu_file.keys()
    pcol.sort(key = lambda line: int(line.split('_')[1]))
    pu_file.save(Columns = pcol)
    lcol = laser_file.keys()
    lcol.sort(key = lambda line: int(line.split('_')[1]))
    laser_file.save(Columns = lcol)
    #ocol = out.keys()
    #ocol.sort(key = lambda line: int(line.split('_')[1]))
    out.save()#(Columns = ocol)
    return out, pu_file, laser_file

###############################################################################
# Creates a indexed array for the x-axis
###############################################################################
def channel(dataNum):
    return np.arange(1, len(dataNum) + 1)

################################################################################
# Plot colors and symbols
################################################################################
def colors(s = 0, c = 0):
    _colors = ['b', 'g', 'r', 'm', 'c', 'y', 'k']
    _shapes = ['o', 's', 'D', '^', 'p', 'H', '8', '>', 'v', '<']
    if s < 0:
        output = (0.1 + 0.9*c/abs(s), 1 - c/abs(s), c/abs(s))
    else:
        output = _shapes[s] + _colors[c]
    return output

        ########################################################################
        # Function: Composite npeaks-gaussians with poissonian envelope
######### f(x, delta_c, coff, s, A, mu) = poisson*sum(gaussians)
        # 
        ########################################################################
def compFun01(x, centr, coff, sigma, amp, mu):
    global _npeaks, _start
    try:
        if isinstance(_npeaks, float):
            _npeaks = int(_npeaks)
        elif isinstance(_npeaks, int):
            pass
        else:
            _npeaks = int(input('Enter number of peaks to use for fit.\n'))
    except(NameError):
        _npeaks = int(input('Enter number of peaks to use for fit.\n'))
    try:
        _start = int(_start)
    except(NameError):
        _start = 0
    centr, sigma, amp, mu = np.abs([centr, sigma, amp, mu])
    amp =  np.sqrt(2*np.pi*sigma**2)*amp
    gauss_cont = sum([compFun02(x, n*centr + coff, sigma, amp*poisson(n, mu))
                      for n in range(_npeaks + _start)])
    return gauss_cont

def compFun01_1(x, p0, p1, p2, sigma, amp, mu):
    global _npeaks, _start
    try:
        if isinstance(_npeaks, float):
            _npeaks = int(_npeaks)
        elif isinstance(_npeaks, int):
            pass
        else:
            _npeaks = int(input('Enter number of peaks to use for fit.\n'))
    except(NameError):
        _npeaks = int(input('Enter number of peaks to use for fit.\n'))
    try:
        _start = int(_start)
    except(NameError):
        _start = 0
    p1, sigma, amp, mu = np.abs([p1, sigma, amp, mu])
    amp =  np.sqrt(2*np.pi*sigma**2)*amp
    gauss_cont = sum([compFun02(x, p0*(n+1)**2 + p1*(n+1) + p2, sigma, amp*poisson(n, mu))
                      for n in range(_npeaks + _start)])
    return gauss_cont

        ########################################################################
        # Function: Composite gaussian function with area A
######### f(x, c, s, A) = constant(x, A)* gaussian(x, c, s)
        # 
        ########################################################################
def compFun02(x, c, s, A):
    c, s, A = np.abs(c), np.abs(s), np.abs(A)
    return A/np.sqrt(2*np.pi*s**2)*gaussian(x, c, s)

        ########################################################################
        # Function: Composite gaussian function with area A and constant
######### offset b
        # f(x, c, s, b, A) = b + A*np.exp(-(x-c)**2/(2*s**2))
        ########################################################################
def compFun03(x, c, s, b, A):
    return b + np.abs(A)/np.sqrt(2*np.pi*s**2)*gaussian(x, c, s)

def compFun03_e(x, c, s, b, A, dc, ds, db, dA):
    terms = np.exp(-(x-c)**2/(2*s**2))
    return np.sqrt(db**2 + terms**2*dA**2 + terms**2*A**2*(x - c)**2/s**4*dc**2 +
                     terms**2*A**2*(x - c)**4/s**6*ds**2)
        ########################################################################
        # Function: Composite poisson function with amplitude A
######### f(n, mu, A) = A*exp(-mu)*mu**n/n!
        # 
        ########################################################################
def compFun04(n, mu, A):
    return A*poisson(n, mu)

        ########################################################################
        # Function: Composite single gaussian, poisson function with amplitude A
######### f(n, mu, A) = A*exp(-mu)*mu**n/n!*np.exp(-(x-(c + coff))**2/(2*s**2))
        # 
        ########################################################################
def compFun05(x, centr, coff, sigma, amp, mu):
    gauss_cont = gaussian(x, centr + coff, sigma)
    poiss_cont = poisson(np.polyval([1/centr, -coff/centr], x), mu)
    #return constant(x, amp)*gauss_cont*poiss_cont
    print('Should not be using this...\n\n\n')
    return 0
        ########################################################################
        # Function: Composite single gaussian with amplitude A and linear offset
######### 
        # 
        ########################################################################
def compFun06(x, c, s, A, b, m):
    return b + m*x + A*np.exp(-(x - c)**2/(2*s**2))
def compFun06_ODR(p, x):
    c, s, A, b, m = p
    return b + m*x + A*np.exp(-(x - c)**2/(2*s**2))


        ########################################################################
        # Function: Pulse fit
######### 
        # 
        ########################################################################
def compFun07(t, t_0, t_r, t_d, A, b):
    return b + A*(1 - np.exp(-(t-t_0)/t_r))*np.exp(-(t-t_0)/t_d)*(t >= t_0)
    
        ########################################################################
        # Function: Pulse fit
######### 
        # 
        ########################################################################
def compFun08(t, t0, tr, td, A, b):
	return b + A*np.exp(-(t-t0)/td)/(1 + np.exp(-(t-t0)/tr))


def compFun09(x, c1, s1, A1, c2, s2, A2):
	return (A1/np.sqrt(2*np.pi*s1**2)*gaussian(x, c1, s1) +
                A2/np.sqrt(2*np.pi*s2**2)*gaussian(x, c2, s2))

def compFun10(x, c, s):
	global mu_10
	global amp_10
	global p_10
	return gaussian(x, c, s)*amp_10*poisson(inv_quad(x, 0, p_10,
                                                         np.zeros(3))[0], mu_10)

def compFun11(p, x):
    """
    This function takes p a set of n-gaussian parameters and computes the sum
    total of the n-gaussians.
    Parameters:
    -----------
    p : list, array of {C_i, S_i, A_i} a list of the parameters of n-gaussians
    x : list, array of x-values
    Return:
    -------
    out : list, array of the sum of the n-gaussians with the parameters set by p
    """
    cent = p[::3]
    sig = p[1::3]
    area = p[2::3]
    out = sum([compFun02(x, *line) for line in reindex([cent, sig, area])])
    return np.asarray(out)


def compFun12(x, c, s, Area, tau, Amp):
    """
    This function is intended to fit a gaussian on an exponential background
    Parameters:
    -----------
    c : float, int of centroid
    s : float, int of sigma
    Area : float, int of gaussian Area
    tau : float, int decay rate of exponential
    Amp : float, int amplitude of decay function the offset is set to zero and absorbed
    x : list, array of x-values
    Return:
    -------
    out : list, array of the gaussian with an exponential decay background
    """  
    return compFun02(x, c, s, Area) + exp_decay(x, 0, tau, Amp)
def compFun12_ODR(p, x):
    return compFun12(x, *p)

def compFun13(x, c, s, Area, taui, Ampi, tauf, Ampf):
    """
    This function is intended to fit a gaussian on an exponential background
    Parameters:
    -----------
    c : float, int of centroid
    s : float, int of sigma
    Area : float, int of gaussian Area
    taui : float, int decay rate of initial exponential
    Ampi : float, int initial amplitude of decay function the offset is set to zero and absorbed
    tauf : float, int decay rate of final exponential
    Ampf : float, int final amplitude of decay function the offset is set to zero and absorbed
    x : list, array of x-values
    Return:
    -------
    out : list, array of the gaussian with an exponential decay background
    """ 
    return compFun02(x, c, s, Area)+ exp_decay(x, 0, taui, Ampi) + exp_decay(x, 0, tauf, Ampf)
def compFun13_ODR(p, x):
    return compFun13(x, *p)

def compFun14(x, c, s, Area, Amp, off):
    return abs(Amp)/(x + off) + compFun02(x, c, s, Area)
def compFun14_ODR(p, x):
    return compFun14(x, *p)

def compFun15(x, c, s, Area, q2, q1, q0):
    return q2*x**2 + q1*x + q0 + compFun02(x, c, s, Area)
def compFun15_ODR(p, x):
    return compFun15(x, *p)


def make_compFun16_ODR(Num):
    global compFun16_ODR
    def compFun16_ODR(p, x):
        slope, offset = abs(p[0]), p[1]
        #dec, Amp = abs(p[2]), abs(p[3])
        #hold = np.asarray(exp_decay(x, 0, dec, Amp))
        #starter = int(round(np.polyval(p[:2], p[2])/3.5, 0))
        hold2 = sum([compFun02(x, slope*(Num + n) + offset, *line)
                     for n, line in enumerate(reindex([p[2::2], p[3::2]]))])
                    #for n, line in enumerate(reindex([p[4::2], p[5::2]]))])
        return hold2# + hold
def compFun16_ODR(p, x):
        slope, offset = abs(p[0]), p[1]
        #dec, Amp = abs(p[2]), abs(p[3])
        #hold = np.asarray(exp_decay(x, 0, dec, Amp))
        #starter = int(round(np.polyval(p[:2], p[2])/3.5, 0))
        hold2 = sum([compFun02(x, slope*(1+ n) + offset, *line)
                     for n, line in enumerate(reindex([p[2::2], p[3::2]]))])
                    #for n, line in enumerate(reindex([p[4::2], p[5::2]]))])
        return hold2# + hold

def compFun17(x, c, s, Area, taui, Ampi, tauf, Ampf):
    return (compFun02(x, c, s, Area) + exp_decay(x, 0, tauf, Ampf) + exp_decay(x, 0, taui, Ampi))

def compFunbg(x, p, select = 0):
    output = compFun02(x, p[1], p[2], p[0])
    if select == 0:
        return output + p[3]/p[4]*exp(-x/p[4])
    elif select == 1:
        return output + p[3]/p[4]*exp(-x/p[4]) + p[5]/p[6]*exp(-x/p[6])
    elif select == 2:
        return output + p[3] + p[4]*x
    elif select == 3:
        return output + p[3] + p[4]*x + p[5]*x**2
    elif select == 4:
        return output + p[3]/(p[4] - x)
    return 0
    

def constant(x, c):
    return c

def data2str(data, Last = False):
    r"""
    This program takes a list of data and converts it to a single tab
    deliminated string.

    Parameters:
    -----------
    data : list of lists or single list of entries to be reduced
    Last : boolean, sets the entry of data as the last entry in the entire set
           so that there is no new line at the end of the file.

    Returns:
    --------
    hold : list of strings that are tab-deliminated if a list of lists was
           initially given or a list with a single string.

    Example:
    --------
    
    """

    hold = []
    if isinstance(data, (list, np.ndarray)):
        if all(isinstance(cluster, (list, np.ndarray)) for cluster in data):
            hold = [''.join(list(np.asarray(list(zip(
                    cluster, ['\t']*(len(cluster) -1) + ['\n']))).flatten()))
                    for cluster in data]
            hold[-1] = hold[-1].replace('\n', '')
        else:
            hold = [''.join(list(np.asarray(list(zip(
                    data, ['\n']*(len(data) -1) + ['']))).flatten()))]
    else:
        if Last != False:
            hold.append(str(data) + '\n')
        else:
            hold.append(str(data))
    return hold

def dataCurveFit(dataX, dataY, param = [1, 1, 1, 1], mask = None, dfspe = 0, file = None, Plot = True, Auto = False, Fun = compFun03):
    r"""
    This program is intended to perform a fit to a given set of data
    Parameters:
    -----------

    Returns:
    --------

    """
    dataX = array(dataX)
    dataY = array(dataY)
    
    if Auto == True:
        peak = param[0]
        counter = 0
        while True:
            interval = findInterval(int(round(param[0], 0)), dataX, dataY)
            hold = quadraticFit(dataX, dataY, interval, Plot = False)
            if counter == 10:
                param = hold[0:-1:2]
                mask = [hold[-1]]
                break
            elif not all(not math.isnan(value) for value in hold[:-1]):
                hold = param
                counter +=1
            else:
                if abs(param[0] - hold[0]) < 0.01*hold[0]:
                    param = hold[0:-1:2]
                    mask = [hold[-1]]
                    break
                else:
                    param = hold[0:-1:2]
                    counter += 1
    elif mask == None:
        mask = [[0, len(dataY) -1]]
        interval = mask[0]
    else:
        interval = mask[0]
    
    delta = int(np.floor(0.1*(interval[-1]- interval[0])))
    if interval[0] - delta < 0 or interval[1] + delta > len(dataY) - 1:
        delta = 0
    gaussX = np.linspace(interval[0]-50,interval[1] + 50,
                         interval[1] - interval[0] + 500)   
    pcov = []
    #(c, s, b, A) = param
    
    ####################
    # Add offset to b 
    ####################
    if Fun != compFun03:
        param = str2num(input('Enter Parameters'))

    
    dataY = dataY + 1
    oldGaussian = Fun(gaussX, *param)
    try:
        for line in mask:
            if line == mask[0]:
                dataFit, rem = removepts(dataY[interval[0]:interval[1]],
                                         ret = True, spe = dfspe)
                xFit = removepts(dataX[interval[0]: interval[1]], rem)
            else:
                special = np.linspace(line[0], line[1],
                                      line[1] - line[0] + 1).astype(int)
                xFit, rem = removepts(xFit, ret = True, spe = special)
                dataFit = removepts(dataFit, rem)
        param, pcov = cfit(Fun, xFit, dataFit, p0 = param,
                           sigma = np.sqrt(dataFit), absolute_sigma = True)
    ####################
    # Remove offset to b
    ####################
        param[2] = param[2] - 1
        
        newGaussian = Fun(gaussX, *param)
#        newerror = compFun03_1(gaussX, *np.insert(param, len(param), np.diag(pcov)))
        if param[0] < 0:
            print('ERROR: y offset is negative!')
        if param[3] < 0:
            print('ERROR: sigma is negative. Run again!')
            param[3] = -param[3]
    except(RuntimeError):
        print('RuntimeError - dataCurveFit', sys.exc_info())
    except(AttributeError):
        print('AttributeError - dataCurveFit', sys.exc_info())
    if Plot:
        try:
            plt.plot(dataX[interval[0] - delta: interval[-1] + delta],
                 dataY[interval[0] - delta: interval[-1] + delta], 'o',
                 markerfacecolor = (1,1,1), markeredgecolor = 'g',
                 markeredgewidth = 1)
            plt.plot(gaussX, oldGaussian, '.r')
            plt.plot(gaussX, newGaussian, '-k', linewidth = 2)
#            plt.errorbar(gaussX, newGaussian, yerr = newerror, fmt = '.-b', linewidth = 2)
            plt.errorbar(xFit, dataFit, yerr = np.sqrt(np.asarray(dataFit)),fmt = 'ok')
            plt.title(file, fontsize = 40)
            plt.tick_params(labelsize = 20)
            plt.xlabel('Channel', fontsize = 30)
            plt.ylabel('Counts', fontsize = 30)
            plt.xlim(interval[0] - delta, interval[1] + delta)
            plt.ylim(ymax = np.floor(1.1*np.max(dataY[interval[0]:interval[1]])))
            plt.show()
        except(NameError):
            print('Plotting error - dataCurveFit', sys.exc_info())
    if pcov == []:
        param = np.ones((1,4))[0].tolist()
        pcov  = np.ones((4,4)).tolist()
    if not isinstance(param, list):
        param = param.tolist()
    if isinstance(pcov, (int, float)):
        pcov  = np.ones((4,4)).tolist()
    elif not isinstance(pcov,list):
        pcov  = pcov.tolist()
    return param, np.sqrt(np.diag(pcov))

def derivate(yval, space, order = 1):
    """
    Calculate the numerical derivative of the yval array using the
    central difference method with correction order space^4
    """
    if order == 1:
        coeff = np.array([1, -8, 0, 8, -1])/(12*space)
    elif order == 2:
        coeff = np.array([-1, 16, -30, 16, -1])/(12*space**2)
    elif order == 3:
        coeff = np.array([1, -8, 13, 0, -13, 8, -1])/(8*space**3)
    elif order == 4:
        coeff = np.array([-1, 12, -39, 56, -39, 12, -1])/(6*space**4)
    elif order == 5:
        coeff = np.array([1, -9, 26, -29, 0, 29, -26, 9, -1])/(6*space**5)
    elif order == 6:
        coeff = np.array([-1, 12, -52, 116, -150, 116, -52, 12, -1])/(4*space**6)
    num = coeff.size
    return np.array(list(np.zeros(int((num - 1)/2)))
                    + list(np.convolve(yval, coeff, 'valid'))
                    + list(np.zeros(int((num - 1)/2))))
    #    return np.insert(np.zeros(num - 1), int((num - 1)/2),
    #                 np.array([sum(yval[val:val + num]*coeff)
    #                           for val in range(len(yval) - num + 1)]))

        ########################################################################
        #  
######### Displays the contents of the input list 
        # 
        ########################################################################
def dis(folder):
    for f, file in enumerate(folder):
        print(index_gen(f, 3), '\t', file)

def e_linear(info):
    (x, dx,  p1, dp1, p0, dp0) = info
    return np.sqrt(x**2*dp1**2 + p1**2*dx**2 + dp0**2)

def e_quadratic(info):
    (x, dx, p2, dp2, p1, dp1, p0, dp0) = info
    return np.sqrt(x**4*dp2**2 + x**2*dp1**2 + dp0**2 +
                   4*x**2*p2**2*dx**2 + p1**2*dx**2)

def exp_rise(t, toff, tau, A):
    return A*(1 - np.exp(-(t-toff)/tau))

def exp_decay(t, toff, tau, A):
    return abs(A)*np.exp(-(t - toff)/abs(tau))



def fileInfo(Folder = None, Path = None, newPath = None, Type = None, New = True):
    r"""
    This program will ask the user for the Path to the folder where the files are
    located and generate a new path / folder to store the files to be generated.
    Parameters:
    -----------
    Folder: list of file names
    Path: str of path to the folder containing the files
    newPath: str of path to save generated files to

    Return:
    -------
    The same information above or default values
    """
    if not Path:
        Path = input('Please enter path to files?\n\'C:\\\n')
    if not Folder:
        Folder = fileSel(Path, Type)
    if not newPath:
        newPath = os.path.join(Path, 'New_text_files')
    if not os.path.exists(newPath) and New:
        os.makedirs(newPath)
    return Folder, Path, newPath


###############################################################################
# Ask the user for path to folder containing the data to be processed and
# allows the user to remove files from the list.
###############################################################################

def fileSel(fileLoc, end = None):
    listFiles = [f for f in os.listdir(fileLoc)
                 if os.path.isfile(os.path.join(fileLoc, f))]
    if end:
        listFiles = [f for f in listFiles if f.endswith(end)]
    listFiles.sort()
    for k, element in enumerate(listFiles):
        print(str(k) + '\t=>' + element)
    print(list(range(len(listFiles))))
    rem = []
    while True:
        remove = input('Enter the numbers of files to remove from processing. Enter e to escape.\n')
        if remove != 'e':
            try:
                rem = str2num(remove)
                temp = [str(k) + '\t=>' + element for k, element in enumerate(listFiles)
                        for index in rem if index == k]
                for k, element in enumerate(listFiles):
                    if all(k!=element for element in rem):
                        print(str(k) + '\t=>' + element)
                    else:
                        print(str(k) + '\t  '+ element)
            except(ValueError, AttributeError):
                print( 'Error in fileSel', sys.exc_info())
        else:
            break
    if rem == []:
        return listFiles
    else:
        return removepts(listFiles, rem)

def findInterval(peak, dataX = np.asarray([]), dataY = np.asarray([]), firDer = np.asarray([]), secDer = np.asarray([])):
    r"""
    This program finds the interval for each peak based on the first and second
    derivative to determine where another peak is starting. 

    Parameters    
    ----------
    dataX: List or array of shape (N,)
        the values are the independent variable
    dataY: List or array of shape (N,)
        the values are the dependent variable
    firDer : List or array of shape (N,)
        the values are the first derivative of the function at each point.
    secDer : List or array of shape (N,)
        the values are the second derivative of the function at each point.
    peak: int
        is the index of where the peak to be considered is located
    Returns
    -------
    interval : List, shape (2)
        of the index corresponding to the start and end of the peak.
    Notes
    -----

    Examples
    --------
	x = np.linspace(0, 10, 201)
	y = 10*np.exp(-(x - 3)**2) + 5 *np.exp(-(x -7)**2)+ np.random.normal(0, 0.05, x.shape)
	yp = firstDer(x, sg_filter(y, 21, 1).tolist())
	ypp = sg_filter(secondDer(x,sg_filter(y, 21, 1)), 21, 1).tolist()
	for peak in findPeak(x, y):
		interval = findInterval(yp, ypp, peak)
		plt.plot(x, y, 'b', linewidth = 3, label = 'Function')
		plt.plot(x, yp, 'r', linewidth = 3, label = 'First Derivative')
		plt.plot(x, ypp, 'k', linewidth = 3, label = 'Second Derivative')
		plt.plot(x[interval[0]:interval[1]], y[interval[0]:interval[1]], 'bD')
		plt.plot(x[interval[0]:interval[1]], yp[interval[0]:interval[1]], 'rD')
		plt.plot(x[interval[0]:interval[1]], ypp[interval[0]:interval[1]], 'kD')
		plt.legend()
		plt.show()
    """
    interval = [0 , 0]
    n = 0
    initial = True
    final = True
    peak = int(peak)
    dataY = array(dataY)
    dataX = array(dataX)
    firDer = array(firDer)
    secDer = array(secDer)
    if not (all([dataY.size, dataX.size]) or all([firDer.size, secDer.size])):
        raise TypeError('findInterval requires two inputs either dataX ' +
                        'and dataY or firDer and secDer')
    if dataY.size and dataX.size:
        if not firDer.size:
            firDer = firstDer(dataX, sg_filter(dataY, 3, 1))
        if not secDer.size:
            secDer = sg_filter(secondDer(dataX,sg_filter(dataY, 3, 1)), 3, 1)
    while True:
        if all(0 >= el for el in firDer[peak - n - 2: peak - n]) and initial:
            initial = False
        if all(0 <= el for el in firDer[peak + n + 1: peak + n + 3]) and final:
            final = False
        if (firDer[peak  - n] > 0 and
            secDer[peak  - n] >= secDer[peak - n + 1]) and initial:
            interval[0] = peak - n
        if (firDer[peak  + n] < 0 and
            secDer[peak  + n] >= secDer[peak + n - 1]) and final:
            interval[1] = peak + n
        if not initial and not final:
            break
        n += 1
    if interval[1] <= interval[0] <= peak:
        interval[1] = peak + 10
    if interval[0] >= interval[1] >= peak or interval[0] == 0:
        interval[0] = peak - 10
    return interval

def findPeak(dataX, dataY, base = 2, per = 5, window = None):
    r"""
    This does a search for peaks given dataY and dataX using the smoothed data
    and smooth fits to the 1st and 2nd derivatives. If the 1st der. changes from +
    to - and the 3 points before it are positive while the 3 points after are
    negative and the two points in the 2nd der. corresponding to the left and right
    of the point are also below zero with the point and the corresponding 5 points
    in the data are above the basenoise as given or calulated from the difference of
    the smooth data to the actual data, then record the point.

    Parameters    
    ----------
    dataX: List or array of shape (N,)
        the values are the independent variable
    dataY: List or array of shape (N,)
        the values are the dependent variable
    firDer : List or array of shape (N,)
        the values are the first derivative of the function at each point.
    secDer : List or array of shape (N,)
        the values are the second derivative of the function at each point.
    baseNoise: float or int
        sets the lower threshold for a peak to be considered
    Returns
    -------
    index : List, shape (N)
        of the index of the peak.
    Notes
    -----

    Examples
    --------
    import numpy as np
    from matplotlib import pyplot as plt
    x = np.linspace(0, 10, 201)
    y = 100*np.exp(-(x - 5)**2) + 5*np.random.random(x.shape)
    peak = findPeak(x, y, base = 1, per = 5)
    A = plt.plot(x, y)
    vertical(peak, y, x)
    plt.show(block = False)

    """
    dataY = array(dataY)
    dataX = array(dataX)
    val = window if window else int(dataY.size*per/100)
    base = base if base else 2
    val = val if val%2 == 1 else val + 1 
    noise = dataY - sg_filter(dataY, val, 0)
    val2 = (val + 1)/2 + ((val + 1)/2)%2
    #nfilt = sg_filter(noise, val2, 0)
    noise = noise - noise[val2:-val2].mean()
    std = noise.std()
    upper = np.where(noise[val2:-val2].mean() + 1.5*base*noise[val2:-val2].std() < noise)[0]
    lower = np.where(noise[val2:-val2].mean() - base*noise[val2:-val2].std() > noise)[0]
    lrange, urange = [], []
    while True:
        try:
            if upper.size:
                ran = np.where(upper.min() > lower)[0]
                lrange.append([lower[0], lower[ran[-1]]])
                lower = np.delete(lower, ran)
            else:
                lrange.append([lower[0], lower[-1]])
                break
            if lower.size:
                ran = np.where(lower.min() > upper)[0]
                urange.append([upper[0], upper[ran[-1]]])
                upper = np.delete(upper, ran)
            else:
                urange.append([upper[0], upper[-1]])
                break
        except(IndexError):
            print('Could not fit peaks!')
            break
    index = []
    for line in urange:
        if line[0] != line[1]:
            part = dataY[line[0]:line[1]]
            index.append(line[0] + np.where(part.max() == part)[0][0])        
    return index

        ########################################################################
        # The first derivative is approximated by fitting a line to 21 
######### consecutive points and taking the slope as the derivative of the 
        # middle point. Equation: m*x + b
        ########################################################################
def firstDer(dataX, dataY, n = None):
    dataY = array(dataY)
    dataX = array(dataX)
    if not n:
        n = int(0.001*dataY.size)
        if n < 3:
            n = 3
    m = dataY.size
    firDer = (dataY[:n] - dataY[1:n+1])/(dataX[:n] - dataX[1:n+1])
    for a in range(n, m - n):
        pfit = np.polyfit(dataX[a:  a + n], dataY[a: a + n], deg = 1)
        firDer = np.insert(firDer, firDer.size, np.asarray([pfit[0]]))
    firDer = np.insert(firDer, firDer.size,
                       (dataY[m - n-1:m-1] - dataY[m-n:])/(dataX[m-n-1:m-1] - dataX[m-n:]))
    return firDer

        ########################################################################
        # This defines the function used by cfit for fitting the gaussian with
######### background to the data
        # 
        ########################################################################
def gaussian(x, c, s):
    return np.exp(-(x-c)**2/(2*s**2))


def gen_fun(peak_num, name = None):
#def gen_fun(peak_num, peak_ind, name = None):
    """
    This function will generate a new function modeled off of compFun02 which is a
    simple gaussian with a given area A centered at c and std s. The generated function
    will be a integar multiples of c with some offset. This is intended for the laser
    calibration of the 3.5 eV pulse laser.
    """
    
    if name:
        pass
    else:
        name = 'Function_{}.py'
        i = 0
        while True:
            if not os.path.exists('C:\\Python33\\{}'.format(name.format(i))):
                name = name.format(i)
                break
            i += 1
    with open(name, 'w') as file:
        nVari = [', C_{}, S_{}, A_{}'.format(int(n), int(n), int(n))
                 for n in peak_num]
        nfun = '(' + 'compFun02(x{}) +\n\t\t'*len(peak_num)
        firstline = 'from modules_V5_7 import compFun02\n'
        secondline = 'def compFunAA(x{}):\n'.format(''.join(nVari))
        thirdline = '\treturn ' + nfun.format(*nVari)[:-5] +')'
        file.write(firstline)
        file.write(secondline)
        file.write(thirdline)

        ########################################################################
        # This will allow the user to import data from a folder into a dataFile
######### class that can then be manipulated by the dataFile class methods.
        # 
        ########################################################################
def impData(Path = None, fileList = None, Verbose = False, Delim = '\t', tline = None):
    fileList, Path, newPath  = fileInfo(fileList, Path, New = False)
    parent = dataFile(fileList, path = Path)
    identical = False
    for k, file in enumerate(fileList):
        parent.set(file, dataFile(path = Path, file = file))
        parent.get(file).set('All', array(loadFile(os.path.join(Path, file))))
        #Remove empty lines
        notemptyline = np.where(parent.get(file).get('All') != '')[0]
        parent.get(file).set('All', parent.get(file).get('All')[notemptyline])
        m = 0
        while identical or tline == 0 or tline:
            if tline == 0 or tline:
                column = tline
            if Verbose:
                print('Using titles from before.\n')
            if column == 'none':
                Ncol = len(parent.get(file).get('All')[start].split(Delim))
                ctitles = ['Column {}'.format(n) for n in range(Ncol)]
                dataStr = reindex([line.split(Delim)
                                   for line in parent.get(file).get('All')[int(column) + 1:3*10**6]
                                   if line != (len(ctitles) - 1)*Delim])
            else:
                ctitles = parent.get(file).get('All')[int(column)].split(Delim)
                dataStr = reindex([line.split(Delim)
                                   for line in parent.get(file).get('All')[int(column) + 1:3*10**6]
                                   if line != (len(ctitles) - 1)*Delim])
            break
        while not identical and not tline and tline != 0:
            for a in range(m + 0, m + 9):
                if a < len(parent.get(file).get('All')):
                    print('{}\t{}'.format(a, parent.get(file).get('All')[a]))
                else:
                    break
            column = input('Which if any of these is the column titles? Enter n for' +
                           ' next, p for previous, e to exit or' +
                           ' none to select data starting point.\n')
            if column == 'e':
                break
            elif column == 'p':
                m = max(0, m - 10)
            elif column == 'n':
                m += 10
            elif column == 'none':
                while True:
                    try:
                        start = int(input('Which row does the data start on?\n'))
                        Ncol = len(parent.get(file).get('All')[start].split(Delim))
                        ctitles = ['Column {}'.format(n) for n in range(Ncol)]
                        print(parent.get(file).get('All')[start].split(Delim))
                        # Limit the loaded data to 1000000 points
                        dataStr = reindex([line.split(Delim)
                                           for line in parent.get(file).get('All')[int(column) + 1:3*10**6]
                                           if line != (len(ctitles) - 1)*Delim])
                        if len(fileList) > 1:
                            action = input('Are these files identical? y/n\n')
                            if action == 'y':
                                identical = True
                        break
                    except(ValueError):#, TypeError):
                        print('Error in No Columns.\n', sys.exc_info())

            else:
                try:
                    ctitles = parent.get(file).get('All')[int(column)].split(Delim)
                    if Verbose:
                        [print(title) for title in ctitles]
                    dataStr = reindex([line.split(Delim)
                                       for line in parent.get(file).get('All')[int(column) + 1:3*10**6]
                                       if line != (len(ctitles) - 1)*Delim])
                    confirm = input('Confirm titles.\n')
                    if confirm == 'y':
                        if len(fileList) > 1:
                            action = input('Are these files identical? y/n\n')
                            if action == 'y':
                                identical = True
                        break
                except(ValueError):
                    print('Entry error', sys.exc_info())
        if Verbose:
            print(file)
        for t, title in enumerate(ctitles):
            cVal = array(dataStr[t])
            try:
                tag = np.where('' != cVal)
                parent.get(file).set(title, cVal[tag].astype(float))
            except(ValueError):
                parent.get(file).set(title, cVal)
                print('Could not convert {} to float.\n'.format(title))
        del m, t
    for key in parent.keys():
        parent.get(key).remove('All')
        if Verbose:
            print('File:\t', key)
    if Verbose:
        print(Path)
    if len(parent.keys()) == 1:
        parent = parent.get(parent.keys()[0])
    return parent

def index_gen(num, length = 3):
    """
    Generates a number of length length filling the left side with zeros.
    """
    num = str(int(num))
    return '0'*(length - len(num)) + num

def initial_poission(dataX, dataY, peaks, LLD, ULD, fit_order, conv = np.asarray([]), file = '', Auto = False):
    """
    This function performs the initial fit for laser_adv to guess
    an appropriate A and mu for the convolution fit.
    Parameters:
    -----------
    dataX : list, array of xaxis
    dataY : list, array of yaxis
    peaks : list, array of peaks from which the poissonian will be fit
    fit_order: int, (1, 2) corresponding to the order of the fit for the peaks
    Returns:
    --------
    LLD : int lower bound of data to be used
    ULD : int upper bound of data to be used
    pos : list [mu, amp] (variance and amplitude)
    """
    global _npeaks, _start, npeaks
    npeaks =peaks
    if fit_order != (1 or 2):
        print('Input error in initial_poissonian!')
        raise
    inverse = [inv_line, inv_quad][fit_order - 1]
    while True:
        _npeaks = len(peaks)
        N = np.arange(_start + 1, _npeaks + _start + 1)
        conv = conv if conv.size else np.polyfit(N, peaks, fit_order)
        Amp = np.asarray([dataY[p] for p in peaks])
        mu = inverse(sum(Amp*np.asarray(peaks))/sum(Amp),
                     0, conv, np.zeros(conv.size))[0]

        if Auto:
            break
        inv_x = inverse(dataX, 0, conv, np.zeros(conv.size))[0]
        while True:
            try:
                rem = np.where(0 == Amp)[0]
                Amp = np.delete(Amp, rem)
                N = np.delete(N, rem)
                pos, posco = cfit(compFun04, N, Amp, p0 = [mu, max(Amp)])
                break
            except:
                print('Error in Initial Poissonian Fit: {}'.format(sys.exc_info()))
                global trouble
                trouble = [Amp, N, mu, conv]
                try:
                    pos = str2num(input('Enter suggested values'))
                    break
                except(ValueError, TypeError):
                    pass
                    
                
        plt.plot(dataX, dataY, 'o', markersize = 10,
                     markerfacecolor = (1, 1, 1))
        plt.plot(dataX[LLD:ULD], dataY[LLD:ULD], 'o', markersize = 10)
        plt.plot(dataX, compFun04(inv_x, *pos), linewidth = 5)
        vertical(peaks[0], dataY)
        plt.title(file, size = 30)
        print('Current LLD, ULD and Offset (marked peak): ',
              [dataX[LLD], dataX[ULD], _start])
        print('Current poissonian fit parameters:\n' +
              '{}'.format(list(np.around(pos))))
        plt.xlim((dataX[LLD], dataX[ULD]))
        plt.ylim((1, 1.1*max(dataY[LLD: ULD])))
        plt.show()
        if not Auto:
            action = input('Enter new settings, [l, u, o] or e to exit.\n')
        else:
            action = 'e'
        if action.startswith('e'):
            del action
            break
        try:
            LLD, ULD, _start = str2num(action)
            LLD = np.where(dataX <= LLD)[0][-1]     #               ADDED on 20160928 so that energy values
            ULD = np.where(ULD <= dataX)[0][0]     #               Could be used...
            if ULD <= LLD:
                print('ULD is less than LLD')
                raise(ValueError)
        except(ValueError, TypeError):
            print('Error in setup: {}'.format(sys.exc_info()))
    if Auto:
        return mu
    else:
        return LLD, ULD, pos

def inv_line(y, dy, p, dp):
        out = y/p[0] - p[1]/p[0]
        sig = np.sqrt(dy**2 + dp[1]**2 + out**2*dp[0]**2)/np.abs(p[0])
        return out, sig

def inv_quad(y, dy, p, dp):
        out = (-p[1] + np.sqrt(p[1]**2 + 4*p[0]*(y - p[2])))/(2*p[0])
        sig = np.sqrt(dy**2 + dp[2]**2 + out**4*dp[0]**2 +
                      out**2*dp[1]**2)/np.abs(2*p[0]*out + p[1])
        return out, sig

def laser_adv(data = None, peaks = np.asarray([]), sort = 0, function = compFun02, dataX = np.asarray([]), usexaxis = None):
    """
    This program is intended to process the calibration files from the laser.
    Parameters:
    -----------
    ran : int value for the initial spread for the interval
    Returns:
    --------
    """
    global _npeaks, _start
    _npeaks = 0
    _start = 0
    
    param = []
    Auto = False
    bypass = True
    if not isinstance(data, dataFile):
        data = impData()
    peaks = np.asarray(peaks)
    out = np.asarray([])
    prior = np.asarray([])
    folder = data.keys()
    dataX = np.asarray(dataX)
    try:
        if sort == 0:
            folder.sort(key = lambda line: int(line.split('_')[1])
                        if len(line.split('_')) > 1 else 0)
        else:
            folder.sort(key = lambda line: int(line.split('_')[-1][:-sort])
                        if len(line.split('_')) > 1 else 0)
    except:
        print('Could not sort the file proceeding unsorted...')
    out_fit = [['00_File', '01_Line (Ch/#)', '02_dL (Ch/#)',
                '03_Const (Ch)', '04_dC (Ch)', '05_Sig (Ch)',
                '06_dS (Ch)',  '07_Amp (Cnts)', '08_dA (Cnts)',
                '09_Mu (#)', '10_dMu (#)']]
    multi_ls = [['00_File', '01_Energy (eV)', '02_Cent (Ch)',
                   '03_dC (Ch)', '04_Sig (Ch)', '05_dS (Ch)',
                   '06_Area (Cnts*Ch)', '07_dA (Cnts*Ch)']]
    multi_odr = [['00_File', '01_Energy (eV)', '02_Cent (Ch)',
                   '03_dC (Ch)', '04_Sig (Ch)', '05_dS (Ch)',
                   '06_Area (Cnts*Ch)', '07_dA (Cnts*Ch)']]
    multi_line_cent = [['00_File', '01_Line (eV/Ch)',
                        '02_dL (eV/Ch)', '03_Const (eV)',
                        '04_dC (eV)']]
    tofullfit = [['00_File', '01_Line (Ch/#)',
                  '02_dL (Ch/#)', '03_Const (Ch)',
                  '04_dC (Ch)','05_Sig (Ch)', '06_dS (Ch)',
                   '07_Area (Cnts)', '08_dA (Cnts)']]

    form = data.file.split('_')[0]
    page1 = data.file.split('_b13-')[-1].replace('.txt', '').replace('_merge', '').replace('_000', '')
    page2 = data.file.split('_b14-')[-1].replace('.txt', '').replace('_merge', '').replace('_000', '')
    page = page1 if len(page1) < len(page2) else page2
    form_page = '_'.join([form, page])
    analysis_path = data.path + '\Analysis_{}'.format(form_page)
    if not os.path.exists(analysis_path):
        os.makedirs(analysis_path)

    multidata = dataFile()
    multidata.path = analysis_path
    multidata.file = 'LS_Peaks_Multi_{}.txt'.format(form_page)

    multiODR = dataFile()
    multiODR.path = analysis_path
    multiODR.file = 'ODR_Peaks_Multi_{}.txt'.format(form_page)
    
    multiline = dataFile()
    multiline.path = analysis_path
    multiline.file = 'Linear_Fit_Multi_{}.txt'.format(form_page)

    linefit = dataFile()
    linefit.path = analysis_path
    linefit.file = 'Linear_Fit_Conv_{}.txt'.format(form_page)

    fullfit = dataFile()
    fullfit.path = analysis_path
    fullfit.file = 'Full_Fit_Conv_{}.txt'.format(form_page)

    energyODR = dataFile()
    energyODR.path = analysis_path
    energyODR.file = 'ODR_Energy_Fits_{}.txt'.format(form_page)
    
    energyLS = dataFile()
    energyLS.path = analysis_path
    energyLS.file = 'LS_Energy_Fits_{}.txt'.format(form_page)

    fitmodel = Model(np.polyval)
    multimodel = Model(compFun11)

    inverse = inv_line
    fun = compFun01
    fitdata = linefit
    fit_order = 1
    
    for f, file in enumerate(folder[:]):
        print('\nStarting to Process: {}'.format(file))
        dataY = data.get(file)
        dataX = np.arange(1, dataY.size + 1) if not dataX.size else dataX
        if isinstance(usexaxis, dataFile):
            dataX = usexaxis.get(usexaxis.keys()[f])
        peaks = np.asarray(findPeak(dataX, dataY)) if not peaks.size else peaks
        peaks = peak_select(dataX, dataY, peaks, file) if not Auto else peaks
        if not param:
            LLD, ULD, _start = 0, dataY.size - 1, 0
        if not Auto and not out.size:
            LLD, ULD, pos = initial_poission(dataX, dataY, peaks, LLD,
                                             ULD, fit_order, file = file)
        ULD = min(ULD, dataY.size)
        action = 'Run'
        if (Auto or out.size) and not isinstance(usexaxis, dataFile):
            temp_peaks = np.polyval(out[-5:-3], range(1, _npeaks + 1))
            print(temp_peaks)
            temp_peaks = np.array([np.where(dataX <= val)[0][-1] for val in temp_peaks])
            print(temp_peaks)
            print(out[-5:-3])
            temp_peaks = np.delete(temp_peaks,  np.where(dataY.size <= temp_peaks)[0])
            mu = initial_poission(dataX, dataY, temp_peaks, LLD, ULD,
                                  fit_order, conv = out[-5:-3], Auto = True)
            print(mu)
            mu = out[-1] # temporary fix...
            print(mu)
        elif (Auto or out.size) and isinstance(usexaxis, dataFile):
            mu = out[-1] # temporary fix...
        if out.size:
            param = list(out[:-1]) + [mu]
        else:
            New = np.arange(_start + 1, len(peaks) + _start + 1)
            conv = np.polyfit(New, dataX[peaks], fit_order)
            sig = np.mean(dataX[peaks[1:]] - dataX[peaks[:-1]])/2
            param = list(conv) + [sig, pos[1], pos[0]]
        attempts = 0
        tally = 0
            
        while True:
            avg_dis = np.mean(dataX[peaks[1:]] - dataX[peaks[:-1]]) if not out.size else out[-5]
            _npeaks = (dataX[ULD] - dataX[LLD])/avg_dis + _start
            if _npeaks < 0 :
                print('Peak locations: ', list(peaks))
                print('dataX[peaks]', list(dataX[peaks]))
                print('OUT values', list(out))
                print('_npeaks: ', _npeaks)
                print('_start: ', _start)
                print('LLD : {}\tULD: {}\tavg_dis: {}'.format(LLD, ULD, avg_dis))
                _npeaks = (dataX[ULD] - dataX[LLD])/avg_dis + _start
            out = np.asarray([])
            if not Auto:
                plt.plot(dataX, dataY, 'bD', markersize = 10)
                plt.plot(dataX, fun(dataX, *param), '--r', linewidth = 5)
            try:
                ydat, rem = removepts(dataY[LLD : ULD], spe = 0, ret = True)
                xdat = removepts(dataX[LLD : ULD], rem = rem)
                out, outerr = cfit(fun, xdat, ydat, p0 = param)
                plt.plot(dataX, fun(dataX, *out), 'k', linewidth = 5)                    

                if abs(out[-4]) > 0.4*abs(out[-5]):
                    if not Auto:
                        print('\nRerunning...')
                    out[-4] = out[-4]%out[-5]
                    param = list(out)
                    attempts += 3
                if Auto:
                    if tally > 2:
                        raise RuntimeError
                    if abs(out[-3]) >= 1.25*abs(param[-3]):
                        out[-3] = param[-3]
                        param = list(out)
                        attempts += 2

                    if attempts != 0:
                        raise Exception('ConvError')
                    attempts += 1
            except (TypeError, RuntimeError):
                Auto = False
                print('Error in initial Convolution ' +
                      'Fit: {}'.format(sys.exc_info()))
            except Exception as element:
                attempts += 1
                print('ConvError,  Attempts: {}'.format(attempts))
                print(element)
            if not Auto and (attempts == 0 or tally > 0):
                print('LLD = {}, ULD = {}, _start = {}'.format(dataX[LLD], dataX[ULD], _start))
                print('Current parameters '+
                      '({}):\n{}'.format(inspect.getargspec(fun)[0][1:],
                                         list(np.around(param, 5))
                                         if not out.size
                                         else list(np.around(out, 5))))
            tally += 1
            if not Auto:
                if attempts > 1 and tally <= 1:
                    action = 'Redo'
                    attempts = 0
                else:
                    plt.xlim((dataX[LLD], dataX[ULD]))
                    plt.ylim((1, 1.1*max(dataY[LLD: ULD])))
                    plt.title(file, size = 30)
                    try:
                        plt.show()
                    except(AttributeError):
                        print('Could not plot data...')
                    action = input('l, u, or p to enter new LLD, ULD or param, ' +
                                   'i to change both l and u, s to skip this file, ' +
                                   '_s to change starting peak or e to escape ' +
                                   '(A to set to Auto from here).\n')
            else:
                print('Tally: {}, Attempts: {}'.format(tally, attempts))
                if tally > 2:
                    Auto = False
                    action = 'Redo'
                elif attempts > 1:
                    action = 'Redo'
                    attempts = 0
                else:
                    print('\nFINAL PARAMETERout = S '+
                          '({}):\n{}'.format(inspect.getargspec(fun)[0][1:],
                                             list(np.around(out, 5))))
                    action = 'e'
            if action.startswith('A'):
                Auto = True
                action = 'e'
            if action.startswith('e'):
                group = [out, np.sqrt(np.diag(outerr))]
                group = np.asarray(reindex(group)).flatten()
                group = list(np.around(group, 5))
                out_fit.append([file] + group)
                break
            if action.startswith('s'):
                out = np.asarray(out_fit[-1][1::2])
                break
            try:
                if action.startswith('l'):
                    LLD = np.where(dataX <= int(action[1:]))[0][-1]
                elif action.startswith('u'):
                    ULD = np.where(int(action[1:]) <= dataX)[0][0]
                elif action.startswith('i'):
                    LLD, ULD = str2num(action[1:])
                    LLD = np.where(dataX <= LLD)[0][-1]
                    ULD = np.where(ULD <= dataX)[0][0]
                elif action.startswith('p'):
                    param = str2num(action[1:])
                    if len(param) != len(inspect.getargspec(fun)[0][1:]):
                        raise TypeError
                    out = np.array([]) # Clear fitted out values
                elif action.startswith('_s'):
                    _start = int(action[2:])
                else:
                    pass
            except(ValueError, TypeError):
                print('Error in edits to Convolusion ' +
                      'Fit: {}'.format(sys.exc_info()))
        
        if action.startswith('s'):
            print('{} is being skipped!'.format(file))
            del action
            
        fit_fun_name = data.file.replace('.txt', '_{}.py'.format(file))
        fit_fun_name = fit_fun_name.replace('-', '_')

        useme = out_fit[-1][1::2]
        #print('useme: ', useme)
        peak_param = np.asarray([1/useme[0], -useme[1]/useme[0]])
        #print('peak_param: ', peak_param)
        new_param = np.asarray([np.polyval(useme[:-3], p)
                                for p in range(1, _start + _npeaks + 1)])
        #print('Starting new_param: ', new_param)
        for_rem = np.where(dataX[LLD] <= np.asarray(new_param))[0]
        new_param = new_param[for_rem]
        #print('After removing peaks below LLD: ', new_param);
        new_ampli = np.asarray([fun(val, *useme) for val in new_param])
        new_area = new_ampli*np.sqrt(2*np.pi*useme[-3]**2)
        #print('Corresponding Areas: ', new_area)
        #use_peaks = np.where(1000 < new_area)[0]
        use_peaks = np.where(0.10*max(new_area) < new_area)[0]
        new_param = np.asarray(new_param)[use_peaks].astype(int)
        new_area = np.asarray(new_area)[use_peaks].astype(int)
        new_sig = useme[-3]*np.ones(new_area.size)
        #print('Final centroids: ', new_param)
        #print('corresponding areas: ', new_area)

        new_param = [np.asarray(line)
                     for line in zip(new_param, new_sig, new_area)]
        new_param = np.asarray(new_param).flatten()
        #print('Final new_param: ', new_param)
        try:
            upp = np.where(int(new_param[-3] + 2*new_param[-2]) <= dataX)[0][0]
            low = np.where(dataX <= int(new_param[0] - 2*new_param[1]))[0][-1]
            upp = min(upp, dataY.size)
            low = max(low, 0)
        except(IndexError):
            upp = dataY.size
            low = 0
        ####        Use a Orthogonal Distance Regression For Multi        ####
        skip = False
        multi_odr_add, low, upp = multi_fit(dataX, dataY, peaks, new_param, upp, low,
                                            fit_fun_name, peak_param, analysis_path,
                                            file, method = 'ODR', Auto = Auto,
                                            ALL = True)
        if multi_odr_add == None:
            skip = True
        else:
            multi_odr = multi_odr + multi_odr_add

            ####        Use a Linear Least Squares Fit For Multi        ####
            new_param = [np.asarray(line)
                         for line in reindex(reindex(multi_odr_add)[2::2])]
            new_param = np.asarray(new_param).flatten()
            #upp = new_param[-3] + 2*new_param[-2]
            #low = new_param[0] - 2*new_param[-2] 
            #gen_fun(new_param[::3], fit_fun_name)
            multi_ls_add = multi_fit(dataX, dataY, peaks,
                                     new_param, upp, low,
                                     fit_fun_name, peak_param,
                                     analysis_path, file, method = 'LS',
                                     Auto = Auto)
            if multi_ls_add == None:
                skip = True
            else:
                multi_ls = multi_ls + multi_ls_add
        
        if not skip:
            try:
                y_energy, x_cent, dx_cent = [np.asarray(line)
                                             for line in reindex(multi_odr_add)[1:4]]
                hc = 1239.84193 # eV*nm
                #wavelength = 355 # nm
                #wavelength_err = 1 # nm
                #wave_energy = hc/wavelength # eV
                #wave_energy = 3.4974739281780685 #eV This is based on the calibration by Hg updated 2015/05/04
                #wave_err = hc/wavelength**2*wavelength_err # eV
                initial_guess = np.polyfit(x_cent, y_energy, 1)
                data_multi = RealData(x_cent, y_energy, sx = dx_cent)
                line_odr = ODR(data_multi, fitmodel, beta0 = initial_guess)
                line_odr.set_job(fit_type = 0)
                line_odr_out = line_odr.run()
                to_mline = np.asarray(reindex([line_odr_out.beta,
                                               line_odr_out.sd_beta])).flatten()
                multi_line_cent.append(['ODR_' + file] + list(to_mline))
                energyODR.set('{}_Energy'.format(index_gen(f)),
                              np.polyval(line_odr_out.beta, dataX))

                multilinefit = [np.asarray(line) for line in reindex(multi_ls_add)[1:4]]
                global trouble
                trouble = [multilinefit, multi_ls_add, initial_guess,
                           dataX, dataY, peaks, new_param, upp, low, fit_fun_name,
                           peak_param, analysis_path, file]
                mline, mlerr = cfit(bin_line, multilinefit[0],
                                    multilinefit[1], p0 = initial_guess,
                                    sigma = multilinefit[2], absolute_sigma = True)
                to_mline = np.asarray(reindex([mline, np.sqrt(np.diag(mlerr))])).flatten()
                multi_line_cent.append(['LS_' + file] + list(to_mline))
                energyLS.set('{}_Energy'.format(index_gen(f)),
                              np.polyval(mline, dataX))
                ###     Perform Advanced Fit    ###
                #global trouble
                #trouble = [new_param[:], mline, multi_ls_add, dataX, dataY, low, upp]
                #raise TypeError
                full_starter = int(round(np.polyval(mline, new_param[0])/3.5, 0))
                new_param = reindex([new_param[1::3], new_param[2::3]])
                new_param = np.asarray([np.asarray(line) for line in new_param])
                new_param = new_param.flatten()
                new_param = np.insert(np.asarray([1/mline[0], -mline[1]/mline[0]]), 2, new_param)
                full_data = RealData(dataX[low: upp], dataY[low:upp], sy = np.sqrt(dataY[low:upp]))
                make_compFun16_ODR(full_starter)
                full_model = Model(compFun16_ODR)
                full_odr = ODR(full_data, full_model, beta0 = new_param)
                full_odr.set_job(fit_type = 0)
                full_odr_out = full_odr.run()
                to_full = np.asarray(reindex([full_odr_out.beta, full_odr_out.sd_beta])).flatten()
                tofullfit.append(['ODR_' + file] + list(to_full))
                

                ###     Plot Residuals          ###
                resid, dres = residual(line_odr_out.beta, line_odr_out.sd_beta,
                                       x_cent, dx_cent, y_energy)
                plt.clf()
                plt.close('all')
                fig, ax = plt.subplots(1)
                ax.errorbar(y_energy, resid, yerr = dres, fmt = 's', markersize = 30,
                            label = 'ODR - Residual')
                resid, dres = residual(mline, np.sqrt(np.diag(mlerr)),
                                       multilinefit[1], multilinefit[2],
                                       multilinefit[0])
                ax.errorbar(multilinefit[0] + 0.25, resid, yerr = dres, # Introduces slight offset
                            fmt = 'D', markersize = 30, label = 'LS - Residual')
                ax.legend(prop = {'size': 50})
                ax.tick_params(labelsize = 75)
                ax.set_ylabel('Residual (eV)', size = 75)
                ax.set_xlabel('Energy (eV)', size = 75)
                figpath = os.path.join(os.path.join(analysis_path, 'Peak_Fits'),
                               file.replace('LS_', '').replace('ODR_', ''))
                #standard_fig_save('Residuals_{}.pdf'.format(file), figpath, Verbose = False)
                plt.close('all')
            except(RuntimeError):
                print('Could not fit the Linear Multi Fit: {}'.format(sys.exc_info()))
     
        try:
            _start = round(inv_line(peaks[0], 0, out[-5:-3],np.zeros(2))[0])
        except(IndexError):
            print('Error in new _start calculation: {}'.format(sys.exc_info()))        
        if bypass and not Auto:
            action = input('Automate the rest?\n')
        if not Auto and action.startswith('y'):
            Auto = True
        elif action.startswith('n'):
            bypass = False
        if skip:
            bypass = False
        if not bypass and 'y' == input('Skip the remaining files?\n'):
            break

        for column in reindex(multi_ls):
            multidata.set(column[0], column[1:])
        multidata.save(Verbose = False)
        for column in reindex(multi_odr):
            multiODR.set(column[0], column[1:])
        multiODR.save(Verbose = False)
        for column in reindex(multi_line_cent):
            multiline.set(column[0], column[1:])
        multiline.save(Verbose = False)
        for column in reindex(out_fit):
            fitdata.set(column[0], column[1:])
        fitdata.save(Verbose = False)
        for column in reindex(tofullfit):
            fullfit.set(column[0], column[1:])
        fullfit.save(Verbose = False)
        energyODR.save(Verbose = False)
        energyLS.save(Verbose = False)

def linear(x, p0, p1):
    r"""
    This function is a simple linear equation
    """
    return p1*x + p0

def loadFile(fileName):
    r"""
    This returns the contents of a file as a list of strings having split the
    file by the deliminating newline '\n'

    Parameters:
    -----------
    fileName: Str, full path and name of file to be open

    Returns:
    --------
    hold : List of str with the contents of the file
    """
    with open(fileName, 'r') as dataStr:
        hold = dataStr.read().split('\n')
        if hold[-1] == '':
            del hold[-1]
    return hold

def modifyData(dataStr, Out = True):
    r"""
    This function locates where the data file starts and end and then removes
    portions of a string list that is not part of the data set and converts the
    output portion to a numeical list.
    """
    dataNum = []
    initial = 0
    final = len(dataStr)
    while True:
        try:
            if all(int(values) >= 0 for values in dataStr[initial: initial + 10]):
                break
        except(ValueError):
            pass
        initial += 1
    while True:
        try:
            if all(int(values) >= 0 for values in dataStr[final-10: final + 1]):
                break
        except(ValueError):
            pass
        final -= 1
    dataNum = np.asarray(dataStr[initial: final]).astype(int)
    if Out:
        return dataNum, initial, final
    else:
        return dataNum

def multi_fit(dataX, dataY, peaks, new_param, upp, low, fit_fun_name, peak_param = np.asarray([]), path = None, file = '', ALL = False, method = 'ODR', forced = False, Auto = True):
    """
    Performs fit to the data based on the inputs...
    """
    global __prior, __range, trouble
    attempt = 0
    #forced = False
    initial_run = 0
    #if method == 'LS':
    #    globals()['func'] = __import__('{}'.format(fit_fun_name.replace('.py', '')))
    multimodel = Model(compFun11)
    fitfun = [o[1] for o in inspect.getmembers(fitFunctions) if inspect.isfunction(o[1])]
    #global trouble
    print('Starting Multi-Fit: {}'.format(method))
    #peak_index = np.arange(peak_index[0], peak_index[0] + peaks.size + 1)######
    #tally = 0
    old_param = new_param[:]
    if not path:
        path = 'C:\\Presentations\\Figures'
    skip = False
    while not skip:
        try:
            for iteration in range(10):
                rem = np.where(0 == dataY[low: upp])[0]
                ydat = np.delete(dataY[low: upp], rem)
                xdat = np.delete(dataX[low: upp], rem)
                if method == 'ODR':
                    trouble = [dataX, dataY, new_param, low, upp]
                    full_data = RealData(xdat, ydat, sy = np.sqrt(ydat))
                    full_odr = ODR(full_data, multimodel, beta0 = new_param)
                    full_odr.set_job(fit_type = 0)
                    full_odr_out = full_odr.run()
                    nout = np.abs(full_odr_out.beta)
                    nerr = full_odr_out.sd_beta
                if method == 'LS':
                    nout, ncov = cfit(fitfun[int(len(new_param)//3) - 1], xdat,
                                      ydat, p0 = new_param,
                                      sigma = np.sqrt(ydat),
                                      absolute_sigma = True)
                    nout = np.abs(nout)
                    nerr = np.sqrt(np.diag(ncov))
                if all(np.abs(new_param - nout) <= 0.001*np.ones(nout.size)):
                    break
                new_param = nout
                avg_peak_dis = np.asarray(new_param[3::3] - new_param[:-3:3]).mean()/2
                if not forced:
                    low = np.where(dataX <= nout[0] - min(nout[-2], avg_peak_dis))[0][-1]
                    upp = np.where(nout[-3] + min(2*nout[-2], avg_peak_dis) <= dataX)[0][0]
                    low = max(low, 0)
                    upp = min(upp, dataX.size)
                ydat, rem = removepts(dataY[low: upp], spe = 0, ret = True)
                xdat = removepts(dataX[low: upp], rem = rem)
            if any([math.isinf(el) or math.isnan(el) for el in nerr]):
                print('\nImproper fit, the error is not finite.\n')
                raise RuntimeError
            if any([el < dataX[low] or el > dataX[upp] for el in new_param[::3]]):
                print('A centroid is outside the bounds.')
                raise RuntimeError
            avg_sigma = new_param[1::3].mean()
            if any([el > 2*avg_sigma or el <avg_sigma/2 for el in new_param[1::3]]):
                print('A sigma varies greatly beyond the other sigmas')
                print('Sigma: {}'.format(new_param[1::3]))
                if Auto and attempt != 0 and 'y' == input('Would you like to continue anyways?\n'):
                    pass
                else:
                    raise RuntimeError
            temp_values = np.asarray(reindex([nout, nerr])).flatten()
            temp_values = [np.asarray(temp_values[a::6]) for a in range(6)]
            file_column = [file]*temp_values[0].size
            peak_index = np.round(np.polyval(peak_param, nout[::3])).astype(int)
            energy_column = 3.5*peak_index if peak_index.size else np.ones(len(file_column))
            temp_values = [file_column, energy_column] + list(np.round(temp_values, 5))
            temp_values = reindex(temp_values)
            lplot = int(max(nout[0] - 3*min(nout[1], avg_peak_dis), 0))
            lplot = np.where(dataX <= lplot)[0][-1]
            uplot = int(min(nout[-3] + 3*min(nout[-2], avg_peak_dis), dataY.size))
            uplot = np.where(dataX <= uplot)[0][0]
            plt.clf()
            """
            plt.bar(dataX[lplot:uplot], dataY[lplot:uplot], width = 1, bottom = 1,
                    color = 'r', label = 'Raw', align = 'center')
            plt.bar(xdat, ydat, width = 1, bottom = 1, color = 'g', label = 'Used',
                    align = 'center')
            plt.plot(dataX[lplot:uplot], compFun11(nout, dataX[lplot:uplot]),
                     'k', linewidth = 5, label = 'Multi-Peak Fits')
            plt.legend(loc = 'upper center', prop = {'size': 40}, ncol = 2)
            save_fig(dataY = dataY, llim = lplot, ulim = uplot, path = path,
                     file = method + '_' + file, title = False, no_clr = True)
            """
            if Auto:
                break
            else:
                print('Allowing for user input')
                raise
        #except(TypeError):
        #    global trouble
        #    trouble = [xdat, ydat, low, upp, old_param, new_param, nout, avg_peak_dis]
        #    raise
        except:
            #raise
            if attempt == 0:
                print('Attempting prior values')
                try:
                    new_param = np.asarray(__prior[:])
                    low, upp = __range
                except NameError:
                    attempt += 1
                    print('No prior values stored...')
            if attempt != 0:
                Auto = False
                hold_on = False
                plt.close('all')
                plt.plot(xdat, ydat, 'o', markersize = 10)
                plt.plot(xdat, compFun11(new_param, xdat), 'k', linewidth = 5)
                print('Possible Error: {}'.format(sys.exc_info()))
                print('new_param:\n{}'.format(list(new_param)))
                print('Interval:\n{}'.format([dataX[low], dataX[upp]]))
                print('Method type is: {}'.format(method))
                plt.show()
                while True:
                    action = input('The multi-gaussian did not converge ' +
                                   'enter s to skip, i to change interval, ' +
                                   'p to change the parameters, ' +
                                   'h to toggle hold_on and change multiple parameters ' +
                                   'o to see the original parameters ' +
                                   'A to switch auto back on.\n')
                    if action.startswith('s'):
                        skip = True
                        break
                    if action.startswith('h'):
                        hold_on = True if not hold_on else False
                        print('hold_on: {}'.format(hold_on))
                    if action.startswith('o'):
                        print('The input parameters were: \n', list(old_param))
                    if action == 'raise':
                        raise
                    try:
                        if action.startswith('i'):
                            if action.endswith('f'):
                                forced = True
                                action = action[:-1]
                            low, upp = str2num(action[1:])
                            low = np.where(dataX <= low)[0][-1]
                            upp = np.where(upp <= dataX)[0][0]
                        elif action.startswith('p'):
                            temp_param = str2num(action[1:])
                            if len(temp_param)%3 != 0:
                                raise TypeError
                            else:
                                new_param = np.asarray(temp_param[:])
                        elif action.startswith('A'):
                            print('Continuing  Auto')
                            Auto = True
                            if hold_on:
                                hold_on = False
                                print('hold_on: {}'.format(hold_on))
                        else:
                            print('Rerunning with prior settings...')
                        if not hold_on:
                            break
                    except(ValueError, TypeError):
                        print('Error in edits to Multi-Fit: {}'.format(sys.exc_info()))
            attempt += 1
    if not skip:
        __prior = nout
        __range = [low, upp]
        if ALL:
            return temp_values, low, upp
        else:
            return temp_values
    else:
        return None

def noisespecta(Folder = None, Path = None, newPath = None, Verbose = False):
    Folder, Path, newPath = fileInfo(Folder, Path, newPath, Type = 'txt')
    for file in Folder:
        hold = loadFile(os.path.join(Path, file))
        for l, line in enumerate(hold):
            hold[l] = line.split(' ')
        hold[0] = ['Frequency (Hz)', 'dBV']
        save(data2str(hold), os.path.join(newPath, file), Verbose = Verbose)

def operation(data = None):
        if not data:
                data = impData()
                while True:
                        action = input('Would you like to add other files? Enter e to exit.\n')
                        if action == 'y':
                                temp0 = impData()
                                for key in temp0.keys():
                                        data.set(key, temp0.get(key))
                        if action == 'e':
                                del action
                                break
        Path = data.path
        folder  = data.keys()
        titles = [['Intensity', 'ID (eV)', 'Peaks (Ch)', 'Cent (Ch)','Ctot (Ch)',
                   'Csig (Ch)', 'Sig (Ch)', 'Stot (Ch)', 'Ssig (Ch)', 'Offset (Cnts)',
                   'Otot (Cnts)', 'Osig (Cnts)', 'Area (Cnts)', 'Atot (Cnts)',
                   'Asig (Cnts)']]
        for f, file in enumerate(folder[:]):
                y = np.asarray(data.get(file).get('Intensity'))
                x = np.asarray(channel(y))
                out = analyze(x, y)
                while True:
                        try:
                                print(file, out[0])
                                ID = str2list(input('Enter peak IDs\n'))
                                break
                        except:
                                print('Error', sys.exc_info())
                out.insert(0, ID)
                out.insert(0, list(y))
                out = reindex(titles + reindex(out, Save = True))
                for column in out:
                        data.get(file).set(column[0], column[1:])
                for key in data.get(file).keys():
                        if not isinstance(data.get(file).get(key), list):
                                data.get(file).list(key)
                data.get(file).save(file, Path = Path,
                                    Columns = titles[0],
                                    Verbose = False)
                print('Sucessfully fitted: {}'.format(folder.pop(0)))

def password(N = 14, special = False):
    options = string.ascii_letters + string.digits 
    if special:
	    options = options + string.punctuation 
    output_index = [np.random.random_integers(0, len(options) - 1) for a in range(N)]
    pw = ''
    for i in output_index:
	    pw = pw + options[i]
    return pw

def peak_select(dataX, dataY, peaks = np.asarray([]), file = ''):
    """
    This function is intended to allow the user to select the peaks from the data
    Parameters:    -----------
    dataX : list, array of the xaxis values
    dataY : list, array of the yaxis values
    peaks : <optional> list, array of suggested peaks
    file  : str name of the data
    Return:
    -------
    peaks : array of selected peaks
    """
    global trouble
    while True:
        trouble = peaks
        peaks = np.array(peaks)
        if peaks.size:
            print('Current peaks: ', list(dataX[peaks]))
        else:
            print('No current peaks: ', peaks)
        A = plt.plot(dataX, dataY, linewidth = 5)
        if peaks.size:
            vertical(peaks, dataY, dataX)
        A = plt.title(file, size = 30)
        plt.show()
        action = input('Enter new peaks or e to escape.\n')
        if action.startswith('e'):
            peaks.sort()
            del action
            break
        try:
            peaks = str2num(action)
            peaks = np.array([np.where(dataX <= val)[0][-1] for val in peaks])
        except(ValueError, TypeError):
            print('Error in peaks input try again: {}'.format(sys.exc_info()))
    return peaks
        ########################################################################
        # Function: Poisson distribution
######### P(n, mu) = exp(-mu)*mu**n/fac(n)
        # 
        ########################################################################
def poisson(n, mu):
    return np.exp(-mu)*mu**n/fac(n)


def pulse(t, t_off, amp, tau_d, tau_r):
    r"""
    This function is to be used for fitting pulse data from the GaGe card on the old Lab
    computer. The amplitude is 
    Parameters:
    -----------
    t     : np.ndarray or float/int value
    t_off : float/int value, the x-axis offset
    Area   : float/int value, the area of the pulse
    tau_d : float/int value, the decay time
    tau_r : float/int value, the rise time

    Return:
    -------
    hold : np.ndarray or float/int value
    Examples:
    
    Notes:
    ---------
    (Area*(tau_d + tau_r)/tau_d**2*(1 - np.exp(-(t-t_off)/tau_r))*np.exp(-(t-t_off)/tau_d))
    """
    t = np.array(t)
    amp = abs(amp)
    tau_d = abs(tau_d)
    tau_r = abs(tau_r)
    hold = amp*(1 - np.exp(-(t-t_off)/tau_r))*np.exp(-(t-t_off)/tau_d)
    hold[t < t_off] = 0
    return hold
    
################################################################################
def pulse_filter(yaxis, gap = 25, length = 55):
    """
        Takes the y values for a trace and determine the location of pulses for a
        given input height thresh:
        Example:
        taxis  = np.arange(16)
        yaxis  = np.array([0]*4 + [1]*8 + [0]*4)
        output = pulse_filter(yaxis, 5, 4)
        output: [0, 0, 0, 0, 0, 0, 0.75, 0.25, -0.25, -0.75, 0, 0, 0, 0, 0, 0]
        """

    if gap%2 == 1: #Gap is odd number of points
        pass
    else:
        print('Gap must be odd')
        raise
    pfilter = np.zeros(2*length + gap)
    pfilter[:length] = 1/length
    pfilter[-length:] = -1/length

    out = np.zeros(len(yaxis))
    if 500 < len(yaxis):
        out[length + int((gap-1)/2):-(length + int((gap-1)/2))] = np.convolve(yaxis, pfilter, 'valid')
    else:
        out[length + int((gap-1)/2):-(length + int((gap-1)/2))] = sig.fftconvolve(yaxis, pfilter, 'valid')
    return out

################################################################################

def pulse_Traces(files = None, Path = None, Record = np.array([]), marker = 0):
    warnings.simplefilter('ignore', np.RankWarning)
    if not Path:
        Path = input('Please enter path to files?\n\'C:\\\n')
    if not files:
        files = fileSel(Path, end = 'rdat')
    if isinstance(files, str):
        files = [files]
    for name in files:
        if not name.endswith('.rdat'):
            print('File {} is not in the rdat format'.format(name))
            raise TypeError
        outfile = dataFile()
        outfile.path = Path
        outfile.file = name.replace('.rdat', '_Traces_{}.txt'.format(index_gen(marker)))
        with open(os.path.join(Path, name), 'rb') as file:
            while True:
                line = file.readline().decode('utf-8')
                if line == '#End of Header\r\n':
                    break
                line = line.split(':')
                if line[0] == 'Timebase':
                    dt = float(line[1].split('\r')[0])*10**6
                elif line[0] == 'Presamples':
                    n  = int(line[1].split('\r')[0])
                elif line[0] == 'Total Samples':
                    N  = int(line[1].split('\r')[0])
                elif line[0] == 'Timestamp offset (s)':
                    t_0  = float(line[1].split('\r')[0])
                elif line[0] == 'Bits':
                    bits = int(line[1].split('\r')[0])
                elif line[0] == 'Range':
                    Range = float(line[1].split('\r')[0])
                elif line[0] == 'Inverted':
                    Polarity = -1 if line[1][0] == 'Y' else 1
                elif line[0] == 'Offset':
                    Offset = float(line[1].split('\r')[0])        
            t_axis = dt*np.linspace(-n, N-n -1, N)
            outfile.set('000_Time', t_axis)
            counter = 0
            while True:
                zero = np.fromfile(file, dtype = 'B', count = 2)
                if not zero.size:
                    break
                thyme = np.fromfile(file, dtype = '>I', count = 1)[0]
                data_or = np.fromfile(file, dtype = '<H', count = 1000)
                v_dat = voltage(data_or, bits, Range, Polarity, Offset)
                if counter in Record:
                    outfile.set('{}_Pulse'.format(index_gen(counter, 5)), v_dat)
                    Record = np.delete(Record, np.where(counter == Record))
                    if not Record.size:
                        break
                counter += 1
            outfile.save()
            
################################################################################
################################################################################
################################################################################

def pulse_Proc(files = None, Path = None, Other = False, Return = False, Record = False, Special = False, **keys):
    r"""
    This program is intended to process the pulse data from the GaGe card.

    Parameters:
    -----------
    files :   Name of file(s) to be reduced, if no file name(s) is(are) given
              then the user will be asked for the file(s) based on the Path
    Path :    Name of the folder where the file to be convereted is located, if
              no folder name is given the user will be asked for a folder
    **keys :  Are keywords for future use
              verbose: boolean sets the return, default is True

    Returns:
    --------
    Location and name of saved file.
    """
    warnings.simplefilter('ignore', np.RankWarning)
    if not Path:
        Path = input('Please enter path to files?\n\'C:\\\n')
    if not files:
        files = fileSel(Path, end = 'rdat')
    if isinstance(files, str):
        files = [files]
    try:
        verbose = keys['verbose']
    except(KeyError):
        verbose = True
    try:
        DEG = keys['DEG']
    except(KeyError):
        DEG = 1
    try:
        A_MAX = keys['A_MAX']
    except(KeyError):
        A_MAX = 0.8
    try:
        A_MIN = keys['A_MIN']
    except(KeyError):
        A_MIN = 0.4
    try:
        END = keys['END']
    except(KeyError):
        END = 0
    try:
        dig_fil = keys['filter']
    except(KeyError):
        dig_fil = False
    try:
        width = keys['width']
    except(KeyError):
        width = 10
    try:
        LFit = keys['LFit']
    except(KeyError):
        LFit = False
        
    for name in files:
        if not name.endswith('.rdat'):
            print('File {} is not in the rdat format'.format(name))
            raise TypeError
        newName = os.path.join(Path, name.replace('.rdat', '.txt'))
        if (Other or Special):
            tofile = [['Pulse #', 'time (s)', 'Amp_raw (V)', 'sig_raw (V)',
                       't_off (us)', 'sig_t (us)', 'Area (V*sec)', 'sig_A (V*sec)',
                       'tau_d (us)', 'sig_d (us)', 'tau_r (us)', 'sig_r (us)',
                       'b_off (V)', 'sig_b (V)']]
            if LFit:
                tofile[0] = tofile[0] + ['m_off (V/s)', 'sig_m (V/s)']
        else:
            tofile = [['Pulse #', 'Amp_raw (V)', 'sig_raw (V)']]
        with open(os.path.join(Path, name), 'rb') as file:
            while True:
                line = file.readline().decode('utf-8')
                if line == '#End of Header\r\n':
                    break
                line = line.split(':')
                if line[0] == 'Timebase':
                    dt = float(line[1].split('\r')[0])*10**6
                elif line[0] == 'Presamples':
                    n  = int(line[1].split('\r')[0])
                elif line[0] == 'Total Samples':
                    N  = int(line[1].split('\r')[0])
                elif line[0] == 'Timestamp offset (s)':
                    t_0  = float(line[1].split('\r')[0])
                elif line[0] == 'Bits':
                    bits = int(line[1].split('\r')[0])
                elif line[0] == 'Range':
                    Range = float(line[1].split('\r')[0])
                elif line[0] == 'Inverted':
                    Polarity = -1 if line[1][0] == 'Y' else 1
                elif line[0] == 'Offset':
                    Offset = float(line[1].split('\r')[0])        
            t_axis = dt*np.linspace(-n, N-n -1, N)
            counter = 0
            if LFit:
                def fun(t, t_off, Area, tau_d, tau_r, b_off, m_off):
                    return b_off + m_off*t + Polarity*pulse(t, t_off, Area, tau_d, tau_r)
            else:
                def fun(t, t_off, Area, tau_d, tau_r, b_off):
                    return b_off + Polarity*pulse(t, t_off, Area, tau_d, tau_r)
            while True:
                counter += 1
                zero = np.fromfile(file, dtype = 'B', count = 2) 
                if not zero.size:
                    break
                if counter > 1:
                    thyme = np.fromfile(file, dtype = '>I', count = 1)[0] - t_0
                else:
                    t_0 = np.fromfile(file, dtype = '>I', count = 1)[0]
                    thyme = 0                    
                data_or = np.fromfile(file, dtype = '<H', count = 1000)
                if dig_fil:
                    d_filter = gaussian(np.linspace(-2, 2, 4/dt + 1), 0, width)/np.sqrt(2*np.pi*width**2)
                    data = np.convolve(d_filter, data_or, 'valid')
                    data_or = data_or[int(len(d_filter)/2): -int(len(d_filter)/2)]
                    if counter == 1:
                        t_axis = t_axis[int(len(d_filter)/2): -int(len(d_filter)/2)]
                    if len(t_axis)!=len(data):
                        print(int(len(d_filter)/2), len(t_axis), len(data))
                        raise
                else:
                    data = data_or
                v_dat = voltage(data, bits, Range, Polarity, Offset)
                v_sig = sp.std(data[:int(0.8*n)])/2**(bits - 1)*Range
                b_off = np.average(v_dat[:int(0.8*n)])
                b_sig = sp.std(v_dat[:int(0.8*n)])
                p_dat = sg_filter(np.abs(v_dat - b_off), 11, 1)
                Ampli = np.max(np.abs(v_dat - b_off))
                if Other:
                    try:
                        Loc_A = np.where(v_dat.min() == v_dat)[0][0]
                        Loc_t = np.where(3*b_sig < p_dat)[0][0]
                        Loc_r = np.where((A_MIN - 0.2)*Ampli < p_dat[:Loc_A])[0][0]
                        Loc_R = np.where(A_MAX*Ampli < p_dat[:Loc_A])[0][0]
                        Loc_d = np.where(A_MAX*Ampli > p_dat[Loc_A:])[0][0] + Loc_A
                        Loc_D = np.where(A_MIN*Ampli > p_dat[Loc_A:])[0][0] + Loc_A

                        t_set = t_axis[Loc_t]
                        t_ris = np.polyfit(p_dat[Loc_r:Loc_R], t_axis[Loc_r:Loc_R], DEG)
                        t_ris = cfit(exp_rise, t_axis[Loc_r:Loc_R], p_dat[Loc_r:Loc_R],
                                     p0 = [t_set, abs(t_ris[0]), Ampli])[0]
                        t_dec = np.polyfit(v_dat[Loc_d:Loc_D], t_axis[Loc_d:Loc_D], DEG)
                        t_dec = cfit(exp_decay, t_axis[Loc_d:Loc_D], p_dat[Loc_d:Loc_D],
                                     p0 = [t_ris[0], abs(t_dec[0]), t_ris[2]])[0]
                        area = (sum(v_dat[Loc_t:Loc_D + END]*dt)
                                - b_off*(t_axis[Loc_D + END] - t_axis[Loc_t]))
                        initial = [t_set, area, abs(t_dec[1]), abs(t_ris[1]), b_off]
                        if LFit:
                            m_off = (v_dat[Loc_t] - v_dat[0])/(t_axis[Loc_t] - t_axis[0])
                            initial.append(m_off)
                        out, outcov= cfit(fun, t_axis[:Loc_D + END], v_dat[:Loc_D + END],
                                          p0 = initial,
                                          sigma = v_sig*np.ones(Loc_D + END),
                                          absolute_sigma = True)
                        vout = np.insert(np.insert(np.abs(out[1:4]), 3, out[4:]),
                                         0, out[0])
                        verr = np.sqrt(np.diag(outcov))
                        output = np.asarray(list(zip(vout, verr))).flatten()
                        # Hard cut on the pulse start to remove some false calculations
                        if abs(out[0]) > 3*abs(t_set):
                            raise
                        if Record:
                            pPath = os.path.join(Path, name.replace('.rdat', ''))
                            if not os.path.exists(pPath):
                                print('Making Folder: {}'.format(pPath))
                                os.makedirs(pPath)
                            psave = os.path.join(pPath,
                                                 name.replace('.rdat',
                                                              '_{}.png'.format(counter)))
                            plt.clf()
                            fig = plt.gcf()
                            fig.set_size_inches((40, 30))
                            plt.plot(t_axis, v_dat, 'o', markersize = 30,
                                     markeredgewidth = '3', markerfacecolor = (1, 1, 1))
                            plt.plot(t_axis[:Loc_D + END], v_dat[:Loc_D + END],
                                     'bo', markersize = 30)
                            plt.plot(t_axis, voltage(data_or, bits, Range, Polarity, Offset)
                                     , linewidth = 7)
                            plt.plot(t_axis, fun(t_axis, *out), linewidth = 10, color = 'g')
                            plt.plot(t_axis, fun(t_axis, *initial), linewidth = 10,
                                     color = 'r')
                            plt.plot(t_axis[Loc_r:Loc_R], v_dat[Loc_r:Loc_R], 'ks',
                                     markersize = 10)
                            plt.plot(t_axis[Loc_r:Loc_R],
                                     b_off - exp_rise(t_axis[Loc_r:Loc_R], *t_ris),
                                     'y-', linewidth = 20)
                            plt.plot(t_axis[Loc_d:Loc_D], v_dat[Loc_d:Loc_D], 'kD',
                                     markersize = 10)
                            plt.plot(t_axis[Loc_d:Loc_D],
                                     b_off - exp_decay(t_axis[Loc_d:Loc_D], *t_dec),
                                     'y-', linewidth = 20)
                            ax = plt.gca()
                            pylab.text(1, 0.75, 'Initial:\nt_off: {}\nb_off: {}\nArea: {}\ndecay: {}\nrise : {}'.format(*initial),
                                       transform = ax.transAxes, fontsize = 50,
                                       horizontalalignment='right',
                                       verticalalignment='top')
                            pylab.text(1, 0, 'Fit:\nt_off: {}\nb_off: {}\nArea: {}\ndecay: {}\nrise : {}'.format(*out),
                                       transform = ax.transAxes, fontsize = 50,
                                       horizontalalignment='right',
                                       verticalalignment='bottom')
                            plt.ylabel('Pulse (V)', size = 100)
                            plt.xlabel(r'Time ($\mu$s)', size = 100)
                            plt.title('Pulse #{}'.format(counter), size = 100)
                            plt.tick_params(labelsize = 100, width = 5, size = 10)
                            fig.savefig(psave, dpi = 100)
                    except:
#                        global A
#                        global B
#                        global C
#                        global D
#                        A = p_dat
#                        B = t_axis
#                        C = b_sig
#                        D = v_dat
#                        raise
                        output = np.zeros(10 + LFit*2)
                    output = np.insert(output, 9, b_sig)
                    output = np.insert(output, 0, Ampli)
                    output = np.insert(output, 0, thyme)
                    output = list(np.insert(output, 0, counter))
                elif Special:
                    try:
                        x = t_axis
                        y = v_dat
                        p1 = np.where(y.min() == y)[0][0]
                        p0 = np.where(y[:np.where(0 == x)[0][0]].max() == y[:np.where(0 == x)[0][0]])[0][-1]
                        if p0 >= 0.5*np.where(0 == x)[0][0]:
                            p = np.polyfit(x[:p0], y[:p0], 1)
                            y_bg = y - np.polyval(p, x)
                            out, outcov = cfit(fun, x[:p1], y_bg[:p1], p0 = [x[np.where(0 == x)[0][0]], 100, 10, 100, 0])
                            vout = np.insert(np.insert(np.abs(out[1:4]), 3, out[4:]),
                                             0, out[0])
                            verr = np.sqrt(np.diag(outcov))
                            output = np.asarray(list(zip(vout, verr))).flatten()
                        b_off = np.average(y_bg[:p0])
                        Ampli = np.max(np.abs(y_bg- b_off))
                    except:
                        output = np.zeros(10 + LFit*2)
                    output = np.insert(output, 9, b_sig)
                    output = np.insert(output, 0, Ampli)
                    output = np.insert(output, 0, thyme)
                    output = list(np.insert(output, 0, counter))                        
                else:
                    output = [counter, Ampli, b_sig]
                    if Record:
                        pPath = os.path.join(Path, name.replace('.rdat', ''))
                        if not os.path.exists(pPath):
                            print('Making Folder: {}'.format(pPath))
                            os.makedirs(pPath)
                        psave = os.path.join(pPath, name.replace('.rdat',
                                                                 '_{}.png'.format(counter)))
                        plt.clf()
                        fig = plt.gcf()
                        fig.set_size_inches((40, 30))
#                        plt.plot(t_axis, data_or, 'o', markersize = 30,
#                                 markeredgewidth = '3', markerfacecolor = (1, 1, 1))
                        plt.plot(t_axis, voltage(data, bits, Range, Polarity, Offset),
                                 'o', markersize = 30, markeredgewidth = '3',
                                 markerfacecolor = (1, 1, 1))
                        ax = plt.gca()
                        plt.ylabel('Pulse (V)', size = 100)
                        plt.xlabel(r'Time ($\mu$s)', size = 100)
                        plt.title('Pulse #{}'.format(counter), size = 100)
                        plt.tick_params(labelsize = 100, width = 5, size = 10)
                        fig.savefig(psave, dpi = 100)
                tofile.append(output)
            save(data2str(tofile), newName, Path, verbose)
    plt.clf()
    if Return:
        return tofile       
################################################################################
################################################################################
################################################################################

def quadratic(x, p0, p1, p2):
    """
    Computes the value of the quadratic:
    y = p0*x**2 + p1*x + p0
    identical to np.polyval(p, x), but is used for computing the calibration.
    Parameters:
    -----------
    x : array whose values are to be calculated
    p0: int or float, coefficent of square term
    p1: int or float, coefficent of linear term
    p2: int or float, constant

    Returns:
    --------
    the value at x for the given quadratic function.
    """
    return p0*x**2 + p1*x + p2
    
################################################################################
################################################################################
################################################################################

def rdat(files = None, Path = None, newPath = None, num_pulse = 0, **keys):
    r"""
    This program takes in an LJHMFF and reduces the size of the file by 50% by
    throwing away the unused channel 0 data.

    Parameters:
    -----------
    files :    Name of file(s) to be reduced, if no file name(s) is(are) given
              then the user will be asked for the file(s) based on the Path
    Path :    Name of the folder where the file to be convereted is located, if
              no folder name is given the user will be asked for a folder
    newPath : Name of new Path where the newFile will be stored, if none is
              given then the newPath will transfer to Data-Processed folder
    **keys :  Are keywords for future use
              verbose: boolean sets the return, default is True

    Returns:
    --------
    Location and name of saved file.
    """
    if not Path:
        Path = input('Please enter path to files?\n\'C:\\\n')
    if not files:
        files = fileSel(Path, 'dat')
    if isinstance(files, str):
        files = [files]
    if not newPath:
        newPath = Path #.replace('Original', 'Processed')
    if not os.path.exists(newPath):
        os.mkdir(newPath)
    if any('verbose' == el for el in keys):
        verbose == keys['verbose']
    else:
        verbose = True
    for name in files:    
        newName = os.path.join(newPath, name.replace('.dat', '.rdat'))
        file = open(os.path.join(Path, name), 'rb')
        generate = True
        current_pulse = 0
        if os.path.exists(os.path.join(newPath, newName)):
            generate = input('File already exists. Enter o to overwrite.\n')
            if generate == 'o':
                os.remove(os.path.join(newPath, newName))
            else:
                generate = False
        if generate:
            with open(newName, 'xb') as newFile:
                while True:
                    line = file.readline()
                    if line.decode('utf-8').split(':')[0] == 'Total Samples':
                        num_samples = 2*int(line.decode('utf-8').split(':')[1])
                    elif line.decode('utf-8') == 'Number of Active Channels: 2\r\n':
                        line = line.decode('utf-8').replace('2', '1').encode('utf-8')
                    elif line.decode('utf-8') == 'Channel: 1.0\r\n':
                        trash = line.decode('utf-8')
                        counter = 0
                        while trash != 'Channel: 1.1\r\n':
                            counter += 1
                            line = file.readline()
                            trash = line.decode('utf-8')
                            if counter == 10:
                                break
                    header = newFile.write(bytes(line))
                    if line.decode('utf-8') == '#End of Header\r\n':
                        break
                while True:
                    trash = file.read(num_samples + 6)
                    pulse = file.read(num_samples + 6)
                    current_pulse += 1
                    if trash == b'' or pulse == b'':
                        print(current_pulse)
                        break
                    pulse = newFile.write(bytes(pulse))
                    if current_pulse == num_pulse:
                        print(current_pulse, num_pulse)
                        break
                    
        if verbose:
            print('Saved: {} in {}'.format(name.replace('.dat', '.rdat'), newPath))
    
def rebin(data, N = 1):
    r"""
    This program will rebin data by combining the adjacent N bins to the left

    Parameters:
    -----------
    data : list or array, to be rebinned into a smaller set
    N    : int, the number of bins to be joined

    Return:
    -------
    out  : list, reduced data set

    Examples:
    --------
    """
    return [sum(data[i*N:(i+1)*N]) for i in range(int(len(data)/N))]

def reindex(data, Save = False):
    """
    This will reindex a list from data_a_b to data_b_a
    Input list should look like:
      a  =      0     1     2     3
    data = [ [ [,] , [,] , [,] , [,], ... [,]] , [...] , ..., [...] ]
      b  =               b_0                      b_1          b_i
    Need to convert this to the form:
      b  =     b_0   b_1   b_2   b_3      b_i
    data = [ [ [,] , [,] , [,] , [,], ... [,]] , [...] , [...] , [...] ]
      a  =                  0                      1       2       3
    """
    if Save:
        maxlen = max([len(data[d]) for d in range(len(data))])
        for c, cluster in enumerate(data):
            data[c] = data[c] + ['']
            while len(data[c]) < maxlen:
                data[c].append('')
    if isinstance(data, list):
        return [list(line) for line in zip(*data)]
    if isinstance(data, np.ndarray):
        return np.array([np.array(list(line)) for line in zip(*data)])


def removepts(data, rem = None, ran = False, ret = False, spe = None):
    r"""Used to remove points from a list or an array.
    Parameters
    ----------
    data : list or array
        the list or array from which entries are to be removed.    
    rem : list or array
        the index of entries to remove from data.
    ran : boolean
        default is False so rem is interpreted as index values. When True the rem
        entries are interpreted as an interval of values to remove.
    ret: boolean
        Return the index values of the removed entires
    spe: int, float, list, array
        Specific values of entries to remove (removes all entries)
    
    Returns
    -------
    data : ndarray, shape (N)
        Returns the data structure with the removed elements.
    rem  :
        List of index of removed entries.

    Notes
    -----

    
    Examples:

    """
    
    if not rem and isinstance(spe, (int, float)):
        spe = str(spe)
        rem = [k for k, item in enumerate(data)
               if item in [ -float('Inf'), float('Inf')]
               or  str(item) == spe or math.isnan(item)]
    elif not rem and isinstance(spe, (list, np.ndarray)):
        if (isinstance(data[0], float) or np.issubdtype(data[0], float)) \
           and (isinstance(spe[0], int) or np.issubdtype(spe[0], int)):
            spe = list(map(str, map(float, spe)))
        elif (isinstance(data[0], int) or np.issubdtype(data[0], int)) \
             and (isinstance(spe[0], float) or np.issubdtype(spe[0], float)):
            spe = list(map(str, map(int, spe)))
        else:
            spe = list(map(str, spe))
        rem = [k for k, item in enumerate(data) if str(item) in spe]
    elif ran:
        rem = np.linspace(rem[0], rem[1],rem[1] - rem[0] + 1)
    rem = list(map(int, rem))
    data = np.delete(data[:], rem)
    if ret == False:
        del rem
        return data
    if ret == True:
        return data , rem

def residual(p = [], dp = [], x_1 = [], dx_1 = [], x_a = [], dx_a = []):
    """
    Calculates the residuals of a line with  known errorbars
    """
    p = np.asarray(p)
    dp = np.asarray(dp)
    x_1 = np.asarray(x_1)
    dx_1 = np.asarray(dx_1)
    x_a = np.asarray(x_a)
    dx_a = np.asarray(dx_a)
    if not dp.size:
        dp = np.zeros(p.size)
    if not dx_1.size:
        dx_1 = np.zeros(x_1.size)
    if not dx_a.size:
        dx_a = np.zeros(x_a.size)

    res = np.asarray(x_a) - np.polyval(p, x_1)
    dres = np.sqrt(dx_a**2 + p[0]**2*dx_1**2 + dp[0]**2*x_1**2 + dp[1]**2)
    return res, dres

def revise_energy(oldX, newX, initial = 5):
    """
    This is intended to generate the matrix that converts from oldX to newX for rebin data
    """
    odX = (oldX[1] - oldX[0])/2
    ndX = (newX[1] - newX[0])/2
    """
    frac = ndX/odX
    if frac >= 1:
        i_th = np.where(newX[0] - ndX < oldX)[0][0]
        ep_sta = [((oldX[i_th] + odX) - (newX[0] - ndX))/(2*odX)]
        [ep_sta.append(1 - (frac - ep_sta[-1])%1) for a in range(len(newX[1:]))]
        compl = [(frac - el)%1 for el in ep_sta]
        alpha = np.asarray([])
        pre = list(np.zeros(i_th))
        
        center = list(np.ones(int(frac - ep_sta)))
        
        post = list(np.zeros(len(newX) - i_th + 2 + len(center)))
    """  
    oldX = np.insert(oldX + odX, 0, oldX[0] - odX)
    newX = np.insert(newX + ndX, 0, newX[0] - ndX)
    alpha = []
    j = 0
    while True:
        alpha.append(np.asarray([alpha_ij(oldX[i], oldX[i+1], newX[j], newX[j+1])
                                 for i in range(oldX.size - 1)]))
        j += 1
        if np.where(0 < alpha[-1])[0].size:
            initial = j
            start = np.where(0 < alpha[-1])[0][0]
            break
        if j >= initial:
            initial += 5
            
    constant = sum(alpha[-1])
    for j in range(initial, newX.size - 1):
        pre_a = np.zeros(start)
        hold = []
        current = constant
        for i in range(start, oldX.size - 1):
            latest = round(alpha_ij(oldX[i], oldX[i +1], newX[j], newX[j+1]), 5)
            hold.append(min(latest, current))
            current = max(current - latest, 0)
            if hold[-1] == 0 and current == 0:
                break
        pre_a = np.insert(pre_a, pre_a.size, np.asarray(hold))
        pre_a = np.insert(pre_a, pre_a.size, np.zeros(oldX.size - pre_a.size - 1))
        alpha = alpha + [pre_a]
        start = np.where(0 < pre_a)[0][0]
    alpha = np.asarray([np.asarray(line) for line in alpha])
    return alpha 

        ########################################################################
######### Saves data to file
        ########################################################################
def save(data, fileName = None, Path = None, Verbose = False):
    if Path == None:
        Path = 'C:\\Data-Processed\\' + datetime.date.isoformat(datetime.date.fromtimestamp(time.time()))
    if not os.path.exists(Path):
        os.makedirs(Path)
    if fileName == None:
        counter = 0
        fileName = 'Processed_data_{}.txt'.format(counter)
        while os.path.exists(os.path.join(Path, fileName)):
            counter += 1
            fileName = 'Processed_data_{}.txt'.format(counter)
    
    with open(os.path.join(Path, fileName), 'w') as saveFile:     
        if isinstance(data, list):
            for element in data:
                saveFile.write(element)
        elif isinstance(data, (int, float, str)):
            saveFile.write(str(data))
    if Verbose:
         print('Saved as: {}'.format(Path, fileName), end = '\r')

def save_fig(dataX = None, dataY = None, yfit = None, energy = None, peak = None, llim = None, ulim = None, lplot = None, uplot = None, path = None, file = None, title = True, no_clr = False):
    """
    Save the figure of data over the selected range
    """
    figpath = os.path.join(os.path.join(path, 'Peak_Fits'),
                           file.replace('LS_', '').replace('ODR_', ''))
    if not os.path.exists(figpath):
        os.makedirs(figpath)
    
    if no_clr:
        ax1 = plt.gca()
        fsize = 50
        pdffile = os.path.join(figpath, 'Multi_Peak_{}.pdf'.format(file))
    else:
        pdffile = os.path.join(figpath, '{}_Peak_{}Ch.pdf'.format(file, peak))
        plt.clf()
        ax1 = plt.gca()
        fsize = 50
        ax1.bar(dataX[llim : ulim], dataY[llim : ulim], width = 1,
                bottom = 1, color = 'r', label = 'Raw')
        ax1.bar(dataX[lplot : uplot], dataY[lplot : uplot], width = 1,
                bottom = 1, color = 'g', label = 'Used')
        ax1.plot(dataX[llim : ulim], yfit[llim : ulim], 'k', linewidth = 7)
    ax1.tick_params(which = 'major', length = 10, width = 3)
    ax1.tick_params(which = 'minor', length = 7)
    ax1.set_ylabel('Counts', size = fsize)
    ax1.set_ylim((max(1, min(dataY[llim : ulim])), 10*max(dataY[llim : ulim])))
    ax1.set_xlim((llim, ulim))
    ax1.set_yscale('log')
    ax1.tick_params(labelsize = fsize)
    ax1.set_xlabel('Channel (Ch)', size = fsize)
    if no_clr:
        ax1.set_title(file, size = fsize)
    else:
        ax1.legend(loc = 'upper center', ncol = 3, prop = {'size': fsize})
        ax1.set_title('Peak @ {}eV'.format(energy), size = fsize)
    fig = plt.gcf()
    fig.set_size_inches((24, 18))
    fig.savefig(pdffile, dpi = 200, format = 'pdf')
    plt.clf()

def secondDer(dataX, dataY, n = None):
    r"""Takes the second derivative of the first n points using the forward
    difference method, the middle is calculated by fitting a quadratic to the 2n+1
    points centered on the middle point, the last n point are calculated using the
    backward difference method.

    """
    dataY = array(dataY)
    dataX = array(dataX)
    if not n:
        n = int(0.001*len(dataY))
        if n < 3:
            n = 3
    m = dataY.size
    secDer = (dataY[:n] - 2*dataY[1:n+1] + dataY[2: n+2])/(dataX[:n] - dataX[1:n+1])**2
    for a in range(n, m-n):
        pfit = np.polyfit(dataX[a:a+n], dataY[a:a+n], deg =2)
        secDer = np.insert(secDer, secDer.size, np.asarray([2*pfit[0]]))
    secDer = np.insert(secDer, secDer.size,
                       (dataY[m-n-2:m-2] - 2*dataY[m-n-1:m-1] + dataY[m-n:]
                        )/(dataX[m-n-1:m-1] - dataX[m-n:])**2)
    return secDer
                       
def sg_filter(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except (ValueError):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * fac(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

        ########################################################################
######### Rounds the input value to n significant figures
        ########################################################################
def sigFigN(x, n = 4):
    xPrime = np.abs(x)-np.floor(np.abs(x))
    if x != 0 and np.floor(np.abs(x)) == 0:
        m = -int(np.floor(np.log10(np.abs(x)/10**(n-1))))
    elif x != 0 and np.floor(np.abs(x)) != 0 and xPrime != 0:
        m = -int(np.floor(np.log10(xPrime/10**(n-2))))        
    else:
        m = 0
    return round(x, m)

        ########################################################################
        # Opens a handle to the file named fileName, returns the data as a 
######### string list, and finally closes the file handle
        # 
        ########################################################################
def spe2txt(Folder = None, Path = None, newPath = None, Type = 'Spe', Reduce = 1024, Verbose = False):

    Folder, Path, newPath = fileInfo(Folder, Path, newPath, Type = Type)
    overwrite = False
    nowrite = False
    identical = False
    print('Reducing bins to: {}'.format(Reduce))
    for file in Folder:
        dataStr = loadFile(os.path.join(Path, file))
        dataNum = modifyData(dataStr, Out = False)
        if Type.endswith('mca'):
            dataNum = dataNum[1:]
        if not identical:
            start = 0
            end = len(dataNum) + 1
            while True:
                ch = channel(dataNum)
                plt.plot(ch, dataNum)
                plt.plot(ch[start:end], dataNum[start:end], 'o', markersize = 5)
                plt.title('Identify start and end of data', size = 30)
                plt.show()
                action = input('Enter data range, e to exit, ey to exit and' +
                               '  set all remaining files as identical.\n')
                if action.startswith('e'):
                    if action == 'ey':
                        identical = True
                    break
                try:
                    start, end = str2num(action)
                except(ValueError, TypeError):
                    pass
        dataRed = dataNum[start: end + Reduce - len(dataNum[start:end])%Reduce]
        if Reduce != 1 and len(dataRed) != Reduce:
            try:
                dataRed = np.asarray(rebin(dataRed, int(len(dataRed)/Reduce)))
            except(ValueError):
                print('Could not reduce {}, reduce value invalid'.format(file))            
        dataRed = np.insert(dataRed.astype('str'), 0, 'Intensity')
        if not os.path.exists(os.path.join(newPath, file.replace(Type, 'txt'))):
            save(data2str(dataRed), file.replace(Type, 'txt'), newPath, Verbose = Verbose)
        else:
            while True:
                if not overwrite and not nowrite:
                    action = input(os.path.join(newPath, file.replace(Type, 'txt')) +
                                ' already exists, y to overwrite this file,' +
                                   ' ya to overwrite all, n to skip,' +
                                   ' or na to skip all copies.\n')
                if overwrite or action == 'ya':
                    save(data2str(dataRed), file.replace(Type, 'txt'),
                         newPath, Verbose = Verbose)
                    overwrite = True
                    break
                elif nowrite or action == 'na':
                    print('Did not convert file: {}'.format(file.replace(Type, 'txt')))
                    nowrite = True
                    break
                elif 'y' == action:
                    save(data2str(dataRed), file.replace(Type, 'txt'),
                         newPath, Verbose = Verbose)
                    break
                elif 'n' == action:
                    print('Did not convert file: {}'.format(file.replace(Type, 'txt')))
                    break

def standard_fig_save(file, path = None, height = 30, width = 30, dpi = 500, ext = 'pdf', fig = None, Verbose = True, binding = 'tight'):
    if not path:
        path = 'C:\\Presentations\\Figures'
    pdffile = os.path.join(path, file)
    if fig == None:
        fig = plt.gcf()
    fig.set_size_inches((width, height))
    fig.savefig(pdffile, dpi = dpi, format = ext, bbox_inches = binding)
    if Verbose:
        print('Saved figure:\n{}'.format(pdffile))

def stj_stat(E, elsig = None, delta = None):
	if not delta:
		delta = 0.7*10**-3
	if elsig == None:
		elsig = 0
	else:
		elsig = np.average(elsig)
	fano = 0.195
	eps = 1.68*delta
	stat = np.sqrt(eps*E*(fano + 1))
	return 2*np.sqrt(2*np.log(2)*(stat**2 + elsig**2))

def str2list(string):
    if not string.startswith('[') or not string.endswith(']'):
        raise(ValueError)
    hold = string[1:-1].split(',')
    return hold    

################################################################################
# This program will convert each string data from file of the form '[#, #, #]' 
# to its numerical counterpart of [#, #, #]
################################################################################
def str2num(string):
    hold = str2list(string)
    try:
        hold = list(map(int, hold))
    except(ValueError):
        hold = list(map(float, hold))
    return hold

def time_stamp():
    folder, path, new = fileInfo(New = False)
    folder.sort()
    hold = []
    for file in folder:
        with open(os.path.join(path, file)) as single:
            temp = single.readlines()[7]
            month, day, year = np.asarray(temp.split(' ')[0].split('/')).astype(int)
            hr, minute, sec = np.asarray(temp.split(' ')[1].split(':')).astype(int)
            thyme = datetime.datetime.timestamp(datetime.datetime(year, month, day, hr, minute, sec, 0))
            hold.append([file, thyme])
    hold.sort(key = lambda line: line[1])
    hold = reindex(hold)
    hold[1] = np.round((np.asarray(hold[1]) - hold[1][0])/3600, 4)
    new = dataFile()
    new.set('File', hold[0])
    new.set('Time (Hrs)', hold[1])
    new.path = path
    new.save()

def vertical(index, dataY, dataX = np.asarray([])):
    r"""
    This program will plot the location of peaks so that the user can see were
    the peaks are located.
    Parameters:
    -----------
    index : int or list of peak locations
    dataY : list or array of data being plotted

    Returns:
    Plotted vertical lines for the data
    """
    if isinstance(index, (list, np.ndarray)):
        if dataX.size:
            [plt.vlines(dataX[peak], 0, dataY[peak], 'k', linewidth = 3) for peak in index]
        else:
            [plt.vlines(peak, 0, dataY[peak], 'k', linewidth = 3) for peak in index]
    else:
        plt.vlines(int(index), 0, dataY[int(index)], 'k', linewidth = 3)

################################################################################
# Convert the byte data of .dat or .rdat file to pulse voltage
################################################################################
def voltage(sample, bits = 12, Range = 2, Polarity = -1, Offset = 0):
    return Offset + (sample/(2**(bits-1)) -1)*Range*Polarity
def voltage_e(dsample, bits = 12, Range = 2, Polarity = -1, Offset = 0):
    return (dsample/2**(bits-1))*Range*Polarity

def weightedvalues(value, weight):
    average = sum(value*weight)/sum(weight)
    std = np.sqrt(value.size*sum(weight*(value - average)**2)/((value.size - 1)*sum(weight)))
    return average, std

def xia_bias(Files = None, Path = None):
    files, path, newPath = fileInfo(Files, Path, newPath = None, Type = 'adc')
    overwrite = False
    nowrite = False
    for name in files:
        if not name.endswith('.adc'):
            print('File {} is not in the adc format'.format(name))
            raise TypeError
        newName = os.path.join(path, name.replace('.adc', '.txt'))
        host = dataFile(['Bias (uV)', 'Current (A)'])
        host.path = newPath
        host.file = newName
        global trouble
        trouble = []
        with open(os.path.join(path, name), 'r') as file:
            counter = 0
            while True:
                counter += 1
                line = file.readline()
                if line == '\n' or line == '':
                    file.readline()
                    break
                line = line.split(' = ')
                if line[0] == 'bias_scan_step_size':
                    step = float(line[1])
                elif line[0] == 'bias_scan_steps':
                    number  = int(line[1])
                elif line[0] == 'bias_scan_start_offset':
                    start = float(line[1])
            leakage = []
            while True:
                line = file.readline()
                if line == '\n' or line == '':
                    break
                else:
                    line = float(line.split('\t')[0])
                leakage.append(line)
            bias = start + np.arange(number)*step
            host.set('Bias (uV)', bias)
            host.set('Current (A)', leakage)
            trouble = host
            host.save(Verbose = False)

def xia_trace(Files = None, Path = None):
    files, path, newPath = fileInfo(Files, Path, newPath = None, Type = 'adc')
    overwrite = False
    nowrite = False
    for name in files:
        if not name.endswith('.adc'):
            print('File {} is not in the adc format'.format(name))
            raise TypeError
        newName = os.path.join(path, name.replace('.adc', '.txt'))
        host = dataFile(['Time (us)', 'Intensity (au)'])
        host.path = newPath
        host.file = newName
        global trouble
        trouble = []
        with open(os.path.join(path, name), 'r') as file:
            counter = 0
            while True:
                counter += 1
                line = file.readline()
                if line == '\n' or line == '':
                    number = int(file.readline())
                    break
                line = line.split(' = ')
                if line[0] == 'Sampling Interval':
                    step = float(line[1])
            intensity = []
            while True:
                line = file.readline()
                if line == '\n' or line == '':
                    break
                else:
                    line = float(line.split('\t')[0])
                intensity.append(line)
            host.set('Time (us)', np.arange(1, number + 1)*step)
            host.set('Intensity (au)', intensity)
            trouble = host
            host.save(Verbose = False)
