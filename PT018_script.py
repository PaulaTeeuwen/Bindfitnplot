'''
This script is used to calculate the maxima in the fluorescence peaks
obtained in a titration experiment for Host-Guest analysis. By applying a
smoothing operator the maximum wavelength and the corresponding
absorbance are obtained for different concentrations of host and guest.
From this a bindingisotherm is obtained. 
This program is applicable for:
* Fluorescence titrations for Jasco FP-8300, saved all spectra as csv files
* The fluorescence spectrum has 2 peaks (adjust yourself if different number)

P.C.P. Teeuwen 
p.teeuwen@student.ru.nl
Last updated: 17-05-2019
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.signal import find_peaks
from scipy.signal import savgol_filter # Used to smooth out the spectra and remove noise, gives back smoothened version with same length array.
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.optimize import root

import seaborn as sns
sns.set()
current_palette = sns.color_palette()


print('')
print('')

os.getcwd()
os.chdir('E:\\Uni\\Bachelor_Internship\\Internship\\Analysis_techniques') # Changes current working directory

##################################################
# Data specifiers

Expnr = 'PT018'
Host = '(+)-H22'
#Host = '(+)-H\u20822'
Guest = '(S,S)-V2'

T = 25 #'\u00b0' +'C '

# Host-stock 1
M_Host = 1435.48 # g/mol   (+)-H23
m_Hstock1 = 1.37 * 10**(-3) #g  The weighed mass of the Host
V_solv_Hstock = 20.00  # mL +- 0.04 The volume of the relevant volumetric flask
m_solv_Hstock1 = 22.76217   # g The weighed mass of the solvent + Host

# Host-stock 2
m_Hstock2 = 1.44 * 10**(-3) #g  The weighed mass of the Host (Held at least 1h under vacuum)
m_solv_Hstock2 = 22.71919   # g The weighed mass of the solvent + Host

# Host-stock 3
m_Hstock3 = 1.54 * 10**(-3) #g  The weighed mass of the Host (Held at least 1h under vacuum, different batch,fibrous)
m_solv_Hstock3 = 22.76069   # g The weighed mass of the solvent + Host

# Guest-stock 1
M_Guest = 728.66 # g/mol   (S,S)-V2
m_Gstock1 = 1.64 * 10**(-3) # g The weighed mass of the Guest
V_solv_Gstock = 10.00 # mL +- 0.025 mL ; The volume of the relevant volumetric flask
m_solv_Gstock1 = 11.36163   # g The weighed mass of the solvent + Guest

# Guest-stock 2
m_Gstock2 = 1.53 * 10**(-3) # g The weighed mass of the Guest
m_solv_Gstock2 = 11.35331  # g The weighed mass of the solvent + Guest

# Guest-stock 3
m_Gstock3 = 1.45 * 10**(-3) # g The weighed mass of the Guest
m_solv_Gstock3 = 11.39119  # g The weighed mass of the solvent + Guest

# Host-measurement
V_solv_Hmeasure = 10.00 # mL +- 0.025 mL ; 
m_Hstock_diluted1 = 0.34025 # g The weighed Host-stock before dilution, with known density calculated earlier
m_Hstock_diluted2 = 0.34260# g
m_Hstock_diluted3 = 0.33671 # g
V_H_titrated = 2.000 # mL The volume of Host measurement solution to which the Guest is titrated. 
# etc..

# Guest-measurement
V_solv_Gmeasure = 5.00 # mL +- 0.025 mL ; 
m_Gstock_diluted1 = 1.28362 # g The weighed Guest-stock before dilution, with known density calculated earlier
m_Gstock_diluted2 = 1.27152 # g
m_Gstock_diluted3 = 1.14782 # g
m_Hstock_diluted_forGmeasurement1 = 0.15911 # g The weighed Host-stock before dilution that is added to the flask
m_Hstock_diluted_forGmeasurement2 = 0.17257 # g
m_Hstock_diluted_forGmeasurement3 = 0.17803 # g
# etc...

# Measurement
n = 3 #triplo
m1 = 32 #amount of measurements
m2 = 32
m3 = 32
#etc.. : m.. = .. 

# Follow around wavelength
lambda1 = 650 #nm
lambda2 = 715 #nm

# What kind of fitting?
F0_free = True # keep this one on True, check is implemented in script
smoothening = True # True if you want to apply a smoothing operator, False if not (keep in mind the error it causes)
formula = 'binding_full' #Choose from binding, binding_alt, binding_full, binding_HGG, binding_HHG
method = 'single_wl' #Choose from 'single_wl' and 'max_wl'


# Guest added in each measurement-----------------------------------------------

for i in np.arange(1,n+1):
    name = 'Guest_add_' + str(i)
    m = 'm'+ str(i)
    vars()[name] = np.array([0,0,10,10,10,10,10,10,10,10,10,10,30,30,30,30,30])  # uL
    endvalues = np.full(vars()[m] - len(vars()[name]),50) # This len = 17 is how many data points until adding 50's
    vars()[name] = np.concatenate(((vars()[name]), endvalues))

# Allow for changes in added guest as compared to planned values----------------
# Format: Guest_add_i[measurement number -1]



##################################################



## 1 Loading data

# 1)Load the data...............................................................
def splitdata(fluorescence, n): #fluorescence is a string
    '''
    This function loads the fluorescence data stored in the current directory.
    The important columns (wavelength, absorbance) are separated by a comma-delimiter.
    '''
    global wavelength # to test output of last country in loop
    global absorbance
    os.chdir('E:\\Uni\\Bachelor_Internship\\Internship\\Analysis_techniques\\Fluorescence\\' + Expnr + '\\set' + str(n))
    wavelength, absorbance = np.genfromtxt(fluorescence + '.csv', unpack = True, delimiter = ',', skip_header = 19, skip_footer = 45, usecols = (0,1))
    return np.array([wavelength , absorbance])
    
# 1a) Define function that will be fitted --------------------------------------

bindingconstant = np.zeros([2,n])
initialF = np.zeros([2,n])
kDHGcollect = np.zeros([2,n])
kHcollect = np.zeros([2,n])
kH0collect = np.zeros([2,n])
kHGcollect = np.zeros([2,n])



def binding(G0, H0, F0, Ka):  # Complex fluorescently silent (quenched) and no dynamic quenching
    #return F0 / (1+Ka*(1/2)*((G0 - H0 - 1/(Ka) + np.sqrt((G0 - H0 - 1/Ka)**2 + 4*G0/Ka))))
    
    #Don't always know whether +/- in the formula: try both. If one doesn't work, try other
    try:
       return F0 / (1+Ka*(1/2)*((G0 - H0 - 1/Ka) + np.sqrt((G0 - H0 - 1/Ka)**2 + 4*G0/Ka)))
    except ZeroDivisionError:
        print('yess')
        return F0 / (1+Ka*(1/2)*((G0 - H0 - 1/Ka) - np.sqrt((G0 - H0 - 1/Ka)**2 + 4*G0/Ka)))
    # based on F0/F = 1+Ka[G]
    #Reference: Thordarson 2010, formula 42

def binding_alt(G0, H0, F0, Ka, kDHG): # No dynamic quenching of host but complex still fluorescently active
    #return F0 + YDHG*(1/2)*(((1/Ka) + H0 + G0) - np.sqrt((G0 + H0 + 1/Ka)**2 - 4*G0*H0))
    
    #Don't always know whether +/- in the formula: try both. If one doesn't work, try other
    try:
       return F0 + kDHG*(1/2)*(((1/Ka) + H0 + G0) - np.sqrt((G0 + H0 + 1/Ka)**2 - 4*G0*H0)) # last + has to be a - 
    except ZeroDivisionError:
        print('yess')
        return F0 + kDHG*(1/2)*(((1/Ka) + H0 + G0) + np.sqrt((G0 + H0 + 1/Ka)**2 - 4*G0*H0))
    
    #Based on DY = Y_DHG([HG])
    #Reference: Thordarson 2010, formula 43
    #Reference: https://pubs.acs.org/doi/suppl/10.1021/jo990112c/suppl_file/jo990112c_s.pdf

def binding_full(G0, H0, F0, Ka, kH, kH0, kHG): #Dynamic quenching possible and complex still fluorescently silent
    G = (1/2)*((G0 - H0 - 1/Ka) + np.sqrt((G0 - H0 - 1/Ka)**2 + 4*G0/Ka))
    return F0 * (kH/kH0 + (kHG/kH0)*Ka*G)/(1+Ka*G)

def Guestformula1(x, G0_f, H0,K1, K2):
    A = K1*K2
    B = K1*(2*K2*H0-K2*G0_f + 1)
    C = K1*(H0 - G0_f + 1)
    D = -1*G0_f
    return A*x**3 + B*x**2 + C*x + D

def binding_HGG(G0, H0, F0, K1, K2, kDHG, kDHG2):
    G = np.array([])
    for i in range(len(G0)):
        G0_f = G0[i]
        Groots = root(Guestformula1, 0, args =  (G0_f, H0,K1, K2) )
        Groots = Groots.x #.x needed to only obtain the coefficients and not all other output
        # Only smallest positive real solution is the only one of relevance. 
        Groots_real = Groots[np.isreal(Groots)].real
        Gmin = min([p for p in Groots_real if p>= 0])
        G = np.concatenate((G,Gmin),axis = None) # None zorgt voor gewoon beide achter elkaar plakken
    return F0 + (kDHG*H0*K1*G + kDHG2*H0*K1*K2*G**2)/(1+K1*G + K1*K2*G**2)
    #Reference: Thordarson 2010

def Guestformula2(x, G0, H0_f, K1, K2):
    A = K1*K2 
    B = K1*(2*K2*G0-K2*H0_f + 1)
    C = K1*(G0-H0_f+1)
    D = -1*H0_f
    return A*x**3 + B*x**2 + C*x + D
    
def binding_HHG(G0, H0, F0, K1, K2, kDHG, kDHG2):
    H = np.array([])
    for i in range(len(G0)):
        H0_f = H0
        Hroots = root(Guestformula2, 0, args = (G0, H0_f, K1, K2))
        Hroots = Hroots.x
        Hroots_real - Hroots[np.isreal(Groots)].real
        Hmin = min([p for p in Hroots_real if p>=0])
        H = np.concatenate((H,Hmin),axis = None)
        coeff = np.array(A[i],B[i],C[i],D[i])
        Hroots = np.roots(coeff) # Might want to try other algorithms too
    return F0 + (kDHG*G0*K1*H + 2*kDHG2*G0*K1*K2*H**2)/(1+K1*H + K1*K2*H**2)
    #Reference: Thordarson 2010

# 2) Expand data array with known concentrations--------------------------------

# Host-stock 1
n_Hstock1 = m_Hstock1 / M_Host # moles
c_Hstock1 = n_Hstock1 / (V_solv_Hstock) # mol/mL
dens_Hstock1 = m_solv_Hstock1 / V_solv_Hstock # g/mL # nodig later
print('Molarity host stock solution => ' + str(round(c_Hstock1*10**(9),2)) + ' \u03BCM')

# Host-stock 2
n_Hstock2 = m_Hstock2 / M_Host # moles
c_Hstock2 = n_Hstock2 / (V_solv_Hstock) # mol/mL
dens_Hstock2 = m_solv_Hstock2 / V_solv_Hstock # g/mL # nodig later
print('Molarity host stock solution => ' + str(round(c_Hstock2*10**(9),2)) + ' \u03BCM')

# Host-stock 3
n_Hstock3 = m_Hstock3 / M_Host # moles
c_Hstock3 = n_Hstock3 / (V_solv_Hstock) # mol/mL
dens_Hstock3 = m_solv_Hstock3 / V_solv_Hstock # g/mL # nodig later
print('Molarity host stock solution => ' + str(round(c_Hstock3*10**(9),2)) + ' \u03BCM')

# Guest-stock 1
n_Gstock1 = m_Gstock1 / M_Guest # moles
c_Gstock1 = n_Gstock1/ (V_solv_Gstock ) # mol/mL
dens_Gstock1 = m_solv_Gstock1 / V_solv_Gstock # g/mL 
print('Guest stock solution => ' + str(round(c_Gstock1*10**(9),2)) + ' \u03BCM')

# Guest-stock 2
n_Gstock2 = m_Gstock2 / M_Guest # moles
c_Gstock2 = n_Gstock2/ (V_solv_Gstock ) # mol/mL
dens_Gstock2 = m_solv_Gstock2 / V_solv_Gstock # g/mL 
print('Guest stock solution => ' + str(round(c_Gstock2*10**(9),2)) + ' \u03BCM')

# Guest-stock 3
n_Gstock3 = m_Gstock3 / M_Guest # moles
c_Gstock3 = n_Gstock3/ (V_solv_Gstock ) # mol/mL
dens_Gstock3 = m_solv_Gstock3 / V_solv_Gstock # g/mL 
print('Guest stock solution => ' + str(round(c_Gstock3*10**(9),2)) + ' \u03BCM')

dens_Hstock = np.array([dens_Hstock1, dens_Hstock2, dens_Hstock3])
c_Hstock = np.array([c_Hstock1, c_Hstock2, c_Hstock3])
dens_Gstock = np.array([dens_Gstock1, dens_Gstock2, dens_Gstock3])
c_Gstock = np.array([c_Gstock1, c_Gstock2, c_Gstock3])

# Host-measurements

for i in range(1, n+1):
    name3 = 'm_Hstock_diluted' + str(i) # Calculate mass of the amount of Host stock that is diluted
    name4 = 'V_Hstock_diluted' + str(i) # Calculate the volume of the amount of Host stock that is diluted (mL)
    name5 = 'Dilution_Hstock' + str(i) # Calculate the dilution factor
    name6 = 'c_Hmeasure' + str(i) # Calculate the resulting concentration of the measurement solution
    vars()[name4] = vars()[name3] / dens_Hstock[i -1] # mL
    vars()[name5] = vars()[name4] / V_solv_Hmeasure
    vars()[name6] = vars()[name5] * c_Hstock[i-1] # mL
    print('Host measurement solution '+ str(i) + '=>' + str(round(vars()[name6]*10**(9),3)) + ' \u03BCM') # From mol/L to uM = umol/ L = 10^9 mol/L

# Guest-measurements
for j in range(1, n+1):
    name7 = 'm_Gstock_diluted' + str(j) # Calculate mass of the amount of Guest stock that is diluted
    name8 = 'V_Gstock_diluted' + str(j) # Calculate the volume of the amount of Guest stock that is diluted (mL)
    name9 = 'Dilution_Gstock' + str(j) # Calculate the dilution factor
    name10 = 'c_Gmeasure' + str(j) # Calculate the resulting concentration of the measurement solution
    vars()[name8] = vars()[name7] / dens_Gstock[j-1]
    vars()[name9] = vars()[name8] / V_solv_Gmeasure
    vars()[name10] = vars()[name9] * c_Gstock[j-1] 
    print('Guest measurement solution '+ str(j) + '=>' + str(round(vars()[name10]*10**(9),3)) + ' \u03BCM Guest')
    
    name11 = 'm_Hstock_diluted_forGmeasurement' + str(j) # Calculate mass of the amount of Host stock that is diluted
    name12 = 'V_Hstock_diluted_forGmeasurement' + str(j) # Calculate the volume of the amount of Host stock that is diluted (mL), using the known density of Host stock
    name13 = 'Dilution_Hstock_forGmeasurement' + str(j) # Calculate the dilution factor
    name14 = 'c_Hmeasure_forGmeasurement' + str(j) # Calculate the resulting concentration of the measurement solution
    vars()[name12] = vars()[name11] / dens_Hstock[j-1]
    vars()[name13] = vars()[name12] / V_solv_Gmeasure
    vars()[name14] = vars()[name13] * c_Hstock[j-1]
    print('Guest measurement solution '+ str(j) + '=>' + str(round(vars()[name14]*10**(9),3)) + ' \u03BCM Host')


# Total cumulative guest in uL--------------------------------------------------
for i in np.arange(1,n+1):
    name = 'Guest_add_' + str(i)
    name2 = 'Guest_total_' + str(i)
    vars()[name2] = np.cumsum(vars()[name]) #uL

# Total cumulative guest in moles ----------------------------------------------
for i in np.arange(1,n+1):
    volumetitration = 'Guest_total_' + str(i) # uL increments
    Gmolestitration = 'Guest_total_moles' + str(i)
    Guestconcentration = 'c_Gmeasure'+ str(i)
    Hmolestitration = 'Host_total_moles' + str(i)
    Hostconcentration = 'c_Hmeasure_forGmeasurement' + str(i)
    Hostconcentration_2000 = 'c_Hmeasure' + str(i)
    vars()[Gmolestitration] = vars()[volumetitration] * 10**(-3) * vars()[Guestconcentration] # Calculate total amount of moles of Guest in eaHmoles_av measurement point, from the additions
    vars()[Hmolestitration] = vars()[volumetitration] * 10**(-3) * vars()[Hostconcentration] # Calculate total amount of moles of Host in eaHmoles_av measurement point, from the additions
   
    vars()[Hmolestitration] = vars()[Hmolestitration] + V_H_titrated * vars()[Hostconcentration_2000] # Add to amount of moles of Host that is already present in the 2000 uL. ###############
    
    ratio = 'Equivalences' + str(i)
    vars()[ratio] = vars()[Gmolestitration]/ vars()[Hmolestitration]
    
    # Let op dat berekend hoeveel mol Host in die 2 mL en hoveel host en guest in ieder stapje

# 3) Plot data
plt.rcParams['xtick.labelsize']=12 # Change font size of all x-axis and y-axis ticks
plt.rcParams['ytick.labelsize']=12

for i in np.arange(1,n+1):
    m = 'm' + str(i)
    name2 = 'Equivalences' + str(i)
    volumetitration = 'Guest_total_' + str(i) # uL increments
    Gmolestitration = 'Guest_total_moles' + str(i)
    Hmolestitration = 'Host_total_moles' + str(i)
    maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
    vars()[maxima_peak1] = np.zeros(vars()[m])
    maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
    vars()[maxima_peak2] = np.zeros(vars()[m])
    maxwl_peak1 = 'Wavelength_max_peak1' + str(i)
    vars()[maxwl_peak1] = np.zeros(vars()[m])
    maxwl_peak2 = 'Wavelength_max_peak2' + str(i)
    vars()[maxwl_peak2] = np.zeros(vars()[m])
    
    # Obtain the baseline curve-------------------------------------------------
    basedata = splitdata(Expnr +'_set' + str(i) + '_1', i)
    baseline = basedata[1]
    
    for j in np.arange(2,vars()[m]+1): # Make this one 1, m+1 if you want to plot the baseline
        file = Expnr+'_set' + str(i) + '_' + str(j)
        data = splitdata(file, i)
        vars()[file] = data
        plt.figure(1)
        plt.subplot(1,n,i)
        graph = 'plot' + str(i)
        absorbance_corrected = absorbance - baseline
        if smoothening == True:
            absorbance_corrected = savgol_filter(absorbance_corrected, 51,3)
        cmap= sns.color_palette("RdGy", vars()[m]) #not equal to m, since the last ones will be white
        vars()[graph] = plt.plot(wavelength, absorbance_corrected, label = str(round(vars()[name2][j-1],2)), color = cmap[j-1])

        lambda1_max1 = np.max
        a = np.where(wavelength == lambda1)
        index_l1 = int(float(a[int(0)])) # trick to obtain the index of the wavelength otherwise arrat([450], dtype = int64)
        b = np.where(wavelength == lambda2)
        index_l2 = int(float(b[int(0)]))
        
        if method == 'max_wl':
            max1 = np.max(absorbance_corrected[index_l1 - 30: index_l1 + 30])
            max1_lambda_index = np.where(absorbance_corrected == max1)
            max1_lambda_index = int(float(max1_lambda_index[int(0)]))
            max1_lambda = wavelength[max1_lambda_index]
        
            max2 = np.max(absorbance_corrected[index_l2 - 30: index_l2 + 30])
            max2_lambda_index = np.where(absorbance_corrected == max2)
            max2_lambda_index = int(float(max2_lambda_index[int(0)]))
            max2_lambda = wavelength[max2_lambda_index]
            
            vars()[maxima_peak1][j-1] = max1 # store maxima
            vars()[maxima_peak2][j-1] = max2
            vars()[maxwl_peak1][j-1] = max1_lambda
            vars()[maxwl_peak2][j-1] = max2_lambda
        if method == 'single_wl':
            max1 = absorbance_corrected[index_l1 ]
            max2 = absorbance_corrected[index_l2]
            vars()[maxima_peak1][j-1] = max1 # store maxima
            vars()[maxima_peak2][j-1] = max2


# 4) Plot data
    
    plt.figure(1)
    plt.title('Replicate' + str(i), fontsize = 16)
    plt.suptitle('Fluorescence spectrum for titration ' + Host + ' + '+ Guest + ', '+ str(T) +' \u00b0' +'C ', fontsize = 20)
    plt.xlabel('Wavenumber (cm^-1)', fontsize = 14)
    plt.ylabel('Fluorescence', fontsize = 14)
    #plt.legend(loc='upper left', title = '[G]0/[H]0', ncol = 4, bbox_to_anchor=(0,-0.15), handlelength = 1.5, columnspacing = 1.0)
    plt.subplots_adjust(bottom=0.31, wspace = 0.23) # This is the slider of figures
    #plt.xlim(wavelength[0], 800)
    plt.xlim(550,800)
    plt.grid
    
    plt.figure(2)
    plt.subplot(1,n,i)
    plt.plot(vars()[name2][1:], vars()[maxima_peak1][1:], label = 'peak 1 at ' + str(lambda1) + ' nm', color = '#ff0000', marker = '.', markerfacecolor = '#ff0000' ) # first index gives background
    
    plt.plot(vars()[name2][1:], vars()[maxima_peak2][1:], label = 'peak 2 at ' + str(lambda2) + ' nm', color = '#e20000', marker = '*', markerfacecolor = '#e20000')
    
    plt.suptitle('Binding isotherm for ' + Host + ' + '+ Guest + ', '+ str(T) +' \u00b0' +'C ', fontsize = 20)
    plt.title('Replicate' + str(i), fontsize = 16)
    plt.subplots_adjust(bottom=0.20, wspace = 0.27) # This is the slider of figures
    plt.xlabel('Equivalences [G]0/[H]0', fontsize = 14)
    plt.ylabel('Fluorescence', fontsize = 14)
    #np.figlegend can be used to only print one legend instead of for each subplot
    plt.grid
    plt.legend(loc = 'upper left', fontsize = 10, handlelength = 1.0, bbox_to_anchor=(0,-0.15))
        
# 5) Fit data
    TotalV = V_H_titrated + vars()[volumetitration] * 10**(-3) # mL
    
    cH0 = vars()[Hmolestitration]/ (TotalV*10**(-3)) # M
    cG0 = vars()[Gmolestitration] / (TotalV*10**(-3))
    G0H0 = cG0/cH0 #ratio on x-axis, this is the exact same as Equivalencesi parameter. n/n = c/c
    guestspace = np.linspace(0, cG0[-1], num = 1000)
    
    Hmoles_av = np.average(cH0) #Takes average of amount of moles Host in sample, since it ideally should be constant (which it isn't and I have in this script the possibility to calculate the exact Host, but for now I keep it this way.
     
    yvalues1 = vars()[maxima_peak1]
    yvalues2 = vars()[maxima_peak2]
    F0_1 = vars()[maxima_peak1][1] # The second index gives the fluorescence of the host solution before added guest
    F0_2 = vars()[maxima_peak2][1] #F0_1 and F0_2 are initial guesses
    print('F0_peak1 = ' + str(F0_1)) # Check if same as fitted values
    print('F0_peak2 = ' + str(F0_2))
    eps = 10**(-15) #Trick to keep H0 not as a variable and only optimize F0 and Ka
    eps2 = 10**(-10)
    
    # Define fitting parameters: Guess and bounds
    if formula == 'binding':
        p0_1 = [Hmoles_av, F0_1, 130000]
        p0_2 = [Hmoles_av, F0_2, 130000]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1),( Hmoles_av +eps, F0_1 + eps2,1000000000))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1),( Hmoles_av +eps, F0_2 + eps2,1000000000))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1),( Hmoles_av +eps, np.inf ,1000000000))
            bounds2 = bounds1
    
    if formula == 'binding_alt': #G0, H0, F0, Ka, kDHG
        p0_1 = [Hmoles_av, F0_1, 130000,-3000000000] # Guess of Ka comes from the formula 'binding'
        p0_2 = [Hmoles_av, F0_2, 130000,-3000000000]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,-np.inf),( Hmoles_av +eps, F0_1 + eps2 ,1000000000,np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,-np.inf),( Hmoles_av +eps, F0_2 + eps2 ,1000000000,np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,-np.inf),( Hmoles_av +eps, np.inf ,1000000000,np.inf))
            bounds2 = bounds1
        
    if formula == 'binding_full': #G0, H0, F0, Ka, kH, kH0, kHG
        p0_1 = [Hmoles_av, F0_1, 130000, 3,   3,  -100]
        p0_2 = [Hmoles_av, F0_2, 130000, 3, 3,  -100]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,-np.inf, -np.inf, -np.inf),( Hmoles_av +eps, F0_1 + eps2 ,1000000000,np.inf, np.inf,np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,-np.inf, -np.inf, -np.inf),( Hmoles_av +eps, F0_2 + eps2 ,1000000000,np.inf,np.inf,np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,-np.inf, -np.inf, -np.inf),( Hmoles_av +eps, np.inf ,1000000000,np.inf, np.inf, np.inf))
            bounds2 = bounds1
            
    if formula == 'binding_HGG': #H0, F0, K1, K2, kDHG, kDHG2
        p0_1 = [Hmoles_av, F0_1, 100000000,100000000,-100000000000000,-1000000]
        p0_2 = [Hmoles_av, F0_2, 100000000,100000000,-100000000000000,-1000000]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,1,-np.inf, -np.inf),( Hmoles_av +eps, F0_1 + eps2 ,1000000000,1000000000, np.inf, np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,1,-np.inf, -np.inf),( Hmoles_av +eps, F0_2 + eps2 ,1000000000,1000000000, np.inf, np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,1,-np.inf, -np.inf),( Hmoles_av +eps, np.inf ,1000000000,1000000000, np.inf, np.inf))
            bounds2 = bounds1
    
    if formula == 'binding_HHG':
        p0_1 = [Hmoles_av, F0_1, 100000,100000,1,1]
        p0_2 = [Hmoles_av, F0_2, 100000,100000,1,1]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,1,-np.inf, -np.inf),( Hmoles_av +eps, F0_1 + eps2 ,1000000000,1000000000, np.inf, np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,1,-np.inf, -np.inf),( Hmoles_av +eps, F0_2 + eps2 ,1000000000,1000000000, np.inf, np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,1,-np.inf, -np.inf),( Hmoles_av +eps, np.inf ,1000000000,1000000000, np.inf, np.inf))
            bounds2 = bounds1

    # Find curve fitting parameters
    popt1, pcov1 = curve_fit(vars()[formula], cG0[1:], yvalues1[1:], p0 = p0_1, bounds = bounds1)
    popt2, pcov2 = curve_fit(vars()[formula], cG0[1:], yvalues2[1:], p0 = p0_2, bounds = bounds2)
    
    
    
    # Store F0 and Ka values ---------------------------------------------------
    
    initialF[0,i-1] = popt1[1]
    initialF[1,i-1] = popt2[1]
    bindingconstant[0,i-1] = popt1[2]
    bindingconstant[1,i-1] = popt2[2]
    if formula == 'binding_alt':
        kDHGcollect[0, i-1] = popt1[3]
        kDHGcollect[1, i-1] = popt2[3]
    if formula == 'binding_full':
        kHcollect[0, i-1] = popt1[3]
        kHcollect[1, i-1] = popt2[3]
        kH0collect[0, i-1] = popt1[4]
        kH0collect[1, i-1] = popt2[4]
        kHGcollect[0, i-1] = popt1[5]
        kHGcollect[1, i-1] = popt2[5]
        
    
    plt.figure(3)
    plt.subplot(1,n,i)
    plt.plot(G0H0[1:], yvalues1[1:], label = 'Peak 1 at ' + str(lambda1) + ' nm' )
    plt.plot(guestspace/Hmoles_av, vars()[formula](guestspace, *popt1),label = 'fit: F_0 = ' + str(round(popt1[1],1)) + u', K\u2090 = ' + str(round(popt1[2],1)))
    plt.plot(G0H0[1:], yvalues2[1:], label = 'Peak 2 at ' + str(lambda2) + ' nm' )
    plt.plot(guestspace/Hmoles_av, vars()[formula](guestspace, *popt2), label = 'fit: F_0 = ' + str(round(popt2[1],1)) + u', K\u2090 = ' + str(round(popt2[2],1)))
    
    plt.suptitle('Fluorescence binding curve '+ formula + ', for ' + Host + ' + '+ Guest +  ', '+ str(T) +' \u00b0' +'C ', fontsize = 20)
    plt.title('Replicate' + str(i), fontsize = 16)
    plt.subplots_adjust(bottom=0.29, wspace = 0.27) # This is the slider of figures
    plt.xlabel('[G]0/[H]0', fontsize = 14)
    plt.ylabel('F', fontsize =14)
    plt.grid
    plt.legend(loc = 'upper left', fontsize = 10, handlelength = 1.0, bbox_to_anchor=(0,-0.15))

    
# Show all plots
plt.show(1)
plt.show(2)
plt.show(3)

T = 25 #'\u00b0' +'C '
TK = T + 273.15
R = 8.31446
# Calculate Delta_r_G
DrG = -1 * R * TK * np.log(bindingconstant)
    
# Print data and statistical analysis
print(' ')
print(' ')
print('The fitted F0 values for peak 1 are:')
print(initialF[0])
print('The fitted F0 values for peak 2 are:')
print(initialF[1])
print('The fitted Ka values for peak 1 are:')
print(bindingconstant[0])
print('The fitted Ka values for peak 2 are:')
print(bindingconstant[1])
print('The DrG values for peak 1 are:')
print(DrG[0])
print('The DrG values for peak 2 are:')
print(DrG[1])
print(' ')
print(' ')
 
 #Alter this to display right values
print('The average F0 for peak 1 is:')
print('{:.2e}'.format(np.average(initialF[0])) + u" \u00B1 " + '{:0=3.2f}e+03'.format(np.std(initialF[0]/10**3)))
print('The average F0 for peak 2 is:')
print('{:.2e}'.format(np.average(initialF[1])) + u" \u00B1 " + '{:0=3.2f}e+03'.format(np.std(initialF[1]/10**3)))
print('The average Ka for peak 1 is:')
print('{:.2e}'.format(np.average(bindingconstant[0])) + u" \u00B1 " + '{:0=3.2f}e+04'.format(np.std(bindingconstant[0]/10**4)))
print('The average Ka for peak 2 is:')
print('{:.2e}'.format(np.average(bindingconstant[1])) + u" \u00B1 " + '{:0=3.2f}e+04'.format(np.std(bindingconstant[1]/10**4)))
print('The average DrG for peak 1 is:')
print('{:.3e}'.format(np.average(DrG[0])) + u" \u00B1 " + '{:0=3.2f}e+03'.format(np.std(DrG[0]/10**3)))
print('The average DrG for peak 2 is:')
print('{:.3e}'.format(np.average(DrG[1])) + u" \u00B1 " + '{:0=3.2f}e+03'.format(np.std(DrG[1]/10**3)))


# TO DO: -----------------------------------------------------------------------
#Option self finding of maximum wavelengths 
#Save spectra
#Clean up script
#Output data nice in a file: guessed and fitted values
#Smoothing operator making more points inbetween, for finding optimum. 
#Make this a function
#Wavelength at the maxima might now just give the first index, however might need to average it and get a more accurate one. Now it gives the lambda of the maxima on 5 nm datapoints instead of more significance. 
#Script finding maxima itself


print(' ')
print(' ')

## Outputfile
os.chdir('E:\\Uni\\Bachelor_Internship\\Internship\\Analysis_techniques\\Fluorescence\\' + Expnr)
if smoothening == True:
    file = open(Expnr+"_output_" + formula + '_' + method+ '_smooth.txt', "w")
if smoothening == False:
    file = open(Expnr+"_output_" + formula + '_' + method+ '.txt', "w")
file.write('Experiment: ' + Expnr + '\n') 
file.write('Outputfile for binding study with fluorescence of ' + Host + ' with ' + Guest + '. \n')
file.write('\nFitting parameters \n')
file.write('Formula: ' + formula + '\n')
file.write('n = \t' + str(n) + '\n')
for i in range(1,n+1):
    m = 'm' + str(i)
    file.write('m' + str(i) + ' = \t' + str(vars()[m]) + '\n')
file.write('lambda =  ' + str(lambda1) + ' nm \n')
file.write('lambda =  ' + str(lambda2) + ' nm \n')
file.write('Method = ' + method + '\n')
file.write('Smoothening = ' + str(smoothening) + '\n')

file.write('\nSolutions\n')
file.write(' \n')
file.write('Host-stock:\n')
file.write('[H]0_1 = ' + str(round(c_Hstock1*10**(9),3)) + ' uM\n')
file.write('[H]0_2 = ' + str(round(c_Hstock2*10**(9),3)) + ' uM\n')
file.write('[H]0_3 = ' + str(round(c_Hstock3*10**(9),3)) + ' uM\n')
file.write('Density_1 = ' + str(round(dens_Hstock1, 2)) + ' g/mL\n')
file.write('Density_2 = ' + str(round(dens_Hstock2, 2)) + ' g/mL\n')
file.write('Density_3 = ' + str(round(dens_Hstock3, 2)) + ' g/mL\n')
file.write('Guest-stock:\n')
file.write('[H]0_1 = ' + str(round(c_Gstock1*10**(9),3)) + ' uM\n')
file.write('[H]0_2 = ' + str(round(c_Gstock2*10**(9),3)) + ' uM\n')
file.write('[H]0_3 = ' + str(round(c_Gstock3*10**(9),3)) + ' uM\n')
file.write('Density_1 = ' + str(round(dens_Gstock1, 2)) + ' g/mL\n')
file.write('Density_2 = ' + str(round(dens_Gstock2, 2)) + ' g/mL\n')
file.write('Density_3 = ' + str(round(dens_Gstock3, 2)) + ' g/mL\n')
for i in range(1,n+1):
    name6 = 'c_Hmeasure' + str(i)
    file.write('Host-measurement ' + str(i) + ':')
    file.write('\t[H]0 = ' + str(round(vars()[name6]*10**(9),3)) + ' uM \n')
for i in range(1,n+1):
    name10 = 'c_Gmeasure' + str(i)
    name14 = 'c_Hmeasure_forGmeasurement' + str(i)
    file.write('Guest-measurement ' + str(i) + ':')
    file.write('\t[H]0 = ' + str(round(vars()[name14]*10**(9),3)) + ' uM \n')
    file.write('\t\t\t[G]0 = ' + str(round(vars()[name10]*10**(9),3)) + ' uM \n')
    
file.write('\nTitration: \n')
file.write(' \n')
for i in range(1,n+1):
    name = 'Guest_add_' + str(i)
    name2 = 'Guest_total_' + str(i)
    ratio = 'Equivalences' + str(i)
    
    file.write('Set '+ str(i) + ':\n')
    file.write('Host measurement volume = ' + str(V_H_titrated * 1000) + ' uL \n')
    file.write('Added guest volume = ')
    file.writelines(str(vars()[name]) + '\n')
    file.write('Cumulative volume = ')
    file.writelines(str(vars()[name2]) + '\n')
    file.write('Equivalences = ' + str(np.around(vars()[ratio],2)) + '\n')
    file.write(' \n')

file.write('\nInitial guesses\n')
file.write('[H]0 = ' + str(Hmoles_av) + '\n')
file.write('F0_1 = ' + str(round(vars()[maxima_peak1][1],1)) + '\n')
file.write('F0_2 = ' + str(round(vars()[maxima_peak2][1],1)) + '\n')
if formula == 'binding':
    file.write('Ka = ' + str(p0_1[2]) + '\n')
if formula == 'binding_alt':
    file.write('Ka = ' + str(p0_1[2]) + '\n')
    file.write('kDHG = ' + str(p0_1[3]) + '\n')
if formula == 'binding_full':
    file.write('Ka = ' + str(p0_1[2]) + '\n') 
    file.write('kH = ' + str(p0_1[3]) + '\n')
    file.write('kH0 = ' + str(p0_1[4]) + '\n')
    file.write('kHG = ' + str(p0_1[5]) + '\n')

file.write('\nOutput data: fitted parameters \n')
if formula == 'binding':
    for i in range(1,n+1): 
        maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
        maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
        file.write('Set ' + str(i) + ':\n')
        file.write('\t F0(fixed) \t F0(fitted) \t Ka \n')
        file.write('Peak 1   ' + str(round(vars()[maxima_peak1][1],1)) + '\t \t ' + str(round(initialF[0,i-1],1)) + '\t \t ' + str(round(bindingconstant[0,i-1],0)) + '\n')
        file.write('Peak 2   ' + str(round(vars()[maxima_peak2][1],1)) + '\t \t ' + str(round(initialF[1,i-1],1)) + '\t \t ' + str(round(bindingconstant[1,i-1],0)) + '\n')
        file.write(' \n')

if formula == 'binding_alt':
    for i in range(1,n+1):
        maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
        maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
        file.write('Set ' + str(i) + ':\n')
        file.write('\t F0(fixed) \t F0(fitted) \t Ka \t \t kDHG (*10^9) \n')
        file.write('Peak 1   ' + str(round(vars()[maxima_peak1][1],1)) + '\t \t ' + str(round(initialF[0,i-1],1)) + '\t \t ' + str(round(bindingconstant[0,i-1],1)) + '\t ' + str(round(kDHGcollect[0,i-1]/10**9,3)) + '\n')
        file.write('Peak 2   ' + str(round(vars()[maxima_peak2][1],1)) + '\t \t ' + str(round(initialF[1,i-1],1)) + '\t \t ' + str(round(bindingconstant[1,i-1],1)) + '\t ' + str(round(kDHGcollect[1,i-1]/10**9,3)) + '\n')
        file.write(' \n')

if formula == 'binding_full':
    for i in range(1,n+1):
        maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
        maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
        kH1 = '{:{align}{width}{prec}f}'.format(kHcollect[0,i-1], align= '<', width = '8', prec= '.3')
        kH2 = '{:{align}{width}{prec}f}'.format(kHcollect[1,i-1], align= '<', width = '8', prec= '.3')
        kH01 = '{:{align}{width}{prec}f}'.format(kH0collect[0,i-1], align= '<', width = '8', prec= '.3')
        kH02 = '{:{align}{width}{prec}f}'.format(kH0collect[1,i-1], align= '<', width = '8', prec= '.3')
        
        file.write('Set ' + str(i) + ':\n')
        file.write('\t F0(fixed) \t F0(fitted) \t Ka \t \t kH \t        kH0 \t \t kHG \n')
        file.write('Peak 1   ' + str(round(vars()[maxima_peak1][1],1)) + '\t \t ' + str(round(initialF[0,i-1],1)) + '\t \t ' + str(round(bindingconstant[0,i-1],1)) + '\t ' + kH1 + '\t' + kH01 + '\t ' + str(round(kHGcollect[0,i-1],3)) +'\n') 
        file.write('Peak 2   ' + str(round(vars()[maxima_peak2][1],1)) + '\t \t ' + str(round(initialF[1,i-1],1)) + '\t \t ' + str(round(bindingconstant[1,i-1],1)) + '\t ' + kH2 + '\t' + kH02 + '\t ' + str(round(kHGcollect[1,i-1],3)) +'\n')
        file.write(' \n') 

file.write('\nAverages \n')
file.write('Peak 1: \n')
file.write('F0 = ' + '{:.2e}'.format(np.average(initialF[0])) + u" \u00B1 " + '{:0=3.2f}e+03'.format(np.std(initialF[0]/10**3)) + '\n')
file.write('Ka = ' + '{:.2e}'.format(np.average(bindingconstant[0])) + u" \u00B1 " + '{:0=3.2f}e+06'.format(np.std(bindingconstant[0]/10**6)) + '\n')
file.write('Peak 2: \n')
file.write('F0 = ' + '{:.2e}'.format(np.average(initialF[0])) + u" \u00B1 " + '{:0=3.2f}e+03'.format(np.std(initialF[0]/10**3)) + '\n')
file.write('Ka = ' + '{:.2e}'.format(np.average(bindingconstant[1])) + u" \u00B1 " + '{:0=3.2f}e+06'.format(np.std(bindingconstant[1]/10**6)) + '\n')

file.close()