'''
This script was used to plot fluorescence data obtained from 
host-guest titrations and it allows the user to easily obtain the 
association constant for the complexes that are researched. It 
allows the user to easily calculate the concentrations of the used 
samples and calculate the equivalences for each titration step. 
The program output is a figure with all the spectra, a figure 
containing the bindingisotherm and a figure with the resulting fit. 
Besides this a separate file is produced and saved containing all 
fitting parameters, sample data and the fitting output. 
This program is applicable for:
* Fluorescence titrations for Jasco FP-8300, saved all spectra as
  csv files
* The fluorescence spectrum has 2 peaks
* 1 stock solution is used for both the host and guest measurement 
  solutions. The script for 2x3 stock solutions is very similar.

The information in the data specifiers section must be filled in by 
hand for each experiment. This script can be used as a template and 
it is possible to adjust it to specific experiments and output 
preferences. Feel free to contact the writer for questions. 

P.C.P. Teeuwen 
p.teeuwen@student.ru.nl
Last updated: 30-05-2019
'''

import os                          # Specify working directory
import numpy as np                 # Work with lists and matrices
import matplotlib.pyplot as plt    # Plot the data in figures
import matplotlib.colors as colors # Specify colors of graphs
import matplotlib.cm as cmx        # Build in colormaps
from scipy.signal import find_peaks    # Find signal peaks
from scipy.signal import savgol_filter # Smoothening spectra
from scipy.optimize import curve_fit   # Fitting curves
from scipy.optimize import fsolve      # Solve equitions
from scipy.optimize import root        # Find roots of equations

import seaborn as sns                  # Plotting lay-out and colors
sns.set()
current_palette = sns.color_palette()

print('')
print('')

############################################################
# Data specifiers

os.getcwd()
os.chdir('E:\\Directory') 
# Change current working directory. Where spectra are stored.

Expnr =                     # Name of experiment, string
Host =                      # Host name, string
Guest =                     # Guest name, string

T = 25                      # Temperature (degrees celcius)
TK = T + 273.15
R = 8.31446

# Host-stock
M_Host =                    # Host molecular weight in g/mol
m_Hstock =                  # Mass of host that is dissolved (g)
V_solv_Hstock =             # Volume of host-stock (mL)
m_solv_Hstock =             # Mass of host-stock solution (g)

# Guest-stock
M_Guest =                   # Guest molecular weight in g/mol
m_Gstock =                  # Mass of guest that is dissolved (g)
V_solv_Gstock =             # Volume  of guest-stock (mL)
m_solv_Gstock =             # Mass of guest-stock solution (g)

# Host-measurement solution
V_solv_Hmeasure =           # Total volume of solution (mL)
m_Hstock_diluted1 =         # Mass of diluted Hstock solution (g)
m_Hstock_diluted2 = 
m_Hstock_diluted3 =
V_H_titrated = 2.0          # Volume of solution in the cuvette (mL)

# Guest-measurement solution
V_solv_Gmeasure =           # Total volume of solution (mL)
m_Gstock_diluted1 =         # Mass of diluted Gstock solution (g)
m_Gstock_diluted2 = 
m_Gstock_diluted3 = 
m_Hstock_diluted_forGmeas1 = # Mass of diluted Hstock (g)
m_Hstock_diluted_forGmeas2 = 
m_Hstock_diluted_forGmeas3 = 
# etc...

# measurement parameters
n = 3                       # Amount of measurement sets 
m1 =                        # Amount of spectra measured per set
m2 =
m3 = 
#etc.. : m.. = .. 

# Follow around wavelength
lambda1 = 650 #nm
lambda2 = 715 #nm

# What kind of fitting?
F0_free = True       # Fitting parameter F0 free, check implemented
smoothening = True   # Apply smoothening of spectra
formula = 'Connors'  # Fitting: SternVolmer, Tsukube or Connors
method = 'single_wl'     # Fitting method: 'single_wl' or 'max_wl'


# Guest added in each measurement (uL)---------------------------

for i in np.arange(1,n+1):
    name = 'Guest_add_' + str(i)
    m = 'm'+ str(i)
    vars()[name] = np.array([0,0,10,10,10,10,10,10,10,10,10,10,30,
                   30,30,30,30])
    endvalues = np.full(vars()[m] - len(vars()[name]),50)
    # fill series with 50
    vars()[name] = np.concatenate(((vars()[name]), endvalues))

# Allow for changes in added guest as compared to planned values
# Format: Guest_add_i[measurement number -1]
Guest_add_2[16] = 35 # Example



##################################################






## Defined functions within script 

# Loading data
def splitdata(fluorescence, n):
    '''
    This function loads the fluorescence data stored in the current 
    directory.The important columns (wavelength, absorbance) are 
    separated by a comma-delimiter. Data must be saved as .csv in 
    the format Expnr_seti_j.csv. Here i is the set number and j 
    is the spectrum number. It obtains each spectra within a submap 
    'seti' in a map 'Expnr'.
    fluorescence = name of data file (string)
    n = set number
    '''
    global wavelength # Need to extract data from this function
    global absorbance
    os.chdir('E:\\Directory\\' + Expnr + '\\set' + str(n))
    wavelength, absorbance = np.genfromtxt(
    fluorescence + '.csv', unpack = True, delimiter = ',', 
    skip_header = 19, skip_footer = 44, usecols = (0,1))
    return np.array([wavelength , absorbance])

# Fitting functions
def SternVolmer(G0, H0, F0, Ka):  
    '''
    Function applicable if complex is fluorescently silent and 
    applies the assumption of no dynamic quenching of the host. 
    Can only be used if quenching is observed.
    '''
    return F0 / (1+Ka*(1/2)*((G0 - H0 - 1/Ka) + np.sqrt((
    G0 - H0 - 1/Ka)**2 + 4*G0/Ka)))
    # based on F0/F = 1+Ka[G]


def Tsukube(G0, H0, F0, Ka, kDHG): 
    '''
    Complex does not need to be fluorescently silent. This function
    applies the assumption of no dynamic quenching of host.
    '''
    return F0 + kDHG*(1/2)*(((1/Ka) + H0 + G0) - np.sqrt(
    (G0 + H0 + 1/Ka)**2 - 4*G0*H0)) 
    #Based on DY = Y_DHG([HG])



def Connors(G0, H0, F0, Ka, kH, kH0, kHG): 
    '''
    Complex does not need to be fluorescently silent. This 
    function does not apply the assumption of no dynamic 
    quenching of the host. 
    '''
    G = (1/2)*((G0 - H0 - 1/Ka) + np.sqrt((G0 - H0 - 1/Ka)**2 
        + 4*G0/Ka))
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
        Groots = Groots.x #.x  to only obtain the coefficients
        # Only smallest positive real solution relevant 
        Groots_real = Groots[np.isreal(Groots)].real
        Gmin = min([p for p in Groots_real if p>= 0])
        G = np.concatenate((G,Gmin),axis = None) 
    return F0 + (kDHG*H0*K1*G + kDHG2*H0*K1*K2*G**2)/(1+K1*G + 
           K1*K2*G**2)
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
        Hroots = np.roots(coeff) # Might want to try other algorithms
    return F0 + (kDHG*G0*K1*H + 2*kDHG2*G0*K1*K2*H**2)/(1+K1*H + 
           K1*K2*H**2)
    #Reference: Thordarson 2010

## Calculate sample data and titration data

# Host-stock -------------------------------------------------------------------

n_Hstock = m_Hstock / M_Host                # Amount of host (n)
c_Hstock = n_Hstock / V_solv_Hstock         # [host] (mol/mL)
dens_Hstock = m_solv_Hstock / V_solv_Hstock # Density stock (g/mL) 
print('Molarity host stock solution => ' 
      +str(round(c_Hstock*10**(9),2)) + ' \u03BCM')

# Guest-stock-------------------------------------------------------------------

n_Gstock = m_Gstock / M_Guest               # Amount of guest
c_Gstock = n_Gstock / V_solv_Gstock         # [guest] (mol/mL)
dens_Gstock = m_solv_Gstock / V_solv_Gstock # Density stock (g/mL) 
print('Guest stock solution => ' 
      + str(round(c_Gstock*10**(9),2)) + ' \u03BCM')

# Host-measurements-------------------------------------------------------------

for i in range(1, n+1):
    # Calculate concentration of each measurement solution 
    mass = 'm_Hstock_diluted' + str(i) 
    volume = 'V_Hstock_diluted' + str(i) 
    dilution = 'Dilution_Hstock' + str(i) 
    concentration = 'c_Hmeasure' + str(i)
    # Volume of stock (mL) 
    vars()[volume] = vars()[mass] / dens_Hstock     
    # Dilution factor
    vars()[dilution] = vars()[volume] / V_solv_Hmeasure  
    # Concentration host
    vars()[concentration] = vars()[dilution] * c_Hstock   
    print('Host measurement solution '+ str(i) + '=>' + str(
    round(vars()[concentration]*10**(9),3)) + ' \u03BCM') 
    # uM = 10^9 mol/L


# Guest-measurements------------------------------------------------------------

for j in range(1, n+1):
    # Calculate the concentration of host and guest of each 
    # measurement solution
    mass2 = 'm_Gstock_diluted' + str(j) 
    volume2 = 'V_Gstock_diluted' + str(j) 
    dilution3 = 'Dilution_Gstock' + str(j) 
    concentrationG = 'c_Gmeasure' + str(j) 
     # Volume of guest-stock (mL)
    vars()[volume2] = vars()[mass2] / dens_Gstock 
    # Dilution factor
    vars()[dilution3] = vars()[volume2] / V_solv_Gmeasure  
    # Concentration guest
    vars()[concentrationG] = vars()[dilution3] * c_Gstock  
    print('Guest measurement solution '+ str(j) + '=>' + str(
    round(vars()[concentrationG]*10**(9),3)) + ' \u03BCM guest')
    
    mass3  = 'm_Hstock_diluted_forGmeas' + str(j) 
    volume3 = 'V_Hstock_diluted_forGmeas' + str(j)
    dilution3 = 'Dilution_Hstock_forGmeas' + str(j) 
    concentrationH = 'c_Hmeasure_forGmeas' + str(j) 
    # Volume of host-stock (mL)
    vars()[volume3] = vars()[mass3] / dens_Hstock    
    # Dilution factor
    vars()[dilution3] = vars()[volume3] / V_solv_Gmeasure 
    # Concentration host
    vars()[concentrationH] = vars()[dilution3] * c_Hstock 
    print('Guest measurement solution '+ str(j) + '=>' + str(
    round(vars()[concentrationH]*10**(9),3)) + ' \u03BCM host')


# Total cumulative guest in uL--------------------------------------------------
for i in np.arange(1,n+1):
    add = 'Guest_add_' + str(i)
    total = 'Guest_total_' + str(i)
    vars()[total] = np.cumsum(vars()[add]) 

# Total cumulative guest in moles ----------------------------------------------
# Take into account the volume change for each addition and that 
# each addition of guest measurement solution also adds extra host. 
for i in np.arange(1,n+1):
    volumetitr = 'Guest_total_' + str(i) # uL increments
    Gmolestitr = 'Guest_total_moles' + str(i)
    Guestcon = 'c_Gmeasure'+ str(i)
    Hmolestitr = 'Host_total_moles' + str(i)
    Hostcon = 'c_Hmeasure_forGmeas' + str(i)
    Hostcon_2000 = 'c_Hmeasure' + str(i)
    # Total amount of moles of guest in cuvette for each point
    vars()[Gmolestitr] = vars()[volumetitr]*10**(-3)*vars()[Guestcon]
    # Total amount of host in cuvette from 2000 uL and additions
    vars()[Hmolestitr] = (vars()[volumetitr]*10**(-3)*vars()[Hostcon] 
                         + V_H_titrated * vars()[Hostcon_2000])
    ratio = 'Equivalences' + str(i)
    vars()[ratio] = vars()[Gmolestitration]/ vars()[Hmolestitration]

## Plot the spectra of each set and fit all parameters
plt.rcParams['xtick.labelsize']=12 # Change font size of axes
plt.rcParams['ytick.labelsize']=12

bindingconstant = np.zeros([2,n]) # Empty matrix to store fitting data
initialF = np.zeros([2,n])
kDHGcollect = np.zeros([2,n])
kHcollect = np.zeros([2,n])
kH0collect = np.zeros([2,n])
kHGcollect = np.zeros([2,n])

for i in np.arange(1,n+1):
    m = 'm' + str(i)
    name2 = 'Equivalences' + str(i)
    volumetitr = 'Guest_total_' + str(i) 
    Gmolestitr = 'Guest_total_moles' + str(i)
    Hmolestitr = 'Host_total_moles' + str(i)
    maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
    vars()[maxima_peak1] = np.zeros(vars()[m])
    maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
    vars()[maxima_peak2] = np.zeros(vars()[m])
    maxwl_peak1 = 'Wavelength_max_peak1' + str(i)
    vars()[maxwl_peak1] = np.zeros(vars()[m])
    maxwl_peak2 = 'Wavelength_max_peak2' + str(i)
    vars()[maxwl_peak2] = np.zeros(vars()[m])
    
    # Obtain the baseline curve
    basedata = splitdata(Expnr +'_set' + str(i) + '_1', i)
    baseline = basedata[1]
    
    # Plot spectra and collect wavelengths and intensities of peaks 
    for j in np.arange(2,vars()[m]+1):
        file = Expnr+'_set' + str(i) + '_' + str(j) # file name
        data = splitdata(file, i)      # Form collumns
        vars()[file] = data            # Data gets name of the file
        plt.figure(1)                  # Open figure 1
        plt.subplot(1,n,i)             # Define subplots
        graph = 'plot' + str(i)
        absorbance_cor = absorbance - baseline # Correct 
        if smoothening == True:                # Smooth spectra
            absorbance_cor = savgol_filter(absorbance_cor, 51,3)
        cmap= sns.color_palette("RdGy", vars()[m]) # Define color
        vars()[graph] = plt.plot(wavelength, absorbance_cor,
        label = str(round(vars()[name2][j-1],2)), color = cmap[j-1])

        lambda1_max1 = np.max
        a = np.where(wavelength == lambda1) # Find index
        index_l1 = int(float(a[int(0)])) 
        # For example extract 450 from arrat([450], dtype = int64)
        b = np.where(wavelength == lambda2)
        index_l2 = int(float(b[int(0)]))
        
        if method == 'max_wl':
            max1 = np.max(
                   absorbance_cor[index_l1 - 30: index_l1 + 30])
            max1_lambda_index = np.where(absorbance_cor == max1)
            max1_lambda_index = int(
                                float(max1_lambda_index[int(0)]))
            max1_lambda = wavelength[max1_lambda_index]
        
            max2 = np.max(
                   absorbance_cor[index_l2 - 30: index_l2 + 30])
            max2_lambda_index = np.where(absorbance_cor == max2)
            max2_lambda_index = int(
                                float(max2_lambda_index[int(0)]))
            max2_lambda = wavelength[max2_lambda_index]
            
            vars()[maxima_peak1][j-1] = max1 # store maxima
            vars()[maxima_peak2][j-1] = max2
            vars()[maxwl_peak1][j-1] = max1_lambda
            vars()[maxwl_peak2][j-1] = max2_lambda
        if method == 'single_wl':
            max1 = absorbance_cor[index_l1 ]
            max2 = absorbance_cor[index_l2]
            vars()[maxima_peak1][j-1] = max1 # store maxima
            vars()[maxima_peak2][j-1] = max2



# 4) Plot data
    
    plt.figure(1)
    plt.title('Duplicate' + str(i), fontsize = 16)
    plt.suptitle('Fluorescence spectrum for titration ' + Host 
    + ' + '+ Guest + ', '+ str(T) +' \u00b0' +'C ', fontsize = 20)
    plt.xlabel('Wavelength (nm)', fontsize = 14)
    plt.ylabel('Fluorescence', fontsize = 14)
    plt.legend(loc='upper left', title = '[G]0/[H]0', ncol = 4, 
               bbox_to_anchor=(0,-0.15), handlelength = 1.5, 
               columnspacing = 1)
    plt.subplots_adjust(bottom=0.31, wspace = 0.23) # Slider 
    plt.xlim(550, wavelength[-1]) # cut off plot
    plt.grid
    
    plt.figure(2)
    plt.subplot(1,n,i)
    plt.plot(vars()[name2][1:], vars()[maxima_peak1][1:], 
             label = 'peak 1 at ' + str(lambda1) + ' nm', 
             color = '#ff0000', marker = '.', 
             markerfacecolor = '#ff0000' ) # skip first index
    
    plt.plot(vars()[name2][1:], vars()[maxima_peak2][1:], 
             label = 'peak 2 at ' + str(lambda2) + ' nm', 
             color = '#e20000', marker = '*', 
             markerfacecolor = '#e20000')
    
    plt.suptitle('Binding isotherm for ' + Host + ' + '+ Guest 
                 + ', '+ str(T) +' \u00b0' +'C ', fontsize = 20)
    plt.title('Replicate' + str(i), fontsize = 16)
    plt.subplots_adjust(bottom=0.20, wspace = 0.27) # Slider 
    plt.xlabel('Equivalences [G]0/[H]0', fontsize = 14)
    plt.ylabel('Fluorescence', fontsize = 14)
    #np.figlegend to print one legend instead of 3
    plt.grid
    plt.legend(loc = 'upper left', fontsize = 10, handlelength = 1.0, 
               bbox_to_anchor=(0,-0.15))
        
# 5) Fit data 
    TotalV = V_H_titrated + vars()[volumetitr] * 10**(-3) # mL
    # Takes into account each addition gives extra volume
    # Starting concentration host
    cH0 = vars()[Hmolestitr]/ (TotalV*10**(-3))  
    # Starting concentration guest
    cG0 = vars()[Gmolestitr] / (TotalV*10**(-3)) 
    G0H0 = cG0/cH0 #ratio on x-axis; n/n = c/c
    guestspace = np.linspace(0, cG0[-1], num = 1000)
    
    Hmoles_av = np.average(cH0) 
    # Takes average of amount of moles Host in sample. 
    # Formula's assume host constant over the titration period. 
     
    yvalues1 = vars()[maxima_peak1]
    yvalues2 = vars()[maxima_peak2]
    # Intensity at peak maxima of Hsolution as initial guess of F0
    F0_1 = vars()[maxima_peak1][1] 
    F0_2 = vars()[maxima_peak2][1] 
    print('F0_peak1 = ' + str(F0_1)) # Check if same as fitted values
    print('F0_peak2 = ' + str(F0_2))
    eps = 10**(-15) #Trick to keep H0 constant
    eps2 = 10**(-10)#Trick to keep F0 constant
    
    # Define fitting parameters: Guess and bounds
    if formula == 'SternVolmer':
        p0_1 = [Hmoles_av, F0_1, 100000]
        p0_2 = [Hmoles_av, F0_2, 100000]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1),
                       ( Hmoles_av +eps, F0_1 + eps2,10**9))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1),
                       ( Hmoles_av +eps, F0_2 + eps2,10**9))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1),
                       ( Hmoles_av +eps, np.inf ,10**9))
            bounds2 = bounds1
    
    if formula == 'Tsukube': #G0, H0, F0, Ka, kDHG
        p0_1 = [Hmoles_av, F0_1, 3500000,-10**8] 
        p0_2 = [Hmoles_av, F0_2, 3500000,-10**8]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,-np.inf),
                       ( Hmoles_av +eps, F0_1 + eps2 ,10**9,np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,-np.inf),
                       ( Hmoles_av +eps, F0_2 + eps2 ,10**9,np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,-np.inf),
                       ( Hmoles_av +eps, np.inf ,10**9,np.inf))
            bounds2 = bounds1

        
    if formula == 'Connors': #G0, H0, F0, Ka, kH, kH0, kHG
        p0_1 = [Hmoles_av, F0_1, 3500000,1,1,0]
        p0_2 = [Hmoles_av, F0_2, 3500000,1,1,0]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 -eps2, 1,-np.inf,
            -np.inf,-np.inf),( Hmoles_av +eps, F0_1 + eps2 ,10**9,
            np.inf, np.inf,np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 -eps2, 1,-np.inf,
            -np.inf,-np.inf),( Hmoles_av +eps, F0_2 + eps2 ,10**9,
            np.inf,np.inf,np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,-np.inf, -np.inf, 
            -np.inf), (Hmoles_av +eps, np.inf ,10**9,np.inf, 
            np.inf, np.inf))
            bounds2 = bounds1
            
    if formula == 'binding_HGG': #H0, F0, K1, K2, kDHG, kDHG2
        p0_1 = [Hmoles_av, F0_1, 10**8,10**8,-10**14,-1000000]
        p0_2 = [Hmoles_av, F0_2, 10**8,10**8,-10**14,-1000000]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,1,-np.inf,
            -np.inf),(Hmoles_av +eps, F0_1 + eps2 ,10**9,10**9, 
            np.inf, np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,1,-np.inf,
            -np.inf),( Hmoles_av +eps, F0_2 + eps2 ,10**9,10**9, 
            np.inf, np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,1,-np.inf, -np.inf),
               ( Hmoles_av +eps, np.inf ,10**9,10**9, np.inf, np.inf))
            bounds2 = bounds1
    
    if formula == 'binding_HHG':
        p0_1 = [Hmoles_av, F0_1, 100000,100000,1,1]
        p0_2 = [Hmoles_av, F0_2, 100000,100000,1,1]
        if F0_free == False:
            bounds1 = ((Hmoles_av - eps, F0_1 - eps2, 1,1,-np.inf, 
            -np.inf),( Hmoles_av +eps, F0_1 + eps2 ,10**9,10**9, 
            np.inf, np.inf))
            bounds2 = ((Hmoles_av - eps, F0_2 - eps2, 1,1,-np.inf, 
            -np.inf),( Hmoles_av +eps, F0_2 + eps2 ,10**9,10**9, 
            np.inf, np.inf))
        if F0_free == True:
            bounds1 = ((Hmoles_av - eps, 0, 1,1,-np.inf, -np.inf),
               ( Hmoles_av +eps, np.inf ,10**9,10**9, np.inf, np.inf))
            bounds2 = bounds1

    # Find curve fitting parameters
    popt1, pcov1 = curve_fit(vars()[formula], cG0[1:], 
                   yvalues1[1:], p0 = p0_1, bounds = bounds1)
    popt2, pcov2 = curve_fit(vars()[formula], cG0[1:], 
                   yvalues2[1:], p0 = p0_2, bounds = bounds2)
    
    # Store F0 and Ka values ---------------------------------------------------
    
    initialF[0,i-1] = popt1[1]
    initialF[1,i-1] = popt2[1]
    bindingconstant[0,i-1] = popt1[2]
    bindingconstant[1,i-1] = popt2[2]
    if formula == 'Tsukube':
        kDHGcollect[0, i-1] = popt1[3]
        kDHGcollect[1, i-1] = popt2[3]
    if formula == 'Connors':
        kHcollect[0, i-1] = popt1[3]
        kHcollect[1, i-1] = popt2[3]
        kH0collect[0, i-1] = popt1[4]
        kH0collect[1, i-1] = popt2[4]
        kHGcollect[0, i-1] = popt1[5]
        kHGcollect[1, i-1] = popt2[5]
        
    plt.figure(3)
    plt.subplot(1,n,i)
    plt.plot(G0H0[1:], yvalues1[1:],label = 'Peak 1 at ' + 
             str(lambda1) + ' nm')
    plt.plot(guestspace/Hmoles_av, vars()[formula](
             guestspace, *popt1), label = 'fit: F_0 = ' + 
             str(round(popt1[1],1)) + u',K\u2090 = ' + 
             str(round(popt1[2],1)))
    plt.plot(G0H0[1:], yvalues2[1:],label = 'Peak 2 at ' + 
             str(lambda2) + ' nm')
    plt.plot(guestspace/Hmoles_av, vars()[formula](
             guestspace, *popt2),label = 'fit: F_0 = ' + 
             str(round(popt2[1],1)) + u', K\u2090 = ' 
             + str(round(popt2[2],1)))
    
    plt.suptitle('Fluorescence binding curve '+ formula + 
    ', for ' + Host + ' + '+ Guest + ', '+ str(T) +' \u00b0' +
    'C ', fontsize = 20)
    plt.title('Replicate' + str(i), fontsize = 16)
    plt.subplots_adjust(bottom=0.29, wspace = 0.27) # Slider
    plt.xlabel('[G]0/[H]0', fontsize = 14)
    plt.ylabel('F', fontsize =14)
    plt.grid
    plt.legend(loc = 'upper left', fontsize = 10, 
    handlelength = 1.0, bbox_to_anchor=(0,-0.15))

# Show all plots
plt.show(1)
plt.show(2)
plt.show(3)

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
print('{:.2e}'.format(np.average(initialF[0])) + u" \u00B1 " 
      + '{:0=3.2f}e+03'.format(np.std(initialF[0]/10**3)))
print('The average F0 for peak 2 is:')
print('{:.2e}'.format(np.average(initialF[1])) + u" \u00B1 " 
      + '{:0=3.2f}e+03'.format(np.std(initialF[1]/10**3)))
print('The average Ka for peak 1 is:')
print('{:.2e}'.format(np.average(bindingconstant[0]))+u" \u00B1 "
      + '{:0=3.2f}e+06'.format(np.std(bindingconstant[0]/10**6)))
print('The average Ka for peak 2 is:')
print('{:.2e}'.format(np.average(bindingconstant[1]))+u" \u00B1 " 
      + '{:0=3.2f}e+06'.format(np.std(bindingconstant[1]/10**6)))
print('The average DrG for peak 1 is:')
print('{:.3e}'.format(np.average(DrG[0])) + u" \u00B1 " + 
      '{:0=3.2f}e+03'.format(np.std(DrG[0]/10**3)))
print('The average DrG for peak 2 is:')
print('{:.3e}'.format(np.average(DrG[1])) + u" \u00B1 " + 
      '{:0=3.2f}e+03'.format(np.std(DrG[1]/10**3)))

## Generate outputfile
os.chdir('Directory' + Expnr)
if smoothening == True:
    file = open(Expnr + "_output_" + formula + '_' + method+ 
           '_smooth.txt', "w")
if smoothening == False:
    file = open(Expnr+ "_output_" + formula + '_' + method+ 
                '.txt', "w")
file.write('Experiment: ' + Expnr + '\n') 
file.write('Outputfile for binding study with fluorescence of ' 
           + Host + ' with ' + Guest + '. \n')
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
file.write('[H]0 = ' + str(round(c_Hstock*10**(9),3)) + ' uM\n')
file.write('Density = ' + str(round(dens_Hstock, 2)) + ' g/mL\n')
file.write('Guest-stock:\n')
file.write('[H]0 = ' + str(round(c_Gstock*10**(9),3)) + ' uM\n')
file.write('Density = ' + str(round(dens_Gstock, 2)) + ' g/mL\n')
for i in range(1,n+1):
    concentration = 'c_Hmeasure' + str(i)
    file.write('Host-measurement ' + str(i) + ':')
    file.write('\t[H]0 = '+str(round(vars()[concentration]*10**(9),3))
               +' uM \n')
for i in range(1,n+1):
    conG = 'c_Gmeasure' + str(i)
    conH = 'c_Hmeasure_forGmeas' + str(i)
    file.write('Guest-measurement ' + str(i) + ':')
    file.write('\t[H]0 = ' + str(round(vars()[conH]*10**(9),3)) 
               + ' uM \n')
    file.write('\t\t\t[G]0 = ' + str(round(vars()[connG]*10**(9),3)) 
               + ' uM \n')
    
file.write('\nTitration: \n')
file.write(' \n')
for i in range(1,n+1):
    add = 'Guest_add_' + str(i)
    total = 'Guest_total_' + str(i)
    ratio = 'Equivalences' + str(i)
    
    file.write('Set '+ str(i) + ':\n')
    file.write('Host measurement volume = ' + str(
                V_H_titrated *1000) +' uL \n')
    file.write('Added guest volume = ')
    file.writelines(str(vars()[add]) + '\n')
    file.write('Cumulative volume = ')
    file.writelines(str(vars()[total]) + '\n')
    file.write('Equivalences = ' + str(np.around(vars()[ratio],2)) 
               + '\n')
    file.write(' \n')

file.write('\nInitial guesses\n')
file.write('[H]0 = ' + str(Hmoles_av) + '\n')
file.write('F0_1 = ' + str(round(vars()[maxima_peak1][1],1)) 
           + '\n')
file.write('F0_2 = ' + str(round(vars()[maxima_peak2][1],1)) 
           + '\n')
if formula == 'SternVolmer':
    file.write('Ka = ' + str(p0_1[2]) + '\n')
if formula == 'Tsukube':
    file.write('Ka = ' + str(p0_1[2]) + '\n')
    file.write('kDHG = ' + str(p0_1[3]) + '\n')
if formula == 'Connors':
    file.write('Ka = ' + str(p0_1[2]) + '\n') 
    file.write('kH = ' + str(p0_1[3]) + '\n')
    file.write('kH0 = ' + str(p0_1[4]) + '\n')
    file.write('kHG = ' + str(p0_1[5]) + '\n')

file.write('\nOutput data: fitted parameters \n')
if formula == 'SternVolmer':
    for i in range(1,n+1): 
        maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
        maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
        file.write('Set ' + str(i) + ':\n')
        file.write('\t F0(fixed) \t F0(fitted) \t Ka \n')
        file.write('Peak 1   ' + str(round(vars()[maxima_peak1][1],1))
         + '\t \t ' + str(round(initialF[0,i-1],1)) + '\t \t ' 
         + str(round(bindingconstant[0,i-1],0)) + '\n')
        file.write('Peak 2   ' + str(round(vars()[maxima_peak2][1],1))
         + '\t \t ' + str(round(initialF[1,i-1],1)) + '\t \t ' 
         + str(round(bindingconstant[1,i-1],0)) + '\n')
        file.write(' \n')

if formula == 'Tsukube':
    for i in range(1,n+1):
        maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
        maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
        file.write('Set ' + str(i) + ':\n')
        file.write('\t F0(fixed) \t F0(fitted) \t Ka \t \t ' +
                   'kDHG (*10^9) \n')
        file.write('Peak 1   ' + str(round(vars()[maxima_peak1][1],1))
         + '\t \t ' + str(round(initialF[0,i-1],1)) + '\t \t ' 
         + str(round(bindingconstant[0,i-1],1)) + '\t ' + 
         str(round(kDHGcollect[0,i-1]/10**9,3)) + '\n')
        file.write('Peak 2   ' + str(round(vars()[maxima_peak2][1],1))
         + '\t \t ' + str(round(initialF[1,i-1],1)) + '\t \t ' + 
        str(round(bindingconstant[1,i-1],1)) + '\t ' + 
        str(round(kDHGcollect[1,i-1]/10**9,3)) + '\n')
        file.write(' \n')

if formula == 'Connors':
    for i in range(1,n+1):
        maxima_peak1 = 'Maximum_absorbance_peak1' + str(i)
        maxima_peak2 = 'Maximum_absorbance_peak2' + str(i)
        kH1 = '{:{align}{width}{prec}f}'.format(kHcollect[0,i-1], 
                align= '<', width = '8', prec= '.3')
        kH2 = '{:{align}{width}{prec}f}'.format(kHcollect[1,i-1], 
                align= '<', width = '8', prec= '.3')
        kH01 = '{:{align}{width}{prec}f}'.format(kH0collect[0,i-1], 
                 align= '<', width = '8', prec= '.3')
        kH02 = '{:{align}{width}{prec}f}'.format(kH0collect[1,i-1], 
                 align= '<', width = '8', prec= '.3')
        
        file.write('Set ' + str(i) + ':\n')
        file.write('\t F0(fixed) \t F0(fitted) \t Ka \t \t kH 
                    \t        kH0 \t \t kHG \n')
        file.write('Peak 1   ' + str(round(vars()[maxima_peak1][1],1))
         + '\t \t ' + str(round(initialF[0,i-1],1)) + '\t \t ' 
         + str(round(bindingconstant[0,i-1],1)) + '\t ' + kH1 + 
         '\t' + kH01 + '\t ' + str(round(kHGcollect[0,i-1],3)) +'\n') 
        file.write('Peak 2   ' + str(round(vars()[maxima_peak2][1],1))
         + '\t \t ' + str(round(initialF[1,i-1],1)) + '\t \t ' 
         + str(round(bindingconstant[1,i-1],1)) + '\t ' + kH2 + 
         '\t' + kH02 + '\t ' + str(round(kHGcollect[1,i-1],3)) +'\n')
        file.write(' \n') 

file.write('\nAverages \n')
file.write('Peak 1: \n')
file.write('F0 = ' + '{:.2e}'.format(np.average(initialF[0])) 
           + u" \u00B1 " + '{:0=3.2f}e+03'.format(
           np.std(initialF[0]/10**3)) + '\n')
file.write('Ka = ' + '{:.2e}'.format(np.average(bindingconstant[0]))
           +u" \u00B1 "+ '{:0=3.2f}e+06'.format(
           np.std(bindingconstant[0]/10**6)) + '\n')
file.write('DrG = ' + '{:.3e}'.format(np.average(DrG[0])) + u" \u00B1 " + 
           '{:0=3.2f}e+03'.format(np.std(DrG[0]/10**3) + 'J/mol')
file.write('Peak 2: \n')
file.write('F0 = ' + '{:.2e}'.format(np.average(initialF[0])) 
           + u" \u00B1 " + '{:0=3.2f}e+03'.format(
           np.std(initialF[0]/10**3)) + '\n')
file.write('Ka = ' + '{:.2e}'.format(np.average(bindingconstant[1]))
           +u" \u00B1 " + '{:0=3.2f}e+06'.format(
           np.std(bindingconstant[1]/10**6)) + '\n')
file.write('DrG = ' + '{:.3e}'.format(np.average(DrG[1])) + u" \u00B1 " + 
           '{:0=3.2f}e+03'.format(np.std(DrG[1]/10**3) + 'J/mol')

file.close()