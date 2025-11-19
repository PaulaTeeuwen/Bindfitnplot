'''
This script summarizes the adjustments in the Bindfitnplot.py
script if 3 stocks for both host and guest are used.
'''

# lines 63 to 73 -------------------------------------------

# Host-stock 1
M_Host = 
m_Hstock1 = 
V_solv_Hstock = 
m_solv_Hstock1 = 

# Host-stock 2
m_Hstock2 = 
m_solv_Hstock2 = 

# Host-stock 3
m_Hstock3 = 
m_solv_Hstock3 = 

# Guest-stock 1
M_Guest = 
m_Gstock1 = 
V_solv_Gstock = 
m_solv_Gstock1 = 

# Guest-stock 2
m_Gstock2 = 
m_solv_Gstock2 = 

# Guest-stock 3
m_Gstock3 = 
m_solv_Gstock3 = 

# lines 229 to 243 ------------------------------------------

# Host-stock 1
n_Hstock1 = m_Hstock1 / M_Host 
c_Hstock1 = n_Hstock1 / V_solv_Hstock
dens_Hstock1 = m_solv_Hstock1 / V_solv_Hstock
print('Molarity host stock solution => ' + 
       str(round(c_Hstock1*10**(9),2)) + ' \u03BCM')

# Host-stock 2
n_Hstock2 = m_Hstock2 / M_Host 
c_Hstock2 = n_Hstock2 / (V_solv_Hstock) 
dens_Hstock2 = m_solv_Hstock2 / V_solv_Hstock 
print('Molarity host stock solution => ' + 
      str(round(c_Hstock2*10**(9),2)) + ' \u03BCM')

# Host-stock 3
n_Hstock3 = m_Hstock3 / M_Host
c_Hstock3 = n_Hstock3 / (V_solv_Hstock) 
dens_Hstock3 = m_solv_Hstock3 / V_solv_Hstock 
print('Molarity host stock solution => ' + 
      str(round(c_Hstock3*10**(9),2)) + ' \u03BCM')

# Guest-stock 1
n_Gstock1 = m_Gstock1 / M_Guest
c_Gstock1 = n_Gstock1/ (V_solv_Gstock )
dens_Gstock1 = m_solv_Gstock1 / V_solv_Gstock
print('Guest stock solution => ' + 
       str(round(c_Gstock1*10**(9),2)) + ' \u03BCM')

# Guest-stock 2
n_Gstock2 = m_Gstock2 / M_Guest 
c_Gstock2 = n_Gstock2/ (V_solv_Gstock ) 
dens_Gstock2 = m_solv_Gstock2 / V_solv_Gstock 
print('Guest stock solution => ' + 
      str(round(c_Gstock2*10**(9),2)) + ' \u03BCM')

# Guest-stock 3
n_Gstock3 = m_Gstock3 / M_Guest 
c_Gstock3 = n_Gstock3/ (V_solv_Gstock ) 
dens_Gstock3 = m_solv_Gstock3 / V_solv_Gstock 
print('Guest stock solution => ' + 
      str(round(c_Gstock3*10**(9),2)) + ' \u03BCM')

dens_Hstock = np.array([dens_Hstock1, dens_Hstock2, dens_Hstock3])
c_Hstock = np.array([c_Hstock1, c_Hstock2, c_Hstock3])
dens_Gstock = np.array([dens_Gstock1, dens_Gstock2, dens_Gstock3])
c_Gstock = np.array([c_Gstock1, c_Gstock2, c_Gstock3])

# line 254 -------------------------------------------------
    vars()[volume] = vars()[mass] / dens_Hstock[j-1]   

# line 258
    vars()[concentration] = vars()[dilution] * c_Hstock [j-1]

# line 273 -------------------------------------------------
    vars()[volume2] = vars()[mass2] / dens_Gstock[j-1] 
    
# line 277 -------------------------------------------------
    vars()[concentrationG] = vars()[dilution3] * c_Gstock[j-1]
    
# line 286 -------------------------------------------------
    vars()[volume3] = vars()[mass3] / dens_Hstock[j-1]

# line 290 -------------------------------------------------  
    vars()[concentrationH] = vars()[dilution3] * c_Hstock[j-1]

