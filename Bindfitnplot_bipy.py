'''
This script summarizes the adjustments in the Bindfitnplot.py
script for complexes which do not quench upon addition of
more guest and instead have an increased and shifted intensity.
'''

# Line 99
# Follow around wavelength
lambda1 = 643 #nm
lambda2 = 708 #nm


# Line 106:
formula = 'Connors'  # Fitting: Tsukube or Connors

# Line 532:
# Find curve fitting parameters
if method == 'single_wl':
    popt1, pcov1 = curve_fit(vars()[formula], cG0[1:],
                   yvalues1[1:], p0 = p0_1, bounds = bounds1)
    popt2, pcov2 = curve_fit(vars()[formula], cG0[1:],
                   yvalues2[1:], p0 = p0_2, bounds = bounds2)
if method == 'max_wl':
    popt1, pcov1 = curve_fit(vars()[formula], cG0[12:],
                   yvalues1[12:], p0 = p0_1, bounds = bounds1)
    popt2, pcov2 = curve_fit(vars()[formula], cG0[12:],
                   yvalues2[12:], p0 = p0_2, bounds = bounds2)

# Line 474
p0_1 = [Hmoles_av, F0_1, 250000,10**9]
p0_2 = [Hmoles_av, F0_2, 250000,10**9]
# The last entry of both p0_1 and p0_2 is positive instead
# of negative. Also these values might need to be altered for
# each experiment, in order to obtain a correct fit, based
# on figure 3.

