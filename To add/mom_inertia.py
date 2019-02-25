# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 09:09:25 2019

@author: Wesley
"""

def mom_inert(stif_z, stif_y, angle):
    """"Enter lists with the z and y co-ords of the stiffeners and the angle\
    of the skin in the straight section."""
    import numpy as np
    
    ail_h   = 205                   # height of the aileron     [mm]
    ail_w   = 605                   # chord of the aileron      [mm]
    
    # we'll start off by finding moment of inertia for a stiffener
    w       = 19            	    # width of stiffener        [mm]
    h       = 16                    # height of stiffener       [mm]
    t       = 1.2                   # thickness of stiffener    [mm]
    
    A_base  = w*t                   # area of the horiz. part   [mm^2]
    A_top   = (h-t)*t               # area of the vert. part    [mm^2]
    A_stif  = A_base + A_top        # area of stiffener         [mm^2]
    
    # here are the z co-ordinates of the stiffeners             [mm]
#    stif_z  = [560.2076385798177, 560.2076385798177, 470.6229157394533, \
#               470.6229157394533, 381.0381928990888, 381.0381928990888, \
#               291.45347005872435, 291.45347005872435, 201.86874721835989, \
#               201.86874721835989, 112.28402437799541, 112.28402437799541, \
#               37.51424233216027, 37.51424233216027, 0]
    
    # same in the y-axis                                        [mm]
#    stif_y  = [8.869774538649953, -8.869774538649953, 26.609323615949847, \
#               -26.609323615949847, 44.348872693249746, -44.348872693249746, \
#               62.08842177054964, -62.08842177054964, 79.82797084784953, \
#               -79.82797084784953, 97.56751992514943, -97.56751992514943, \
#               78.07387583997614, -78.07387583997614, 0]
    
#    angle   = 0.2012196329874848    # angle of skin             [rad]
    c_y     = 0                     # centroid y co-ord (symm)  [mm]
    c_z     = 263.5208490384224     # centroid z co-ord         [mm]
    stif_zz = 0                     # initial value mom inert   [mm^4]
    stif_yy = 0                     # initial value mom inert   [mm^4]
    
    # now calculate the Steiner terms for all stiffeners
    for i in range(len(stif_z)):
        stein_z = A_stif*(stif_z[i]-c_y)**2
        stif_zz = stif_zz + stein_z
    
    for i in range(len(stif_y)):
        stein_y = A_stif*(stif_y[i]-c_z)**2
        stif_yy = stif_yy + stein_y
    
    print ('---------------------------------------------------------------------')
    print ('Moment of inertia due to the Steiner term of stiffeners:')
    print ('I_u =', '%e' % stif_zz, 'mm^4, and I_v =', '%e' % stif_yy, 'mm^4')
    print ('---------------------------------------------------------------------')
    
    skin_z  = (ail_w - (ail_h/2))/2                             # skin centroid         [mm]
    skin_t  = 1.1                                               # skin thickness        [mm]
    skin_l  = np.sqrt((ail_w - (ail_h/2))**2 + (ail_h/2)**2)    # skin length           [mm]
    
    # add the moment of inertia for the straight part of the skin
    skin_zz     = ((skin_t*skin_l**3)*np.cos(angle)**2)/12
    skin_yy     = ((skin_t*skin_l**3)*np.sin(angle)**2)/12
    skin_stein  = skin_t*skin_l*(skin_z-c_z)**2
    
    print ('Moment of inertia due to the straight part of the skin:')
    print ('I_u =', '%e' % skin_zz, 'mm^4, and I_v =', '%e' % (skin_yy+skin_stein), 'mm^4')
    print ('---------------------------------------------------------------------')
    
    z_out       = (ail_h/2) - (4*ail_h/2)/(3*np.pi)             # centroid outer        [mm]
    z_in        = (ail_h/2) - (4*(ail_h-skin_t)/2)/(3*np.pi)    # centroid inner        [mm]
    A_out       = (np.pi*(ail_h/2)**2)/2                        # area outer circ       [mm^2]
    A_in        = (np.pi*((ail_h-skin_t)/2)**2)/2               # area inner circ       [mm^2]
    circ_z      = (A_out*z_out + A_in*z_in)/(A_out+A_in)        # centroid z            [mm]
    
    # add the moment of inertia for the circular part of the skin
    circ_zz     = A_out*(ail_h/2)**2/4 - A_in*((ail_h-skin_t)/2)**2/4
    circ_yy     = A_out*(ail_h/2)**2/4 - A_in*((ail_h-skin_t)/2)**2/4
    circ_stein  = (A_in+A_out)*(circ_z-c_z)**2
    
    print ('Moment of inertia due to the circular part of the skin:')
    print ('I_u =', '%e' % circ_zz, 'mm^4, and I_v =', '%e' % (circ_yy+circ_stein), 'mm^4')
    print ('---------------------------------------------------------------------')
    
    spar_z      = ail_h/2           # centroid of the spar      [mm]
    spar_h      = ail_h             # height of the spar        [mm]
    spar_t      = 2.8               # thickness of spar         [mm]
    
    # finally, add the moment of inertia for the spar
    spar_zz     = (spar_h*spar_t**3)/12
    spar_yy     = (spar_t*spar_h**3)/12
    spar_stein  = spar_h*spar_t*(spar_z-c_z)**2
    
    print ('Moment of inertia due to the spar:')
    print ('I_u =', '%e' % spar_zz, 'mm^4, and I_v =', '%e' % (spar_yy+spar_stein), 'mm^4')
    print ('---------------------------------------------------------------------')
    
    # then add all the moment of inertia together
    I_u         = stif_zz + skin_zz + circ_zz + spar_zz
    I_v         = stif_yy + skin_yy + skin_stein + circ_yy + circ_stein + spar_yy + spar_stein
    
    print ('Total moment of inertia:')
    print ('I_u =', '%e' % I_u, 'mm^4, and I_v =', '%e' % I_v, 'mm^4')
    print ('---------------------------------------------------------------------')
    
    # now switch to the new axis system
    phi         = np.radians(28)    # aileron rot angle         [rad]
    I_zz        = (I_u + I_v)/2 + ((I_u - I_v)/2)*np.cos(2*phi)
    I_yy        = (I_u + I_v)/2 - ((I_u - I_v)/2)*np.cos(2*phi)
    
    print ('Total moment of inertia in the rotated axis system:')
    print ('I_zz =', '%e' % I_zz, 'mm^4, and I_yy =', '%e' % I_yy, 'mm^4')
    print ('---------------------------------------------------------------------')
    
    return (I_zz, I_yy)

###############################################################################
############################# IGNORE THE NEXT BIT #############################
####################### ONLY NEEDED FOR TESTING PURPOSES ######################
###############################################################################

stif_z  = [560.2076385798177, 560.2076385798177, 470.6229157394533, \
           470.6229157394533, 381.0381928990888, 381.0381928990888, \
           291.45347005872435, 291.45347005872435, 201.86874721835989, \
           201.86874721835989, 112.28402437799541, 112.28402437799541, \
           37.51424233216027, 37.51424233216027, 0]

stif_y  = [8.869774538649953, -8.869774538649953, 26.609323615949847, \
           -26.609323615949847, 44.348872693249746, -44.348872693249746, \
           62.08842177054964, -62.08842177054964, 79.82797084784953, \
           -79.82797084784953, 97.56751992514943, -97.56751992514943, \
           78.07387583997614, -78.07387583997614, 0]

angle   = 0.2012196329874848