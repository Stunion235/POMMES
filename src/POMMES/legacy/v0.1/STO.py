import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

fig, ax1 = plt.subplots()
all_int_val = []

def calc_structure_factor(L1, L2, ratio_sampledata,color_line):
    
    # step size can be made smaller or larger
    step_size = 0.005
    # plots and calculates for occupancy values of 0.5 to 1.0
    occupancy_Sr1 = 0.5
    end_occ_Sr1 = 1.0
    occ = []
    ratios = []

    while (occupancy_Sr1 < (end_occ_Sr1 + step_size)): 
        def intensities_func(h,k,l):
            # wavelength for xrd; constant
            lamb = 1.5406
            
            # lattice constants 
            # change depending on film
            a = 3.92574
            b = 3.92574
            c_1 = 7.92937
            
            # calculate d
            d = 1 / (math.sqrt((h/a)**2 + (k/b)**2 + (l/c_1)**2))
        
            # calculate atomic form factors
            def approx(a1,b1,a2,b2,a3,b3,a4,b4,c):
                Q = (2*math.pi) / d
                sum1 = a1 * math.e**(-1 * b1 * (Q/(4*math.pi))**2) 
                sum2 = a2 * math.e**(-1 * b2 * (Q/(4*math.pi))**2)
                sum3 = a3 * math.e**(-1 * b3 * (Q/(4*math.pi))**2)
                sum4 = a4 * math.e**(-1 * b4 * (Q/(4*math.pi))**2)
            
                f_approx = sum1+sum2+sum3+sum4+c
                return f_approx
            
            f_Sr_approx = approx(18.0874,1.4907,8.1373,12.6963,2.5654,24.5651,-34.193,-0.0138,41.4025)
            f_Ta_approx = approx(29.1587,1.50711,18.8407,0.116741,12.8268,6.31524,5.38695,12.4244,1.78555)
            f_O_approx = approx(3.0485,13.2771,2.2868,5.7011,1.5463,0.3239,0.867,32.9089,0.2508)

            # take Ta and O displacements into account
            # do this for each distinct atom in the atomic basis
            def Ta_disp(occSr1):
                Ta_z_initial = 0.25
                # 0.0095 is film-dependent difference between Ta position at occ = 0.5 and occ = 1.0
                Ta_z = Ta_z_initial + (0.0095/0.50)*(occupancy_Sr1-0.5)
                return Ta_z
            
            def O1_disp(occSr1):
                O1_z_initial = 0.25
                # -0.02008 is film-dependent difference between O1 position at occ = 0.5 and occ = 1.0
                O1_z = O1_z_initial + (-0.02008/0.50)*(occupancy_Sr1-0.5)
                return O1_z
            def O3_disp(occSr1):
                O3_z_initial = 0.25
                # -0.02001 is film-dependent difference between O3 position at occ = 0.5 and occ = 1.0
                O3_z = O3_z_initial + (-0.02001/0.50)*(occupancy_Sr1-0.5)
                return O3_z
            
            # calculate structure factor
            def structure_factor(fSr,fTa,fO,occSr1):
                # all are film-dependent, so this should change if not calculating ETO film values
                # the following numbers make up the ETO atomic position basis
                
                Sr1_xyz = [1,1,0]
                Sr2_xyz = [1,1,0.5]
                
                Ta1_xyz = [0.5,0.5,Ta_disp(occSr1)]
                Ta2_xyz = [0.5,0.5,1-Ta_disp(occSr1)]
                
                O1_xyz = [0.5,1,O1_disp(occSr1)]
                O2_xyz = [0.5,0.0,1-O1_disp(occSr1)]
                O3_xyz = [0.0,0.5,O3_disp(occSr1)]
                O4_xyz = [1,0.5,1-O3_disp(occSr1)]
                O5_xyz = [0.50000,0.50000,1.00000]
                O6_xyz = [0.50000,0.50000,0.50000]
        
                factor = (-1j*2*math.pi)
                
                # calculate the terms that will sum to the structure factor

                Sr1_S_term = occSr1 * fSr * math.e ** (factor * (h * Sr1_xyz[0] + k * Sr1_xyz[1] + l * Sr1_xyz[2]))
                Sr2_S_term = (1-occSr1) * fSr * math.e ** (factor * (h * Sr2_xyz[0] + k * Sr2_xyz[1] + l * Sr2_xyz[2]))
                
                Ta1_S_term = fTa*math.e**(factor*(h*Ta1_xyz[0]+k*Ta1_xyz[1]+l*Ta1_xyz[2]))
                Ta2_S_term = fTa*math.e**(factor*(h*Ta2_xyz[0]+k*Ta2_xyz[1]+l*Ta2_xyz[2]))
        
                O1_S_term = fO*math.e**(factor*(h*O1_xyz[0]+k*O1_xyz[1]+l*O1_xyz[2]))
                O2_S_term = fO*math.e**(factor*(h*O2_xyz[0]+k*O2_xyz[1]+l*O2_xyz[2]))
                O3_S_term = fO*math.e**(factor*(h*O3_xyz[0]+k*O3_xyz[1]+l*O3_xyz[2]))
                O4_S_term = fO*math.e**(factor*(h*O4_xyz[0]+k*O4_xyz[1]+l*O4_xyz[2]))
                O5_S_term = fO*math.e**(factor*(h*O5_xyz[0]+k*O5_xyz[1]+l*O5_xyz[2]))
                O6_S_term = fO*math.e**(factor*(h*O6_xyz[0]+k*O6_xyz[1]+l*O6_xyz[2]))
                
                # put all terms in a list
                all_func = [Sr1_S_term,Sr2_S_term,Ta1_S_term,Ta2_S_term,O1_S_term,O2_S_term,O3_S_term,O4_S_term,O5_S_term,O6_S_term]
                
                # calculate structure factor ("sum")
                sum = 0
                for i in all_func:
                    sum+=i
                return sum
        
            S = structure_factor(f_Sr_approx,f_Ta_approx,f_O_approx,occupancy_Sr1)
            
            # lorentz polarization factor
            def ang_dep():
                theta_rad = math.asin(1.5406/(2*d))
                twotheta_rad = theta_rad*2
                Lp = (1 + (math.cos(twotheta_rad) ** 2)) / ((math.sin(theta_rad)) ** 2 * math.cos(theta_rad))
                return(Lp)
            
            # multiplicity
            # the method for calculating this CHANGES depending on symmetry / type of crystal
            # so this may be different for non-ETO films
            def mult_orth():
                count = 0
                h_k_l = [h, k, l]
                for i in h_k_l:
                    if i > 0:
                        count += 1
                M = 2 ** count
                return M
        
            # S can be imaginary number so magnitude is found by using abs() function
            S_mag = abs(S)
            
            # temperature factor has been ignored in this code
            
            theta_rad = math.asin(1.5406/(2*d))
            twotheta_rad = theta_rad*2
            theta_M_rad = 22.65*(math.pi/180)
            Int = (S_mag**2) * (1+ (math.cos(2*theta_M_rad)**4)*(math.cos(twotheta_rad)**2))/(math.sin(theta_rad)*math.cos(theta_rad))
            return Int
            
            '''
            Int = (S_mag ** 2) * mult_orth() * ang_dep()
            return Int
            '''
            
            
        Int_hkl1 = intensities_func(0,0,L1)
        Int_hkl2 = intensities_func(0,0,L2)
        
        # calculate ratio
        ratio = Int_hkl1/Int_hkl2
        
        occ.append(occupancy_Sr1)
        ratios.append(ratio)
        occupancy_Sr1 += step_size
        
        
        
        
            
    plt.rcParams['figure.dpi'] = 1200
    ratios_sampledata1 = [ratio_sampledata]
    inter_func1 = interp1d(ratios,occ, kind='nearest')
    new = inter_func1(ratios_sampledata1)
    newarr = np.array(new)
    
    '''
    with open('file_path_here' + str(L1) + '_2_to_L_v2.xy', 'a') as file:
        # Iterate through the data list and write each item to the file
        for i in occ:
            file.write(str(i) + '\n')
        file.write('\n')
        for i in ratios:
            file.write(str(i) + '\n')
        file.write('\n')
        for i in new:
            file.write(str(i) + '\n')
        file.write('\n')
        for i in ratios_sampledata1:
            file.write(str(i) + '\n')
        file.write('\n')
    '''
    # different ylim possibilities - play around based on specific L1 and film
    
    plt.ylim(-0.1,3)
    #plt.gca().set_yticklabels([])
    #plt.gca().set_xticklabels([])
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')

    # Optionally, you can customize the tick direction and length
    ax1.tick_params(axis='both', direction='in')

    #plt.ylim(-0.2,6)
    
    #order_param = 2*occ - 1;
   
    ax1.set_ylabel('Ratio of 00' + str(L1) + '/00L Peak Intensities')
    ax1.set_xlabel('Occupancy of Sr1 Layer')
    ax1.plot(occ, ratios, color=color_line, label='00' + str(L2))
    ax1.set_xlim(0.48,1.02)
   
    # Create ax2 as a twin of ax1
    #ax2 = ax1.twiny()
    #ax2.set_ylabel('Order Parameter')  # Use 'set_ylabel' for the y-axis label on ax2
    #ax2.set_xlim(0, 1.0)  # Set the x-axis limits for ax2
    #ax2.set_xscale('linear')
    
    ax1.scatter(new, ratios_sampledata1, color='black')
    
    print('interpolated occupancy=' + str(new) + ', ratio=' + str(ratios_sampledata1))
    plt.legend()

    # find average interpolated occupancy
    all_int_val.append(new[0])

    
    #return int(ratios_sampledata1[0])
  

# first parameter is L1 in hk L1 : hk L2 ratio 
# second parameter is L2 in hk L1 : hk L2 ratio
# third parameter is sample data intensity ratio (hk L1 : hk L2) for interpolation
# fourth parameter is color of structure factor line to be plotted

# 1 for 1/2, 3 for 3/2 order peaks

#Usage: calc_structure_factor(L1,L2,ratio,color) where ratio is the 00(L1)/00(L2) intensity ratio
L1=1
calc_structure_factor(L1, 2, 64/2914,'magenta')
calc_structure_factor(L1, 4, 64/1874,'purple')
calc_structure_factor(L1, 6, 64/81,'green')

def Average(lst): 
    return sum(lst) / len(lst) 

print("Average: "+str(Average(all_int_val)))

plt.title('Occupancy Based on 00' + str(L1) +'/00L Intensities')
#plt.savefig('1_2_ETO0031_v4.svg')
plt.show()
