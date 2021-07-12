#%%

import math

# CONSTANTS (Table 1)

#matboard
global sigmaC_matboard
global sigmaT_matboard
global mass_matboard
global dimensions_matboard
global shear_strength_matboard
global youngs_mod
global poissons_ratio

#contact cement
global shear_strength_cement

#misc
global distance_joint_AL
global distance_AL
global P
global t
global E
global b1,b2,b3
global b

#matboard
sigmaC_matboard = 6.0 #MPa
sigmaT_matboard = 30.0 #MPa
mass_matboard = 0.75 #kg
dimensions_matboard = 813.0 * 1016.0 * 1.27 #mm**3
shear_strength_matboard = 4.0 #MPa
poissons_ratio = 0.2

#contact cement
shear_strength_cement = 2.0 #MPa

#misc
distance_RS_AL = 280.0 #mm       NOTE /!\: AL = applied load(s)
distance_AL = 390.0 #mm
P = 1000.0 #N
t = 1.27 #mm #thickness of the cross section
E = 4000.0 #MPa
b1 = 100.0 #mm
b2 = t #mm
b3 = 10.0 #mm
b = 60.0 #mm
length_of_bridge = 950.0 #mm

#variables
h1 = float(input("Enter h1."))
h2 = float(input("Enter h2."))


#=============================================================================================================================================#
###############################################################################################################################################
#=============================================================================================================================================#


#SPEC FINDER: gather the cross section specs to then be processed in the EXECUTION

def FindCrossSectionSpecs(h1,h2):
	global t
	global E
	global b1,b2,b3
	global b
	global P
	global distance_joint_AL
	global distance_AL

	#find total area of cross section

	area1 = b1*h1
	area2 = 2.0*b2*h2 
	area3 = 2.0*t*b3 #constant
	A = area1 + area2 + area3

	#find y_bar
	y_bar = (area1*(h2+h1/2) + area2*h2/2 + area3*(h2-t/2))/A

	#find I
	I1 = (b1*h1**3)/12
	I2 = (b2*h2**3)/12
	I3 = (b3*t**3)/12
	I = I1 + area1*(y_bar-(h2+h1/2))**2 + I2 + area2*(y_bar-h2/2)**2 + I3 + area3*(y_bar-(h2-t/2))**2

	# find y_c
	h_total = h1+h2
	y_c = h_total - y_bar

	#find y_t
	y_t = y_bar

	#find Q
	Q = t*y_bar**2 #simplified equation

	#FIND BEAM DEFLECTION
	deltaB = (P*distance_RS_AL**3)/(6*E*I) + (P*distance_RS_AL**2 * distance_AL)/(4*E*I) + (P*distance_RS_AL * distance_AL**2)/(16*E*I)
	
	return A, y_bar, I, y_c, y_t, Q, deltaB


#=============================================================================================================================================#
###############################################################################################################################################
#=============================================================================================================================================#


#DIAGRAMS FINDER: FIND SFD, BMD to find MMAX  	
#/!\ this is constant: it does not depend on h1, h2 and the SPECS found above

			# 			  |	 P/2		|  P/2
			# 			  v             v
			# ________________________________________
			# |										 |
			# |	 									 |
			# |______________________________________|
			# /\				/\					ooo
			#  ^				 ^					 ^
			#  |P_RS  			 |P_support		 	 |P_RS

def FindSupportReactionForces():
	global P
	global distance_RS_AL
	global distance_AL

	P_support = 3/(distance_RS_AL+distance_AL/2)**3 * ((P*distance_RS_AL**3)/3 + (P*distance_RS_AL**2 * \
		distance_AL)/2 + (P*distance_RS_AL * distance_AL**2)/8)
	P_RS = (P-P_support)/2

	return P_support,P_RS

#SFD
def SFD_left_P_RS_to_left_AL():
	global P_RS

	P_support,P_RS = FindSupportReactionForces()

	SFD_left_P_RS_to_left_AL_VAR = P_RS
	return SFD_left_P_RS_to_left_AL_VAR

def SFD_left_AL_to_P_support():
	global P

	SFD_left_P_RS_to_left_AL_VAR = SFD_left_P_RS_to_left_AL()
	SFD_left_AL_to_P_support_VAR = SFD_left_P_RS_to_left_AL_VAR - P/2
	return SFD_left_AL_to_P_support_VAR

def SFD_P_support_to_right_AL(): 
	P_support,P_RS = FindSupportReactionForces()

	SFD_left_AL_to_P_support_VAR = SFD_left_AL_to_P_support()
	SFD_P_support_to_right_AL_VAR = SFD_left_AL_to_P_support_VAR + P_support
	return SFD_P_support_to_right_AL_VAR

def SFD_right_AL_to_right_RS():
	global P

	SFD_P_support_to_right_AL_VAR = SFD_P_support_to_right_AL()
	SFD_right_AL_to_right_RS_VAR = SFD_P_support_to_right_AL_VAR - P/2
	return SFD_right_AL_to_right_RS_VAR

#BMD
def BMD_left_AL():
	global distance_RS_AL

	SFD_left_P_RS_to_left_AL_VAR = SFD_left_P_RS_to_left_AL()
	BMD_left_AL_VAR = SFD_left_P_RS_to_left_AL_VAR * distance_RS_AL
	return BMD_left_AL_VAR

def BMD_P_support():
	global distance_AL

	SFD_left_AL_to_P_support_VAR = SFD_left_AL_to_P_support()
	BMD_left_AL_VAR = BMD_left_AL()
	BMD_P_support_VAR = SFD_left_AL_to_P_support_VAR * distance_AL/2 + BMD_left_AL_VAR
	return BMD_P_support_VAR

def BMD_right_AL():
	global distance_AL

	SFD_P_support_to_right_AL_VAR = SFD_P_support_to_right_AL()
	BMD_P_support_VAR = BMD_P_support()
	BMD_right_AL_VAR = SFD_P_support_to_right_AL_VAR * distance_AL/2 + BMD_P_support_VAR
	return BMD_right_AL_VAR

def FindMmax():
	a = abs(BMD_left_AL())
	b = abs(BMD_P_support())
	c = abs(BMD_right_AL())
	Mmax = max(a,b,c)
	return Mmax


#=============================================================================================================================================#
################################################################################################################################################
#=============================================================================================================================================#


#EXECUTION: Use specs found in SPEC FINDER using the 11 summary points AND FIND VALUES OF P_FAIL for each type of failure

def FlexuralFailure(): #Navier's Equation for Flexural Stresses
	global sigmaC_matboard
	global sigmaT_matboard

	A, y_bar, I, y_c, y_t, Q, deltaB = FindCrossSectionSpecs(h1,h2)
	Mmax = FindMmax()
	#compression
	P_fail_C = (I*sigmaC_matboard)/((Mmax/P)*y_c)

	#tension
	P_fail_T = (I*sigmaT_matboard)/((Mmax/P)*y_t)

	P_fail_FF = min(P_fail_C,P_fail_T)
	print("Due to flexural stresses, the bridge will fail at:",P_fail_FF,"N.")
	return P_fail_FF

def ShearFailure(): #Jourawski’s Equation for shear stress			#instantiate b in the constants 
	global P_fail
	global shear_strength_matboard
	global t

	A, y_bar, I, y_c, y_t, Q, deltaB = FindCrossSectionSpecs(h1,h2)

	P_fail_zone = (2*shear_strength_matboard*I*2*t)/Q

	#find different b and Q /!\
	glue_shear_P_failure = (2*shear_strength_matboard*I*2*t)/Q

	P_fail_SF = min(P_fail_zone,glue_shear_P_failure) #find P_fail
	print("Due to shear forces, the bridge will fail at:",P_fail_SF,"N.")
	return P_fail_SF



def BucklingFailure(): #failure due to buckling
    ######################################supposed to work, but didn't#############################################
    global E
    global poissons_ratio
    global t
    global b
    global b1,b2,b3
    global distance_AL

    A, y_bar, I, y_c, y_t, Q, deltaB = FindCrossSectionSpecs(h1,h2)
    Mmax = FindMmax()

    plate_compression = (Mmax*y_t)/I

    comp_top_flange = ((4 * math.pi**2 * E) * ((h1/b)**2))/(12 * (1 - poissons_ratio**2)) #calculate this by hand
    P_fail_BC1 = comp_top_flange/plate_compression

    tips_flange = ((0.425 * math.pi**2 * E)*((h1/20)**2))/(12 * (1 - poissons_ratio**2)) #calculate this by hand
    P_fail_BC2 = tips_flange/plate_compression

    webs_flexural = ((6 * math.pi**2 * E) * ((t/y_c)**2))/(12 * (1 - poissons_ratio**2)) #calculate this by hand
    P_fail_webs_flex = webs_flexural / plate_compression

    webs_shear = ((5 * math.pi**2 * E) * ((t/h2)**2 + (t/(distance_AL/2))**2))/(12 * (1 - poissons_ratio**2)) #calculate this by hand
    P_fail_webs_shear = (2 * webs_shear * I * 2*t)/Q

    P_fail_BF = min(P_fail_BC1,P_fail_BC2,P_fail_webs_flex,P_fail_webs_shear)
    ########################################supposed to work, but didn't######################################

    print("Due to buckling, the bridge will fail at:",P_fail_BF,"N.")
    return P_fail_BF

#=============================================================================================================================================#
###############################################################################################################################################
#=============================================================================================================================================#


#OUTPUT: USING THE VALUES FROM THE EXECUTION, find lowest P_fail

def LowestP_fail():
	a = FlexuralFailure()
	b = ShearFailure()
	c = BucklingFailure() #insert buckling here

	lowestPfail = min(a,b,c)
	print("The bridge will fail at:",lowestPfail,"N.")
	return lowestPfail

#=============================================================================================================================================#
###############################################################################################################################################
#=============================================================================================================================================#


# ALGORITHM TRAINING: IF THE LOWEST VALUE ISN'T >= 1000, ITERATE THRU INPUT COMBINATIONS \
# AND RECORD THE INPUTS FOR WHICH LOWESTPFAIL IS >= 1000


def Optimize():	
	A_O, y_bar_O, I_O, y_c_O, y_t_O, Q_O, deltaB_O, h1, h2, L = 0,0,0,0,0,0,0,0,0,0

	#constant values, regardless of the values of h1,h2
	FindSupportReactionForces()
	SFD_left_P_RS_to_left_AL()
	SFD_left_AL_to_P_support()
	SFD_P_support_to_right_AL()
	SFD_right_AL_to_right_RS()
	BMD_left_AL()
	BMD_P_support()
	BMD_right_AL()
	FindMmax()

	for i in range(2,20,1):
		for j in range(50,180,1):
			h1,h2 = i,j
			(A, y_bar, I, y_c, y_t, Q, deltaB) = FindCrossSectionSpecs(h1,h2)
			FF = FlexuralFailure()
			SF = ShearFailure()
			BF = BucklingFailure()
			L = LowestP_fail()
			if L >= 1000:
				A_O, y_bar_O, I_O, y_c_O, y_t_O, Q_O, deltaB_O = FindCrossSectionSpecs(h1,h2)
				break
	return A_O, y_bar_O, I_O, y_c_O, y_t_O, Q_O, deltaB_O, h1, h2, L
	return True

# ALGORITHM FURTHER TRAINING: ITERATE THROUGH INPUT COMBINATIONS AND RECORD THE INPUTS \
# FOR WHICH LOWESTPFAIL IS MAXIMAL

def DeepOptimize(): #update
	A_DO, y_bar_DO, I_DO, h_total_DO, y_c_DO, y_t_DO, Q_DO, deltaB_DO = 0,0,0,0,0,0,0,0

	o = Optimize()
	if o == True: #if the basic optimization is successful, i.e. lowestPfail >= 1000
	
		A_DO, y_bar_DO, I_DO, h_total_DO, y_c_DO, y_t_DO, Q_DO, deltaB_DO, h1, h2, L = 0,0,0,0,0,0,0,0,0,0 #set h1,h2 to 0

		#constant values, regardless of the values of h1,h2
		FindSupportReactionForces()
		SFD_left_P_RS_to_left_AL()
		SFD_left_AL_to_P_support()
		SFD_P_support_to_right_AL()
		SFD_right_AL_to_right_RS()
		BMD_left_AL()
		BMD_P_support()
		BMD_right_AL()
		FindMmax()

		for i in range(2,20,1):
			for j in range(50,180,1):
				h1,h2 = i,j
				(A, y_bar, I, y_c, y_t, Q, deltaB) = FindCrossSectionSpecs(h1,h2)
				FF = FlexuralFailure()
				SF = ShearFailure()
				BF = BucklingFailure()
				L = LowestP_fail()
				if FindCrossSectionSpecs(h1-1,h2-1) < FindCrossSectionSpecs(h1,h2) and \
						FindCrossSectionSpecs(h1+1,h2+1) < FindCrossSectionSpecs(h1,h2):   #save local max values and update with iteration
						A_DO, y_bar_DO, I_DO, h_total_DO, y_c_DO, y_t_DO, Q_DO, deltaB_DO = FindCrossSectionSpecs(h1,h2)
						break
		return A_DO, y_bar_DO, I_DO, h_total_DO, y_c_DO, y_t_DO, Q_DO, deltaB_DO, h1, h2, L
		return True

#=============================================================================================================================================#
###############################################################################################################################################
#=============================================================================================================================================#


def StrengthWeightRatio(): #calculate the strengh-weight ratio based on the failure load and the weight of the bridge
	pass


if __name__ == "__main__":

	print("=============FindCrossSectionSpecs(h1,h2)=============")
	print(FindCrossSectionSpecs(h1,h2))
	print("	A 			y_bar 				I 					y_c 				y_t 			Q 				deltaB")

	print("\n=============THE FOLLOWING VALUES ARE CONSTANT FOR ANY h1,h2=============\n")
	print("=============FindSupportReactionForces()=============")
	print(FindSupportReactionForces())
	print("=============SFD_left_P_RS_to_left_AL()=============")
	print(SFD_left_P_RS_to_left_AL())
	print("=============SFD_left_AL_to_P_support()=============")
	print(SFD_left_AL_to_P_support())
	print("=============SFD_P_support_to_right_AL()=============")
	print(SFD_P_support_to_right_AL())
	print("=============SFD_right_AL_to_right_RS()=============")
	print(SFD_right_AL_to_right_RS())
	print("=============BMD_left_AL()=============")
	print(BMD_left_AL())
	print("=============BMD_P_support()=============")
	print(BMD_P_support())
	print("=============BMD_right_AL()=============")
	print(BMD_right_AL())
	print("=============FindMmax()=============")
	print(FindMmax())
	print("\n=============END OF THE CONSTANT VALUES=============\n")

	print("\n=============FlexuralFailure()=============")
	print(FlexuralFailure())
	print("=============ShearFailure()=============")
	print(ShearFailure())
	print("=============BucklingFailure()=============")
	print(BucklingFailure(),"\n")

	print("\n=============LowestP_fail()=============")
	print(LowestP_fail())

	# print("=============Optimize()=============")
	# print(Optimize())

	#print("=============DeepOptimize()=============")
	#print(DeepOptimize())

	#print("=============StrengthWeightRatio()=============")
	#print(StrengthWeightRatio())

	
# %%
