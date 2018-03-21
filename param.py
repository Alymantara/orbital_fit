#Orbital parameters file for orbital.py

object 			= 'PNV Aql-Halfa ajustadas por Guassianas 750 kms'		# Name of the object
file 			= 'gauss_750_kms.txt'			# Name of the input file
psname 			= 'PNV_gauss_750kms'	# Name of output plot file
output 			= 'pdf'			# Can choose between: pdf, eps or png
data 			= False			# If True then exit data will put in file *.txt
plot 			= True			# Plot in Python window
plotlim	  		= 3			# Plot limits. 1 = close fit.
lim_res	  		= 1.2			# Plot limits for residuals. 1 = close fit.
#----------- Inital Values -----------
porb 			= 0.05698			# Period in days
hjd0 			= 57566.6		# HJD Zero point for fit
gama 			= -21.0			# Systemic Velocity in km/s
k1			   = -50.0			# Semi-amplitude velocity in km/s
#----------- Deal with errors -----------
sigma		= 5		# Default sigma if errors are unknown. If errors are known, use None
scale_errors	= True			# Scales errors so chi^2 = 1
do_bootstrap 	= False			# To perform bootstrap
boot_iter   	= 10000			# Number of bootstrap copies
#----------- Fix any parameter -----------
fix_porb		= False	        # Change to False to fix parameter
fix_hjd0		= True			# Change to False to fix parameter
fix_gama		= True			# Change to False to fix parameter
fix_k1      	= True			# Change to False to fix parameter
#----------- For Looping Program -----------
variable 		= 'porb'		# Variable to loop: porb, hjd0, gama, k1
initial 		= 0.056		# Initial value
delta 			= 0.000005			# Step value for loop
num 			= 1000			# Number of loops
errors          = True
