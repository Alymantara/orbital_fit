#Orbital parameters file

object='V630 Cas'    # Name of the object
file='ha45-6.txt'             # Name of the input file
psname='halfa45'              # Name of output plot file
errors          = False       # Does the input file has errors?

output 			= 'pdf'			# Can choose between: pdf, eps or png
data 			= False			# If True then exit data will put in file *.txt
plot 			= True			# Plot in Python window
plotlim	  		= 2			# Plot limits. 1 = close fit.
lim_res	  		= 1.2			# Plot limits for residuals. 1 = close fit.
#----------- Inital Values -----------
porb 			= 0.256387			# Period in days
hjd0 			= 58025.96		# HJD Zero point for fit
gama 			=  0.0			# Systemic Velocity in km/s
k1			   = 150.0			# Semi-amplitude velocity in km/s
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
initial 		= 0.2550		# Initial value
delta 			= 0.000005			# Step value for loop
num 			= 1000			# Number of loops

### 2.56387
