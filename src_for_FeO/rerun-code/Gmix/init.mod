# USER Parameter 

variable mass1 equal 55.845               # mass
variable up equal 1.0e-2                  # strain  

# metal units, elastic constants in GPa
units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

# Define MD parameters
variable timestep equal 0.001             # timestep
variable tdamp equal 0.01                 # time constant for thermostat
variable nequil equal 1

read_data ${atomdir}/read.in
