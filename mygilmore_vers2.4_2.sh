#!/bin/bash

scriptpath=$(dirname $(realpath $0))
thisdir=$(pwd)


   Rstart=29.5e-6 #5.45126e-6
   RdotStart=1470. #m/s
   Rn=73e-6     #m
   deltaT=1e-11  #s
   Tstart=0.    #s
   Tend=115e-6  #s
   patm=101315  #Pa
   pac=0.0       #Pa
   f=26120.0        # Hz
   bVan=0.0000364 # in m³/mol: air according to http://de.wikipedia.org/wiki/Van-der-Waals-Gleichung
   #mu=0.001     # Pa s (I think)
   #mu=0.00089008220776922 # Pa s # this is correct for 25°C. taken from: http://www.peacesoftware.de/einigewerte/wasser_dampf_e.html
   mu=1.002e-3  # water: 1.002e-3  # was: 0.0000186 # Pa s. correct value at 20°C:1.002mPa*s taken from: http://en.wikipedia.org/wiki/Viscosity
   sigma=0.0725  #water: 0.0725 # Pa m (I think)
   BTait=3046e5 # Pa
   nTait=7.15  #water: 7.15   # no unit
   
   # kappa: fractions don't work!!! Only floating point!
   # down below for gnuplot: put kappa in manually!
   #kappa=1.3333333333333333333333333333333333333333333333333333333333333333333333333333333333
   kappa=1.4
   pv=0 #2337.     # Pa
   #rho0=998.    # density of liquid in kg/m³
   rho0=998.20608789369  # water: 998.20608789369 #<-- is correct for 20°C. This one is correct for 25°C: 997.04802746938
   #rho0v=0.46638524217467 # density of gas in kg/m³
   #rho0v=1.2042092351 #0.08737912 # // kg/m^3 of vapour or gas at starting radius!!!!! Has to be adapted for each Rn!!!!!!
   SpecGasConst=287. #287 for dry air, 462 for H2O vapour. Both in J/molK
   TempRef=20. # in deg Celsius.
   T0=1e-2
   epsilon=1e-2
   phase=-107.0 #in deg

if [[ ! -e ${scriptpath}/mygilmore_vers2.4_2 ]]
then 
   echo "g++ -O3 mygilmore_vers2.4_2.c -o mygilmore_vers2.4_2"
   g++ -O3 ${scriptpath}/mygilmore_vers2.4_2.c -o ${scriptpath}/mygilmore_vers2.4_2
fi

${scriptpath}/mygilmore_vers2.4_2 $Rstart $RdotStart $Rn $deltaT $Tstart $Tend $patm $pac $f $bVan $mu $sigma $BTait $nTait $kappa $pv $rho0 $SpecGasConst $TempRef $T0 $epsilon $phase

if [ $thisdir != $scriptpath ]
then
 mv gilmore2.dat $scriptpath/
fi

# echo "reset
# set term postscript eps color enhanced solid
# set encoding utf8
# set output \"gilmore.eps\"
# set grid
# set key above
# set xlabel \"t [{/Symbol m}s]\"
# set ylabel \"R [{/Symbol m}m]\"
# set y2label \"p_{ac}+p_{stat} [bar]\"
# set y2tics
# set title \"Comparing PA75\% experiments with gilmore model at f=${f}\n\\
# r_{start}=${Rstart}m, {/Symbol g}=${kappa}, R_{n}=${Rn}m, p_{stat}=${patm}Pa, p_{ac}=${pac}Pa, {/Symbol f}=${phase}^o\"
# 
# p \"${fitfile}\" u ((\$1)):((\$2)) w l t \"${name}\",\\
#   \"gilmore2.dat\"    u ((\$1)*1e6):((\$2)*1e6) w l t \"Gilmore fit\",\\
#   \"\"                u ((\$1)*1e6):((\$6)/1e5) w l axes x1y2 t \"ac. pressure\",\\
#   1                axes x1y2 t \"\"" > gilmore_2.gnuplot
# 
# gnuplot gilmore_2.gnuplot
# # 
# epstopdf ${fitfile}.eps
# # 
# rm ${fitfile}.eps
