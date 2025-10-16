#!/bin/bash

# An equilibration workflow.
# Designed for "standard" protein or protein/ligand systems
# in explicit solvent. Not optimised for membrane protein systems.

#### Edit the variables below ####
prmtop_file="abl_ligand.prmtop"

inpcrd_file="abl_ligand.inpcrd"

solute=265 #number of residues in solute
T="310.0"

#### Stop editing ####
PMEMD="pmemd.cuda"

# 1K step energy minimization with strong restraints on heavy atoms, no shake
cat > step1.in <<EOF
Min explicit solvent heavy atom rest no shake
&cntrl
  imin = 1, ncyc = 50, maxcyc = 1000,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 5.0,
&end
EOF

# NTV MD with strong restraints on heavy atoms, shake, dt=.001, 15 ps
cat > step2.in <<EOF
MD explicit solvent heavy atom rest shake dt 0.001
&cntrl
  imin = 0, nstlim = 30000, dt=0.001,
  ntpr = 50, ntwr = 500,
  iwrap = 1, 
  ntc = 2, ntf = 2,
  ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
  ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 5.0,
&end
EOF

# Energy minimization with relaxed restraints on heavy atoms, no shake
cat > step3.in <<EOF
Min explicit solvent relaxed heavy atom rest no shake
&cntrl
  imin = 1, ncyc = 50, maxcyc = 1000,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 2.0,
&end
EOF

# Energy minimization with minimal restraints on heavy atoms, no shake
cat > step4.in <<EOF
Min explicit solvent minimal heavy atom rest no shake
&cntrl
  imin = 1, ncyc = 50, maxcyc = 1000,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 0.1,
&end
EOF

# Energy minimization with no restraints, no shake
cat > step5.in <<EOF
Min explicit solvent no heavy atom res no shake
&cntrl
  imin = 1, ncyc = 50, maxcyc = 1000,
  ntwx = 1000, ntpr = 50, ntwr = 500,
&end
EOF

# NTP MD with shake and low restraints on heavy atoms, 20 ps dt=.002
cat > step6.in <<EOF
MD explicit solvent heavy atom low rest shake dt 0.002
&cntrl
  imin = 0, nstlim = 10000, dt = 0.002,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  iwrap = 1, 
  ntc = 2, ntf = 2, ntb = 2,
  ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
  ntp = 1, taup = 1.0,
  ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 1.0,
&end
EOF

# NTP MD with shake and minimal restraints on heavy atoms
cat > step7.in <<EOF
MD explicit solvent heavy atom minimal rest shake 20 ps, dt=.002
&cntrl
  imin = 0, nstlim = 10000, dt=0.002,
  ntx = 5, irest = 1,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  iwrap = 1, 
  ntc = 2, ntf = 2, ntb = 2,
  ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
  ntp = 1, taup = 1.0,
  ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 0.5,
&end
EOF

# NTP MD with shake and minimal restraints on backbone atoms, dt=0.002, 20 ps
cat > step8.in <<EOF
MD explicit solvent heavy atom minimal BB rest shake dt 0.002
&cntrl
  imin = 0, nstlim = 10000, dt=0.002,
  ntx = 5, irest = 1,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  iwrap = 1, 
  ntc = 2, ntf = 2, ntb = 2,
  ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
  ntp = 1, taup = 1.0,
  ntr = 1, restraintmask = ":1-$solute@H,N,CA,HA,C,O", restraint_wt = 0.5,
&end
EOF

# NTP MD with shake and no restraints, dt=0.002, 200 ps
cat > step9.in <<EOF
MD explicit solvent heavy atom no rest shake dt 0.002
&cntrl
  imin = 0, nstlim = 10000, dt=0.002,
  ntx = 5, irest = 1,
  ntwx = 1000, ntpr = 50, ntwr = 500,
  iwrap = 1, 
  ntc = 2, ntf = 2, ntb = 2,
  ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
  ntp = 1, taup = 1.0,
&end
EOF

START="`date +%s.%N`"

# Minimization Phase - reference coords are updated each run
for RUN in step1 step2 step3 step4 step5 ; do
 echo "------------------------"
 echo "Minimization phase: $RUN"
 echo "------------------------"
 if [[ ! -f $RUN.rst7 ]]; then
     echo "File -- $RUN.rst7 -- does not exists. Running job..."
     $PMEMD -O -i $RUN.in -p $prmtop_file -c $inpcrd_file -ref $inpcrd_file -o $RUN.out -x $RUN.nc -r $RUN.rst7 -inf $RUN.mdinfpo
 else
     echo "File -- $RUN.rst7 -- exists.  Checking the next step."
 fi
 echo ""
 inpcrd_file="$RUN.rst7"
done

# Equilibration phase - reference coords are last coords from minimize phase
REF=$inpcrd_file
for RUN in step6 step7 step8 step9 ; do

 echo "------------------------"
 echo "Equilibration phase: $RUN"
 echo "------------------------"
 if [[ ! -f $RUN.rst7 ]]; then
     echo "File -- $RUN.rst7 -- does not exists. Running job..."
     $PMEMD -O -i $RUN.in -p $prmtop_file -c $inpcrd_file -ref $REF -o $RUN.out -x $RUN.nc -r $RUN.rst7 -inf $RUN.mdinfo
 else
     echo "File -- $RUN.rst7 -- exists.  Checking the next step."
 fi
  echo ""
  inpcrd_file="$RUN.rst7"

done

# Reset the time in the restrat file to zero:
sed -i 's/0.3000000E+02/0.0000000E+00/g' step9.rst7

STOP="`date +%s.%N`"
TIMING=`echo "scale=4; $STOP - $START;" | bc`
echo "$TIMING seconds."
echo ""

exit 0
