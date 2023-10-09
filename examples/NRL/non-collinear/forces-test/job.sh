#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`
TBPARAM_DIR=/home/rcardias/TBKOSTER/tb_parameters
BIN_DIR=/home/rcardias/TBKOSTER/linux/bin


# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "Damped Spin Dynamics of a Cr dimer."

# set the needed environment variables
#. ../../../environment_variables

rm -f tempo tempo2
rm -f Etot_vs_force*

cat > Etot_vs_forces.dat << EOF
@# a f Etot
EOF

for a in 2.12 2.14 2.16 2.18 2.20 2.22 2.24 2.26 2.28 2.30 2.32 2.34 2.36 2.38; do

  echo "a= ${a}"


cat > in_master.txt<<EOF
&units
 energy = 'eV'
 length = 'ang'
 time = 'fs'
 /
&calculation
 processing = 'scf'
 post_processing = 'forces'
 post_processing_dir = 'scf'
 /
&element
 ne = 1
 symbol(1) = 'Fe'
 q(1) = 8.0
 q_d(1) = 7.0
 u_lcn(1) = 20.0000000000000000
 i_stoner_d(1) = 0.95
 xi_so_d(1) = 0.0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = ${a}
 v(1,:) = 1.000000000000000, 0.0000000000000000, 0.0000000000000000
 v(2,:) = 0.0000000000000000, 1.000000000000000, 0.0000000000000000
 v(3,:) = 0.0000000000000000, 0.0000000000000000, 1.0000000000000000
 /
&atom
 ns = 4
 na = 2
 ntag = 1
 stag(1) = 2
 tag(1) = 'Fe_dimer'
 pbc = 0, 0, 0
 r(1,:) = 0.0000000000000000, 0.0000000000000000, 0.0
 r(2,:) = 0.0000000000000000, 0.0000000000000000, 1.0
 r_coord = 'direct'
 m(1,:) = 3.0, 0.0, 0.0
 m(2,:) = 3.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 1, 1, 1
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.05
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.00001
 delta_q  = 0.00001
 verbose = .false.
 ni_max = 1000
 /
&forces
 computed=.true.
 /
&sd
 integrator = 'st_1'
 t_i = 0
 t_f = 100
 dt = 0.01
 fixed_time_step = .true.
 alpha = 0.0
 temp = 0.0
 verbose = .true.
 /
EOF

# Set TBKOSTER root directory in in_master.txt
#sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
#mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cat > alat << EOF
a= $a
EOF

grep Fz ./scf/out_log.txt | awk '{print $8}' | head -1 > fz

grep -e 'a=' alat | awk '{print $2}' > a

grep -e 'en =' out_energy.txt | awk '{print $3}' > en

paste a fz en >> Etot_vs_forces.dat

cp ./scf/out_log.txt ./scf/out_log${a}.txt
done

