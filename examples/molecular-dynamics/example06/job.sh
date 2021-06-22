#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`


# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "Forces for the Cr triangular trimer to a non-collinear magnetic calculation."

# set the needed environment variables
. ../../environment_variables

rm -f out*
rm -f scf/*
mkdir scf
rm -f tempo tempo2
rm -f Etot_vs_force*

cat > Etot_vs_forces.dat << EOF
@# a f Etot
EOF

for dx in  -0.04 -0.035 -0.03 -0.025 -0.02 -0.01 0 0.01 0.02 0.025 0.03 0.035 0.04  ; do

  echo "dx= ${dx}"




cat > in_master.txt<<EOF
&units
 energy = 'eV'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&calculation
 processing = 'scf'
 post_processing='forces'
 post_processing_dir = 'scf'
 /
&element
 ne = 1
 symbol(1) = 'Cr'
 q(1) = 6.0
 q_d(1) = 5.0
 u_lcn(1) = 20.0000000000000000
 i_stoner_d(1) = 0.82
 xi_so_d(1) = 0.0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/cr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = 2.
 v(1,:) = 1.000000000000000, 0.0000000000000000, 0.0000000000000000
 v(2,:) = 0.0000000000000000, 1.000000000000000, 0.0000000000000000
 v(3,:) = 0.0000000000000000, 0.0000000000000000, 1.0000000000000000
 /
&atom
 ns = 4
 na = 3
 ntag = 1
 stag(1) = 3
 tag(1) = 'Cr_trimer'
 pbc = 0, 0, 0
 r(1,:) = $dx, 0.0000000000000000, 0.0
 r(2,:) = 0.86602540378443864676,-0.5, 0.0
 r(3,:) = 0.86602540378443864676, 0.5, 0.0
 m_listing = 'by_atom'
 m_coord = 'spherical'
 m(1,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(2,:) = 3.1800000000000000, 90.000000000000000, 90.0000000000000000
 m(3,:) = 3.1800000000000000, 90.0000000000000000, 180.0000000000000000
 lambda_pen_listing= 'by_tag'
 lambda_pen(1) = 0.0
 m_coord = 'spherical'
 /
&mesh
 type = 'mp'
 gx = 1, 1, 1
 dx = 0, 0, 0
 /
&hamiltonian_tb
 m_penalization = 'r,theta,phi'
 /
&energy
 smearing = 'mv'
 degauss = 0.2
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.00001
 delta_q  = 0.00001
 verbose = .true.
 ni_max = 1500
 /
&forces
 computed=.true.
 /
EOF

# Set DyNaMol root directory in in_master.txt
#sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
#mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x 

cat > alat << EOF
a= $(echo "$dx*2.00" | bc -l)
EOF

grep 'Atom 1' ./scf/out_log.txt | awk '{print $4}' | head -1 > fx

grep -e 'a=' alat | awk '{print $2}' > a

grep -e 'en =' out_energy.txt | awk '{print $3}' > en

paste a fx en >> Etot_vs_forces.dat

cp ./scf/out_log.txt ./scf/out_log$(echo "$dx*2.03" | bc -l).txt
done

rm -f a alat en fx
