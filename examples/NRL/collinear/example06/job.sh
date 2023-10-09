#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the total energy vs magnetization of Febcc"
$ECHO "penalization technique"

# set the needed environment variables
. ../../environment_variables
rm -f tempo tempo2
rm -f out*

cat > Etot_vs_mag.dat << EOF
@# mag  Etot
EOF

a=2.87
for mag in 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.3 2.35 2.4 2.45 2.5 2.6 2.8; do

  echo "mag=' $mag"
cat > Etot_vs_mag << EOF
@#  mag(muB)  Etot(eV)
EOF

cat > in_master.txt<<EOF
&calculation
 processing='scf'
 post_processing='txt2xyz'
 post_processing_dir='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Fe'
 q(1)   = 8.0
 q_d(1) = 7.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.95
 xi_so_d(1) = 0.06
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) =  1.0  0.0,  0.0
 v(2,:) =  0.0, 1.0,  0.0
 v(3,:) =  0.0, 0.0,  1.0
 /
&atom
 ns = 2
 na = 2
 ntag = 1
 stag(1) = 2
 tag(1) = 'Fe_bulk'
 r(1,:) = 0.0, 0.0, 0.0
 r(2,:) = 0.5, 0.5, 0.5 
 m(1,:) = $mag, 0.0, 0.0
  lambda_pen_listing = 'by_atom'
  lambda_pen(1) = 10.0
  lambda_pen(2)=0.0
 /
&mesh
 type = 'mp'
 gx = 10, 10, 10
 dx = 0, 0, 0
 /
&hamiltonian_tb
m_penalization = 'r'
/
&energy
 smearing = 'mv'
 degauss = 0.1
 /
&mixing
 alpha = 0.05
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 ni_max=200
 verbose=.true.
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cat > tempo << EOF
mag= $mag
EOF

cat tempo out_energy.txt>>tempo2

cat out_log.txt >> out_log_all.txt


done
grep -e 'mag=' -e 'en =' tempo2 | awk '/mag/{mag = $(NF)}/en/{print mag, $(NF)}' >> Etot_vs_mag.dat

rm -f tempo tempo2
