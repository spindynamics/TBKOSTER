#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to the total energy and magnetization versus a for Cofcc"

# set the needed environment variables
. ../../environment_variables

rm -f tempo tempo2 tempo3 tempo4
rm -f out*
mkdir scf

cat > Etot_vs_a.dat << EOF
@#  a(A)  Etot(eV)
EOF

cat > M_vs_a.dat << EOF
@#  a(A)  Mtot(muB)
EOF

for a in 3.0 3.05 3.15 3.20 3.30 3.40 3.50 3.60 3.70 3.80 4.0; do

$ECHO "a= $a"



cat > in_master.txt<<EOF
&calculation
 processing='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Co'
 q(1)   = 9.0
 q_d(1) = 8.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 1.10
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/co_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 0.0 0.5 0.5
 v(2,:) = 0.5 0.0 0.5
 v(3,:) = 0.5 0.5 0.0
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Co'
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 10, 10, 10
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.2
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .true.
 ni_max = 100
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cat > tempo << EOF
a= $a
EOF

grep m_r_tot out_log.txt | tail -1 >tempo3

cat tempo out_energy.txt>>tempo2
cat tempo tempo3>>tempo4



done
grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> Etot_vs_a.dat
grep -e 'a=' -e 'm_r_tot =' tempo4 | awk '/a/{a = $(NF)}/m_r_tot/{print a, $(NF-1)}' >> M_vs_a.dat

rm -f tempo tempo2 tempo3 tempo4
