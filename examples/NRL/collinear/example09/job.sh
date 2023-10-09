#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to evaluate the influence of number of kpoints"


# set the needed environment variables
. ../../../environment_variables

rm -f tempo tempo2
rm -f out*

rm -rf results
mkdir results

a=2.87

cat > results/Etot_vs_nK.dat << EOF
@#  nK  Etot
EOF

for nK in 5 10 15 20 25 ; do

$ECHO "nK=" $nK"x"$nK"x"$nK

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
 symbol(1) = 'Fe'
 q(1)   = 8.0
 q_d(1) = 7.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.95
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Fe_bulk'
 r(1,:) = 0.0, 0.0, 0.0
 m(1,:) = 2, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = $nK, $nK, $nK
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
 delta_en = 0.0001
 delta_q  = 0.0001
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cat > tempo << EOF
nK= $nK
EOF

cat tempo out_energy.txt>>tempo2

done
grep -e 'nK=' -e 'en =' tempo2 | awk '/nK/{nK = $(NF)}/en/{print nK, $(NF)}' > results/Etot_vs_nK.dat

rm -f tempo tempo2

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to evaluate the influence of the smearing method"

for smearing in fd g mp mv; do
  echo "smearing= $smearing"
  
cat > results/Etot_vs_degauss.$smearing.dat << EOF
@#  degauss(eV)  Etot(eV)
EOF

for degauss in 0.005 0.01 0.05 0.1 0.15  0.2 0.3 ; do

 echo "degauss= $degauss"

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
 symbol(1) = 'Fe'
 q(1)   = 8.0
 q_d(1) = 7.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.95
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Fe_bulk'
 r(1,:) = 0.0, 0.0, 0.0
 m(1,:) = 2, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 10, 10, 10
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = '$smearing'
 degauss = $degauss
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cat > tempo << EOF
degauss= $degauss
EOF

cat tempo out_energy.txt>>tempo2


done
grep -e 'degauss=' -e 'en =' tempo2 | awk '/degauss/{degauss = $(NF)}/en/{print degauss, $(NF)}' > results/Etot_vs_degauss.$smearing.dat

rm -f tempo tempo2

done



