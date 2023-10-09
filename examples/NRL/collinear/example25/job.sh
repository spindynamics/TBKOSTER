#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the total energy of zirconium fcc, bcc and hcp(c/a)"

# set the needed environment variables
. ../../environment_variables

rm -f out*
mkdir scf
rm -f results/Etot_vs_a_hcp.dat


#covera=$(echo "sqrt(8/3)" |bc -l)


for covera in  1.58 1.59 1.60 1.61 1.62 1.63 1.64 1.65 1.66  ; do

rm -f tempo tempo2

cat >> results/Etot_vs_a_hcp.dat << EOF
@# a  Etot(eV)   covera=$covera
EOF

$ECHO "hcp"
$ECHO "c/a= $covera"

for a in 3.00 3.05 3.10 3.11 3.12 3.13 3.14 3.15 3.16 3.17 3.18 3.19 3.20 3.25 3.30 3.35 ; do

$ECHO "a= $a"


cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 /
&units
 length='ang'
 energy='ev'
 time='fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Zr'
 q(1)  = 4.0
 q_d(1) = 3.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = .0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/zr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = $(echo "sqrt(3)/2" |bc -l) 0.5 0.0
 v(2,:) = $(echo "-sqrt(3)/2" |bc -l) 0.5 0.0
 v(3,:) = 0.0 0.0 $covera 
 /
&atom
 ns = 2
 na = 2
 ntag = 1
 tag(1) = 'Zr' 
 stag(1)=2
 pbc = 5, 5, 5
 r_coord='direct'
 r(1,:) =      $(echo "1.0/3.0" |bc -l)     $(echo "2.0/3.0" |bc -l)      0.25
 r(2,:) =      $(echo "2.0/3.0" |bc -l)     $(echo "1.0/3.0" |bc -l)      0.75
 /
&mesh
 type = 'mp'
 gx = 15, 15, 15
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
 delta_en=0.0001
 delta_q=0.0001
 verbose=.true.
 ni_max=200
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

cat tempo out_energy.txt>>tempo2

done

grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> results/Etot_vs_a_hcp.dat

rm -f tempo tempo2 
done


rm -f out*
mkdir scf
rm -f results/Etot_vs_a_fcc.dat

cat > results/Etot_vs_a_fcc.dat << EOF
@# a  Etot(eV)
EOF


for a in 3.00 3.05 3.10 3.11 3.12 3.13 3.14 3.15 3.16 3.17 3.18 3.19 3.20 3.25 3.30 3.35 ; do

$ECHO "fcc"

$ECHO "a= $a"

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 /
&units
 length='ang'
 energy='ev'
 time='fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Zr'
 q(1)  = 4.0
 q_d(1) = 3.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = .0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/zr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $(echo "$a*sqrt(2)" |bc -l)
 v(1,:) = 0.0 0.5 0.5
 v(2,:) = 0.5 0.0 0.5
 v(3,:) = 0.5 0.5 0.0
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 tag(1) = 'Zr' 
 stag(1)=1
 pbc = 5, 5, 5
 r_coord='direct'
 r(1,:) =      0.0,  0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 15, 15, 15
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
 delta_en=0.0001
 delta_q=0.0001
 verbose=.true.
 ni_max=200
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

cat tempo out_energy.txt>>tempo2

done

grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> results/Etot_vs_a_fcc.dat

rm -f tempo tempo2 


rm -f out*
mkdir scf
rm -f results/Etot_vs_a_bcc.dat

cat > results/Etot_vs_a_bcc.dat << EOF
@# a  Etot(eV)
EOF


for a in 3.00 3.05 3.10 3.11 3.12 3.13 3.14 3.15 3.16 3.17 3.18 3.19 3.20 3.25 3.30 3.35 ; do

$ECHO "bcc"

$ECHO "a= $a"

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 /
&units
 length='ang'
 energy='ev'
 time='fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Zr'
 q(1)  = 4.0
 q_d(1) = 3.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = .0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/zr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $(echo "$a*2/sqrt(3)" |bc -l)
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 tag(1) = 'Zr' 
 stag(1)=1
 pbc = 5, 5, 5
 r_coord='direct'
 r(1,:) =      0.0,  0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 15, 15, 15
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
 delta_en=0.0001
 delta_q=0.0001
 verbose=.true.
 ni_max=200
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

cat tempo out_energy.txt>>tempo2

done

grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> results/Etot_vs_a_bcc.dat

rm -f tempo tempo2 
