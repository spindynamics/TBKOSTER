#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate band structure and DOS with SOC for Febcc"

# set the needed environment variables
. ../../environment_variables

rm -f out*
rm -rf band
mkdir band

a=2.87



cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 post_processing = 'band'
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
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 4
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Fe_bulk'
 r(1,:) = 0.0, 0.0, 0.0
 m(1,:) = 1.0, 0.0, 0.0
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
 /
EOF

cat > band/in_mesh.txt<<EOF
&mesh
 type = 'path'
 gxs = 20
 nxs = 6
 xs_label(1) = 'G'
 xs_label(2) = 'H'
 xs_label(3) = 'N'
 xs_label(4) = 'G'
 xs_label(5) = 'P'
 xs_label(6) = 'H'
 xs(1,:) =  0  ,  0  ,  0
 xs(2,:) =  0.5, -0.5,  0.5
 xs(3,:) =  0  ,  0  ,  0.5
 xs(4,:) =  0  ,  0  ,  0
 xs(5,:) =  0.25, 0.25, 0.25
 xs(6,:) =  0.5, -0.5,  0.5
 /
EOF

cat > band/in_energy.txt<<EOF
&energy
 smearing = 'mv'
 degauss = 0.2
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > band/in_dos.txt<<EOF
&dos
 nen=100
 na_dos=1
 ia= 1
 en_min=-10
 en_max=10
 /
EOF



# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 


rm -rf scf
mkdir scf
rm -rf dos
mkdir dos

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 post_processing = 'dos'
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
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 4
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Fe_bulk'
 r(1,:) = 0.0, 0.0, 0.0
 m(1,:) = 1.0, 0.0, 0.0
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
 /
EOF

cat > dos/in_energy.txt<<EOF
&energy
 smearing = 'g'
 degauss = 0.1
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > dos/in_dos.txt<<EOF
&dos
 nen=200
 na_dos=1
 ia= 1
 en_min=-10
 en_max=10
 /
EOF

cat > dos/in_mesh.txt<<EOF
&mesh
 type = 'mp'
 gx = 10, 10, 10
 dx = 0, 0, 0
 /
EOF




# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 







