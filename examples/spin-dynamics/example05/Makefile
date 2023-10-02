SHELL:=bash
.ONESHELL:
.NOTPARALLEL:

TBKOSTER_ROOT_DIR=/home/barreto/github/TBKOSTER
EXE = TBKOSTER.x
TBKOSTER = ${TBKOSTER_ROOT_DIR}/linux/bin/${EXE}
TBPARAM_DIR=${TBKOSTER_ROOT_DIR}/tb_parameters


define generate_file =
 cat > in_master.orig << EOF
&units
 energy = 'eV'
 length = 'ang'
 time = 'fs'
 /
&calculation
 pre_processing = 'none'
 pre_processing_dir = '.'
 processing = 'scf'
 post_processing = 'txt2xyz'
 post_processing_dir = '.'
 /
&element
 ne = 1
 symbol(1) = 'Fe'
 q(1) = 8.0
 q_d(1) = 7.0
 u_lcn(1) = 20.0000000000000000
 i_stoner_d(1) = 0.95
 xi_so_d(1) = 5.9999999999999998E-002
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = 2.27
 v(1,:) = 10.000000000000000, 0.0000000000000000, 0.0000000000000000
 v(2,:) = 0.0000000000000000, 10.000000000000000, 0.0000000000000000
 v(3,:) = 0.0000000000000000, 0.0000000000000000, 5.0000000000000000
 /
&atom
 ns = 4
 na = 5
 ntag = 1
 stag(1) = 5
 tag(1) = 'Fe_wire'
 pbc = 0, 0, 1
 r(1,:) = 0.0000000000000000, 0.0000000000000000, 0.0
 r(2,:) = 0.0000000000000000, 0.0000000000000000, 0.2
 r(3,:) = 0.0000000000000000, 0.0000000000000000, 0.4
 r(4,:) = 0.0000000000000000, 0.0000000000000000, 0.6
 r(5,:) = 0.0000000000000000, 0.0000000000000000, 0.8
 m_listing = 'by_atom'
 m(1,:) = 3.1800000000000000, 1.0000000000000000, 0.0000000000000000
 m(2,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(3,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(4,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(5,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 lambda_pen_listing= 'by_tag'
 lambda_pen(1) = 30.0
 m_coord = 'spherical'
 /
&mesh
 type = 'mp'
 gx = 1, 1, 100
 dx = 0, 0, 0
 /
&hamiltonian_tb
 m_penalization = 'theta,phi'
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
 /
&sd
 integrator = 'euler'
 t_i = 0
 t_f = 30
 dt = 0.05
 fixed_time_step = .true.
 alpha = 1.0
 temp = 0.0
 verbose = .true.
 /
EOF
endef

define content =
	for i in 2.05 2.10 2.15 2.20 2.25 2.30; do
	rm -rf ${i}	
	mkdir ${i}
        sed "s|processing = 'scf'|processing = 'sd'|g" in_master.orig > temp
        sed "s|v_factor = 2.27|v_factor = ${i}|g" temp > ${i}/in_master_${i}.txt
        rm temp
        mv ${i}/in_master_${i}.txt ${i}/in_master.txt
        cd ${i}
	/home/barreto/github/TBKOSTER/linux/bin/TBKOSTER.x >job_out_${i}.txt
	cd ..
	done
endef


in_master: ; @$(value generate_file)

.ONESHELL:

run: in_master
	@echo 'SCF collinear spin calculation of a 5-atom Fe wire with magnetic penalization on atom 1'
	@echo ${TBPARAM_DIR}
	@cp in_master.orig in_master.txt
	@${TBKOSTER} >job_out.txt

loop: in_master
	@echo "Perform a spindynamics on a loop when the space between the atoms is varying"
	@rm -f in_master.txt
	@echo ${TBPARAM_DIR}	
	@$(value content)

clean:
	@rm -fr out* sd job_out.txt in_master* 2.*
