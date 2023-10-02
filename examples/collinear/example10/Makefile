.ONESHELL:

TBKOSTER_ROOT_DIRECTORY=/home/barreto/github.com/TBKOSTER
TBKOSTER_COMMAND=$TBKOSTER_ROOT_DIRECTORY/bin/TBKOSTER.x
define generate_file =
cat > in_master.orig << EOF
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 /
&calculation
 processing = 'scf'
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
 filename(1) = 'TBKOSTER_ROOT_DIRECTORY/tb_parameters/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = 2.27
 v(1,:) = 10.0,  0.0,  0.0
 v(2,:) =  0.0, 10.0,  0.0
 v(3,:) =  0.0,  0.0,  5.0
 /
&atom
 ns = 2
 na = 5
 ntag = 1
 stag(1) = 5
 tag(1) = 'Fe_wire'
 pbc = 0, 0, 1
 r(1,:) = 0.0, 0.0, 0.0
 r(2,:) = 0.0, 0.0, 0.2
 r(3,:) = 0.0, 0.0, 0.4
 r(4,:) = 0.0, 0.0, 0.6
 r(5,:) = 0.0, 0.0, 0.8
 m_listing = 'by_tag'
 m_coord = 'spherical'
 m(1,:) = 3.0, 0.0, 0.0
 lambda_pen_listing = 'by_tag'
 lambda_pen(1) = 30.0
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
 degauss = 0.2
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.00001
 delta_q  = 0.00001
 verbose = .true.
 /
EOF
endef

run: help in_master
	@sed "s|TBKOSTER_ROOT_DIRECTORY|$(TBKOSTER_ROOT_DIRECTORY)|g" in_master.orig > in_master.txt
	@$(IN_TBKOSTER_COMMAND) >job_out.txt

verify:
	@$(info $$TBKOSTER_ROOT_DIRECTORY is [$(TBKOSTER_ROOT_DIRECTORY)])
	@$(info $$TBKOSTER_COMMAND is [$(TBKOSTER_COMMAND)])

in_master:
ifeq ($(TBKOSTER_ROOT_DIRECTORY),)
@$(error TBKOSTER_ROOT_DIRECTORY is not set)
endif
ifeq ($(TBKOSTER_COMMAND),)
@$(error TBKOSTER_COMMAND is not set)
endif
	@$(value generate_file)

clean:
	@rm -f out* job_out.txt in_master*

help:
	@$(info example 10:)
	@$(info SCF collinear spin calculation of a 5-atom Fe wire with magnetic penalization on atom 1)
