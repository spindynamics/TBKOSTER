SHELL:=bash
.ONESHELL:
.NOTPARALLEL:

TBKOSTER_ROOT_DIRECTORY=/home/barreto/github/TBKOSTER/
TBKOSTER_COMMAND=$(TBKOSTER_ROOT_DIRECTORY)/linux/bin/TBKOSTER.x
TBPARAM_DIR=$(TBKOSTER_ROOT_DIRECTORY)/tb_parameters
EXAMPLE_DIR=$(TBKOSTER_ROOT_DIRECTORY)/linux/exampleXX

define get_x =
$(shell grep 'r(1,:) =' $(EXAMPLE_DIR)/results/force_log_$(1).txt | cut -d"=" -f2 | cut -d"," -f1)
endef

define diag_x =
$(shell echo "$(shell grep 'r(1,:) =' $(EXAMPLE_DIR)/results/force_log_$(1).txt | cut -d"=" -f2 | cut -d"," -f1)*2.03*sqrt(3)"| bc -l)
endef

define get_e =
$(shell grep 'en =' $(EXAMPLE_DIR)/results/out_energy_$(1).txt | cut -d"=" -f2)
endef

define get_fx =
$(shell grep 'Atom 1' $(EXAMPLE_DIR)/results/force_log_$(1).txt |cut -d"=" -f2 |cut -d"F" -f1)
endef

define get_variable =
@$(info $(1) $(call get_x,$(1)) $(call diag_x,$(1)) $(call get_e,$(1)) $(call get_fx,$(1)))
endef

define corpse =
$(call generate_file,$(TBPARAM_DIR),$(1),1)
@mv in_master.orig in_master.txt	
@$(TBKOSTER_COMMAND) >/dev/null
@mv scf/out_log.txt results/force_log_$(1).txt
@mv out_energy.txt results/out_energy_$(1).txt
endef

define generate_file =
cat > in_master.orig << EOF
&calculation
processing = 'scf'
post_processing= 'forces'
post_processing_dir= 'scf'
/
&units
length='ang'
energy='ev'
time='fs'
mass='hau'
/
&element
ne=1
symbol(1)='Co'
q(1)=9.0
q_d(1)=8.0
u_lcn(1)=20
i_stoner_d(1) = 1.10
xi_so_d(1) = 0.0
/
&element_tb
filename(1) = '$(1)/Co_no_onsite.par'
/
&lattice
v_factor = 2.03
v(1,:) = 1, 0, 0
v(2,:) = 0, 1, 0
v(3,:) = 0, 0, 1
/
&atom
ns = 1
na = 3
ntag = 1
tag(1)='Co'
stag(1)= 3
pbc = 0, 0, 0
r_coord='direct'
r(1,:) = $(2),0.0,0.0
r(2,:) = 0.86602540378443864676,-0.5,0.0
r(3,:) = 0.86602540378443864676,0.5,0.0
m_listing = 'by_atom'
m_coord = 'spherical'
m(1,:) = 3.0000000000000000, 0.0000000000000000, 0.0000000000000000
m(2,:) = 3.0000000000000000, 0.0000000000000000, 0.0000000000000000
m(3,:) = 3.0000000000000000, 0.0000000000000000, 0.0000000000000000
lambda_pen_listing= 'by_tag'
lambda_pen(1) = 5.0
/
&mesh
type = 'mp'
gx = $(3), $(3), $(3)
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
delta_en=0.0001
delta_q=0.0001
verbose=.false.
ni_max= 1500
/
&forces
computed=.true.
/
EOF
endef

run:
	@mkdir -p scf
	$(call corpse,-0.05)
	@$(call corpse,-0.04)
	@$(call corpse,-0.03)
	@$(call corpse,-0.02)
	@$(call corpse,-0.01)
	@$(call corpse,0.00)
	@$(call corpse,0.01)
	@$(call corpse,0.02)
	@$(call corpse,0.03)
	@$(call corpse,0.04)
	@$(call corpse,0.05)
	
post_proc:
	@$(info # dx r(1,1) r(1,1)*2.03*sqrt(3) Energy force(1,1))
	@$(call get_variable,-0.05)
	@$(call get_variable,-0.04)
	@$(call get_variable,-0.03)
	@$(call get_variable,-0.02)
	@$(call get_variable,-0.01)
	@$(call get_variable,0.00)
	@$(call get_variable,0.01)
	@$(call get_variable,0.02)
	@$(call get_variable,0.03)
	@$(call get_variable,0.04)
	@$(call get_variable,0.05)
	
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

clean:
	@rm -f out* job_out.txt in_master* job_out* force*
	@rm -fr scf

help:
	@$(info example XX:)
	@$(info Co trimer with no onsite term)
