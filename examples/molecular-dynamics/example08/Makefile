SHELL:=bash
.ONESHELL:
.NOTPARALLEL:

DYNAMOL_ROOT_DIRECTORY=/home/barreto/github/DyNaMol
DYNAMOL_COMMAND=$(DYNAMOL_ROOT_DIRECTORY)/linux/bin/DyNaMol.x
TBPARAM_DIR=$(DYNAMOL_ROOT_DIRECTORY)/tb_parameters
EXAMPLE_DIR=$(DYNAMOL_ROOT_DIRECTORY)/build/examples

define get_x =
$(shell grep 'r(1,:) =' $(EXAMPLE_DIR)/force_log_$(1).txt | cut -d"=" -f2 | cut -d"," -f1)
endef

define diag_x =
$(shell echo "$(shell grep 'r(1,:) =' $(EXAMPLE_DIR)/force_log_$(1).txt | cut -d"=" -f2 | cut -d"," -f1)*2.03*sqrt(3)"| bc -l)
endef

define get_e =
$(shell grep 'en =' $(EXAMPLE_DIR)/out_energy_$(1).txt | cut -d"=" -f2)
endef

define get_fx =
$(shell grep 'Atom 1' $(EXAMPLE_DIR)/force_log_$(1).txt |cut -d"=" -f2 |cut -d"F" -f1)
endef

define get_variable =
@$(info $(1) $(call get_x,$(1)) $(call diag_x,$(1)) $(call get_e,$(1)) $(call get_fx,$(1)))
endef

define corpse =
$(call generate_file,$(TBPARAM_DIR),$(1),1)
@mv in_master.orig in_master.txt	
@$(DYNAMOL_COMMAND) 
endef

define generate_file =
cat > in_master.orig << EOF
&calculation
processing = 'md'
post_processing= 'none'
post_processing_dir= ''
/
&units
length='ang'
energy='ev'
time='fs'
mass='g/mol'
/
&element
ne=1
symbol(1)='Pt'
/
&element_tb
filename(1) = '$(1)/pt_par_fcc_bcc_sc_lda_fl'
/
&lattice
v_factor = 2.0
v(1,:) = 1, 0, 0
v(2,:) = 0, 1, 0
v(3,:) = 0, 0, 1
/
&atom
ns = 1
na = 2
ntag = 1
tag(1)='Pt'
stag(1)= 2
pbc = 0, 0, 0
r_coord='direct'
 r(1,:) =      0.0000000     0.0000000    -0.50000000
 r(2,:) =      0.0000000     0.0000000     0.5000000
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
delta_en=0.0001
delta_q=0.0001
verbose=.false.
ni_max= 1500
/
&forces
computed=.true.
/
&md
ensemble='nvt'
t_i=0
t_f=10
dt=0.01
verbose=.true.
gamma=10.0
temperature=0.0
/
EOF
endef

run:
	$(call corpse,-0.01)
	
verify:
	@$(info $$DYNAMOL_ROOT_DIRECTORY is [$(DYNAMOL_ROOT_DIRECTORY)])
	@$(info $$DYNAMOL_COMMAND is [$(DYNAMOL_COMMAND)])

in_master:
ifeq ($(DYNAMOL_ROOT_DIRECTORY),)
@$(error DYNAMOL_ROOT_DIRECTORY is not set)
endif
ifeq ($(DYNAMOL_COMMAND),)
@$(error DYNAMOL_COMMAND is not set)
endif

clean:
	@rm -f out* in_master*

help:
	@$(info example 07:)
	@$(info Molecular Dynamics for Co trimer)
