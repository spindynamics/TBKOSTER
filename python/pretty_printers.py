import gdb
import numpy
import re
import subprocess

# GDB version
gdb_version = re.search("GNU gdb \([^)]*\) (.*)",
	str(subprocess.check_output(["gdb","--version"]))).group(1)

def register_pretty_printers():
	gdb.pretty_printers.append(lookup_function)

def lookup_function(val):
	# If val is a pointer, replace val by its target
	if val.type.code == gdb.TYPE_CODE_PTR:
		val = val.dereference()

	# Get the tag
	tag = val.type.tag

	if tag == None:
		return None

	# Match tag with pretty printer
	if tag == 'atom_tb':
		return pretty_printer_atom_tb(val)
	elif tag == 'atom':
		return pretty_printer_atom(val)
	elif tag == "calculation":
		return pretty_printer_calculation(val)
	elif tag == "charge":
		return pretty_printer_charge(val)
	#elif tag == "density_of_states":
	#	return pretty_printer_density_of_states(val)
	elif tag == "element_tb":
		return pretty_printer_element_tb(val)
	elif tag == "element":
		return pretty_printer_element(val)
	elif tag == "energy":
		return pretty_printer_energy(val)
	elif tag == "hamiltonian_tb":
		return pretty_printer_hamiltonian_tb(val)
	elif tag == "lattice":
		return pretty_printer_lattice(val)
	elif tag == "mesh":
		return pretty_printer_mesh(val)
	elif tag == "mixing":
		return pretty_printer_mixing(val)
	elif tag == "self_consistent_field":
		return pretty_printer_self_consistent_field(val)
	#elif tag == "spin_dynamics":
	#	return pretty_printer_spin_dynamics(val)
	elif tag == 'units':
		return pretty_printer_units(val)

	# Cannot find a pretty printer, return None
	return None

def get_gdb_pointer(var,var_type):
	if int(gdb_version[0]) >= 8:
		return var.cast(gdb.lookup_type(var_type)).address
	else:
		return var.cast(gdb.lookup_type(var_type).pointer())

class pretty_printer_atom_tb(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		# parent type
		parent = pretty_printer_atom(self.val["atom"])
		string = parent.to_string()
		return string

class pretty_printer_atom(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		na = int(self.val["na"])
		ntag = int(self.val["ntag"])
		stag = numpy.zeros(ntag,numpy.dtype(int))
		# header
		string = "\n&atom"
		# ns
		string += "\n ns = " + str(self.val["ns"])
		# na
		string += "\n na = " + str(self.val["na"])
		# ntag
		string += "\n ntag = " + str(self.val["ntag"])
		# stag
		ptr = get_gdb_pointer(self.val["stag"],"integer(kind=4)")
		for itag in range(0,ntag):
			string += "\n stag(" + str(itag+1) + ") = " + str((ptr+itag).dereference())
		# tag
		ptr = get_gdb_pointer(self.val["tag"],"character")
		for itag in range(0,ntag):
			var = ""
			for i in range(0,132):
				var += chr((ptr+itag*132+i).dereference())
			string += "\n tag(" + str(itag+1) + ") = '" + var.rstrip() + "'"
		# ia2ie
		ptr = get_gdb_pointer(self.val["ia2ie"],"integer(kind=4)")
		for ia in range(0,na):
			string += "\n ia2ie(" + str(ia+1) + ") = " + str((ptr+ia).dereference())
		# pbc
		string += "\n pbc = " \
		 + str(self.val["pbc"][1]) + ", " \
		 + str(self.val["pbc"][2]) + ", " \
		 + str(self.val["pbc"][3])
		# k_spiral
		string += "\n k_spiral = " \
		 + str(self.val["k_spiral"][1]) + ", " \
		 + str(self.val["k_spiral"][2]) + ", " \
		 + str(self.val["k_spiral"][3])
		# r_coord
		string += "\n r_coord = " + str(self.val["r_coord"]).replace(" ","")
		# r
		ptr = get_gdb_pointer(self.val["r"],"real(kind=8)")
		for ia in range(0,na):
			string += "\n r(" + str(ia+1) + ",:) = " \
			 + str((ptr+ia).dereference()) + ", " \
			 + str((ptr+ia+na).dereference()) + ", " \
			 + str((ptr+ia+2*na).dereference())
		# p
		ptr = get_gdb_pointer(self.val["p"],"real(kind=8)")
		for ia in range(0,na):
			string += "\n p(" + str(ia+1) + ",:) = " \
			 + str((ptr+ia).dereference()) + ", " \
			 + str((ptr+ia+na).dereference()) + ", " \
			 + str((ptr+ia+2*na).dereference())
		# m
		ptr = get_gdb_pointer(self.val["m"],"real(kind=8)")
		for ia in range(0,na):
			string += "\n m(" + str(ia+1) + ",:) = " \
			 + str((ptr+ia).dereference()) + ", " \
			 + str((ptr+ia+na).dereference()) + ", " \
			 + str((ptr+ia+2*na).dereference())
		# lambda_pen
		ptr = get_gdb_pointer(self.val["lambda_pen"],"real(kind=8)")
		for ia in range(0,na):
			string += "\n lambda_pen(" + str(ia+1) + ") = " + str((ptr+ia).dereference())
		# b_pen
		ptr = get_gdb_pointer(self.val["b_pen"],"real(kind=8)")
		for ia in range(0,na):
			string += "\n b_pen(" + str(ia+1) + ",:) = " \
			 + str((ptr+ia).dereference()) + ", " \
			 + str((ptr+ia+na).dereference()) + ", " \
			 + str((ptr+ia+2*na).dereference())
		# footer
		string += "\n /"
		return string

class pretty_printer_calculation(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		# header
		string = "\n&calculation"
		# pre_processing
		string += "\n pre_processing = " \
		 + str(self.val["pre_processing"]).replace(" ","")
		# pre_processing_dir
		ptr = get_gdb_pointer(self.val["pre_processing_dir"],"character")
		var = ""
		for i in range(0,132):
			var += chr((ptr+i).dereference())
		string += "\n pre_processing_dir = '" + var.rstrip() + "'"
		# processing
		string += "\n processing = " \
		 + str(self.val["processing"]).replace(" ","")
		# post_processing
		string += "\n post_processing = " \
		 + str(self.val["post_processing"]).replace(" ","")
		# post_processing_dir
		ptr = get_gdb_pointer(self.val["post_processing_dir"],"character")
		var = ""
		for i in range(0,132):
			var += chr((ptr+i).dereference())
		string += "\n post_processing_dir = '" + var.rstrip() + "'"
		# footer
		string += "\n /"
		return string

class pretty_printer_charge(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		e = self.val["e"]["_data"]
		no_max = int(e["no_max"])
		a = self.val["a"]["_data"]
		ns = int(a["ns"])
		na = int(a["na"])
		# header
		string = "\n&charge"
		# q_mul_in
		ptr = get_gdb_pointer(self.val["q_mul_in"],"real(kind=8)")
		for ia in range(0,na):
		 for l in range(0,3):
		  for iis in range(0,ns):
		   string += "\n q_mul_in(" + str(ia+1) + "," + str(l+1) + "," + str(iis) + ") = " \
		    + str((ptr+ia+l*na+iis*na*3).dereference())
		# q_mul_out
		ptr = get_gdb_pointer(self.val["q_mul_out"],"real(kind=8)")
		for ia in range(0,na):
		 for l in range(0,3):
		  for iis in range(0,ns):
		   string += "\n q_mul_out(" + str(ia+1) + "," + str(l+1) + "," + str(iis) + ") = " \
		    + str((ptr+ia+l*na+iis*na*3).dereference())
		# delta_q_mul
		string += "\n delta_q_mul = " + str(self.val["delta_q_mul"])
		# rho_net_in
		ptr = get_gdb_pointer(self.val["rho_net_in"],"complex(kind=8)")
		for ia in range(0,na):
		 for io1 in range(0,no_max):
		  for io2 in range(0,no_max):
		   for iis in range(0,ns):
		    string += "\n rho_net_in(" + str(ia+1) + "," + str(io1+1) + "," \
		     + str(io2+1) + "," + str(iis+1) + ") = " \
		     + str((ptr+ia+io1*na+io2*na*no_max+iis*na*no_max*no_max).dereference())
		# rho_net_out
		ptr = get_gdb_pointer(self.val["rho_net_out"],"complex(kind=8)")
		for ia in range(0,na):
		 for io1 in range(0,no_max):
		  for io2 in range(0,no_max):
		   for iis in range(0,ns):
		    string += "\n rho_net_out(" + str(ia+1) + "," + str(io1+1) + "," \
		     + str(io2+1) + "," + str(iis+1) + ") = " \
		     + str((ptr+ia+io1*na+io2*na*no_max+iis*na*no_max*no_max).dereference())
		# delta_rho_net
		string += "\n delta_rho_net = " + str(self.val["delta_rho_net"])
		# rho_net_out_diagonal
		string += "\n rho_net_out_diagonal = " + str(self.val["rho_net_out_diagonal"])
		# header
		string += "\n /"
		return string

class pretty_printer_element_tb(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		ne = int(self.val["element"]["ne"])
		# parent type
		parent = pretty_printer_element(self.val["element"])
		string = parent.to_string()
		# header
		string += "\n&element_tb"
		# filename
		ptr = get_gdb_pointer(self.val["filename"],"character")
		for ie in range(0,ne):
			var = ""
			for i in range(0,132):
				var += chr((ptr+ie*132+i).dereference())
			string += "\n filename(" + str(ie+1) + ") = '" + var.rstrip() + "'"
		# footer
		string += "\n /"
		return string

class pretty_printer_element(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		ne = int(self.val["ne"])
		no = numpy.zeros(ne,numpy.dtype(int))
		no_max = int(self.val["no_max"])
		# header
		string = "\n&element"
		# ne
		string += "\n ne = " + str(self.val["ne"])
		# symbol
		ptr = get_gdb_pointer(self.val["symbol"],"character")
		for ie in range(0,ne):
			var = ""
			for i in range(0,2):
				var += chr((ptr+ie*2+i).dereference())
			string += "\n symbol(" + str(ie+1) + ") = '" + var.rstrip() + "'"
		# no
		ptr = get_gdb_pointer(self.val["no"],"integer(kind=4)")
		for ie in range(0,ne):
			no[ie] = int((ptr+ie).dereference())
			string += "\n no(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# o
		ptr = get_gdb_pointer(self.val["o"],"integer(kind=4)")
		for ie in range(0,ne):
			string += "\n o(" + str(ie+1) + ",:) = " + str((ptr+ie).dereference())
			for io in range(1,no[ie]):
				string += ", " + str((ptr+ie+io*ne).dereference())
		# q
		ptr = get_gdb_pointer(self.val["q"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n q(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# q_s
		ptr = get_gdb_pointer(self.val["q_s"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n q_s(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# q_p
		ptr = get_gdb_pointer(self.val["q_p"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n q_p(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# q_d
		ptr = get_gdb_pointer(self.val["q_d"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n q_d(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# u_lcn
		ptr = get_gdb_pointer(self.val["u_lcn"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n u_lcn(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# u_lcn_d
		ptr = get_gdb_pointer(self.val["u_lcn_d"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n u_lcn_d(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# i_stoner_d
		ptr = get_gdb_pointer(self.val["i_stoner_d"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n i_stoner_d(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# b
		ptr = get_gdb_pointer(self.val["b"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n b(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# j_dd
		ptr = get_gdb_pointer(self.val["j_dd"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n j_dd(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# u_dd
		ptr = get_gdb_pointer(self.val["u_dd"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n u_dd(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# xi_so_p
		ptr = get_gdb_pointer(self.val["xi_so_p"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n xi_so_p(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# xi_so_d
		ptr = get_gdb_pointer(self.val["xi_so_d"],"real(kind=8)")
		for ie in range(0,ne):
			string += "\n xi_so_d(" + str(ie+1) + ") = " + str((ptr+ie).dereference())
		# footer
		string += "\n /"
		return string

class pretty_printer_energy(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		nen_k = int(self.val["nen_k"])
		# header
		string = "\n&energy"
		# smearing
		string += "\n smearing = " + str(self.val["smearing"])
		# degauss
		string += "\n degauss = " + str(self.val["degauss"])
		# fixed_fermi_level
		string += "\n fixed_fermi_level = " + str(self.val["fixed_fermi_level"])
		# en_f_ffl
		string += "\n en_f_ffl = " + str(self.val["en_f_ffl"])
		# fixed_spin_moment
		string += "\n fixed_spin_moment = " + str(self.val["fixed_spin_moment"])
		# m_fsm
		string += "\n m_fsm = " + str(self.val["m_fsm"])
		# en_min
		string += "\n en_min = " + str(self.val["en_min"])
		# en_max
		string += "\n en_max = " + str(self.val["en_max"])
		# nen_k
		string += "\n nen_k = " + str(self.val["nen_k"])
		# en_k
		ptr = get_gdb_pointer(self.val["en_k"],"real(kind=8)")
		for ien_k in range(0,nen_k):
			string += "\n en_k(" + str(ien_k+1) + ") = " + str((ptr+ien_k).dereference())
		# en_k_fsm
		if(self.val["fixed_spin_moment"]):
			ptr = get_gdb_pointer(self.val["en_k_fsm"],"real(kind=8)")
			for ien_k in range(0,nen_k):
				for iis in range(0,2):
					string += "\n en_k_fsm(" + str(ien_k+1) + "," + str(iis+1) + ") = " \
					 + str((ptr+ien_k+iis*nen_k).dereference())
		# indx
		if(not self.val["fixed_spin_moment"]):
			ptr = get_gdb_pointer(self.val["indx"],"integer(kind=4)")
			for ien_k in range(0,nen_k):
				string += "\n indx(" + str(ien_k+1) + ") = " + str((ptr+ien_k).dereference())
		# indx_fsm
		if(self.val["fixed_spin_moment"]):
			ptr = get_gdb_pointer(self.val["indx_fsm"],"integer(kind=4)")
			for ien_k in range(0,nen_k):
				for iis in range(0,2):
					string += "\n indx_fsm(" + str(ien_k+1) + "," + str(iis+1) + ") = " \
					 + str((ptr+ien_k+iis*nen_k).dereference())
		# f_k
		ptr = get_gdb_pointer(self.val["f_k"],"real(kind=8)")
		for ien_k in range(0,nen_k):
			string += "\n f_k(" + str(ien_k+1) + ") = " + str((ptr+ien_k).dereference())
		# f_k_fsm
		if(self.val["fixed_spin_moment"]):
			ptr = get_gdb_pointer(self.val["f_k_fsm"],"real(kind=8)")
			for ien_k in range(0,nen_k):
				for iis in range(0,2):
					string += "\n f_k_fsm(" + str(ien_k+1) + "," + str(iis+1) + ") = " \
					 + str((ptr+ien_k+iis*nen_k).dereference())
		# en_f
		string += "\n en_f = " + str(self.val["en_f"])
		# en_band
		string += "\n en_band = " + str(self.val["en_band"])
		# en_band_f
		string += "\n en_band_f = " + str(self.val["en_band_f"])
		# en_band_local

		# en_band_local_f

		# s_t
		string += "\n s_t = " + str(self.val["s_t"])
		# en_dc_eei
		string += "\n en_dc_eei = " + str(self.val["en_dc_eei"])
		# en_dc_lcn
		string += "\n en_dc_lcn = " + str(self.val["en_dc_lcn"])
		# en_dc_pen
		string += "\n en_dc_pen = " + str(self.val["en_dc_pen"])
		# en_in
		string += "\n en_in = " + str(self.val["en_in"])
		# en_out
		string += "\n en_out = " + str(self.val["en_out"])
		# delta_en
		string += "\n delta_en = " + str(self.val["delta_en"])
		# footer
		string += "\n /"
		return string

class pretty_printer_hamiltonian_tb(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		e = self.val["e_tb"]["_data"]["element"]
		ne = int(e["ne"])
		no = numpy.zeros(ne,numpy.dtype(int))
		ptr = get_gdb_pointer(e["no"],"integer(kind=4)")
		for ie in range(0,ne):
			no[ie] = int((ptr+ie).dereference())
		no_max = int(e["no_max"])
		a = self.val["a_tb"]["_data"]["atom"]
		ns = int(a["ns"])
		nsp = int(a["nsp"])
		na = int(a["na"])
		nn_max = int(a["nn_max"])
		ia2ie = numpy.zeros(na,numpy.dtype(int))
		ptr = get_gdb_pointer(a["ia2ie"],"integer(kind=4)")
		for ia in range(0,na):
			ia2ie[ia] = int((ptr+ia).dereference())
		# header
		string = "\n&hamiltonian_tb"
		# e_e_interaction
		string += "\n e_e_interaction = " + str(self.val["e_e_interaction"]).replace(" ","")
		# m_penalization
		string += "\n m_penalization = " + str(self.val["m_penalization"]).replace(" ","")
		# nh
		string += "\n nh = " + str(self.val["nh"])
		# iaos2ih
		ptr = get_gdb_pointer(self.val["iaos2ih"],"integer(kind=4)")
		for ia in range(0,na):
		 ie = ia2ie[ia]-1
		 for io in range(0,no[ie]):
		  for isp in range(0,nsp):
		   string += "\n iaos2ih(" + str(ia+1) + "," + str(io+1) + "," + str(isp+1) + ") = " \
		    + str((ptr+ia+io*na+isp*na*no_max).dereference())
		# en_intra
		ptr = get_gdb_pointer(self.val["en_intra"],"real(kind=8)")
		for ia in range(0,na):
		 ie = ia2ie[ia]-1
		 for io in range(0,no_max):
		  string += "\n en_intra(" + str(ia+1) + "," + str(io+1) + ") = " \
		   + str((ptr+ia+io*na).dereference())
		# h_r
		ptr = get_gdb_pointer(self.val["h_r"],"real(kind=8)")
		for ia in range(0,na):
		 for iin in range(0,nn_max+1):
		  for io1 in range(0,no_max):
		   for io2 in range(0,no_max):
		    string += "\n h_r(" + str(ia+1) + "," + str(iin+1) \
		     + "," + str(io1+1) + "," + str(io2+1) + ") = " \
		     + str((ptr+ia+iin*na+io1*na*nn_max+io2*na*nn_max*no_max).dereference())
		# s_r
		ptr = get_gdb_pointer(self.val["s_r"],"real(kind=8)")
		for ia in range(0,na):
		 for iin in range(0,nn_max+1):
		  for io1 in range(0,no_max):
		   for io2 in range(0,no_max):
		    string += "\n s_r(" + str(ia+1) + "," + str(iin+1) \
		     + "," + str(io1+1) + "," + str(io2+1) + ") = " \
		     + str((ptr+ia+iin*na+io1*na*nn_max+io2*na*nn_max*no_max).dereference())
		# delta_h_eei
		ptr = get_gdb_pointer(self.val["delta_h_eei"],"complex(kind=8)")
		for ia in range(0,na):
		 for io1 in range(0,no_max):
		  for io2 in range(0,no_max):
		   for iis in range(0,ns):
		    string += "\n delta_h_eei(" + str(ia+1) + "," + str(io1+1) \
		     + "," + str(io2+1) + "," + str(iis+1) + ") = " \
		     + str((ptr+ia+io1*na+io2*na*no_max+iis*na*no_max*no_max).dereference())
		# delta_v_lcn
		ptr = get_gdb_pointer(self.val["delta_v_lcn"],"complex(kind=8)")
		for ia in range(0,na):
		 for l in range(0,3):
		  for iis in range(0,ns):
		   string += "\n delta_v_lcn(" + str(ia+1) + "," + str(l+1) + "," + str(iis+1) + ") = " \
		    + str((ptr+ia+l*na+iis*na*3).dereference())
		# delta_v_pen
		ptr = get_gdb_pointer(self.val["delta_v_pen"],"complex(kind=8)")
		for ia in range(0,na):
		 for l in range(0,3):
		  for iis in range(0,ns):
		   string += "\n delta_v_pen(" + str(ia+1) + "," + str(l+1) + "," + str(iis+1) + ") = " \
		    + str((ptr+ia+l*na+iis*na*3).dereference())
		# footer
		string += "\n /"
		return string

class pretty_printer_lattice(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		# header
		string = "\n&lattice"
		# v
		for i in range(1,4):
			string += "\n v(" + str(i) + ",:) = " \
			 + str(self.val["v"][i][1]) + ", " \
			 + str(self.val["v"][i][2]) + ", " \
			 + str(self.val["v"][i][3])
		# footer
		string += "\n /"
		return string

class pretty_printer_mesh(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		nx = int(self.val["nx"])
		nxs = int(self.val["nxs"])
		# header
		string = "\n&mesh"
		# type
		string += "\n type = " + str(self.val["type"]).replace(" ","")
		# gx
		string += "\n gx = " \
		 + str(self.val["gx"][1]) + ", " \
		 + str(self.val["gx"][2]) + ", " \
		 + str(self.val["gx"][3])
		# dx
		string += "\n dx = " \
		 + str(self.val["dx"][1]) + ", " \
		 + str(self.val["dx"][2]) + ", " \
		 + str(self.val["dx"][3])
		# gxs
		string += "\n gxs = " + str(self.val["gxs"])
		# nxs
		string += "\n nxs = " + str(self.val["nxs"])
		# xs_label
		ptr = get_gdb_pointer(self.val["xs_label"],"character")
		for ixs in range(0,nxs):
			var = ""
			for i in range(0,2):
				var += chr((ptr+ixs*2+i).dereference())
			string += "\n filename(" + str(ixs+1) + ") = '" + var.rstrip() + "'"
		# xs_coord
		if(nxs > 0):
			string += "\n xs_coord = " + str(self.val["xs_coord"]).replace(" ","")
		# xs
		ptr = get_gdb_pointer(self.val["xs"],"real(kind=8)")
		for ixs in range(0,nxs):
			string += "\n xs(" + str(ixs+1) + ",:) = " \
			 + str((ptr+ixs).dereference()) + ", " \
			 + str((ptr+ixs+nxs).dereference()) + ", " \
			 + str((ptr+ixs+2*nxs).dereference())
		# nx
		string += "\n nx = " + str(self.val["nx"])
		# x_coord
		string += "\n x_coord = " + str(self.val["x_coord"]).replace(" ","")
		# x
		ptr = get_gdb_pointer(self.val["x"],"real(kind=8)")
		for ix in range(0,nx):
			string += "\n x(" + str(ix+1) + ",:) = " \
			 + str((ptr+ix).dereference()) + ", " \
			 + str((ptr+ix+nx).dereference()) + ", " \
			 + str((ptr+ix+2*nx).dereference())
		# w
		ptr = get_gdb_pointer(self.val["w"],"real(kind=8)")
		for ix in range(0,nx):
			string += "\n w(" + str(ix+1) + ") = " + str((ptr+ix).dereference())
		# footer
		string += "\n /"
		return string

class pretty_printer_mixing(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		# header
		string = "\n&mixing"
		# type
		string += "\n type = " + str(self.val["type"])
		# alpha
		string += "\n alpha = " + str(self.val["alpha"])
		# n_init
		string += "\n n_init = " + str(self.val["n_init"])
		# n_hist
		string += "\n n_hist = " + str(self.val["n_hist"])
		# footer
		string += "\n /"
		return string

class pretty_printer_self_consistent_field(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		# header
		string = "\n&scf"
		# ni_min
		string += "\n ni_min = " + str(self.val["ni_min"])
		# ni_max
		string += "\n ni_max = " + str(self.val["ni_max"])
		# delta_en
		string += "\n delta_en = " + str(self.val["delta_en"])
		# delta_q
		string += "\n delta_q = " + str(self.val["delta_q"])
		# footer
		string += "\n /"
		return string

class pretty_printer_units(object):
	def __init__(self, val):
		self.val = val

	def to_string(self):
		# header
		string = "\n&units"
		# energy
		string += "\n energy = " + str(self.val["energy"]).replace(" ","")
		# length
		string += "\n length = " + str(self.val["length"]).replace(" ","")
		# time
		string += "\n time = " + str(self.val["time"]).replace(" ","")
		# footer
		string += "\n /"
		return string
