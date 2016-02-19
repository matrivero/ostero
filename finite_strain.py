#!/usr/bin/env python

#Main of Ostero. Read input and boundary files, applies boundary conditions
#and calls external functions. Solver is called here. Numpy linalg is used. 
#Copyright (C) 2016 Matias Rivero
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#matias.rivero@bsc.es

import sys
sys.path.insert(0,"/home/matute/Desktop/mpi4py-1.3.1/Execs/lib/python/")
from mpi4py import MPI 

import math
import numpy as np
import writeout
from collections import defaultdict
import external

#=================================================================| COMMDOM |===#
sys.path.append("/home/matute/Desktop/PLE_2016ENE28/LIBPLEPP/Wrappers/Python")
sys.path.append("/home/matute/Desktop/code_saturne-3.3.4/Execs/lib/")
import Commdomm

print (MPI.get_vendor())
# Open MPI all handles are pointers
# MPICH2 they are plain 'int'.

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()
#===============================================================================#

#PRINT USAGE IF THERE ARE LESS THAN TWO ARGUMENTS
if (len(sys.argv) < 3):
	print "Usage ./finite_strain.py input_file boundary_file"
	sys.exit()

#READ INPUT FILE
input_file = open(sys.argv[1],'r')
case_name = 'NONAME' #default
for line in input_file:
	line = line.strip()
	if line.startswith('case_name'):
		line = line.split()
		case_name = line[1]
	elif line.startswith('mesh_path'):
		line = line.split()
		mesh_path = line[1]
	elif line.startswith('geometrical_treatment'):
		line = line.split()
		geom_treatment = line[1]
	elif line.startswith('constitutive_model'):
		line = line.split()
		model = line[1]
	elif line.startswith('sub_model'):
		line = line.split()
		submodel = line[1]
	elif line.startswith('time_step_size'):
		line = line.split()
		time_step_size = line[1]
	elif line.startswith('total_steps'):
		line = line.split()
		total_steps = line[1]
input_file.close()

#================================================================| COMMDOM |===#
app_type = "ostero_FEM_code"
app_name = case_name
CD = Commdomm.CommDom() 
CD.init() 
CD.set_app_type(app_type);
CD.set_app_name(app_name);
CD.set_world_comm(world_comm)

local_comm = MPI.COMM_NULL
local_comm = CD.set_mpi_comms()
local_comm.Barrier()
local_rank = local_comm.Get_rank()

namei = 'IDENTER'
namej = 'BLOCK'
#==============================================================================#

try:
	mesh_path
except NameError:
	print "You must define the mesh path... bye!"
	sys.exit()

try:
	geom_treatment
except NameError:
	print "You must define the geometrical treatment: LINEAR or NONLINEAR... bye!"
	sys.exit()

try:
	time_step_size
except NameError:
	print "You must define the time step size... bye!"
	sys.exit()

try:
	total_steps
except NameError:
	print "You must define the total number of calculation steps... bye!"
	sys.exit()

#READ MESH FILE
physical_names = defaultdict(list)
element_groups = defaultdict(list)
nodes = []
mshfile = open(mesh_path,'r')
for line in mshfile:
	line = line.strip()
	if line.startswith('$'):
		if line == '$PhysicalNames':
			readmode = 1
		elif line == '$Nodes':
			readmode = 2
		elif line == '$Elements':			
			readmode = 3
		else:
			readmode = 0
	elif readmode:
		if readmode == 1:
			if len(line.split(' ')) == 1:
				num_pn = eval(line)
			else:
				(dim, tag, name) = line.split(' ')
				name = name.replace('"','')
				physical_names[name].append(int(dim))
				physical_names[name].append(int(tag))
		if readmode == 2:
			if len(line.split(' ')) == 1:
				num_nodes = eval(line)
			else:
				line = line.split(' ') 
				line = [eval(val) for val in line]
				nodes.append(line)
		if readmode == 3:
			if len(line.split(' ')) == 1:
				num_elements = eval(line)
			else:
				line = line.split(' ') 
				line = [eval(val) for val in line]
				element_groups[line[3]].append(line)
mshfile.close()

#PRINT SOME INFORMATION
if geom_treatment == 'NONLINEAR':
	try:
		model
	except NameError:
		print "Select a valid model for NONLINEAR geometrical treatment: ISOL, BELY, ZIEN or LAUR... bye!"
		sys.exit()
	if model == 'ISOL':
		try:
			submodel
		except NameError:
			print "Select a valid submodel for NONLINEAR geometrical treatment, ISOLIN model: PLANE_STRESS or PLANE_STRAIN... bye!"
			sys.exit()
		if submodel == 'PLANE_STRESS':
			header_output = case_name+'_NONLINEAR_ISOLIN_PLANE_STRESS'
			print "ISOLINEAL MATERIAL MODEL / NONLINEAR FORMULATION / PLANE STRESS APPROXIMATION"
		elif submodel == 'PLANE_STRAIN':
			header_output = case_name+'_NONLINEAR_ISOLIN_PLANE_STRAIN'
			print "ISOLINEAL MATERIAL MODEL / NONLINEAR FORMULATION / PLANE STRAIN APPROXIMATION"  
		else:
			print "You need a submodel for ISOL model: PLANE_STRESS or PLANE_STRAIN... bye!"
			sys.exit()
	elif model == 'BELY':
		header_output = case_name+'_NONLINEAR_BELYTSCHKO'
		print "BELYTSCHKO's BOOK / NEO-HOOKEAN MATERIAL MODEL"
	elif model == 'ZIEN':
		header_output = case_name+'_NONLINEAR_ZIENKIEWICZ'
		print "ZIENKIEWICZ's BOOK / NEO-HOOKEAN MATERIAL MODEL" 
	elif model == 'LAUR': 
		header_output = case_name+'_NONLINEAR_LAURSEN'
		print "LAURSEN's BOOK / NEO-HOOKEAN MATERIAL MODEL"
	else:
		print "This model is not valid for NONLINEAR geometrical treatment... bye!"
		sys.exit()
elif geom_treatment == 'LINEAR':
	try:
		model
	except NameError:
		model = 'DUMMY' #add something in order to define the variable.

	try:
		submodel
	except NameError:
		print "Select a valid submodel for LINEAR geometrical treatment: PLANE_STRESS or PLANE_STRAIN... bye!"
		sys.exit()
	if submodel == 'PLANE_STRESS':
		header_output = case_name+'_LINEAR_PLANE_STRESS'
		print "PLANE STRESS / LINEAR FORMULATION"
	elif submodel == 'PLANE_STRAIN':
		header_output = case_name+'_LINEAR_PLANE_STRAIN'
		print "PLANE STRAIN / LINEAR FORMULATION"
	else:
		print "Select a valid submodel for LINEAR geometrical treatment: PLANE_STRESS or PLANE_STRAIN... bye!"
		sys.exit()

#CHECK IF IS A 2D or 3D CASE
accum_x = 0.0
accum_y = 0.0
accum_z = 0.0
for node in range(num_nodes):
	accum_x = accum_x + nodes[node][1]
	accum_y = accum_y + nodes[node][2]
	accum_z = accum_z + nodes[node][3]
if ((accum_x == 0.0) or (accum_y == 0.0) or (accum_z == 0.0)):
	ndime = 2
else:
	ndime = 3
	print "3D cases: available soon"
	sys.exit()

#READ BOUNDARY FILE
volume_conditions = defaultdict(list)
boundary_condition_disp = defaultdict(list)
boundary_condition_press = defaultdict(list)
boundary_condition_contact = defaultdict(list)
link_boundary_volume_elem = defaultdict(list)
boundary_file = open(sys.argv[2],'r')
for line in boundary_file:
	line = line.strip()
	if line.startswith('$VolumeDefinition'):
		readmode = 1
	elif line.startswith('$BoundaryConditionsDisplacement'):
		readmode = 2
	elif line.startswith('$BoundaryConditionsPressure'):
		readmode = 3
	elif line.startswith('$BoundaryConditionsContact'):
		readmode = 4
	elif readmode:
		if readmode == 1: #BULK DEFINITIONS
			if len(line.split()) == 1:
				num_vo = eval(line)
			else:
				(name, Y, nu) = line.split()
				name = name.replace('"','')
				volume_conditions[name].append(eval(Y))
				volume_conditions[name].append(eval(nu))
		elif readmode == 2: #DISPLACEMENT BOUNDARY DEFINITIONS
			if len(line.split()) == 1:
				num_bc_disp = eval(line)
			else:
				(name, fix_x, fix_y, dx, dy) = line.split()
				name = name.replace('"','')
				boundary_condition_disp[name].append(bool(eval(fix_x)))
				boundary_condition_disp[name].append(bool(eval(fix_y)))
				boundary_condition_disp[name].append(eval(dx))
				boundary_condition_disp[name].append(eval(dy))
		elif readmode == 3: #PRESSURE BOUNDARY DEFINITIONS
			if len(line.split()) == 1:
				num_bc_press = eval(line)
			else:
				(name, pressure) = line.split()
				name = name.replace('"','')
				boundary_condition_press[name].append(eval(pressure))
				# LINK BETWEEN BOUNDARY ELEMENTS AND VOLUME ELEMENTS FOR 2D MESHES (2D to 1D)
				for eg in element_groups[physical_names[name][1]]: #forall pressure boundary elements 
					for vo in volume_conditions: #forall volume domains
						for eg2 in element_groups[physical_names[vo][1]]: #forall volume elements
							if eg2[1] == 2: #TRIANGLE ELEMENT
								if ((eg[5] == eg2[5]) or (eg[5] == eg2[6]) or (eg[5] == eg2[7])) and \
									((eg[6] == eg2[5]) or (eg[6] == eg2[6]) or (eg[6] == eg2[7])):
									link_boundary_volume_elem[name].append(int(eg2[0]))
							elif eg2[1] == 3: #QUAD ELEMENT
								if ((eg[5] == eg2[5]) or (eg[5] == eg2[6]) or (eg[5] == eg2[7]) or (eg[5] == eg2[8])) and \
									((eg[6] == eg2[5]) or (eg[6] == eg2[6]) or (eg[6] == eg2[7]) or (eg[6] == eg2[8])):
									link_boundary_volume_elem[name].append(int(eg2[0]))
		elif readmode == 4: #CONTACT BOUNDARY DEFINITIONS
			if len(line.split()) == 1:
				name = line
				boundary_condition_contact[name].append(name)
				# LINK BETWEEN BOUNDARY ELEMENTS AND VOLUME ELEMENTS FOR 2D MESHES (2D to 1D)
				for eg in element_groups[physical_names[name][1]]: #forall contact boundary elements 
					for vo in volume_conditions: #forall volume domains
						for eg2 in element_groups[physical_names[vo][1]]: #forall volume elements
							if eg2[1] == 2: #TRIANGLE ELEMENT
								if ((eg[5] == eg2[5]) or (eg[5] == eg2[6]) or (eg[5] == eg2[7])) and \
									((eg[6] == eg2[5]) or (eg[6] == eg2[6]) or (eg[6] == eg2[7])):
									link_boundary_volume_elem[name].append(int(eg2[0]))
							elif eg2[1] == 3: #QUAD ELEMENT
								if ((eg[5] == eg2[5]) or (eg[5] == eg2[6]) or (eg[5] == eg2[7]) or (eg[5] == eg2[8])) and \
									((eg[6] == eg2[5]) or (eg[6] == eg2[6]) or (eg[6] == eg2[7]) or (eg[6] == eg2[8])):
									link_boundary_volume_elem[name].append(int(eg2[0]))
			

boundary_file.close()

nodes = np.array(nodes)

#ASSIGN MATERIALS TO VOLUMES
num_elements_bulk = 0
elements_bulk = []
for vo in volume_conditions:
	num_elements_bulk += len(element_groups[physical_names[vo][1]])
	for eg in element_groups[physical_names[vo][1]]:
		elements_bulk.append(eg)	
elements_bulk = np.array(elements_bulk)
young = np.zeros((num_elements_bulk)) 
poisson = np.zeros((num_elements_bulk))
ielem = 0
for vo in volume_conditions:
	for eg in element_groups[physical_names[vo][1]]:
		young[ielem] = volume_conditions[vo][0]
		poisson[ielem] = volume_conditions[vo][1]
		ielem += 1

#================================================================| COMMDOM |===#
local_comm.Barrier()

commij = MPI.COMM_NULL  
if( (CD.__get_app_name__() == namei) and (CD.__get_friends__(namej) == 1) ):
  commij = CD.get_mpi_commij(namej)
if( (CD.__get_app_name__() == namej) and (CD.__get_friends__(namei) == 1) ):
  commij = CD.get_mpi_commij(namei)

n_vertices_i = num_nodes
n_elements_i = num_elements_bulk

vertex_num = []
vertex_type = []
for el in range(num_elements_bulk):
	if (elements_bulk[el][1] == 2): #TRIANGLE ELEMENT
		vertex_type.append(int(10)) #---> AS ALYA TYPE
		vertex_num.append(int(elements_bulk[el][5]))
		vertex_num.append(int(elements_bulk[el][6]))
		vertex_num.append(int(elements_bulk[el][7]))
	if (elements_bulk[el][1] == 3): #QUAD ELEMENT
		vertex_type.append(int(12)) #---> AS ALYA TYPE
		vertex_num.append(int(elements_bulk[el][5]))
		vertex_num.append(int(elements_bulk[el][6]))
		vertex_num.append(int(elements_bulk[el][7]))
		vertex_num.append(int(elements_bulk[el][8]))
vertex_num_i = Commdomm.iarray(len(vertex_num))
for i in range( len(vertex_num) ): vertex_num_i[i] = vertex_num[i]
vertex_type_i = Commdomm.iarray(num_elements_bulk)
for i in range( num_elements_bulk ): vertex_type_i[i] = vertex_type[i]

vertex_coords_i = Commdomm.darray(num_nodes*ndime)
vertex_coords_j = Commdomm.darray(num_nodes*ndime)
#==============================================================================#

weight_bounda = np.zeros((2))
weight_bounda[0] = 1.0
weight_bounda[1] = 1.0
#shape1_bound = 0.5*(1-s)
#shape2_bound = 0.5*(1+s)
shape_bound = np.zeros((2,2))
shape_bound[0][0] = 0.788675 #node 1, gauss point 1 
shape_bound[1][0] = 0.211324 #node 2, gauss point 1
shape_bound[0][1] = 0.211324 #node 1, gauss point 2
shape_bound[1][1] = 0.788675 #node 2, gauss point 2
deriv_bound = np.zeros((2))
deriv_bound[0] = -0.5
deriv_bound[1] = 0.5

displ = np.zeros((num_nodes*2))
strain = np.zeros((3,num_nodes))
stress = np.zeros((3,num_nodes))
writeout.writeoutput(header_output,-1,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress)

external.mod_fortran.num_nodes = num_nodes
external.mod_fortran.num_elements_bulk = num_elements_bulk
if geom_treatment == 'NONLINEAR':
	external.mod_fortran.model = model
if ((geom_treatment == 'LINEAR') | (model == 'ISOL')):
	external.mod_fortran.submodel = submodel
external.mod_fortran.init()
external.mod_fortran.nodes = nodes
external.mod_fortran.elements_bulk = elements_bulk
external.mod_fortran.young = young
external.mod_fortran.poisson = poisson

###################################################################JACOBIAN AND DERIVATIVES CALCULATION###################################################################
external.mod_fortran.deriv_and_detjac_calc()
##########################################################################################################################################################################

#############################################################################VMASS CALCULATION############################################################################
external.mod_fortran.vmass_calc()
##########################################################################################################################################################################

it_counter_global = 0

#================================================================| COMMDOM |===#
num_max_pnode = 4
def get_element_coords_j(vert_type_i,vert_coords_i,vert_num_i,n_dist_j,d_locations_i,d_coords_j,shfunc_j):
	
	vertices_i = Commdomm.iarray(num_max_pnode)
	vol_coords_j = Commdomm.darray(num_max_pnode)
	point_j = np.zeros((ndime))
	if (n_dist_j > 0): 
		for ii in range(n_dist_j):
			for jj in range(ndime): point_j[jj] = d_coords_j[ndime*ii + jj]

			ielem = d_locations_i[ii] - 1
			if vert_type_i[ielem] == 10:
				pnode = 3
			elif vert_type_i[ielem] == 12:
				pnode = 4
			
			for jj in range(pnode): vertices_i[jj] = vert_num_i[pnode*ielem + jj]
			
			element_interpolation(vert_coords_i,vertices_i,point_j,vol_coords_j,pnode)
			
			for jj in range(pnode): shfunc_j[pnode*ii + jj] = vol_coords_j[jj]


def element_interpolation(coords_all_i,vert_i,point,vol_coords,pnode):

	elem_coord = np.zeros((pnode,ndime))

	for ii in range(pnode):
		coord_idx = vert_i[ii] - 1
		for jj in range(ndime): elem_coord[ii][jj] = coords_all_i[ndime*coord_idx + jj]
	#-----> SO FAR, SO GOOD (AT LEAST FOR TRIANGLES)
	#-----> AT THIS POINT I HAVE THE COORDINATES OF EACH VERTEX THAT SURROUNDS THE POINT_J (THE "TACHITA" POINT)
	#-----> NOW I HAVE TO FIND THE BARYCENTRIC COORDINATES OF THAT POINT_J INSIDE THE ELEMENT DEFINED BY elem_coord

	if pnode == 3:
		v0x = elem_coord[1][0] - elem_coord[0][0] # vertex b_x - vertex a_x 
		v0y = elem_coord[1][1] - elem_coord[0][1] # vertex b_y - vertex a_y
		v0 = [v0x,v0y]
		v1x = elem_coord[2][0] - elem_coord[0][0] # vertex c_x - vertex a_x
		v1y = elem_coord[2][1] - elem_coord[0][1] # vertex c_y - vertex a_y
		v1 = [v1x,v1y]
		v2x = point[0] - elem_coord[0][0] # point_j_x - vertex a_x 
		v2y = point[1] - elem_coord[0][1] # point_j_y - vertex_a_y
		v2 = [v2x,v2y]

		d00 = np.dot(v0,v0)
		d01 = np.dot(v0,v1)
		d11 = np.dot(v1,v1)
		d20 = np.dot(v2,v0)
		d21 = np.dot(v2,v1)

		denom = (d00*d11) - (d01*d01)
		vol_coords[1] = (d11*d20 - d01*d21)/denom
		vol_coords[2] = (d00*d21 - d01*d20)/denom
		vol_coords[0] = 1.0 - vol_coords[1] - vol_coords[2] 

	#elif pnode = 4:
	
def propi2propj(propi1,propi2,propi3,vert_type_i,vert_num_i,n_dist_j,d_locations_i,shpf_j,propj1,propj2,propj3):

	vertices_i = Commdomm.iarray(num_max_pnode)
	baryc_coord = np.zeros((num_max_pnode))
	prop1 = np.zeros((num_max_pnode*2))
	prop2 = np.zeros((num_max_pnode*3))
	prop3 = np.zeros((num_max_pnode*3))
	
	if (n_dist_j > 0): 
		for ii in range(n_dist_j):
			
			ielem = d_locations_i[ii] - 1
			if vert_type_i[ielem] == 10:
				pnode = 3
			elif vert_type_i[ielem] == 12:
				pnode = 4
			
			for jj in range(pnode): 
				vertices_i[jj] = vert_num_i[pnode*ielem + jj]
				baryc_coord[jj] = shpf_j[pnode*ii + jj]
				
				#displ
				prop1[2*jj] = propi1[2*(vertices_i[jj]-1)]
				prop1[2*jj+1] = propi1[2*(vertices_i[jj]-1)+1]
				#strain
				prop2[3*jj] = propi2[0][vertices_i[jj]-1]
				prop2[3*jj+1] = propi2[1][vertices_i[jj]-1]
				prop2[3*jj+2] = propi2[2][vertices_i[jj]-1]
				#stress
				prop3[3*jj] = propi3[0][vertices_i[jj]-1]
				prop3[3*jj+1] = propi3[1][vertices_i[jj]-1]
				prop3[3*jj+2] = propi3[2][vertices_i[jj]-1]

				#displ
				propj1[2*ii] = propj1[2*ii] + baryc_coord[jj]*prop1[2*jj] 
				propj1[2*ii+1] = propj1[2*ii+1] + baryc_coord[jj]*prop1[2*jj+1] 
				#strain
				propj2[3*ii] = propj2[3*ii] + baryc_coord[jj]*prop2[3*jj] 
				propj2[3*ii+1] = propj2[3*ii+1] + baryc_coord[jj]*prop2[3*jj+1] 
				propj2[3*ii+2] = propj2[3*ii+2] + baryc_coord[jj]*prop2[3*jj+2]  
				#stress
				propj3[3*ii] = propj3[3*ii] + baryc_coord[jj]*prop3[3*jj] 
				propj3[3*ii+1] = propj3[3*ii+1] + baryc_coord[jj]*prop3[3*jj+1] 
				propj3[3*ii+2] = propj3[3*ii+2] + baryc_coord[jj]*prop3[3*jj+2] 
	
				#print app_name,jj,vertices_i[jj],prop1[2*jj],prop1[2*jj+1],baryc_coord[jj]
				#print app_name,vertices_i[jj],prop3[3*jj],prop3[3*jj+1],prop3[3*jj+2]
				#print app_name,jj,vertices_i[jj],propj1[2*ii],propj1[2*ii+1],baryc_coord[jj]
	
#==============================================================================#

for z in range(int(total_steps)):

	print ' '
	print 'Solving time step',z+1,'...'
	it_counter = 0
	eps = 0.0001
	norm_ddispl = 100*eps

	if geom_treatment == 'NONLINEAR':
		#IMPOSE DISPLACEMENTS 
		for bc in boundary_condition_disp:
			for eg in element_groups[physical_names[bc][1]]:
				#X COORD OF FIRST AND SECOND NODE (1D SURFACE ELEMENT)
				displ[(eg[5]-1)*2] = (boundary_condition_disp[bc][2]/int(total_steps))*(z+1)
				displ[(eg[6]-1)*2] = (boundary_condition_disp[bc][2]/int(total_steps))*(z+1)
				#Y COORD OF FIRST AND SECOND NODE (1D SURFACE ELEMENT)
				displ[(eg[5]-1)*2+1] = (boundary_condition_disp[bc][3]/int(total_steps))*(z+1)
				displ[(eg[6]-1)*2+1] = (boundary_condition_disp[bc][3]/int(total_steps))*(z+1)
	
	while (norm_ddispl > eps):
	
		if geom_treatment == 'NONLINEAR':	
			external.mod_fortran.displ = displ
			external.mod_fortran.stress_calc_on = False
			external.mod_fortran.assembly_nonlinear()

		if geom_treatment == 'LINEAR':	
			external.mod_fortran.stress_calc_on = False
			external.mod_fortran.assembly_linear()

		k_tot = external.mod_fortran.k_tot
		r_tot = external.mod_fortran.r_tot
		k_tot_ORIG = k_tot
		r_tot_ORIG = r_tot

		#IMPOSE PRESSURE BOUNDARY CONDITIONS
		for bc in boundary_condition_press:
			cont = 0
			for eg in element_groups[physical_names[bc][1]]:
				#(eg[5]-1)*2 component X of the first node
				#(eg[5]-1)*2+1 component Y of the first node
				#(eg[6]-1)*2 component X of the second node
				#(eg[6]-1)*2+1 component Y of the second node
				coord_nodes = np.zeros((2,2))
				coord_nodes[0][0] = nodes[eg[5]-1][1] 
				coord_nodes[0][1] = nodes[eg[5]-1][2] 
				coord_nodes[1][0] = nodes[eg[6]-1][1] 
				coord_nodes[1][1] = nodes[eg[6]-1][2]
				tangent_x = coord_nodes[0][0]*deriv_bound[0] + coord_nodes[1][0]*deriv_bound[1]
				tangent_y = coord_nodes[0][1]*deriv_bound[0] + coord_nodes[1][1]*deriv_bound[1]
				normal_x = tangent_y 
				normal_y = -tangent_x
				length_bound_elem = np.sqrt(tangent_x*tangent_x + tangent_y*tangent_y)
				norm_normal = np.sqrt(normal_x*normal_x + normal_y*normal_y)
				tangent_x = tangent_x/length_bound_elem
				tangent_y = tangent_y/length_bound_elem
				normal_x = normal_x/norm_normal
				normal_y = normal_y/norm_normal
				# CHECK NORMAL DIRECTION IF IT IS INWARDS (NOT OK) OR OUTWARDS (OK)
				offset_elem_bulk = elements_bulk[0][0]
				center_of_mass_x = 0
				center_of_mass_y = 0
				if elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][1] == 2: #TRIANGLE ELEMENT
					pnode = 3
				elif elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][1] == 3: #QUAD ELEMENT
					pnode = 4
				for node in range(pnode):
					center_of_mass_x += nodes[elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][5+node] - 1][1]
					center_of_mass_y += nodes[elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][5+node] - 1][2]
				center_of_mass_x = center_of_mass_x/pnode
				center_of_mass_y = center_of_mass_y/pnode
				dot_prod = (center_of_mass_x - coord_nodes[0][0])*normal_x + (center_of_mass_y - coord_nodes[0][1])*normal_y 
				if dot_prod > 0:
					normal_x = -normal_x
					normal_y = -normal_y
					tangent_x = -tangent_x
					tangent_y = -tangent_y
				cont += 1
				#########################
				tract = np.zeros((2))
				tract[0] = normal_x * (boundary_condition_press[bc][0]/int(total_steps))*(z+1)
				tract[1] = normal_y * (boundary_condition_press[bc][0]/int(total_steps))*(z+1)
				r_elem = np.zeros((4))
				for gauss_point in range(2):
					for inode in range(2):
						idofn = inode*2 -1
						for idime in range(2):
							idofn = idofn + 1
							r_elem[idofn] = r_elem[idofn] + tract[idime]*weight_bounda[gauss_point]*length_bound_elem*shape_bound[inode][gauss_point]
				for inode in range(2):
					r_tot[(eg[5+inode]-1)*2] = r_tot[(eg[5+inode]-1)*2] + r_elem[(2*inode)] 
					r_tot[(eg[5+inode]-1)*2+1] = r_tot[(eg[5+inode]-1)*2+1] + r_elem[(2*inode)+1]

		#IMPOSE DISPLACEMENT BOUNDARY CONDITIONS
		if geom_treatment == 'NONLINEAR':
			for bc in boundary_condition_disp:
				for eg in element_groups[physical_names[bc][1]]:
					#(eg[5]-1)*2 component X of the first node
					#(eg[5]-1)*2+1 component Y of the first node
					#(eg[6]-1)*2 component X of the second node
					#(eg[6]-1)*2+1 component Y of the second node
					for no in range(5,7):
						if(boundary_condition_disp[bc][0]): #ASK FOR FIX_X
							for nn in range(num_nodes*2):
								k_tot[(eg[no]-1)*2][nn] = 0.0 #put zero on the row 	
								k_tot[nn][(eg[no]-1)*2] = 0.0 	
							k_tot[(eg[no]-1)*2][(eg[no]-1)*2] = 1.0
							r_tot[(eg[no]-1)*2] = 0.0	
						
						if(boundary_condition_disp[bc][1]): #ASK FOR FIX_Y
							for nn in range(num_nodes*2):
								k_tot[(eg[no]-1)*2+1][nn] = 0.0	
								k_tot[nn][(eg[no]-1)*2+1] = 0.0
							k_tot[(eg[no]-1)*2+1][(eg[no]-1)*2+1] = 1.0	
							r_tot[(eg[no]-1)*2+1] = 0.0

		elif geom_treatment == 'LINEAR':
			for bc in boundary_condition_disp:
				for eg in element_groups[physical_names[bc][1]]:
					#(eg[5]-1)*2 component X of the first node
					#(eg[5]-1)*2+1 component Y of the first node
					#(eg[6]-1)*2 component X of the second node
					#(eg[6]-1)*2+1 component Y of the second node
					for no in range(5,7):
						if(boundary_condition_disp[bc][0]): #ASK FOR FIX_X
							k_tot[(eg[no]-1)*2][(eg[no]-1)*2] = 1.0
							r_tot[(eg[no]-1)*2] = (boundary_condition_disp[bc][2]/int(total_steps))*(z+1)
							for nn in range(num_nodes*2):
								if(nn != (eg[no]-1)*2):
									r_tot[nn] = r_tot[nn] - k_tot[nn][(eg[no]-1)*2]*(boundary_condition_disp[bc][2]/int(total_steps))*(z+1)
									k_tot[nn][(eg[no]-1)*2] = 0.0
									k_tot[(eg[no]-1)*2][nn] = 0.0

						if(boundary_condition_disp[bc][1]): #ASK FOR FIX_Y
							k_tot[(eg[no]-1)*2+1][(eg[no]-1)*2+1] = 1.0
							r_tot[(eg[no]-1)*2+1] = (boundary_condition_disp[bc][3]/int(total_steps))*(z+1)
							for nn in range(num_nodes*2):
								if(nn != (eg[no]-1)*2+1):
									r_tot[nn] = r_tot[nn] - k_tot[nn][(eg[no]-1)*2+1]*(boundary_condition_disp[bc][3]/int(total_steps))*(z+1)
									k_tot[nn][(eg[no]-1)*2+1] = 0.0
									k_tot[(eg[no]-1)*2+1][nn] = 0.0

	
		ddispl = np.linalg.solve(k_tot,r_tot)
		#external.mod_fortran.dealloca_global_matrices()
		if geom_treatment == 'NONLINEAR':
			norm_ddispl = eps
			displ_ant = displ
			displ = displ_ant + ddispl
			norm_ddispl = np.linalg.norm(displ-displ_ant)/np.linalg.norm(displ)
			it_counter = it_counter + 1
			it_counter_global = it_counter_global + 1
			print "Newton-Raphson iteration:",it_counter
			print "Displacement increment error:",norm_ddispl
		elif geom_treatment == 'LINEAR':
			displ = ddispl
			norm_ddispl = eps
			print "Linear solution ok!"

	#external.mod_fortran.displ = displ
	#external.mod_fortran.stress_calc_on = True
	#if geom_treatment == 'NONLINEAR':
	#	external.mod_fortran.assembly_nonlinear()
	#elif geom_treatment == 'LINEAR':
	#	external.mod_fortran.assembly_linear()
	#strain = external.mod_fortran.strain
	#stress = external.mod_fortran.stress

	#================================================================| COMMDOM |===#
	vertex_coords = []
	for no in range(num_nodes):
		vertex_coords.append(nodes[no][1] + displ[2*no])
		vertex_coords.append(nodes[no][2] + displ[2*no+1]) 
	for i in range( len(vertex_coords) ): vertex_coords_i[i] = vertex_coords[i]

	#for i in range( len(vertex_coords) ): vertex_coords_j[i] = vertex_coords[i]
	#n_vertices_j = num_nodes
	
	#if app_name == 'BLOCK':
	count_tmp = 0
	id_tmp = np.zeros((num_nodes),dtype=np.int)
	map_bound_intli = np.zeros((num_nodes),dtype=np.int)
	for bc in boundary_condition_contact:
		for eg in element_groups[physical_names[bc][1]]:
			id_tmp[eg[5]-1] = 1
			id_tmp[eg[6]-1] = 1
	for no in range(num_nodes): 
		if (id_tmp[no] != 0):
			vertex_coords_j[2*int(count_tmp)] = nodes[no][1] + displ[2*no]
			vertex_coords_j[2*int(count_tmp)+1] = nodes[no][2] + displ[2*no+1]
			map_bound_intli[int(count_tmp)] = no+1	
			count_tmp = count_tmp + id_tmp[no]
	n_vertices_j = int(count_tmp)
				

	CD.locator_create2(local_comm, commij, 0.001)
	CD.locator_set_cs_mesh(n_vertices_i,n_elements_i,vertex_coords_i,vertex_num_i,vertex_type_i,n_vertices_j,vertex_coords_j,ndime) 
	#CD.save_dist_coords(z+1, local_comm)
	n_recv = CD.get_n_interior()
	n_send = CD.get_n_dist_points() # NSEND ARE THE NUMBER OF MY ELEMENTS THAT I WILL SEND TO THE OTHER CODE
					# MY ELEMENTS WHICH HAVE AT LEAST ONE NODE OF THE OTHER CODE
	dist_locations_i = Commdomm.iarray(n_send)
	dist_coords_j = Commdomm.darray(n_send*ndime)
	CD.__locator_get_dist_locations__( dist_locations_i )  
	CD.__locator_get_dist_coords__(    dist_coords_j    )
        #locations = []
        #for ii in range(n_send): locations.append( dist_locations_i[ii] ) 
        #print app_name,n_send,locations # FOR EACH DOMAIN, THIS WILL PRINT THE ELEMENTS WHICH HAVE INSIDE NODES FROM THE OTHER DOMAIN

	shapef_j = np.zeros((num_max_pnode*n_send)) # BARYCENTRIC COORDINATES - INTERPOLATION POINTS
	get_element_coords_j(vertex_type_i,vertex_coords_i,vertex_num_i,n_send,dist_locations_i,dist_coords_j,shapef_j)

	# FOR EACH APP_NAME I HAVE TO INTERPOLATE THE PROPERTIES OF THE NODES TO THE "TACHITA" POINT 
	#displ = np.zeros((num_nodes*2))
	#strain = np.zeros((3,num_nodes))
	#stress = np.zeros((3,num_nodes))
	
	displ_ij = Commdomm.darray(n_send*ndime)
	strain_ij = Commdomm.darray(n_send*3)
	stress_ij = Commdomm.darray(n_send*3)
	for ii in range(n_send*2): displ_ij[ii] = 0.0
	for ii in range(n_send*3): strain_ij[ii] = 0.0
	for ii in range(n_send*3): stress_ij[ii] = 0.0

	displcoord = np.zeros((num_nodes*2))
	for node in range(num_nodes):
		displcoord[2*node] = displ[2*node] + nodes[node][1]
		displcoord[2*node+1] = displ[2*node+1] + nodes[node][2]
	
	propi2propj(displcoord,strain,stress,vertex_type_i,vertex_num_i,n_send,dist_locations_i,shapef_j,displ_ij,strain_ij,stress_ij)

	displ_ji = Commdomm.darray(n_recv*ndime)
	strain_ji = Commdomm.darray(n_recv*3) 
	stress_ji = Commdomm.darray(n_recv*3)
	for ii in range(n_recv*2): displ_ji[ii] = 0.0
	for ii in range(n_recv*3): strain_ji[ii] = 0.0
	for ii in range(n_recv*3): stress_ji[ii] = 0.0
	
	CD.__locator_exchange_double_scalar__(displ_ij,displ_ji,2)
	CD.__locator_exchange_double_scalar__(strain_ij,strain_ji,3)
	CD.__locator_exchange_double_scalar__(stress_ij,stress_ji,3)
	
	interior_list_j = Commdomm.iarray(n_recv)
	CD.__locator_get_interior_list__(interior_list_j)

	#==============================================================================#

	proje_points_ij = Commdomm.darray(n_send*ndime)
	for ii in range(n_send*ndime): proje_points_ij[ii] = 0.0
	proje_points_ji = Commdomm.darray(n_recv*ndime)
	for ii in range(n_recv*ndime): proje_points_ji[ii] = 0.0

	slope_boun_ij = Commdomm.darray(n_send)
	for ii in range(n_send): slope_boun_ij[ii] = 0.0
	slope_boun_ji = Commdomm.darray(n_recv)
	for ii in range(n_recv): slope_boun_ji[ii] = 0.0

	normal_ij = Commdomm.darray(n_send*ndime)
	for ii in range(n_send*ndime): normal_ij[ii] = 0.0
	normal_ji = Commdomm.darray(n_recv*ndime)
	for ii in range(n_recv*ndime): normal_ji[ii] = 0.0

	if app_name == 'IDENTER': #---> IMPOSE BOUNDARY CONDITIONS ON NODES, ONLY Y DISPLACEMENT
		#dist_coords_j ---> coords of the other code that are inside of my geometry. need to check if are boundary nodes
		proje_tmp = np.zeros((2))
		acc_proje_ponts = 0
		for isend in range(n_send):
			for bc in boundary_condition_contact:
				cont = 0
				for eg in element_groups[physical_names[bc][1]]:
					#x_coord: dist_coords_j[2*isend]
					#y_coord: dist_coords_j[2*isend+1]
					p1_x = nodes[eg[6]-1][1] + displ[(eg[6]-1)*2]
					p1_y = nodes[eg[6]-1][2] + displ[(eg[6]-1)*2+1]
					p2_x = nodes[eg[5]-1][1] + displ[(eg[5]-1)*2]
					p2_y = nodes[eg[5]-1][2] + displ[(eg[5]-1)*2+1]
					slope_boun = (p1_y - p2_y)/(p1_x - p2_x)
					oo_boun = p1_y - slope_boun*p1_x
					proje_tmp[0] = dist_coords_j[2*isend]
					proje_tmp[1] = slope_boun*dist_coords_j[2*isend] + oo_boun
					tangent_x = p1_x - p2_x
					tangent_y = p1_y - p2_y
					normal_x = tangent_y 
					normal_y = -tangent_x
					# CHECK NORMAL DIRECTION IF IT IS INWARDS (NOT OK) OR OUTWARDS (OK)
					offset_elem_bulk = elements_bulk[0][0]
					center_of_mass_x = 0
					center_of_mass_y = 0
					if elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][1] == 2: #TRIANGLE ELEMENT
						pnode = 3
					elif elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][1] == 3: #QUAD ELEMENT
						pnode = 4
					for node in range(pnode):
						center_of_mass_x += nodes[elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][5+node] - 1][1]
						center_of_mass_y += nodes[elements_bulk[link_boundary_volume_elem[bc][cont] - offset_elem_bulk][5+node] - 1][2]
					center_of_mass_x = center_of_mass_x/pnode
					center_of_mass_y = center_of_mass_y/pnode
					dot_prod = (center_of_mass_x - nodes[eg[5]-1][1])*normal_x + (center_of_mass_y - nodes[eg[5]-1][2])*normal_y 
					if dot_prod > 0:
						normal_x = -normal_x
						normal_y = -normal_y
						tangent_x = -tangent_x
						tangent_y = -tangent_y
					cont += 1
					#########################
					##########################################################
					dxc = proje_tmp[0] - p2_x
					dyc = proje_tmp[1] - p2_y
					dxl = p1_x - p2_x
					dyl = p1_y - p2_y
					cross = dxc*dyl - dyc*dxl
					if (cross < 1e-6): #check if the projected point is inside the boundary	
						d1x = proje_tmp[0] - p2_x
						d1y = proje_tmp[1] - p2_y
						d1 = np.sqrt(d1x*d1x + d1y*d1y)
						d2x = proje_tmp[0] - p1_x
						d2y = proje_tmp[1] - p1_y
						d2 = np.sqrt(d2x*d2x + d2y*d2y)
						d3x = p1_x - p2_x
						d3y = p1_y - p2_y
						d3 = np.sqrt(d3x*d3x + d3y*d3y)
						if ((d1 <= d3) & (d2 <= d3)):
							proje_points_ij[2*acc_proje_ponts] = proje_tmp[0]
							proje_points_ij[2*acc_proje_ponts+1] = proje_tmp[1]
							slope_boun_ij[acc_proje_ponts] = slope_boun
							normal_ij[2*acc_proje_ponts] = normal_x
							normal_ij[2*acc_proje_ponts+1] = normal_y
							acc_proje_ponts = acc_proje_ponts + 1
					##########################################################
	
	CD.__locator_exchange_double_scalar__(proje_points_ij,proje_points_ji,2)
	CD.__locator_exchange_double_scalar__(slope_boun_ij,slope_boun_ji,1)
	CD.__locator_exchange_double_scalar__(normal_ij,normal_ji,2)

	
	#===============================================================================#

	#if app_name == 'BLOCK':
	#	for ii in range(n_recv):
	#		k_tot[(point_intli-1)*2][(point_intli-1)*2] = 1.0
	#		r_tot[(point_intli-1)*2] = 0.0
	#		k_tot[(point_intli-1)*2+1][(point_intli-1)*2+1] = 1.0
	#		r_tot[(point_intli-1)*2+1] = proje_points_ji[2*ii+1]-nodes[point_intli-1][2]
	#		for nn in range(num_nodes*2):
	#			if(nn != (point_intli-1)*2):
	#				r_tot[nn] = r_tot[nn] - 0.0
	#				k_tot[nn][(point_intli-1)*2] = 0.0
	#				k_tot[(point_intli-1)*2][nn] = 0.0
	#		for nn in range(num_nodes*2):
	#			if(nn != (point_intli-1)*2+1):
	#				r_tot[nn] = r_tot[nn] - k_tot[nn][(point_intli-1)*2+1]*(proje_points_ji[2*ii+1]-nodes[point_intli-1][2])
	#				k_tot[nn][(point_intli-1)*2+1] = 0.0
	#				k_tot[(point_intli-1)*2+1][nn] = 0.0
	
        #############################################################################################################

	#if app_name == 'BLOCK':
	#	BB = np.zeros((n_recv,num_nodes*2))
	#	hh = np.zeros((n_recv))
	#	ZZ = np.zeros((n_recv,n_recv))
	#	for ii in range(n_recv):
	#		point_intli = map_bound_intli[interior_list_j[ii]-1]
	#		x_pos = (point_intli-1)*2
	#		y_pos = (point_intli-1)*2+1
	#		hh[ii] = proje_points_ji[2*ii+1]-nodes[point_intli-1][2]
	#		BB[ii][x_pos] = -slope_boun_ji[ii]
	#		BB[ii][y_pos] = 1
	#	k_tot_block = np.bmat( [[k_tot, BB.transpose()], [BB, ZZ]] )
	#	r_tot_block = np.concatenate( [r_tot, hh] )
	#	
	#if app_name == 'BLOCK':
	#	if(n_recv == 0):
	#		displ = np.linalg.solve(k_tot,r_tot)
	#	else:
	#		displ = np.linalg.solve(k_tot_block,r_tot_block)
	#if app_name == 'IDENTER':
	#	displ = np.linalg.solve(k_tot,r_tot)
	
	external.mod_fortran.displ = displ
	external.mod_fortran.stress_calc_on = True
	if geom_treatment == 'NONLINEAR':
		external.mod_fortran.assembly_nonlinear()
	elif geom_treatment == 'LINEAR':
		external.mod_fortran.assembly_linear()
	strain = external.mod_fortran.strain
	stress = external.mod_fortran.stress

	if app_name == 'BLOCK': 
		release_nodes = np.zeros((num_nodes),dtype=np.int)
		nodes_in_contact = np.ones((n_recv),dtype=np.int)

	acc_release = 0
	while_counter = 0
        #############################################################################################################
	while True:
		if app_name == 'IDENTER':
			break
		if app_name == 'BLOCK':
			while_counter += 1
			acc_release_old = acc_release
			acc_release = 0
			for ii in range(n_recv):
				point_intli = map_bound_intli[interior_list_j[ii]-1]
				reac_x = normal_ji[2*ii]*stress[0][point_intli-1] + normal_ji[2*ii+1]*stress[2][point_intli-1]
				reac_y = normal_ji[2*ii]*stress[2][point_intli-1] + normal_ji[2*ii+1]*stress[1][point_intli-1]
				Rn = reac_x*normal_ji[2*ii] + reac_y*normal_ji[2*ii+1]
				if (Rn > 0): #Rn > 0 must be released
					release_nodes[point_intli-1] = 1
					nodes_in_contact[ii] = 0

			for ii in range(num_nodes):
				acc_release += release_nodes[ii]

			if ((acc_release_old >= acc_release) & (while_counter >= 2)):
				break

			BB = np.zeros((n_recv - acc_release,num_nodes*2))
			hh = np.zeros((n_recv - acc_release))
			ZZ = np.zeros((n_recv - acc_release,n_recv - acc_release))
			count_tmp = 0
			for ii in range(n_recv):
				point_intli = map_bound_intli[interior_list_j[ii]-1]
				if (release_nodes[point_intli-1] == 0): #apply boundary conditions
					x_pos = (point_intli-1)*2
					y_pos = (point_intli-1)*2+1
					hh[count_tmp] = proje_points_ji[2*ii+1]-nodes[point_intli-1][2]
					BB[count_tmp][x_pos] = -slope_boun_ji[ii]
					BB[count_tmp][y_pos] = 1
					count_tmp += 1
			k_tot_block = np.bmat( [[k_tot, BB.transpose()], [BB, ZZ]] )
			r_tot_block = np.concatenate( [r_tot, hh] )

		if app_name == 'BLOCK':
			if(n_recv == 0):
				displ = np.linalg.solve(k_tot,r_tot)
			else:
				displ = np.linalg.solve(k_tot_block,r_tot_block)
		if app_name == 'IDENTER':
			displ = np.linalg.solve(k_tot,r_tot)

		external.mod_fortran.dealloca_stress_strain_matrices()
		external.mod_fortran.displ = displ
		external.mod_fortran.stress_calc_on = True
		if geom_treatment == 'NONLINEAR':
			external.mod_fortran.assembly_nonlinear()
		elif geom_treatment == 'LINEAR':
			external.mod_fortran.assembly_linear()
		strain = external.mod_fortran.strain
		stress = external.mod_fortran.stress
        #############################################################################################################

	
	if app_name == 'BLOCK':
		nsend = 1
		nrecv = 0
		send  = Commdomm.iarray(nsend)
		recv  = Commdomm.iarray(nrecv)
		send[0] = int(nodes_in_contact.sum(dtype=np.int))
		CD.__mpi_sendrecv_int__(send, nsend, recv, nrecv, local_comm, commij)

		RESIDUAL = np.dot(k_tot_ORIG,displ[0:num_nodes*2]) - r_tot_ORIG
		nsend = 4*int(nodes_in_contact.sum(dtype=np.int))
		nrecv = 1
		send = Commdomm.darray(nsend)
		recv = Commdomm.darray(nrecv)
		count_tmp = 0
		for ii in range(n_recv):
			if (nodes_in_contact[ii] == 1):
				point_intli = map_bound_intli[interior_list_j[ii]-1]
				send[int(count_tmp*4)] = nodes[point_intli-1][1] + displ[(point_intli-1)*2]
				send[int(count_tmp*4+1)] = nodes[point_intli-1][2] + displ[(point_intli-1)*2+1]
				send[int(count_tmp*4+2)] = RESIDUAL[(point_intli-1)*2]
				send[int(count_tmp*4+3)] = RESIDUAL[(point_intli-1)*2+1]
				count_tmp += 1
		CD.__mpi_sendrecv_real__(send, nsend, recv, nrecv, local_comm, commij)


	if app_name == 'IDENTER':
		nsend = 0
		nrecv = 1
		send  = Commdomm.iarray(nsend)
		recv  = Commdomm.iarray(nrecv)
		CD.__mpi_sendrecv_int__(send, nsend, recv, nrecv, local_comm, commij)

		nsend = 0
		nrecv = recv[0]*4
		send  = Commdomm.darray(nsend)
		recv  = Commdomm.darray(nrecv)
		CD.__mpi_sendrecv_real__(send, nsend, recv, nrecv, local_comm, commij)
	
		#EN ESTE PUNTO EL IDENTER YA TIENE LA INFO DE LAS COORDENADAS DE LOS NODOS DE CONTACTO DEL BLOCK Y SUS RESIDUOS
		#AHORA HAY QUE PASAR ESOS RESIDUOS A LOS NODOS DEL IDENTER, HACIENDO UNA INTERPOLACION
		#PRIMERO HAY QUE VER EN QUE ELEMENTO DE BOUNDARY DEL IDENTER CAE CADA NODO DEL BLOCK Y UNA VEZ CONSEGUIDO ESO
		#HAY QUE OBTENER LA COORDENADA PARAMETRICA PARA HACER LA DISTRIBUCION DEL RESIDUO A LOS NODOS DEL IDENTER
		#Y POR ULTIMO HAY QUE CALCULAR EL EQUILIBRIO DEL IDENTER


        #############################################################################################################
	#if app_name == 'BLOCK':
	#	acc_release = 0
	#	for ii in range(n_recv):
	#		point_intli = map_bound_intli[interior_list_j[ii]-1]
	#		reac_x = normal_ji[2*ii]*stress[0][point_intli-1] + normal_ji[2*ii+1]*stress[2][point_intli-1]
	#		reac_y = normal_ji[2*ii]*stress[2][point_intli-1] + normal_ji[2*ii+1]*stress[1][point_intli-1]
	#		Rn = reac_x*normal_ji[2*ii] + reac_y*normal_ji[2*ii+1]
	#		if (Rn > 0): #Rn > 0 must be released
	#			release_nodes[point_intli-1] = 1
	#
	#	for ii in range(num_nodes):
	#		acc_release += release_nodes[ii]
	#
	#	BB = np.zeros((n_recv - acc_release,num_nodes*2))
	#	hh = np.zeros((n_recv - acc_release))
	#	ZZ = np.zeros((n_recv - acc_release,n_recv - acc_release))
	#	count_tmp = 0
	#	for ii in range(n_recv):
	#		point_intli = map_bound_intli[interior_list_j[ii]-1]
	#		if (release_nodes[point_intli-1] == 0): #apply boundary conditions
	#			x_pos = (point_intli-1)*2
	#			y_pos = (point_intli-1)*2+1
	#			hh[count_tmp] = proje_points_ji[2*ii+1]-nodes[point_intli-1][2]
	#			BB[count_tmp][x_pos] = -slope_boun_ji[ii]
	#			BB[count_tmp][y_pos] = 1
	#			count_tmp += 1
	#	k_tot_block = np.bmat( [[k_tot, BB.transpose()], [BB, ZZ]] )
	#	r_tot_block = np.concatenate( [r_tot, hh] )
	#
	#if app_name == 'BLOCK':
	#	if(n_recv == 0):
	#		displ = np.linalg.solve(k_tot,r_tot)
	#	else:
	#		displ = np.linalg.solve(k_tot_block,r_tot_block)
	#if app_name == 'IDENTER':
	#	displ = np.linalg.solve(k_tot,r_tot)
	#
	#external.mod_fortran.dealloca_stress_strain_matrices()
	#external.mod_fortran.displ = displ
	#external.mod_fortran.stress_calc_on = True
	#if geom_treatment == 'NONLINEAR':
	#	external.mod_fortran.assembly_nonlinear()
	#elif geom_treatment == 'LINEAR':
	#	external.mod_fortran.assembly_linear()
	#strain = external.mod_fortran.strain
	#stress = external.mod_fortran.stress
        #############################################################################################################

	external.mod_fortran.dealloca_global_matrices()
	writeout.writeoutput(header_output,z,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress)
	external.mod_fortran.dealloca_stress_strain_matrices()
	
	CD.locator_destroy()

external.mod_fortran.dealloca_init()
