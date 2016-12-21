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
import math
import numpy as np
import writeout
from collections import defaultdict
import external

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
		print "Select a valid constitutive model for NONLINEAR geometrical treatment: ISOL, BELY, ZIEN or LAUR... bye!"
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
link_boundary_volume_elem = defaultdict(list)
gravity_present = False
transient_problem = False
boundary_file = open(sys.argv[2],'r')
for line in boundary_file:
	line = line.strip()
	if line.startswith('$VolumeDefinition'):
		readmode = 1
	elif line.startswith('$BoundaryConditionsDisplacement'):
		readmode = 2
	elif line.startswith('$BoundaryConditionsPressure'):
		readmode = 3
	elif line.startswith('$Gravity'):
		readmode = 4
	elif line.startswith('$Transient'):
		readmode = 5
	elif readmode:
		if readmode == 1: #BULK DEFINITIONS
			if line.startswith('#') or line.startswith('!'): continue
			try:
				(name, Y, nu, density) = line.split()
			except:
				print "Despite density will only matters (and will be used) if Gravity is activated, I need the input in the Volume Definition"
				print "If Gravity is not activated, just add some random number for each material and re-run. Bye :("
				sys.exit()
			name = name.replace('"','')
			volume_conditions[name].append(eval(Y))
			volume_conditions[name].append(eval(nu))
			volume_conditions[name].append(eval(density))
		elif readmode == 2: #DISPLACEMENT BOUNDARY DEFINITIONS
			if line.startswith('#') or line.startswith('!'): continue
			(name, fix_x, fix_y, dx, dy, start, end) = line.split()
			name = name.replace('"','')
			boundary_condition_disp[name].append(bool(eval(fix_x)))
			boundary_condition_disp[name].append(bool(eval(fix_y)))
			boundary_condition_disp[name].append(eval(dx))
			boundary_condition_disp[name].append(eval(dy))
			boundary_condition_disp[name].append(eval(start))
			boundary_condition_disp[name].append(eval(end))
		elif readmode == 3: #PRESSURE BOUNDARY DEFINITIONS
			if line.startswith('#') or line.startswith('!'): continue
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
		elif readmode == 4: #GRAVITY DEFINITION
			if line.startswith('#') or line.startswith('!'): continue
			gravity_present = True
			(grav_magnitude, direction_x, direction_y) = line.split()
			n_x = eval(direction_x)/np.sqrt(eval(direction_x)*eval(direction_x) + eval(direction_y)*eval(direction_y))
			n_y = eval(direction_y)/np.sqrt(eval(direction_x)*eval(direction_x) + eval(direction_y)*eval(direction_y))
			grav_direction = np.zeros((2))
			grav_direction[0] = n_x
			grav_direction[1] = n_y
		elif readmode == 5: #NEWMARK - BETA AND GAMMA VALUES
			if line.startswith('#') or line.startswith('!'): continue
			transient_problem = True
			(beta_new, gamma_new) = line.split()
			beta_newmark = eval(beta_new)
			gamma_newmark = eval(gamma_new)

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
density = np.zeros((num_elements_bulk))
ielem = 0
for vo in volume_conditions:
	for eg in element_groups[physical_names[vo][1]]:
		young[ielem] = volume_conditions[vo][0]
		poisson[ielem] = volume_conditions[vo][1]
		density[ielem] = volume_conditions[vo][2]
		ielem += 1

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
veloc = np.zeros((num_nodes*2)) #---> initial condition for velocity = 0
accel = np.zeros((num_nodes*2)) #---> only initializing the acceleration array
dt = float(time_step_size)
dt2 = float(time_step_size)*float(time_step_size)
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
external.mod_fortran.density = density
#IMPOSE GRAVITY
external.mod_fortran.gravity_calc_on = gravity_present
if (gravity_present):
	external.mod_fortran.grav_magnitude = eval(grav_magnitude)
	external.mod_fortran.grav_direction = grav_direction


###################################################################JACOBIAN AND DERIVATIVES CALCULATION###################################################################
external.mod_fortran.deriv_and_detjac_calc()
##########################################################################################################################################################################

#############################################################################VMASS CALCULATION############################################################################
external.mod_fortran.vmass_calc()
##########################################################################################################################################################################

it_counter_global = 0

for z in range(int(total_steps)):

	print ' '
	print 'Solving time step',z+1,'...'
	it_counter = 0
	eps = 0.0001
	norm_ddispl = 100*eps

	if geom_treatment == 'NONLINEAR':

		#NEWMARK - ESTIMATE NEXT SOLUTION FOR NONLINEAR MODEL (ASSUMING A NULL ACCELERATION FOR THE FIRST TIME STEP)
		if transient_problem:
			displ_mono = displ + dt*veloc + (dt2/2)*(1-2*beta_newmark)*accel 
			displ = displ_mono
			veloc_mono = veloc + (1-gamma_newmark)*dt*accel

		#IMPOSE DISPLACEMENTS 
		for bc in boundary_condition_disp:
			for eg in element_groups[physical_names[bc][1]]:
				for ii in range(len(boundary_condition_disp[bc])/6):
					if (boundary_condition_disp[bc][6*ii+4] <= (z+1) <= boundary_condition_disp[bc][6*ii+5]):
						diff_tstep = boundary_condition_disp[bc][6*ii+5] - boundary_condition_disp[bc][6*ii+4] + 1
						if (ii == 0):
							bcx_ant = 0.0
							bcy_ant = 0.0
							step_ant = 0
						else:
							bcx_ant = boundary_condition_disp[bc][6*(ii-1)+2]
							bcy_ant = boundary_condition_disp[bc][6*(ii-1)+3]
							step_ant = boundary_condition_disp[bc][6*(ii-1)+5]
						#X COORD OF FIRST AND SECOND NODE (1D SURFACE ELEMENT)
						displ[(eg[5]-1)*2] = (((boundary_condition_disp[bc][6*ii+2]-bcx_ant)/int(diff_tstep))*(z+1-step_ant)+bcx_ant)
						if (physical_names[bc][0] == 1): #0D is a node, 1D element is a line with two nodes. 
							displ[(eg[6]-1)*2] = (((boundary_condition_disp[bc][6*ii+2]-bcx_ant)/int(diff_tstep))*(z+1-step_ant)+bcx_ant)
						#Y COORD OF FIRST AND SECOND NODE (1D SURFACE ELEMENT)
						displ[(eg[5]-1)*2+1] = (((boundary_condition_disp[bc][6*ii+3]-bcy_ant)/int(diff_tstep))*(z+1-step_ant)+bcy_ant)
						if (physical_names[bc][0] == 1): #0D is a node, 1D element is a line with two nodes.
							displ[(eg[6]-1)*2+1] = (((boundary_condition_disp[bc][6*ii+3]-bcy_ant)/int(diff_tstep))*(z+1-step_ant)+bcy_ant)

	while (norm_ddispl > eps):
	
		external.mod_fortran.dealloca_global_matrices() 
		if geom_treatment == 'NONLINEAR':
			external.mod_fortran.displ = displ
			external.mod_fortran.stress_calc_on = False
			external.mod_fortran.assembly_nonlinear()

		if geom_treatment == 'LINEAR':	
			external.mod_fortran.stress_calc_on = False
			external.mod_fortran.assembly_linear()

		k_tot = external.mod_fortran.k_tot
		m_tot = external.mod_fortran.m_tot
		r_tot = external.mod_fortran.r_tot

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

		#NEWMARK STUFF
		if transient_problem:
			#INERTIAL CONTRIBUTION TO JACOBIAN
			k_tot = np.add(m_tot*(1/(beta_newmark*dt2)),k_tot)
			#INERTIAL CONTRIBUTION TO RHS
			if geom_treatment == 'NONLINEAR':
				accel = (1/(beta_newmark*dt2))*(displ-displ_mono)
				veloc = veloc_mono + gamma_newmark*dt*accel
				r_tot = np.add(np.dot(-m_tot,accel),r_tot)
			if geom_treatment == 'LINEAR':
				if it_counter_global == 0: #FIRST GLOBAL ITERATION IN LINEAR MODEL - SOLVE FOR ACCELERATION
					accel = np.linalg.solve(m_tot,r_tot - np.dot(k_tot,displ))
				displ_upd = displ + dt*veloc + (0.5-beta_newmark)*dt2*accel
				r_tot = r_tot + np.dot(m_tot*(1/(beta_newmark*dt2)),displ_upd)

		#IMPOSE DISPLACEMENT BOUNDARY CONDITIONS
		if geom_treatment == 'NONLINEAR':
			for bc in boundary_condition_disp:
				for eg in element_groups[physical_names[bc][1]]:
					#(eg[5]-1)*2 component X of the first node
					#(eg[5]-1)*2+1 component Y of the first node
					#(eg[6]-1)*2 component X of the second node
					#(eg[6]-1)*2+1 component Y of the second node
					for ii in range(len(boundary_condition_disp[bc])/6):
						if (boundary_condition_disp[bc][6*ii+4] <= (z+1) <= boundary_condition_disp[bc][6*ii+5]):
							if (physical_names[bc][0] == 1): #1D boundary element (line)
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
							else:  #0D boundary element (point)
								for no in range(5,6):  
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
					for ii in range(len(boundary_condition_disp[bc])/6):
						if (boundary_condition_disp[bc][6*ii+4] <= (z+1) <= boundary_condition_disp[bc][6*ii+5]):
							diff_tstep = boundary_condition_disp[bc][6*ii+5] - boundary_condition_disp[bc][6*ii+4] + 1 
							if (ii == 0):
								bcx_ant = 0.0
								bcy_ant = 0.0
								step_ant = 0
							else:
								bcx_ant = boundary_condition_disp[bc][6*(ii-1)+2]
								bcy_ant = boundary_condition_disp[bc][6*(ii-1)+3]
								step_ant = boundary_condition_disp[bc][6*(ii-1)+5]
							if (physical_names[bc][0] == 1): #1D boundary element (line)
								for no in range(5,7):
									if(boundary_condition_disp[bc][6*ii]): #ASK FOR FIX_X
										k_tot[(eg[no]-1)*2][(eg[no]-1)*2] = 1.0
										r_tot[(eg[no]-1)*2] = (((boundary_condition_disp[bc][6*ii+2]-bcx_ant)/int(diff_tstep))*(z+1-step_ant)+bcx_ant)
										for nn in range(num_nodes*2):
											if(nn != (eg[no]-1)*2):
												r_tot[nn] = r_tot[nn] - \
												k_tot[nn][(eg[no]-1)*2]*\
												(((boundary_condition_disp[bc][6*ii+2]-bcx_ant)/int(diff_tstep))*(z+1-step_ant)+bcx_ant)
												k_tot[nn][(eg[no]-1)*2] = 0.0
												k_tot[(eg[no]-1)*2][nn] = 0.0
		
									if(boundary_condition_disp[bc][6*ii+1]): #ASK FOR FIX_Y
										k_tot[(eg[no]-1)*2+1][(eg[no]-1)*2+1] = 1.0
										r_tot[(eg[no]-1)*2+1] = (((boundary_condition_disp[bc][6*ii+3]-bcy_ant)/int(diff_tstep))*(z+1-step_ant)+bcy_ant)
										for nn in range(num_nodes*2):
											if(nn != (eg[no]-1)*2+1):
												r_tot[nn] = r_tot[nn] - \
												k_tot[nn][(eg[no]-1)*2+1]*\
												(((boundary_condition_disp[bc][6*ii+3]-bcy_ant)/int(diff_tstep))*(z+1-step_ant)+bcy_ant)
												k_tot[nn][(eg[no]-1)*2+1] = 0.0
												k_tot[(eg[no]-1)*2+1][nn] = 0.0
							else: #0D boundary element (point)
								for no in range(5,6):
									if(boundary_condition_disp[bc][6*ii]): #ASK FOR FIX_X
										k_tot[(eg[no]-1)*2][(eg[no]-1)*2] = 1.0
										r_tot[(eg[no]-1)*2] = (((boundary_condition_disp[bc][6*ii+2]-bcx_ant)/int(diff_tstep))*(z+1-step_ant)+bcx_ant)
										for nn in range(num_nodes*2):
											if(nn != (eg[no]-1)*2):
												r_tot[nn] = r_tot[nn] - \
												k_tot[nn][(eg[no]-1)*2]*\
												(((boundary_condition_disp[bc][6*ii+2]-bcx_ant)/int(diff_tstep))*(z+1-step_ant)+bcx_ant)
												k_tot[nn][(eg[no]-1)*2] = 0.0
												k_tot[(eg[no]-1)*2][nn] = 0.0
		
									if(boundary_condition_disp[bc][6*ii+1]): #ASK FOR FIX_Y
										k_tot[(eg[no]-1)*2+1][(eg[no]-1)*2+1] = 1.0
										r_tot[(eg[no]-1)*2+1] = (((boundary_condition_disp[bc][6*ii+3]-bcy_ant)/int(diff_tstep))*(z+1-step_ant)+bcy_ant)
										for nn in range(num_nodes*2):
											if(nn != (eg[no]-1)*2+1):
												r_tot[nn] = r_tot[nn] - \
												k_tot[nn][(eg[no]-1)*2+1]*\
												(((boundary_condition_disp[bc][6*ii+3]-bcy_ant)/int(diff_tstep))*(z+1-step_ant)+bcy_ant)
												k_tot[nn][(eg[no]-1)*2+1] = 0.0
												k_tot[(eg[no]-1)*2+1][nn] = 0.0

	
		ddispl = np.linalg.solve(k_tot,r_tot)
		external.mod_fortran.dealloca_global_matrices()
		if geom_treatment == 'NONLINEAR':
			norm_ddispl = eps
			displ_ant = displ
			displ = displ_ant + ddispl
			norm_ddispl = np.linalg.norm(displ-displ_ant)/np.linalg.norm(displ)
			it_counter = it_counter + 1
			print "Newton-Raphson iteration:",it_counter
			print "Displacement increment error:",norm_ddispl
		elif geom_treatment == 'LINEAR':
			displ = ddispl
			norm_ddispl = eps
			if transient_problem:
				accel_ant = accel
				accel = (1/(beta_newmark*dt2))*(displ - displ_upd)
				veloc = veloc + dt*((1-gamma_newmark)*accel_ant + gamma_newmark*accel) 
			print "Linear solution ok!"
		it_counter_global = it_counter_global + 1

	external.mod_fortran.displ = displ
	external.mod_fortran.stress_calc_on = True
	if geom_treatment == 'NONLINEAR':
		external.mod_fortran.assembly_nonlinear()
	elif geom_treatment == 'LINEAR':
		external.mod_fortran.assembly_linear()
	strain = external.mod_fortran.strain
	stress = external.mod_fortran.stress

	external.mod_fortran.dealloca_global_matrices()
	writeout.writeoutput(header_output,z,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress)
	external.mod_fortran.dealloca_stress_strain_matrices()

external.mod_fortran.dealloca_init()

print "\nEXECUTION FINISHED SUCCESSFULLY!\n"
