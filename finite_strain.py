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
for line in input_file:
	line = line.strip()
	if line.startswith('mesh_path'):
		line = line.split()
		mesh_path = line[1]
	elif line.startswith('constitutive_model'):
		line = line.split()
		model = line[1]
	elif line.startswith('time_step_size'):
		line = line.split()
		time_step_size = line[1]
	elif line.startswith('total_steps'):
		line = line.split()
		total_steps = line[1]
	elif line.startswith('space_dimension'):
		line = line.split()
		space_dim = line[1]
input_file.close()

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

#READ BOUNDARY FILE
volume_conditions = defaultdict(list)
boundary_condition_disp = defaultdict(list)
boundary_condition_press = defaultdict(list)
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
				# LINK BETWEEN BOUNDARY ELEMENTS AND VOLUME ELEMENTS
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
lame1_mu = np.zeros((num_elements_bulk))
lame2_lambda = np.zeros((num_elements_bulk))
plane = np.zeros((num_elements_bulk))
ielem = 0
for vo in volume_conditions:
	for eg in element_groups[physical_names[vo][1]]:
		Y = volume_conditions[vo][0]
		nu = volume_conditions[vo][1]
		lame2_lambda[ielem] = nu*Y/((1+nu)*(1-2*nu)) 
		lame1_mu[ielem] = Y/(2*(1+nu))
		plane[ielem] = (1.0-2.0*nu)/(1.0-nu)
		ielem += 1
#TODO: check if it is plane stress or plane deformation

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
writeout.writeoutput(model,-1,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress)

external.mod_fortran.num_nodes = num_nodes
external.mod_fortran.num_elements_bulk = num_elements_bulk
external.mod_fortran.model = model
external.mod_fortran.init()
external.mod_fortran.nodes = nodes
external.mod_fortran.elements_bulk = elements_bulk
external.mod_fortran.plane = plane
external.mod_fortran.lame1_mu = lame1_mu
external.mod_fortran.lame2_lambda = lame2_lambda

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
		
		external.mod_fortran.displ = displ
		external.mod_fortran.stress_calc_on = False
		external.mod_fortran.assembly()
		k_tot = external.mod_fortran.k_tot
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
				dot_prod = (center_of_mass_x - coord_nodes[0][0])*normal_x + (center_of_mass_y -coord_nodes[0][1])*normal_y 
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
		for bc in boundary_condition_disp:
			for eg in element_groups[physical_names[bc][1]]:
				#(eg[5]-1)*2 component X of the first node
				#(eg[5]-1)*2+1 component Y of the first node
				#(eg[6]-1)*2 component X of the second node
				#(eg[6]-1)*2+1 component Y of the second node
				for no in range(5,7):
					for nn in range(num_nodes*2):
						if(boundary_condition_disp[bc][0]): #ASK FOR FIX_X
							k_tot[(eg[no]-1)*2][nn] = 0.0	
							k_tot[nn][(eg[no]-1)*2] = 0.0	
						
						if(boundary_condition_disp[bc][1]): #ASK FOR FIX_Y
							k_tot[(eg[no]-1)*2+1][nn] = 0.0	
							k_tot[nn][(eg[no]-1)*2+1] = 0.0
	
					if(boundary_condition_disp[bc][0]): #ASK FOR FIX_X
						k_tot[(eg[no]-1)*2][(eg[no]-1)*2] = 1.0	
						r_tot[(eg[no]-1)*2] = 0.0	

					if(boundary_condition_disp[bc][1]): #ASK FOR FIX_Y
						k_tot[(eg[no]-1)*2+1][(eg[no]-1)*2+1] = 1.0	
						r_tot[(eg[no]-1)*2+1] = 0.0
					
		ddispl = np.linalg.solve(k_tot,r_tot)
		external.mod_fortran.dealloca_global_matrices()
		displ_ant = displ
		displ = displ_ant + ddispl
		norm_ddispl = np.linalg.norm(displ-displ_ant)/np.linalg.norm(displ)
		it_counter = it_counter + 1
		it_counter_global = it_counter_global + 1
		print "Newton-Raphson iteration:",it_counter
		print "Displacement increment error:",norm_ddispl
	
	external.mod_fortran.stress_calc_on = True
	external.mod_fortran.assembly()
	strain = external.mod_fortran.strain
	stress = external.mod_fortran.stress
	writeout.writeoutput(model,z,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress)
	external.mod_fortran.dealloca_stress_strain_matrices()
external.mod_fortran.dealloca_init()
