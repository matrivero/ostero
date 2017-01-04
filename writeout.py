#Function that generates the output in vtk format.
#Paraview can be used to postprocess the results.
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
from numpy import *

def writeoutput(header,z,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress,damage,tau,tau_history,force):
	
	vtkfile_ini = open(header+'_out.'+str(z+1)+'.vtk','w')
	vtkfile_ini.write('# vtk DataFile Version 2.0\n')
	vtkfile_ini.write('Generated with Solid Solver\n')
	vtkfile_ini.write('ASCII\n')
	vtkfile_ini.write('DATASET UNSTRUCTURED_GRID\n')
	vtkfile_ini.write('POINTS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write(str(num_nodes))
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float \n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(nodes[i][1]))
		vtkfile_ini.write(' ')
		vtkfile_ini.write(str(nodes[i][2]))
		vtkfile_ini.write(' ')
		vtkfile_ini.write(str(nodes[i][3]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('CELLS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write(str(num_elements_bulk))
	vtkfile_ini.write(' ')
	num_tets = 0
	num_quads = 0
	for i in range(num_elements_bulk):
		if elements_bulk[i][1] == 2: #TRIANGLE ELEMENT
			num_tets = num_tets + 1
		elif elements_bulk[i][1] == 3: #QUAD ELEMENT
			num_quads = num_quads + 1
	vtkfile_ini.write(str(num_tets*4 + num_quads*5)) #CUIDADO ACA!!!! ARREGLAR ESTO
	vtkfile_ini.write('\n')
	for i in range(num_elements_bulk):
		if elements_bulk[i][1] == 2: #TRIANGLE ELEMENT
			vtkfile_ini.write('3')
			vtkfile_ini.write(' ')
			vtkfile_ini.write(str(int(elements_bulk[i][5])-1))
			vtkfile_ini.write(' ')
			vtkfile_ini.write(str(int(elements_bulk[i][6])-1))
			vtkfile_ini.write(' ')	
			vtkfile_ini.write(str(int(elements_bulk[i][7])-1))
			vtkfile_ini.write('\n')	
		elif elements_bulk[i][1] == 3: #QUAD ELEMENT
			vtkfile_ini.write('4')
			vtkfile_ini.write(' ')
			vtkfile_ini.write(str(int(elements_bulk[i][5])-1))
			vtkfile_ini.write(' ')
			vtkfile_ini.write(str(int(elements_bulk[i][6])-1))
			vtkfile_ini.write(' ')	
			vtkfile_ini.write(str(int(elements_bulk[i][7])-1))
			vtkfile_ini.write(' ')	
			vtkfile_ini.write(str(int(elements_bulk[i][8])-1))
			vtkfile_ini.write('\n')	
	vtkfile_ini.write('CELL_TYPES')
	vtkfile_ini.write(' ')
	vtkfile_ini.write(str(num_elements_bulk))
	vtkfile_ini.write('\n')
	for i in range(num_elements_bulk):
		if elements_bulk[i][1] == 2: #TRIANGLE ELEMENT
			vtkfile_ini.write('5 \n')
		elif elements_bulk[i][1] == 3: #QUAD ELEMENT
			vtkfile_ini.write('9 \n')
	vtkfile_ini.write('POINT_DATA')
	vtkfile_ini.write(' ')
	vtkfile_ini.write(str(num_nodes))
	vtkfile_ini.write('\n')
	vtkfile_ini.write('VECTORS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('Displacement')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(displ[2*i]))
		vtkfile_ini.write(' ')
		vtkfile_ini.write(str(displ[(2*i)+1]))
		vtkfile_ini.write(' ')
		vtkfile_ini.write('0')
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('EPSXX')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(strain[0][i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('EPSYY')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(strain[1][i]))
		vtkfile_ini.write('\n')
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('EPSXY')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(strain[2][i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('SIGXX')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(stress[0][i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('SIGYY')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(stress[1][i]))
		vtkfile_ini.write('\n')
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('SIGXY')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(stress[2][i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('VECTORS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('int_force')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	for i in range(num_nodes):
		vtkfile_ini.write(str(force[2*i]))
		vtkfile_ini.write(' ')
		vtkfile_ini.write(str(force[(2*i)+1]))
		vtkfile_ini.write(' ')
		vtkfile_ini.write('0')
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('CELL_DATA')
	vtkfile_ini.write(' ')
	vtkfile_ini.write(str(num_elements_bulk))
	vtkfile_ini.write('\n')
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('Damage')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_elements_bulk):
		vtkfile_ini.write(str(damage[i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('Tau')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_elements_bulk):
		vtkfile_ini.write(str(tau[i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.write('SCALARS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('tau_history')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	vtkfile_ini.write('LOOKUP_TABLE')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('default\n')
	for i in range(num_elements_bulk):
		vtkfile_ini.write(str(tau_history[i]))
		vtkfile_ini.write('\n')	
	vtkfile_ini.close()
