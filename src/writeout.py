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

def writeoutput(header,z,num_nodes,nodes,num_elements_bulk,elements_bulk,displ,strain,stress,damage,tau,tau_history,force,ndime):
	
	nvert    = { 'POINT': 1, 'SEGMENT': 2, 'TRIANGLE': 3, 'QUADRANGLE': 4, 'TETRAHEDRON': 4,  'HEXAHEDRON': 8 ,  '6NPRISM': 6 }
	vtkkind  = { 'POINT': 1, 'SEGMENT': 3, 'TRIANGLE': 5, 'QUADRANGLE': 9, 'TETRAHEDRON': 10, 'HEXAHEDRON': 12 , '6NPRISM':13}
	gmshkind = { 15:'POINT', 1:'SEGMENT', 2:'TRIANGLE', 3:'QUADRANGLE', 4:'TETRAHEDRON', 5:'HEXAHEDRON', 6:'6NPRISM'}

	
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
	numero = 0
	for i in range(num_elements_bulk):
		numero += nvert[gmshkind[elements_bulk[i][1]]] + 1

	vtkfile_ini.write(str(numero)) #CUIDADO ACA!!!! ARREGLAR ESTO (ESTa BIEN ASi ?)
	vtkfile_ini.write('\n')
	for i in range(num_elements_bulk):
		vtkfile_ini.write(str(nvert[gmshkind[elements_bulk[i][1]]]))
		vtkfile_ini.write(' ')
		for n in range(nvert[gmshkind[elements_bulk[i][1]]]):
			vtkfile_ini.write(str(int(elements_bulk[i][5+n])-1))
			vtkfile_ini.write(' ')
		vtkfile_ini.write('\n')	

	vtkfile_ini.write('CELL_TYPES')
	vtkfile_ini.write(' ')
	vtkfile_ini.write(str(num_elements_bulk))
	vtkfile_ini.write('\n')
	for i in range(num_elements_bulk):
		vtkfile_ini.write(str(vtkkind[gmshkind[elements_bulk[i][1]]])+'\n')
		
	if ndime == 2:
		strain_names=['EPSXX','EPSYY','EPSXY']
		stress_names=['SIGXX','SIGYY','SIGXY']

	elif ndime == 3:
		strain_names=['EPSXX','EPSYY','EPSZZ','EPSYZ','EPSXY','EPSXZ']
		stress_names=['SIGXX','SIGYY','SIGZZ','SIGYZ','SIGXY','SIGXZ']

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
		for d in range(ndime):
			vtkfile_ini.write(' ')
			vtkfile_ini.write(str(displ[(ndime*i)+d]))
		for d in range(3-ndime):
			vtkfile_ini.write(' ')
			vtkfile_ini.write('0')
		vtkfile_ini.write('\n')	
	
	for voigt_index in range(len(strain_names)):
		vtkfile_ini.write('SCALARS')
		vtkfile_ini.write(' ')
		vtkfile_ini.write(strain_names[voigt_index])
		vtkfile_ini.write(' ')
		vtkfile_ini.write('float\n')
		vtkfile_ini.write('LOOKUP_TABLE')
		vtkfile_ini.write(' ')
		vtkfile_ini.write('default\n')
		for i in range(num_nodes):
			vtkfile_ini.write(str(strain[voigt_index][i]))
			vtkfile_ini.write('\n')
			
		vtkfile_ini.write('SCALARS')
		vtkfile_ini.write(' ')
		vtkfile_ini.write(stress_names[voigt_index])
		vtkfile_ini.write(' ')
		vtkfile_ini.write('float\n')
		vtkfile_ini.write('LOOKUP_TABLE')
		vtkfile_ini.write(' ')
		vtkfile_ini.write('default\n')
		for i in range(num_nodes):
			vtkfile_ini.write(str(stress[voigt_index][i]))
			vtkfile_ini.write('\n')
			
	vtkfile_ini.write('VECTORS')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('int_force')
	vtkfile_ini.write(' ')
	vtkfile_ini.write('float\n')
	for i in range(num_nodes):
		for d in range(ndime):
			vtkfile_ini.write(str(force[ndime*i+d]))
			vtkfile_ini.write(' ')
		for d in range(3-ndime):
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
