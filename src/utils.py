#Function to calculate the net force and displacement on the boundary

import sys
import numpy as np


def calc_force_and_disp(physical_entity_name, element_groups, physical_names, nodes, displ, force):

	force_elem = [0, 0]
	force_aver_normal = 0
	force_aver_tangent = 0
	u_averange = [0, 0]
	length_tot = 0

	for eg in element_groups[physical_names[physical_entity_name][1]]:
		tangent = [nodes[eg[6]-1][1] - nodes[eg[5]-1][1], nodes[eg[6]-1][2] - nodes[eg[5]-1][2]]
		length = np.linalg.norm(tangent)
		length_tot += length
		tangent = tangent / length
		normal = [tangent[1], -tangent[0]]
		u_averange[0] += (displ[(eg[5]-1)*2+0] + displ[(eg[6]-1)*2+0])/2 * length
		u_averange[1] += (displ[(eg[5]-1)*2+1] + displ[(eg[6]-1)*2+1])/2 * length
		force_elem[0] = (force[(eg[5]-1)*2+0] + force[(eg[6]-1)*2+0])/2
		force_elem[1] = (force[(eg[5]-1)*2+1] + force[(eg[6]-1)*2+1])/2
		force_aver_normal += np.dot(force_elem, normal) * length
		force_aver_tangent += np.dot(force_elem, tangent) * length

	return [np.dot(u_averange, normal)/length_tot, np.dot(u_averange, tangent)/length_tot, force_aver_normal/length_tot, force_aver_tangent/length_tot]
