
#Function to calculate the net force and displacement on the boundary

import sys
import numpy as np
    
def calc_force_and_disp(physical_entity_name,element_groups,physical_names,nodes,displ,stress):
        fn = 0
        ft = 0
        u_averange = [0,0]
        s_averange = [[0,0],[0,0]]
        length_tot = 0
	for eg in element_groups[physical_names[physical_entity_name][1]]:
	    tangent = [ nodes[eg[6]-1][1] - nodes[eg[5]-1][1] , nodes[eg[6]-1][2] - nodes[eg[5]-1][2] ]
	    length = np.linalg.norm(tangent)
	    length_tot += length 
	    tangent = tangent / length
	    normal =  [ tangent[1], -tangent[0] ]
	    u_averange[0] += (displ[(eg[5]-1)*2+0] + displ[(eg[6]-1)*2+0] )/2 * length
	    u_averange[1] += (displ[(eg[5]-1)*2+1] + displ[(eg[6]-1)*2+1] )/2 * length
	    s_averange += [[(stress[0][(eg[5]-1)*2] + stress[0][(eg[6]-1)*2])/2 , (stress[2][(eg[5]-1)*2] + stress[2][(eg[6]-1)*2])/2], [(stress[2][(eg[5]-1)*2] + stress[2][(eg[6]-1)*2])/2 , (stress[1][(eg[5]-1)*2] + stress[1][(eg[6]-1)*2])/2]]
	    fn += np.linalg.norm(np.dot(s_averange,normal)) * length
	    ft += np.linalg.norm(np.dot(s_averange,tangent)) * length
	    
	return [np.dot(u_averange,normal)/length_tot, np.dot(u_averange,tangent)/length_tot, fn/length_tot , ft/length_tot ]