import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

from optparse import OptionParser

############################################################
__version__ = "1.041"
# Add line mode for Oxygen to line between Nb
# No abs
############################################################

#---------------

def ReadPoscr2Cartesian(filename):
	"""
	Filename: VASP, POSCAR TYPE
	Function: Read POSCAR type 
	"""
	poscar=[]; unitcell=[]; compound=[]; position=[];scales=[];num_atoms=0
	with open(filename, 'r') as f:
		i = 1; j = 1; k =0; Selective=False; kn = 0
		for line in f:
			if line == None: break
			if i > 2 and i <6:
				unitcell.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
			elif i > 5 and i < 8:
				compound.append(line.split())
			elif i == 8 :
				if j == 1 :
					if line.split()[0][0] == 'S':
						# print 'Selective Dynamics are applied'
						Selective= True
						i = i - 1
					j = 2
				if j == 2:
					if line.split()[0][0] == 'C':
						scales=[[1.,0,0],[0,1.,0],[0,0,1.]]
					elif line.split()[0][0] == 'D':
						scales=[[ float(unitcell[0][0]), float(unitcell[0][1]), float(unitcell[0][2])],\
								[ float(unitcell[1][0]), float(unitcell[1][1]), float(unitcell[1][2])],\
								[ float(unitcell[2][0]), float(unitcell[2][1]), float(unitcell[2][2])]]
				if num_atoms == 0:
					for temp in compound[1]: num_atoms=num_atoms+int(temp)

			elif i > 8 :
				if i <= 9 + num_atoms:
					x = scales[0][0] * float(line.split()[0]) + scales[1][0] * float(line.split()[1]) + scales[2][0] * float(line.split()[2])
					y = scales[0][1] * float(line.split()[0]) + scales[1][1] * float(line.split()[1]) + scales[2][1] * float(line.split()[2])
					z = scales[0][2] * float(line.split()[0]) + scales[1][2] * float(line.split()[1]) + scales[2][2] * float(line.split()[2])
					if k <= int(compound[1][kn]):
						if Selective:
							position.append([compound[0][kn], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
						else:
							position.append([compound[0][kn], x, y, z])
					k= k+1

					if k == int(compound[1][kn]):
						kn = kn + 1
						k = 0


			if i == 8 + num_atoms:
				return unitcell, compound, position
			else:
				i = i + 1

#-----------

def ExpandPosition(unitcell, position, dimension):
	"""
	To make a supercell, expand the cell and add atoms properly
	output datas are new_unitcell, new_compound, new_position
	"""
	new_position=[]; new_unitcell=[[0,0,0],[0,0,0],[0,0,0]]; new_compound=[]; temp=[]
	for temp1 in range(len(position)):
		## expanding to a
		for temp2 in range(dimension[0]): 
			## expanding to b
			for temp3 in range(dimension[1]):
				## expanding to c
				for temp4 in range(dimension[2]):
					if len(position[temp1]) == 4: 
						new_position.append(\
							[position[temp1][0],\
							 position[temp1][1] + temp2 * float(unitcell[0][0]) + temp3 * float(unitcell[1][0]) + temp4 * float(unitcell[2][0]) , \
							 position[temp1][2] + temp2 * float(unitcell[0][1]) + temp3 * float(unitcell[1][1]) + temp4 * float(unitcell[2][1]) , \
							 position[temp1][3] + temp2 * float(unitcell[0][2]) + temp3 * float(unitcell[1][2]) + temp4 * float(unitcell[2][2]) , ])
					if len(position[temp1]) == 7:
						new_position.append(\
							[position[temp1][0],\
							 position[temp1][1] + temp2 * float(unitcell[0][0]) + temp3 * float(unitcell[1][0]) + temp4 * float(unitcell[2][0]) , \
							 position[temp1][2] + temp2 * float(unitcell[0][1]) + temp3 * float(unitcell[1][1]) + temp4 * float(unitcell[2][1]) , \
							 position[temp1][3] + temp2 * float(unitcell[0][2]) + temp3 * float(unitcell[1][2]) + temp4 * float(unitcell[2][2]) , \
							 position[temp1][4], position[temp1][5], position[temp1][6]])

	for a in unitcell: temp.append(a)
	for temp1 in range(len(temp)):
		# each axis
		for temp2 in range(len(temp[temp1])):
			#each component
			new_unitcell[temp1][temp2] = float(temp[temp1][temp2]) * dimension[temp1]

	return new_unitcell, new_position

#-----------


def vector_size(array):
	array=np.array(array)
	return (array[0]**2+array[1]**2+array[2]**2)**0.5

#-----------

def GetPolyhedronPosition(filename, opts):
	latt_vec, [atom_list,atom_num_list],cartesian_coor = ReadPoscr2Cartesian(filename)
	original_center_positions = []

	for atom_position in cartesian_coor:
		if atom_position[0] == opts.center:
			original_center_positions.append([\
				atom_position[1] +  1*(float(latt_vec[0][0]) +  float(latt_vec[1][0]) +  float(latt_vec[2][0])) , \
				atom_position[2] +  1*(float(latt_vec[0][1]) +  float(latt_vec[1][1]) +  float(latt_vec[2][1])) , \
				atom_position[3] +  1*(float(latt_vec[0][2]) +  float(latt_vec[1][2]) +  float(latt_vec[2][2])) ])			 

	super_vertex=[]
	expaned_latt_vec, expaned_cartesian_coor = ExpandPosition(latt_vec, cartesian_coor, [3,3,3] )
	for atom_position in expaned_cartesian_coor:
		if atom_position[0] == opts.vertex:
			super_vertex.append(atom_position[1:])

	center_vertex=[]
	for original_center_position in original_center_positions:
		vertex_atoms=[]
		temp_df = pd.DataFrame()
		for atom_v in super_vertex:
			temp_df = temp_df.append({'orig':original_center_position, 
									  'vertex':atom_v, 
									  'dist': vector_size(np.array(original_center_position) - np.array(atom_v))},
									 ignore_index=True)
	
		temp_df=temp_df.sort_values(['dist']).head(opts.NumberPolygon)
		vertex_atoms = [x_df[0] for x_df in temp_df.loc[:,['vertex']].values.tolist()]
		center_vertex.append([original_center_position,vertex_atoms])
	return center_vertex
#-----------

def AnisotropicQuadraticElongation(center_vertex):
	[lmx, lmy, lmz ]=[0, 0, 0]
	[vlmx, vlmy, vlmz]=[0, 0, 0]
	for i, x in enumerate(center_vertex):
		[lmx0, lmy0, lmz0]=[0, 0, 0]
		[vlmx0, vlmy0, vlmz0]=[0, 0, 0]
		center, outer = x
		'''
		l_0 calculation for different polyhedrons
		additional polyhedron formulas should be added
		formulas depand on the number of Polygonal Vertex
		'''
		if opts.NumberPolygon > 3:
			PolyVolume = ConvexHull(np.array(outer), qhull_options='Pp').volume
			if opts.NumberPolygon == 6:
				l_0 = (3.0*PolyVolume/4.0) ** (1.0/3.0)
			elif opts.NumberPolygon == 12:
				l_0 = 0.9510565163*((12/5)*PolyVolume*(1/(3+(5)**0.5)))**(1/3)
		elif opts.NumberPolygon == 2:
			l_0 = vector_size(np.array(outer[0])-np.array(outer[1]))/2
		for y in outer:
			vlmx0 = ( y[0] - center[0])/ l_0 + vlmx0
			vlmy0 = ( y[1] - center[1])/ l_0 + vlmy0
			vlmz0 = ( y[2] - center[2])/ l_0 + vlmz0
			lmx0 =( y[0] - center[0])**2/ l_0**2 + lmx0
			lmy0 =( y[1] - center[1])**2/ l_0**2 + lmy0
			lmz0 =( y[2] - center[2])**2/ l_0**2 + lmz0
		# Each QE are calculated (but later... current mode: average)
		# print(lmx0/6, lmy0/6 , lmz0/6, lmx0/6 + lmy0/6 +lmz0/6)
		vlmx = vlmx0/opts.NumberPolygon + vlmx
		vlmy = vlmy0/opts.NumberPolygon + vlmy
		vlmz = vlmz0/opts.NumberPolygon + vlmz
		lmx = lmx0/opts.NumberPolygon + lmx
		lmy = lmy0/opts.NumberPolygon + lmy
		lmz = lmz0/opts.NumberPolygon + lmz
	vlmx = vlmx / (i + 1)
	vlmy = vlmy / (i + 1)
	vlmz = vlmz / (i + 1) 

	lmx = lmx / (i + 1)
	lmy = lmy / (i + 1)
	lmz = lmz / (i + 1)
	return vlmx, vlmy, vlmz, lmx, lmy, lmz,i 

############################################################
def command_line_arg():
	usage = "usage: %prog [options] arg1 arg2"  
	par = OptionParser(usage=usage, version= __version__)

	par.add_option("-i", '--input', 
			action='append', type="string", dest='inputs',
			default=[],
			help='location of the POSCAR')

	par.add_option("-c", '--center', 
			action='store', type="string", dest='center',
			default='Nb',
			help='Polyhedra center atom name')

	par.add_option("-v", '--vertex', 
			action='store', type="string", dest='vertex',
			default='O',
			help='Polyhedra vertex atom name')

	par.add_option("--np","--numberofpolygon",
			action='store', type="int", dest='NumberPolygon',
			default=6,
			help='Set number of atoms to define Polyhedron')

	par.add_option("--detail","-d",
			action='store_true', dest='detail',
			default=False,
			help='If there is many polyhedrons, print all and averaged QE')

	par.add_option("--tocsv","-s",
			action='store', type="string", dest='saveCSV',
			default=False,
			help='Store result as csv file (exel type)')

	par.add_option("--multi","-m",
			action='store_true', dest='multi',
			default=False,
			help='Useful only one folder tree system (Crystal_structure_KNO), using with -i folderpath')
	return  par.parse_args( )

############################################################
if __name__ == '__main__':
	from time import time
	t0 = time()
	opts, args = command_line_arg()
	if opts.saveCSV: df = pd.DataFrame()

	if opts.multi:
		import os
		path= opts.inputs[0]
		center_vertex_list = []
		structure_list=[]
		for structure in os.listdir(path):
			if structure == '.DS_Store' or '.csv' in structure or '.png' in structure: pass
			else:
				for strains in os.listdir('%s/%s'%(path,structure)):
					if strains == '.DS_Store': pass
					else:
						for temp in os.listdir('%s/%s/%s/' %(path,structure,strains)):
							if temp == 'POSCAR':   filename = '%s/%s/%s/POSCAR' %(path,structure,strains)
							elif '.vasp' in temp : filename = '%s/%s/%s/%s' %(path,structure,strains, temp)

						print('[CODE] Start to reading POSCAR in ( %s )'%filename)
						center_vertex = GetPolyhedronPosition(filename,opts)
						center_vertex_list.append(center_vertex)
						structure_list.append('%s %s' %(structure.replace('Structure_','').replace('structure_',''),strains.split('_')[1]))

		print('[CODE] Finish for reading structures (case: %i ) '%len(center_vertex_list))		
		print('+-----------+--------+--------+-------+-------+-------+---------+')
		print('| Structure | Strain | N_atom |  QE_x |  QE_y |  QE_z |   QE_s  |')
		print('+-----------+--------+--------+-------+-------+-------+---------+')
		for structure_strain, center_vertex in zip(structure_list, center_vertex_list):
			[structure, strain] = structure_strain.split()
            ##dir XXXX_m1.vasp sorting
			if list(strain)[0] == 'm':
				strain = -float(list(strain)[1])
				print(strain)
			else :
				strain = float(strain)/100
			vlmx, vlmy, vlmz, lmx, lmy, lmz,i =AnisotropicQuadraticElongation(center_vertex)
			
			if opts.detail:
				for j, center_vertex_target in enumerate(center_vertex):
					vlmx,vlmy,vlmz,lmx,lmy,lmz,i  = AnisotropicQuadraticElongation([center_vertex_target])
					if opts.saveCSV:
						infodic = {'Structure': structure, 'Strain':strain, 'N_oct': j+1, 
						'$\Lambda_x$':vlmx, '$\Lambda_y$':vlmy, '$\Lambda_z$':vlmz, '$\Lambda$':sum([lmx,lmy,lmz]) }
						df = df.append(infodic,ignore_index=True)

					print('| %9.8s | %6.3f | %6.0i | %4.3f | %4.3f | %4.3f | %4.5f |'
						%( structure, strain, j+1, vlmx,vlmy,vlmz, sum([lmx,lmy,lmz])))
				print('+-----------+--------+--------+-------+-------+-------+---------+')
				print('| %9.8s | %6.3f | %6.6s | %4.3f | %4.3f | %4.3f | %4.5f |'
					%( structure, strain,'mean', vlmx,vlmy,vlmz, sum([lmx,lmy,lmz])))
				print('+-----------+--------+--------+-------+-------+-------+---------+')
			else:
				print('| %9.8s | %6.3f | %6.0i | %4.3f | %4.3f | %4.3f | %4.5f |'
					%( structure, strain, i+1, vlmx,vlmy,vlmz, sum([lmx,lmy,lmz])))
				print('+-----------+--------+--------+-------+-------+-------+---------+')
				if opts.saveCSV:
					infodic = {'Structure': structure, 'Strain':strain, 'N_oct': 0, 
					'$\Lambda_x$':vlmx, '$\Lambda_y$':vlmy, '$\Lambda_z$':vlmz, '$\Lambda$':sum([lmx,lmy,lmz]) }
					df = df.append(infodic,ignore_index=True)

		if opts.saveCSV:
			df = df.sort_values(['Structure','Strain'])
			if '.csv' in opts.saveCSV:
				df.to_csv('%s' %opts.saveCSV)
			else:
				print('[CODE] the outputfilename (%s) is not defined well (no .csv)' %opts.saveCSV)
				print('[CODE] Save file with default name: temp.csv')
				df.to_csv('temp.csv')

			import matplotlib.pyplot as plt
			import matplotlib
			df = df.sort_values(['Structure','Strain'])
			df.to_csv('temp.csv')

			matplotlib.rcParams.update({'font.size': 14})
			fig,axes = plt.subplots (2,2, figsize=(12,12))


			for phase, ax  in zip(df.Structure.unique(), [axes[0,0],axes[1,0],axes[0,1],axes[1,1]]):
			#     display(df[df.Structure == phase])
			    temp=df[df.Structure == phase]
			    ax1 =ax.twinx()
			    ax1.plot(temp.Strain, temp.iloc[:,1],'-x',c = 'red', label = temp.columns[1], markersize = 10)
			    ax1.scatter(temp.Strain, temp.iloc[:,2],color= "none", edgecolor="red", label = temp.columns[2], s=100)
			    ax1.plot(temp.Strain, temp.iloc[:,3],'-s',c = 'red', label = temp.columns[3], markersize = 10)
			    ax1.set(ylabel='Anisotropic Quadratic Distortion $\\langle\\Lambda_i\\rangle$',ylim=[-0.01,0.15],)

			    ax1.legend(loc=1)
			#     ax1.spines["right"].set_edgecolor('red')
			    ax1.tick_params(axis='y', colors='red')
			    ax1.yaxis.label.set_color('red')
			    
			    ax.plot(temp.Strain, temp.iloc[:,0],'-o',c = 'blue', label = temp.columns[0], markersize = 10)
			    ax.set(title = 'Phase:%s'%phase, xlim=[-0.055,0.055], ylim=[0.98,1.06],
			           ylabel='Isotropic Quadratic Distortion $\\langle\\Lambda\\rangle$',
			           xlabel='Strain')
			    ax.legend(loc=2)
			    
			#     ax.spines["left"].set_edgecolor('blue')
			    ax.tick_params(axis='y', colors='blue')
			    ax.yaxis.label.set_color('blue')
			    ax.grid()
			fig.tight_layout()
			fig.savefig('out.png',dpi=500)

	else:

		print('+--------+-------+-------+-------+---------+')
		print('| N_Poly |  QE_x |  QE_y |  QE_z |   QE_s  |')
		print('+--------+-------+-------+-------+---------+')
		for filename in opts.inputs:
			center_vertex = GetPolyhedronPosition(filename, opts)
			vlmx, vlmy, vlmz,lmx,lmy,lmz,i  = AnisotropicQuadraticElongation(center_vertex)	
			if opts.detail:
				for j, center_vertex_target in enumerate(center_vertex):
					vlmx, vlmy, vlmz, lmx,lmy,lmz,i  = AnisotropicQuadraticElongation([center_vertex_target])
					print('| %6.0i | %4.3f | %4.3f | %4.3f | %4.5f |'
						%( j+1, vlmx,vlmy,vlmz, sum([lmx,lmy,lmz])))
				print('+--------+-------+-------+-------+---------+')
				print('| %6.6s | %4.3f | %4.3f | %4.3f | %4.5f |'
					%( 'mean', vlmx,vlmy,vlmz, sum([lmx,lmy,lmz])))
				print('+--------+-------+-------+-------+---------+')
			else:
				print('| %6.0i | %4.3f | %4.3f | %4.3f | %4.5f |'
					%( i+1, vlmx,vlmy,vlmz, sum([lmx,lmy,lmz])))
				print('+--------+-------+-------+-------+---------+')
				print('*Average %i number of polyhedrons '%(i+1))



	t1 = time()
	print ('\nDOS plot completed! Time Used: %.2f [sec]\n' % (t1 - t0))
#-----------