# 4/21/2017 Lindsay Bassman bassman@usc.edu
import getopt, sys, os
import shutil
import numpy as np
from pymatgen import Structure, Lattice
from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.transformations.site_transformations import RemoveSitesTransformation, ReplaceSiteSpeciesTransformation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.matproj.snl import StructureNL
from pymatgen import __version__ as pmgversion

#elements in TMDCs
elems = ["Mo", "W", "S", "Se", "Te"]
elidx = {"Mo": 0, "W": 1, "S": 2, "Se": 3, "Te": 4}
molec = ["MoS2", "MoSe2", "MoTe2", "WS2", "WSe2", "WTe2"]
molidx = {"MoS2": 0, "MoSe2": 1, "MoTe2": 2, "WS2": 3, "WSe2": 4, "WTe2": 5}

def invalid(reason):
	# Prints error message and exits.
	print("Invalid "+reason+".")
	sys.exit(2)

# define the history node
def history_node(t):
    pmg_github_url = 'https://github.com/materialsproject/pymatgen/tree/{}'
    return {'name': t.__module__ + "." + t.__class__.__name__,
            'url': pmg_github_url.format(pmgversion),
            'description': str(t)}

#creat SNl for submission to MP database
def create_SNL(dirbase, molecules, atoms, spc_present, num_each_spc, struct, s):
	layers = len(molecules)
	with MPRester("sm5RbuEp83T9Wo7P") as m:	
		first_mol = struct[0]
		mono_or_homo = 0
		#if system is a monolayer or homogeneous use its proper .cif file, else use generic WTe2 for heterostructures
		if (layers == 1) or all(x==first_mol for x in struct):
			mono_or_homo = 1
			if (first_mol == molec[0]):
				structure = m.get_structure_by_material_id("mp-2815") #MoS2
				ref = m.get_materials_id_references("mp-2815")
				r1 = np.array([0,2,4])
			elif (first_mol == molec[1]):
				structure = m.get_structure_by_material_id("mp-1634") #MoSe2
				ref = m.get_materials_id_references("mp-1634")
				r1 = np.array([0,2,4])
			elif (first_mol == molec[2]):
				structure = m.get_structure_by_material_id("mp-602")  #MoTe2
				ref = m.get_materials_id_references("mp-602")
				r1 = np.array([1,2,5])
			elif (first_mol == molec[3]):
				structure = m.get_structure_by_material_id("mp-224") #WS2
				ref = m.get_materials_id_references("mp-224")
				r1 = np.array([0,3,5])
			elif (first_mol == molec[4]):
				structure = m.get_structure_by_material_id("mp-1821") #WSe2
				ref = m.get_materials_id_references("mp-1821")
				r1 = np.array([0,2,4])
			elif (first_mol == molec[5]):
				structure = m.get_structure_by_material_id("mp-1019322") #WTe2	
				ref = m.get_materials_id_references("mp-1019322")
				r1 = np.array([0,3,5])
		else:
			structure = m.get_structure_by_material_id("mp-1019322") #WTe2
			ref = m.get_materials_id_references("mp-1019322")
			r1 = np.array([0,3,5])

		# initialize history
		history = []
		
		#half the height of original unit cell...to be used for vacuum length calculation later
		halfz = (structure.lattice.c)/2
		
		#make supercell if necessary
		levels = layers
		if (levels%2 == 1): levels = levels+1
		tsuper = SupercellTransformation([[1,0,0],[0,1,0],[0,0,(levels)/2]])
		history.append(history_node(tsuper))
		supercell = tsuper.apply_transformation(structure)

		#make species replacements for heterostructures with more than one layer
		levels = layers 
		if (levels%2 == 1): levels = levels+1
		#if heterostructure has more than one layer:
		if (mono_or_homo == 0):
			for i in range(0,len(molecules)):
				if (molecules[i] == 5): 
					continue
				else:			
					TMspc = elems[atoms[2*i]]
					TMloc = (levels*2) + (i%2)*(levels/2) + int(np.floor((i)/2))
					DCspc = elems[atoms[2*i+1]]
					DCloc1 = (levels - (levels/2)) - i%2*(levels/2) + int(np.floor((i)/2))
					DCloc2 = levels + i%2*(levels/2) + int(np.floor((i)/2))
					t1 = ReplaceSiteSpeciesTransformation({TMloc:TMspc})
					t2 = ReplaceSiteSpeciesTransformation({DCloc1:DCspc})
					t3 = ReplaceSiteSpeciesTransformation({DCloc2:DCspc})
					history.append(history_node(t1))
					history.append(history_node(t2))
					history.append(history_node(t3))
					supercell = t1.apply_transformation(supercell)
					supercell = t2.apply_transformation(supercell)
					supercell = t3.apply_transformation(supercell)
                
		#remove top layer of atom if necessary
		mult_factor = (layers+1)/2 -1
		r = r1 + (r1+1)*mult_factor
		tremove = RemoveSitesTransformation(r)
		if (layers%2 == 1):
			supercell = tremove.apply_transformation(supercell)
			history.append(history_node(tremove))

		#sort structure
		supercell = supercell.get_sorted_structure()		
		
		#extend z-axis cell vector to add vaccuum to supercell
		vacuum = 10.0
		old_lattice = supercell.lattice
		if (layers%2 == 1):
			new_c = old_lattice.c - halfz + vacuum
		else:
			new_c = old_lattice.c + vacuum
		new_lattice = Lattice.from_parameters(old_lattice.a, old_lattice.b, new_c, old_lattice.alpha, old_lattice.beta, old_lattice.gamma)
		final_structure = Structure(new_lattice,supercell.species,supercell.frac_coords*np.array([1., 1., (old_lattice.c/new_lattice.c)]), coords_are_cartesian=False)
		hnode = {'name': 'add vaccuum', 'url': '','description': 'increase z-direction cell vector by 10 angstroms'}
		history.append(hnode)

		#creat final SNL
		authors = [{"name": "Lindsay Bassman", "email": "bassman@usc.edu"}]
		projects = ["TMDC-Heterostructures"]
		remarks = ["MAGICS calculation of band structures of 2D TMDC stacked heterostructures"] 
		final_snl = StructureNL(final_structure, authors, projects=projects, remarks=remarks, references=ref, history=history)
		
		#optionally write POSCAR file
		poscar = Poscar(final_structure, s)
		poscar.write_file(dirbase+"POSCAR", direct=False)

		#submit snl
	#with MPRester("sm5RbuEp83T9Wo7P",endpoint="https://www.materialsproject.org/rest/v1") as m2:
	#	m2.submit_snl(final_snl)
	
def main():
	# Parse command line arguments.
	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:d:")
	except getopt.GetoptError as err:
		invalid(err)
	force = True
	for o, a in opts:
		if (o == "-s"):
			s = a
		elif (o == "-d"):
			dirbase = a
	struct = s.split("-")
	if len(struct) == 0:
		print("Invalid structure.")
		sys.exit(2)
	molecules = [0]* (len(struct))
	atoms = [0]*(2*len(struct))
	for i, elem in enumerate(struct):
		if elem not in molec:
			print("Invalid structure.")
		for j, mol in enumerate(molec):
			if mol == elem:
				molecules[i] = j
				atoms[2*i] = int(np.floor(j/3))
				atoms[2*i+1] = j%3 + 2
	num_each_spc = [0]*5
	spc_present = []
	for i, idx in enumerate(atoms):
		# Double count for every other element.
		num_each_spc[idx] += 1+i%2
	for idx, contains in enumerate(num_each_spc):
		if contains > 0:
			spc_present.append(idx)
	# Create folders.
	dirbase = s+"/"
	if not os.path.exists(dirbase):
		os.makedirs(dirbase)
	#create structure SNL
	create_SNL(dirbase, molecules, atoms, spc_present, num_each_spc, struct, s)

if __name__ == "__main__":
	main()
