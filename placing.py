import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis import transformations
from scipy.spatial.transform import Rotation as R
import random

# Running this code will place one chain at a time

u = mda.Universe("2_molec.pdb")

def rotate_and_position():
    og_sugar_chain = u.select_atoms("resname 3VA or resname 3LB or resname 0SA or resname ROH")
    protein = u.select_atoms("not (resname 3VA or resname 3LB or resname 0SA or resname ROH)")

    residues = [21] # Enter resid here, was originally meant to fill the list with multiple resids, but that feature was scratched.

    extra_angle = -30

    sugar_chains = mda.AtomGroup([], u)
    for residue in residues:
        sugar_chain = og_sugar_chain.copy()

        # amino acid hydrogen atom
        if "SER" in str(protein.residues[residue]):
            for atom in protein.residues[residue].atoms:
                if atom.name == 'HG':
                    target_coords = atom.position
                    break
        elif "THR" in str(protein.residues[residue]):
            for atom in protein.residues[residue].atoms:
                if atom.name == 'HG1':
                    target_coords = atom.position
                    break

        # Rotate with variable angles until the chains are in the expected position none of the chains are overlapping. If the chains ended up in a wrong position,
        # manually delete them from the chain.pdb file

        r = R.from_euler('X', 270, degrees=True)
        r2 = R.from_euler('Y', 90 + extra_angle, degrees=True)
        r3 = R.from_euler('Z', 45, degrees=True)

        for atom in sugar_chain.atoms:
            if "SER" in str(protein.residues[residue]):
                #atom.position = r.apply(atom.position)
                #atom.position = r2.apply(atom.position)
                atom.position = r3.apply(atom.position)
            elif "THR" in str(protein.residues[residue]):
                atom.position = r.apply(atom.position)
                atom.position = r2.apply(atom.position)
        extra_angle += 30 # Was meant to be used in the case of applying the modificatinos to multiple residues at one go, currently I recommend applying the chains one at a time and
        # visually examining the PDB after each application
        
        u.atoms.write("rotated_structure.pdb")

        # Translate
        base_sugar_atom_position = sugar_chain.residues[0].atoms[1].position

        translate_cords = [target_coords[0] - base_sugar_atom_position[0], target_coords[1] - base_sugar_atom_position[1], target_coords[2] - base_sugar_atom_position[2]]
    
        for atom in sugar_chain.atoms:
            return_to_origin = []
            for mol_cord, target_c in zip(atom.position, translate_cords):
                return_to_origin.append(mol_cord+target_c)
            atom.position = return_to_origin

        u.atoms.write("rotated_translated_structure.pdb")
        print(sugar_chains)
        sugar_chains = sugar_chains + sugar_chain
        print(sugar_chains)
    
        sug = mda.Merge(sugar_chains)
        sug.atoms.write('sug.pdb')
        final = mda.Merge(protein, sugar_chains)
        final.atoms.write("test_final.pdb")
        get_chain_coords("test_final.pdb")

def get_chain_coords(filename):
    chains = open("chain.pdb", "a")
    with open(filename, "r") as f:
        for line in f:
            if "3VA" in line or "3LB" in line or "0SA" in line:
                chains.write(line)
    chains.close()
                

rotate_and_position()