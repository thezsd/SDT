import os

def read_mol2(file_path):
    atom_data = []
    with open(file_path, 'r') as file:
        read_atoms = False
        for line in file:
            if line.startswith('@<TRIPOS>ATOM'):
                read_atoms = True
                continue
            elif line.startswith('@<TRIPOS>BOND'):
                read_atoms = False
                break
            elif read_atoms and line.strip():
                atom_info = line.split()
                atom_data.append([float(atom_info[2]), float(atom_info[3]), float(atom_info[4])])
    return atom_data

def get_pocket_dimensions(atom_data):
    min_coords = [min(coords[i] for coords in atom_data) for i in range(3)]
    max_coords = [max(coords[i] for coords in atom_data) for i in range(3)]
    center = [(max_coords[i] + min_coords[i]) / 2 for i in range(3)]
    size = [max_coords[i] - min_coords[i] for i in range(3)]
    return center, size

if __name__ == "__main":
    
    