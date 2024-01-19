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

def create_cfg_file(folder_name, center, size):
    cfg_template = f'''
receptor = /home/lbbe06/SDT/input/all/{folder_name}/receptor.pdb
input_folder = /mnt/e/pdb/pdb_without_IC50_peptide_icon_water/{folder_name}/{folder_name}
center_x = {center[0]}
center_y = {center[1]}
center_z = {center[2]}
size_x = {size[0]}
size_y = {size[1]}
size_z = {size[2]}
max_conformations =  9
'''
    with open(f'/mnt/e/pdb/cfg/{folder_name}.cfg', 'w') as cfg_file:
        cfg_file.write(cfg_template.strip())

def process_folders(base_folder):
    for folder_name in os.listdir(base_folder):
        mol2_path = os.path.join(base_folder, folder_name, f'{folder_name}_ligand.mol2')
        if os.path.isfile(mol2_path):
            atom_data = read_mol2(mol2_path)
            center, size = get_pocket_dimensions(atom_data)
            create_cfg_file(folder_name, center, size)

# Run the script
if __name__ == "__main__":
    base_folder = '/home/lbbe06/SDT/input/all'
    process_folders(base_folder)