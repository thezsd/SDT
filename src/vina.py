import os
import math
import subprocess
import sys


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

def read_pdb(file_path):
    atom_data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):  # Only read lines that start with ATOM
                try:
                    # Extracting the x, y, z coordinates from fixed positions
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_data.append([x, y, z])
                except ValueError:
                    # In case of a conversion issue, we skip the line
                    continue
    return atom_data

def get_pocket_dimensions(atom_data):
    min_coords = [min(coords[i] for coords in atom_data) for i in range(3)]
    max_coords = [max(coords[i] for coords in atom_data) for i in range(3)]
    center = [math.ceil((max_coords[i] + min_coords[i]) / 2) for i in range(3)]
    size = [math.ceil(max_coords[i] - min_coords[i] + 5) for i in range(3)]
    return center, size

def vina_dock(center,size,folder_name,root):
    vina_dock_command=f'''
    /public/home/wujunjun/vina_1.2.5_linux_x86_64 --receptor {root}{folder_name}/receptor.pdbqt --batch {root}{folder_name}/decoys_final_ligand* --center_x {center[0]}  --center_y {center[1]} --center_z {center[2]}  --size_x {size[0]} --size_y {size[1]} --size_z {size[2]} --dir {root}{folder_name}/decoys_final --exhaustiveness 32 --num_modes 1 > {root}{folder_name}/decoys_final/result.txt
    '''
    mkdir_command=f'''mkdir -p {root}{folder_name}/decoys_final '''

    # print(mkdir_command)
    # print(vina_dock_command)
    subprocess.run(mkdir_command, shell=True ,check=True)

    try:
        subprocess.run(vina_dock_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vina_dock in {folder_name}: {e}")


def vina_split(folder_name,root):
    vina_split_command=f'''
    /public/home/wujunjun/vina_split_1.2.5_linux_x86_64 --input {root}{folder_name}/decoys_final.pdbqt
    '''
    # print(vina_split_command)
    try:
        subprocess.run(vina_split_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vina_splite in {folder_name}: {e}")



# Run the script
if __name__ == "__main__":
     # 检查确保至少有一个参数被传递
    if len(sys.argv) < 2:
        print("Usage: script.py <root_path>")
        sys.exit(1)

    root=sys.argv[1]
    for folder_name in os.listdir(root):
        pdb_path=os.path.join(root,folder_name,'receptor.pdb')
        atom_data=read_pdb(pdb_path)
        center,size=get_pocket_dimensions(atom_data)
        vina_split(folder_name,root)
        vina_dock(center,size,folder_name,root)


    
