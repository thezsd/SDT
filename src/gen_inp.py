import os
import sys

def find_pdbqt_files(base_path,save_path):
    """
    This function searches for '.pdbqt' files in the specified base_path.
    It looks for files in the pattern 'actives_final_ligand_[number]_out.pdbqt' 
    within each 'actives_final' folder in the subdirectories of base_path.

    :param base_path: Path to the base directory containing protein folders.
    :return: A list of paths to the '.pdbqt' files found.
    """
    pdbqt_files = []
    for protein_folder in os.listdir(base_path):
        protein_path= os.path.join(base_path,protein_folder,'receptor.pdb')
        actives_final_path = os.path.join(base_path, protein_folder, 'actives_final')
        if os.path.isdir(actives_final_path):
            for file in os.listdir(actives_final_path):
                if file.startswith('actives_final_ligand_') and file.endswith('_out.pdbqt'):
                    ligand_path= os.path.join(actives_final_path, file)

                    # print(f"{protein_path} {ligand_path}")
                    # pdbqt_files.append(f"{protein_path} {ligand_path}")
                    with open (save_path,'a')as f :
                        print(f"{protein_path} {ligand_path}",file=f)

    

# Base directory path
if __name__ == '__main__':

    base_path = '/mnt/e/SDT/input/all'
    save_path= '/mnt/e/SDT/src/inp.dat'
    find_pdbqt_files(base_path,save_path)

# Finding and printing the pdbqt files
# pdbqt_file_paths = find_pdbqt_files(base_path)
# for file_path in pdbqt_file_paths:
#     print(file_path)
