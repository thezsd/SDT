import os
import concurrent.futures
from Bio import PDB
from Bio.PDB.Polypeptide import three_to_one

# 设置PDB文件所在的目录
pdb_directory = '/home/lbbe06/SDT/input/all/'
output_directory = '/home/lbbe06/SDT/input/all'

# 确保输出目录存在
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# 解析PDB文件并提取序列的函数
def process_pdb_file(filename):
    pdb_path = os.path.join(pdb_directory, filename, 'receptor.pdb')
    if os.path.exists(pdb_path):
        parser = PDB.PDBParser()
        structure = parser.get_structure(filename, pdb_path)
        sequences = []
        for chain in structure.get_chains():
            chain_id = chain.id
            sequence = ''
            for residue in chain.get_residues():
                if PDB.is_aa(residue, standard=True):
                    try:
                        sequence += three_to_one(residue.get_resname())
                    except KeyError:
                        sequence += 'X'
            if sequence:
                sequences.append(f'>{filename}_{chain_id}\n{sequence}\n')
        return sequences
    return []

# 使用线程池处理所有PDB文件
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
    futures = [executor.submit(process_pdb_file, filename) for filename in os.listdir(pdb_directory)]
    with open(os.path.join(output_directory, 'output.fasta'), 'w') as fasta_file:
        for future in concurrent.futures.as_completed(futures):
            fasta_file.writelines(future.result())

print("处理完成！")
