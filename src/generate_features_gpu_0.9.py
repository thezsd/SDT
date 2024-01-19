#!/usr/bin/env python
"""
使用了gpu来实现接触特征的产生，但是在原子类型的识别上，较为落后，没有细化原子类别，同时修改了默认的层数为60
同时提供了将ligand转化为pdb格式的能力
"""



import numpy as np
import pandas as pd
import mdtraj as mt
import itertools
import multiprocessing
import sys, os
from biopandas.mol2 import PandasMol2
import argparse
from argparse import RawDescriptionHelpFormatter
from rdkit import Chem
import subprocess as sp
import torch
import itertools
import time



# Define required atomic element type
elements_ligand = ["H", "C", "O", "N", "P", "S", "HAX", "DU"]
elements_protein = ["H", "C", "O", "N", "P", "S", "HAX", "DU"]

ALL_ELEMENTS = ["H", "C", "O", "N", "P", "S", "HAX", "DU"]


class Molecule(object):
    """Small molecule parser object with Rdkit package.

    Parameters
    ----------
    in_format : str, default = 'smile'
        Input information (file) format. 
        Options: smile, pdb, sdf, mol2, mol

    Attributes
    ----------
    molecule_ : rdkit.Chem.Molecule object
    mol_file : str
        The input file name or Smile string
    converter_ : dict, dict of rdkit.Chem.MolFrom** methods
        The file loading method dictionary. The keys are:
        pdb, sdf, mol2, mol, smile


    """

    def __init__(self, in_format="smile"):

        self.format = in_format
        self.molecule_ = None
        self.mol_file = None
        self.converter_ = None
        self.mol_converter()

    def mol_converter(self):
        """The converter methods are stored in a dictionary.

        Returns
        -------
        self : return an instance of itself

        """
        self.converter_ = {
            "pdb": Chem.MolFromPDBFile,
            "mol2": Chem.MolFromMol2File,
            "mol": Chem.MolFromMolFile,
            "smile": Chem.MolFromSmiles,
            "sdf": Chem.MolFromMolBlock,
            "pdbqt": self.babel_converter,
        }

        return self

    def babel_converter(self, mol_file, output):
        if os.path.exists(mol_file):
            try:
                cmd = 'obabel %s -O %s > /dev/null' % (mol_file, output)
                job = sp.Popen(cmd, shell=True)
                job.communicate()

                self.molecule_ = self.converter_['pdb']()
                return self.molecule_
            except:
                return None

    def load_molecule(self, mol_file):
        """Load a molecule to have a rdkit.Chem.Molecule object

        Parameters
        ----------
        mol_file : str
            The input file name or SMILE string

        Returns
        -------
        molecule : rdkit.Chem.Molecule object
            The molecule object

        """

        self.mol_file = mol_file
        if not os.path.exists(self.mol_file):
            print("Molecule file not exists. ")
            return None

        if self.format not in ["mol2", "mol", "pdb", "sdf", "pdbqt"]:
            print("File format is not correct. ")
            return None
        else:
            try:
                self.molecule_ = self.converter_[self.format](self.mol_file)
            except RuntimeError:
                return None

            return self.molecule_

class ProteinParser(object):
    """Featurization of Protein-Ligand Complex based on
    onion-shape distance counts of atom-types.

    Parameters
    ----------
    pdb_fn : str
        The input pdb file name. The file must be in PDB format.

    Attributes
    ----------
    pdb : mdtraj.Trajectory
        The mdtraj.trajectory object containing the pdb.
    receptor_indices : np.ndarray
        The receptor (protein) atom indices in mdtraj.Trajectory
    rec_ele : np.ndarray
        The element types of each of the atoms in the receptor
    pdb_parsed_ : bool
        Whether the pdb file has been parsed.
    distance_computed : bool
        Whether the distances between atoms in receptor and ligand has been computed.

    Examples
    --------
    >>> pdb = ProteinParser("input.pdb")
    >>> pdb.parsePDB('protein and chainid 0')
    >>> pdb.coordinates_
    >>> print(pdb.rec_ele)

    """

    def __init__(self, pdb_fn):
        self.pdb = mt.load_pdb(pdb_fn)

        self.receptor_indices = np.array([])
        self.rec_ele = np.array([])

        self.pdb_parsed_ = False
        self.coordinates_ = None

    def get_coordinates(self):
        """
        Get the coordinates in the pdb file given the receptor indices.

        Returns
        -------
        self : an instance of itself

        """
        self.coordinates_ = self.pdb.xyz[0][self.receptor_indices]

        return self

    def parsePDB(self, rec_sele="protein"):
        """
        Parse the pdb file and get the detail information of the protein.

        Parameters
        ----------
        rec_sele : str,
            The string for protein selection. Please refer to the following link.

        References
        ----------
        Mdtraj atom selection language: http://mdtraj.org/development/atom_selection.html

        Returns
        -------

        """
        top = self.pdb.topology
        # obtain the atom indices of the protein
        self.receptor_indices = top.select(rec_sele)
        _table, _bond = top.to_dataframe()

        # fetch the element type of each one of the protein atom
        self.rec_ele = _table['element'][self.receptor_indices].values
        # fetch the coordinates of each one of the protein atom
        self.get_coordinates()

        self.pdb_parsed_ = True

        return self


class LigandParser(object):
    """Parse the ligand with biopanda to obtain coordinates and elements.

    Parameters
    ----------
    ligand_fn : str,
        The input ligand file name.

    Methods
    -------

    Attributes
    ----------
    lig : a biopandas mol2 read object
    lig_data : a panda data object holding the atom information
    coordinates : np.ndarray, shape = [ N, 3]
        The coordinates of the atoms in the ligand, N is the number of atoms.

    """

    def __init__(self, ligand_fn):
        self.lig_file = ligand_fn
        self.lig = None
        self.lig_data = None

        self.lig_ele = None
        self.coordinates_ = None
        self.mol2_parsed_ = False

    def _format_convert(self, input, output):
        mol = Molecule(in_format=input.split(".")[-1])
        mol.babel_converter(input, output)
        return self

    def get_element(self):
        ele = list(self.lig_data["atom_type"].values)
        self.lig_ele = list(map(get_ligand_elementtype, ele))
        return self

    def get_coordinates(self):
        """
        Get the coordinates in the pdb file given the ligand indices.

        Returns
        -------
        self : an instance of itself

        """
        self.coordinates_ = self.lig_data[['x', 'y', 'z']].values
        return self

    def parseMol2(self):
        if not self.mol2_parsed_:
            if self.lig_file.split(".")[-1] != "mol2":
                out_file = self.lig_file + ".mol2"
                self._format_convert(self.lig_file, out_file)
                self.lig_file = out_file

            if os.path.exists(self.lig_file):
                try:
                    self.lig = PandasMol2().read_mol2(self.lig_file)
                except ValueError:
                    templ_ligfile = self.lig_file+"templ.pdb"
                    self._format_convert(self.lig_file, templ_ligfile)
                    if os.path.exists(templ_ligfile):
                        self.lig = mt.load_pdb(templ_ligfile)
                        top = self.lig.topolgy
                        table, bond = top.to_dataframe()
                        self.lig_ele = list(table['element'])
                        self.coordinates_ = self.lig.xyz[0] * 10.0
                        self.lig_data = table
                        self.lig_data['x'] = self.coordinates_[:, 0]
                        self.lig_data['y'] = self.coordinates_[:, 1]
                        self.lig_data['z'] = self.coordinates_[:, 2]
                        self.mol2_parsed_ = True
                        os.remove(templ_ligfile)
                        return self
            else:
                return None

            self.lig_data = self.lig.df
            self.get_element()
            self.get_coordinates()
            self.mol2_parsed_ = True

        return self


def get_protein_elementtype(e):
    # if e in elements_protein:
    #     return e
    # else:
    #     return "DU"
    
    if e in elements_protein:
        return e
    elif e in ['Cl', 'Br', 'I', 'F']:
        return 'HAX'
    else:
        return "DU"


def get_ligand_elementtype(e):
    if e == "C.ar":
        return "CAR"
    elif e.split(".")[0] in elements_ligand:
        return e.split(".")[0]
    elif e in ['Cl', 'Br', 'I', 'F']:
        return 'HAX'
    else:
        return "DU"

# Assuming torch is installed and a GPU is available
# We need to check if a GPU is available and move computations to GPU if it is.
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Convert the numpy array distance computation to use PyTorch for GPU acceleration
def atomic_distance_torch(dat):
    return torch.sqrt(torch.sum((dat[0] - dat[1]) ** 2))

def distance_pairs_torch(coord_pro, coord_lig):

    # 扩展 coord_pro 和 coord_lig 以生成所有可能的坐标对
    coord_pro_expanded = coord_pro.unsqueeze(1).expand(-1, coord_lig.size(0), -1)
    coord_lig_expanded = coord_lig.unsqueeze(0).expand(coord_pro.size(0), -1, -1)

    # 计算所有坐标对之间的距离
    distances = torch.sqrt(torch.sum((coord_pro_expanded - coord_lig_expanded) ** 2, dim=2))

    return distances.flatten()

def distance2counts_torch(megadata):
    d = torch.tensor(megadata[0], device=device)
    c = torch.tensor(megadata[1], device=device)
    return torch.sum((d <= c) * 1.0)

def generate_features(args):
    pro_fn, lig_fn, n_cutoffs = args
    print("INFO: Processing %s and %s ..." % (pro_fn, lig_fn))

    pro = ProteinParser(pro_fn)
    pro.parsePDB()
    protein_data = pd.DataFrame([])
    protein_data["element"] = pro.rec_ele
    # print(pro.rec_ele)
    for i, d in enumerate(['x', 'y', 'z']):
        # coordinates by mdtraj in unit nanometer
        protein_data[d] = pro.coordinates_[:, i]

    lig = LigandParser(lig_fn)
    lig.parseMol2()
    ligand_data = pd.DataFrame()
    ligand_data['element'] = lig.lig_ele
    for i, d in enumerate(['x', 'y', 'z']):
        # the coordinates in ligand are in angstrom
        ligand_data[d] = lig.coordinates_[:, i] * 0.1
        
    # 将 DataFrame 转换为 PyTorch 张量
    protein_tensor = torch.tensor(protein_data[['x', 'y', 'z']].values, device=device)
    ligand_tensor = torch.tensor(ligand_data[['x', 'y', 'z']].values, device=device)

    onionnet_counts = pd.DataFrame()

    # 使用张量运算而非显式循环
    for el in elements_ligand:
        for ep in elements_protein:
            protein_mask = (protein_data['element'] == ep).to_numpy()
            ligand_mask = (ligand_data['element'] == el).to_numpy()

            protein_xyz = protein_tensor[protein_mask]
            ligand_xyz = ligand_tensor[ligand_mask]

            counts = torch.zeros(len(n_cutoffs), device=device)

            if protein_xyz.shape[0] and ligand_xyz.shape[0]:
                distances = distance_pairs_torch(protein_xyz, ligand_xyz)

                for i, c in enumerate(n_cutoffs):
                    single_count = distance2counts_torch((distances, c))
                    if i > 0:
                        single_count -= torch.sum(counts[:i])
                    counts[i] = single_count

            feature_id = "%s_%s" % (el, ep)
            onionnet_counts[feature_id] =counts.cpu().numpy()

    # 提取第四个单词和最后一个单词来构成新的标签
    fourth_word = pro_fn.split('/')[3]
    last_word = lig_fn.split('/')[-1].split('.')[0]
    new_label = fourth_word + "_" + last_word
    
    return list(onionnet_counts.values.ravel()), new_label


def main():
    d = """
       Predicting protein-ligand binding affinities (pKa) with OnionNet model.
       
       Citation: Zheng L, Fan J, Mu Y. arXiv preprint arXiv:1906.02418, 2019.
       Author: Liangzhen Zheng (zhenglz@outlook.com)
       This script is used to generate inter-molecular element-type specific
       contact features. Installation instructions should be refered to
       https://github.com/zhenglz/onionnet-v2
       
       Examples:
       Show help information
       python generate_features.py -h
       Run the script
       python generate_features.py -inp input_samples.dat -out features_samples.csv
       # tutorial example
       cd tuttorials/PDB_samples
       python ../../generate_features.py -inp input_PDB_testing.dat -out 
       PDB_testing_features.csv
    """
    parser = argparse.ArgumentParser(description=d,
                                     formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-inp", type=str, default="input.dat",
                        help="Input. The input file containg the file path of each \n"
                             "of the protein-ligand complexes files (in pdb format.)\n"
                             "There should be only 1 column, each row or line containing\n"
                             "the input file path, relative or absolute path.")
    parser.add_argument("-out", type=str, default="output.csv",
                        help="Output. Default is output.csv \n"
                             "The output file name containing the features, each sample\n"
                             "per row. ")
    parser.add_argument("-nt", type=int, default=1,
                        help="Input, optional. Default is 1. "
                             "Use how many of cpu cores.")
    parser.add_argument("-upbound", type=float, default=3.1,
                        help="Input, optional. Default is 3.1 nm. "
                             "The largest distance cutoff.")
    parser.add_argument("-lowbound", type=float, default=0.1,
                        help="Input, optional. Default is 0.1 nm. "
                             "The lowest distance cutoff.")
    parser.add_argument("-nbins", type=int, default=60,
                        help="Input, optional. Default is 60. "
                             "The number of distance cutoffs.")

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    num_threads = args.nt

    n_cutoffs = np.linspace(args.lowbound,
                            args.upbound, args.nbins)

    with open(args.inp) as lines:
        codes = [x.split() for x in lines if
                 ("#" not in x and len(x.split()))]

    input_arguments = []
    for c in codes:
        r_fn, l_fn = c[0], c[1]
        input_arguments.append((r_fn, l_fn, n_cutoffs))

    multiprocessing.set_start_method('spawn', force=True)

    if num_threads <= 1:
        dat = []
        labels = []
        for item in input_arguments:
            ocontacts, fn = generate_features(item)
            dat.append(ocontacts)
            labels.append(fn)
    else:
        pool = multiprocessing.Pool(num_threads)
        results = pool.map(generate_features, input_arguments)

        pool.close()
        pool.join()
        dat = [x[0] for x in results]
        labels = [x[1] for x in results]

    elecombin = list(itertools.product(elements_ligand, elements_protein))
    elecombin = ["_".join(x) for x in elecombin]
    columns = []
    for i in range(n_cutoffs.shape[0]):
        columns += [e+"_"+str(i+1) for e in elecombin]

    features = pd.DataFrame(dat, index=labels, columns=columns)
    features.to_csv(args.out, header=True, index=True)

    print("INFO: Feature extraction completed.")


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()  # 记录程序结束时间
    duration = end_time - start_time  # 计算运行持续时间

    print("Program completed in {:.2f} seconds".format(duration))