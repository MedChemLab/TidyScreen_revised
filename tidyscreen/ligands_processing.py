from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import pandas as pd
from pandarallel import pandarallel
import tarfile
import os
from tidyscreen import general_functions as gen_func
import sqlite3

def process_ligand_table_to_pdbqt(conn,table_name,pdbqt_table,chemspace_db):
    sql = f"""SELECT SMILES, Name FROM '{table_name}';"""
    df = pd.read_sql_query(sql,conn)
    pandarallel.initialize(progress_bar=True)
    df.parallel_apply(lambda row: append_ligand_pdqbt_blob_object_to_table(row,pdbqt_table,chemspace_db), axis=1)

def append_ligand_pdqbt_blob_object_to_table(row,pdqbt_table,chemspace_db):
    pdb_tar_file, pdbqt_tar_file = compute_pdqbt_file(row["SMILES"])
    pdb_tar_blob_obj = gen_func.create_blob_object(pdb_tar_file)
    pdbqt_tar_blob_obj = gen_func.create_blob_object(pdbqt_tar_file)
    inchi_key = compute_inchi_key(row['SMILES'])
    
    try:
        conn = sqlite3.connect(chemspace_db)
        cursor = conn.cursor()
        sql = f"""INSERT INTO '{pdqbt_table}' (SMILES,Name,inchi_key,stored_pdbqt_file, pdbqt_blob_item, stored_pdb_file,pdb_blob_item)
              VALUES (?,?,?,?,?,?,?);"""
        data_tuple = (row["SMILES"],row["Name"],inchi_key,inchi_key+".pdbqt.tar.gz",pdbqt_tar_blob_obj,inchi_key+".pdb.tar.gz",pdb_tar_blob_obj)
        cursor.execute(sql,data_tuple)
        conn.commit()
        cursor.close()

        # Clean /tmp files
        os.remove(pdb_tar_file)
        os.remove(pdbqt_tar_file)
        
    except:
        print("ERROR")
    
def compute_pdqbt_file(smiles):
    mol_hs = mol_from_smiles(smiles)
    confs, mol_hs, ps = generate_ligand_conformers(mol_hs)
    selected_mol = get_conformer_rank(confs, mol_hs, ps)
    pdb_file, inchi_key = save_ligand_pdb_file(selected_mol)    
    atoms_dict = create_meeko_atoms_dict()
    pdb_tar_file, pdbqt_tar_file = mol_meeko_processing(pdb_file, inchi_key,atoms_dict)
    return pdb_tar_file, pdbqt_tar_file

def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol_hs = Chem.AddHs(mol)
    return mol_hs

def generate_ligand_conformers(mol_hs,nbr_confs=50,mmff='MMFF94',maxIters=10, conformer_rank=0):
    props = AllChem.ETKDG()
    props.pruneRmsThresh = 0.25
    props.useRandomCoords = True
    props.numThreads = 1
    
    # Create the object containing the conformers
    confs = AllChem.EmbedMultipleConfs(mol_hs,nbr_confs,props)
    ps = AllChem.MMFFGetMoleculeProperties(mol_hs)
    AllChem.MMFFOptimizeMoleculeConfs(mol_hs, mmffVariant=mmff, maxIters=maxIters)
    
    return confs, mol_hs, ps

def get_conformer_rank(confs, mol_hs, ps,rank=0):
    
    conformers_energies_dict={} # Empty dictionary to store conformers

    for conf in confs:
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol_hs,ps,confId=conf)
            ff.Minimize()
            energy_value = ff.CalcEnergy()
            conformers_energies_dict[conf] = energy_value
        except:
            continue

    # Sort the dictionary 
    conformers_energies_dict_sorted=sorted(conformers_energies_dict.items(), key=lambda x: x[1])

    # The following will store the lowest energy conformer as .pdbqt
    selected_conformer = Chem.MolToMolBlock(mol_hs,confId=conformers_energies_dict_sorted[rank][0])
    selected_mol = Chem.MolFromMolBlock(selected_conformer, removeHs=False)

    return selected_mol

def save_ligand_pdb_file(selected_mol):
    inchi_key = Chem.MolToInchiKey(selected_mol)
    pdb_file = f'/tmp/{inchi_key}.pdb'
    Chem.MolToPDBFile(selected_mol,pdb_file)
    return pdb_file, inchi_key

def compute_inchi_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    inchi_key = Chem.MolToInchiKey(mol)
    return inchi_key

def mol_meeko_processing(pdb_file, inchi_key,atoms_dict):
    mol = Chem.MolFromPDBFile(f'/tmp/{inchi_key}.pdb',removeHs=False)
    mk_prep = MoleculePreparation(add_atom_types=atoms_dict)
    mol_setup_list = mk_prep(mol)
    molsetup = mol_setup_list[0]

    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)
    
    with open(f'/tmp/{inchi_key}.pdbqt','w') as pdbqt_file:
        pdbqt_file.write(pdbqt_string[0])
        
    # This will output the .pdbqt file 
    pdbqt_file = f'/tmp/{inchi_key}.pdbqt'
    blob_pdbqt_file = f'{pdbqt_file}.tar.gz'
    tar = tarfile.open(blob_pdbqt_file, "w:gz")
    tar.add(pdbqt_file,arcname=f"{inchi_key}.pdbqt")
    tar.close()
    os.remove(pdbqt_file)
    
    # This will output the .pdb file 
    pdb_file = f'/tmp/{inchi_key}.pdb'
    blob_pdb_file = f'{pdb_file}.tar.gz'
    tar = tarfile.open(blob_pdb_file, "w:gz")
    tar.add(pdb_file,arcname=f"{inchi_key}.pdb")
    tar.close()    
    os.remove(pdb_file)

    return blob_pdb_file, blob_pdbqt_file


def create_meeko_atoms_dict():
    atoms_dict = [{"smarts": "[#1]", "atype": "H",},
    {"smarts": "[#1][#7,#8,#9,#15,#16]","atype": "HD"},
    {"smarts": "[#5]", "atype": "B"},
    {"smarts": "[C]", "atype": "C"},
    {"smarts": "[c]", "atype": "A"},
    {"smarts": "[#7]", "atype": "NA"},
    {"smarts": "[#8]", "atype": "OA"},
    {"smarts": "[#9]", "atype": "F"},
    {"smarts": "[#12]", "atype": "Mg"},
    {"smarts": "[#14]", "atype": "Si"},
    {"smarts": "[#15]", "atype": "P"},
    {"smarts": "[#16]", "atype": "S"},
    {"smarts": "[#17]", "atype": "Cl"},
    {"smarts": "[#20]", "atype": "Ca"},
    {"smarts": "[#25]", "atype": "Mn"},
    {"smarts": "[#26]", "atype": "Fe"},
    {"smarts": "[#30]", "atype": "Zn"},
    {"smarts": "[#35]", "atype": "Br"},
    {"smarts": "[#53]", "atype": "I"},
    {"smarts": "[#7X3v3][a]", "atype": "N", "comment": "pyrrole, aniline"},
    {"smarts": "[#7X3v3][#6X3v4]", "atype": "N", "comment": "amide"},
    {"smarts": "[#7+1]", "atype": "N", "comment": "ammonium, pyridinium"},
    {"smarts": "[SX2]", "atype": "SA", "comment": "sulfur acceptor"},
    {"smarts": "[#1][#6X3]:[#6X3]([#6X4])[#7]:[#7][#7][#6X4]", "atype": "HD", "comment": "4,5-H in 1,2,3-triazole"},
    ]

    return atoms_dict

def enumerate_stereoisomers_info(smiles,name):
    """
    This function will compute the possible stereoisomers combinations given a SMILES string.

    The function will calculate and return a dictionary, containing ordered lists with the corresponding stereoisomer SMILES, overal stereoconfig and number of total stereoisomers for the specific compound.

    ------
    Parameters:
    ------
    - smiles: a SMILES string that is the input to the computation.
    """
    
    try:
        mol_hs = mol_from_smiles(smiles)
        isomers = list(EnumerateStereoisomers(mol_hs))
        isomers_list = [] # The stereoisomer SMILES will be appended here.
        isomers_name_list = [] # The stereoisomer SMILES will be appended here.
        stereo_config_list = [] # The stereo config will be appended here.
        # Compute each stereoisomer property
        for isomer in isomers:
            isomer_no_hs = Chem.RemoveHs(isomer)
            isomer_smiles = Chem.MolToSmiles(isomer_no_hs)
            stereo_config = Chem.FindMolChiralCenters(isomer,force=True,includeUnassigned=True,useLegacyImplementation=True)
            isomers_list.append(isomer_smiles)
            stereo_config_list.append(stereo_config)
            isomers_name_list.append(name)
   
        number_of_isomers = len(isomers_list)
        number_of_optical_isomers = count_chiral_centers(smiles)
        number_of_isomers_list = []

        # This will create a list with n-times the number of isomers for further exploding the list in a dataframe.
        for nbr in range(number_of_isomers):
            number_of_isomers_list.append(number_of_optical_isomers)

        # This will clean the stereo_config_list_1 to retain only the config in 'SS' (or whatever config) formate.
        stereo_config_list_2 = []
    
        for isomer_nbr in range(number_of_isomers):
            stereo_config = ''        
            chiral_centers_nbr = len(stereo_config_list[isomer_nbr])
            for chiral_center_nbr in range(chiral_centers_nbr):
                center_config = (stereo_config_list[isomer_nbr][chiral_center_nbr][1])    
                stereo_config = stereo_config + center_config

                # Construct the final config list.
                stereo_config_list_2.append(stereo_config)

        results_dict = {'isomers_smiles':isomers_list,'isomers_name':isomers_name_list ,'isomers_config':stereo_config_list_2,'chiral_ctrs_nbrs':number_of_isomers_list}

        return results_dict
    
    except: # Will proceed to the next function call 
        pass
    
def count_chiral_centers(smiles):
    """
    This function will count and return the number of chiral centers for a given SMILES string:

    ------
    Parameters:
    ------
    - smiles: the smiles strings of the molecule to be processed
    
    ------
    Returns:
    ------
    An integer value corresponding to the number of stereocenters.
    """

    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_hs = Chem.AddHs(mol)
        # This will compute all stereoisomers: optical plus spatial
        isomers = list(EnumerateStereoisomers(mol_hs))

        # This will isolate the count of optical isomers:
        for isomer in isomers:
            stereo_config = Chem.FindMolChiralCenters(isomer,force=True,includeUnassigned=True,useLegacyImplementation=True)
            chiral_centers_nbr = len(stereo_config) 

        return chiral_centers_nbr
    
    except Exception as error:
        print(error)
        
def explode_stereo_column(df):
    
    # Check if the 'stereo_info' column exists in the dataframe, and in that case continue with the data explosion
    if "stereo_info" in df.columns:
        df_stereo = df['stereo_info'].apply(pd.Series)
        df_stereo_2 = df_stereo.explode(['isomers_smiles','isomers_name', 'isomers_config', 'chiral_ctrs_nbrs'],).reset_index().rename(columns={"isomers_smiles":"SMILES",'isomers_name':"Name","isomers_config":"iso_config","chiral_ctrs_nbrs":"chrl_ctrs_nbrs"})
        df_stereo_2['inchi_key'] = df_stereo_2['SMILES'].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
        
        df_stereo_2 = df_stereo_2[["SMILES","Name","inchi_key","iso_config","chrl_ctrs_nbrs"]]

        # Drop duplicated columns in case some kind of enumeration bug happens on the SMILES column
        df_stereo_2.drop_duplicates(subset=['inchi_key','iso_config','chrl_ctrs_nbrs'],inplace=True)
        
        # Create a mol_id index and make it the first column
        df_stereo_2["mol_id"] = df_stereo_2.index
        first_column = df_stereo_2.pop('mol_id')
        df_stereo_2.insert(0, 'mol_id', first_column)
        
        return df_stereo_2
        
    else:
        print("No column named 'stereo_info' located in the dataframe to be processed")