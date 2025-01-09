from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D
from rdkit.Chem import rdmolops
import pandas as pd
from pandarallel import pandarallel
import tarfile
import os
from tidyscreen import general_functions as gen_func
import sqlite3
import json
import sys

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
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        return inchi_key
    except:
        pass

def compute_inchi_key_for_whole_df(df):
    pandarallel.initialize(progress_bar=False)
    df["inchi_key"] = df["SMILES"].parallel_apply(lambda smiles: compute_inchi_key(smiles))
    return df

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

def filter_ligands_table_with_smarts(conn,table_name,filters_instances_list, filters_criteria_list, filters_SMARTS_list):

    sql = f"SELECT * FROM {table_name}"
    ligands_df = pd.read_sql(sql,conn)

    counter = 1
    for index, smarts_filter in enumerate(filters_SMARTS_list):
        match_column = f"match_{str(counter)}"
        instances_column = f"instances_{str(counter)}"
        # This will assign the filter criteria to apply to each filter, i.e. 'I' (include) or 'E' (exclude)
        instances_nbr = filters_instances_list[index]
        selection_criteria = filters_criteria_list[index]
        pandarallel.initialize(progress_bar=False)
        ligands_df[[match_column,instances_column]] = ligands_df["SMILES"].parallel_apply(lambda x: pd.Series(match_smarts_filter(x,smarts_filter,instances_nbr,selection_criteria)))
        counter +=1

    ligands_df_filtered = filter_tagged_molecules_dataframe(ligands_df)
    
    print(f'A total of {len(ligands_df_filtered)} matched the filter')

    return ligands_df_filtered
        
def match_smarts_filter(smiles,smarts_filter,instances_nbr,selection_criteria):
    
    mol = Chem.MolFromSmiles(smiles)
    smarts_mol_filter = Chem.MolFromSmarts(smarts_filter)
    if mol.HasSubstructMatch(smarts_mol_filter):
        if selection_criteria == 'I':
            # Get the number of smarts matches
            smarts_matches = len(mol.GetSubstructMatches(smarts_mol_filter))
            
            # Chech is the number of matches is higher than the allowed instances. In that case exclude
            if smarts_matches > instances_nbr:
                #print(f"Matches: {smarts_matches} - Instance: {instances_nbr}")
                return 0, int(smarts_matches) # The molecule is excluded because the number of matches is above the limit
            else:
                return 1, int(smarts_matches) # First value: matchig flag; second value: nbr of instances of matches
        else:
            return 0, 0 # First value: matchig flag; second value: nbr of instances of matches
    else:
        if selection_criteria == 'I':
            return 0, 0 # First value: matchig flag; second value: nbr of instances of matches
        else:
            return 1, 0 # First value: matchig flag; second value: nbr of instances of matches
        
def filter_tagged_molecules_dataframe(ligands_df):
    # Get the column names corresponding to matching tags
    columns_with_tags = [col for col in ligands_df.columns if "match" in col or "instance" in col]

    # Filter sequentially the tags in the dataframe
    counter = 1
    for tag_column in columns_with_tags:
        if counter == 1: # This will hit for the first cleaning iteration
            ligands_df_filtered = ligands_df[ligands_df[tag_column] == 1]
            # Drop the column that has been processed
            ligands_df_filtered = ligands_df_filtered.drop(columns=[tag_column])
            counter +=1
        else:
            ligands_df_filtered = ligands_df_filtered[ligands_df_filtered[tag_column] == 1]
            # Drop the column that has been processed
            ligands_df_filtered = ligands_df_filtered.drop(columns=[tag_column])
            counter +=1
            
    return ligands_df_filtered.reset_index(drop=True)

def apply_unimolecular_reaction(molecules_df, table_name, reaction_smarts, scheme_index):
    
    # Apply in a parallelized scheme the corresponding reaction
    pandarallel.initialize(progress_bar=False)
    product_column = f"product_step_{scheme_index}"
    
    # If the reaction is the first one in the schem (i.e. index = 0), the start reactions with the given table name:
    if scheme_index == 0:
        molecules_df[f"SMILES_{product_column}"] = molecules_df[f"SMILES_{table_name}"].parallel_apply(lambda smiles:react_molecule_1_component(smiles,reaction_smarts))

    # Id the reaction is a continuation of a serial scheme, then the reactants smiles are stored in index-1 location
    else:
        # Check if the table name has been correctly identified as a unimolecules synthetic workflow based on serial reactions (i.e. '->'). If not, inform and exit.add(element)
        if table_name != "->":
            print(f"The table named '{table_name}' is not part of a synthetic pipeline involved in a unimolecular reaction step. It should be indicated as '->'.")
            sys.exit()

        # Apply the subsequent step of the reaction pipeline.
        molecules_df[f"SMILES_{product_column}"] = molecules_df[f"SMILES_product_step_{scheme_index-1}"].parallel_apply(lambda smiles:react_molecule_1_component(smiles,reaction_smarts))

    return molecules_df

def apply_bimolecular_reaction(conn,molecules_df, tables_names, reaction_smarts, scheme_index):

    # Apply in a parallelized scheme the corresponding reaction
    pandarallel.initialize(progress_bar=False)
    product_column = f"product_step_{scheme_index}"
    
    # If the reaction is the first one in the schem (i.e. index = 0), the start reactions with the given table name:
    if scheme_index == 0:
    
        reactants_1_df = pd.DataFrame(molecules_df.iloc[:, 0]).dropna()
        reactants_2_df = pd.DataFrame(molecules_df.iloc[:, 1]).dropna()
        colname_1 = f'SMILES_{tables_names[0]}'
        colname_2 = f'SMILES_{tables_names[1]}'

        # This dataframe will store all the products generated throughout the iterations
        df_all_products = pd.DataFrame()

        # In order to combinatorialy prepare the products, the loop on the dataframe should be performed on the one containing the lowest number of compounds, otherwise it will fail. So here the length of both dataframes are compared, and actioned in accordance

        if len(reactants_1_df) < len(reactants_2_df): # The first dataframe contains less reactants, so loop on it;

            for index, row in reactants_1_df.iterrows():
                current_smiles_1 = row[colname_1]

                molecules_df[f"SMILES_{product_column}"] = reactants_2_df[colname_2].parallel_apply(lambda smiles_2:react_molecule_2_component(current_smiles_1,smiles_2,reaction_smarts))

                #molecules_df["smiles1"] = reactants_2_df[colname_2].parallel_apply(lambda smiles_2:react_molecule_2_component(current_smiles_1,smiles_2,reaction_smarts))

                # Append the dataframe generated in this iteration to the overall dataframe
                df_all_products = pd.concat([df_all_products,molecules_df],axis=0)

        else: # The df2 contains less reactants that df1, so loop over if;

            for index, row in reactants_2_df.iterrows():
                current_smiles_2 = row[colname_2]

                molecules_df[f"SMILES_{product_column}"] = reactants_1_df[colname_1].parallel_apply(lambda smiles_1:react_molecule_2_component(smiles_1,current_smiles_2,reaction_smarts))

                # Append the dataframe generated in this iteration to the overall dataframe
                df_all_products = pd.concat([df_all_products,molecules_df],axis=0)

    else:
        # Check if the table name has been correctly identified in the synthetic workflow based on serial reactions (i.e. '->'). If not, inform and exit.add(element)
        if tables_names[0] != "->":
            print(f"The table named '{tables_names[0]}' is not part of a synthetic pipeline involved in a next reaction step. It should be indicated as '->'.")
            sys.exit()
        
        reactants_1_df = pd.DataFrame(molecules_df[f"SMILES_product_step_{scheme_index-1}"])
        reactants_2_df = retrieve_reactants_df(conn,tables_names[1])

        colname_1 = f"SMILES_product_step_{scheme_index-1}"
        colname_2 = f'SMILES_{tables_names[1]}'

        # This dataframe will store all the products generated throughout the iterations
        df_all_products = pd.DataFrame()

        if len(reactants_1_df) < len(reactants_2_df): # The first dataframe contains less reactants, so loop on it;

            for index, row in reactants_1_df.iterrows():
                current_smiles_1 = row.iloc[0]
                molecules_df[f"SMILES_{product_column}"] = reactants_2_df[colname_2].parallel_apply(lambda smiles_2:react_molecule_2_component(current_smiles_1,smiles_2,reaction_smarts))

                # Append the dataframe generated in this iteration to the overall dataframe
                df_all_products = pd.concat([df_all_products,molecules_df],axis=0)

        else: # The df2 contains less reactants that df1, so loop over if;
        
            for index, row in reactants_2_df.iterrows():
                current_smiles_2 = row[colname_2]

                molecules_df[f"SMILES_{product_column}"] = reactants_1_df[colname_1].parallel_apply(lambda smiles_1:react_molecule_2_component(smiles_1,current_smiles_2,reaction_smarts))

                # Append the dataframe generated in this iteration to the overall dataframe
                df_all_products = pd.concat([df_all_products,molecules_df],axis=0)

    return df_all_products # Return the overall df


def react_molecule_1_component(smiles,smarts):
    try:
        mol = Chem.MolFromSmiles(smiles)
        reaction_mol = AllChem.ReactionFromSmarts(smarts)
        product_mol = reaction_mol.RunReactants((mol,))[0][0]
        product_smiles = Chem.MolToSmiles(product_mol)
        
        return product_smiles
    except Exception as error:
        print(f"The reaction outputed the following error: \n {error} \n for compound: \n {smiles}")

def react_molecule_2_component(smiles_1,smiles_2,smarts):
    try:
        mol_1 = Chem.MolFromSmiles(smiles_1)
        mol_2 = Chem.MolFromSmiles(smiles_2)
        reaction_mol = AllChem.ReactionFromSmarts(smarts)
        product_mol = reaction_mol.RunReactants((mol_1,mol_2))[0][0]
        
        product_smiles = Chem.MolToSmiles(product_mol)

        #return product_smiles
        return smiles_1
    
    except Exception as error:
        print(f"The reaction outputed the following error: \n {error} \n for compound: \n {smiles_1}")

def retrieve_initial_reactants(conn,table_names):
    
    # Search the tables corresponding to the first reaction
    for table in table_names[0]:
        sql = f"SELECT * FROM {table}"
        molecules_df = pd.read_sql(sql,conn)
        molecules_df.rename(columns={'SMILES':f'SMILES_{table}'},inplace=True)

    return molecules_df

def retrieve_initial_reactants_updated(conn,table_names):
    # Search the tables corresponding to the first reaction
    
    molecules_df = pd.DataFrame()
    
    for table in table_names[0]:
        sql = f"SELECT SMILES FROM {table}"
        current_df = pd.read_sql(sql,conn)
        current_df.rename(columns={'SMILES':f'SMILES_{table}'},inplace=True)
        molecules_df = pd.concat([molecules_df,current_df],axis=1)
        
    return molecules_df

def retrieve_reactants_df(conn,table_name):
    sql = f"SELECT SMILES FROM {table_name}"
    molecules_df = pd.read_sql(sql,conn)
    molecules_df.rename(columns={'SMILES':f'SMILES_{table_name}'},inplace=True)

    return molecules_df

def retrieve_initial_reactants_as_lists(conn,table_names):
    # This will generate a list lists corresponding to the reactants tables involved in the reaction scheme
    reactants_lists_of_lists = []
    for table in table_names[0]: # Will get the table/s associated to the first reaction step
        sql = f"SELECT SMILES FROM {table}"
        molecules_df = pd.read_sql(sql,conn)
        currrent_list = molecules_df["SMILES"].to_list()
        reactants_lists_of_lists.append(currrent_list)

    return reactants_lists_of_lists

def depict_ligands_table(conn,table_name,output_path,max_mols_ppage):
    sql=f"SELECT SMILES, inchi_key, Name FROM {table_name};"
    molecules_df = pd.read_sql_query(sql,conn)

    # Since the molecules might be processed in chunks, it is usefull to create lists with the information in order to process it.add
    smiles_list, inchi_key_list, names_list = molecules_df["SMILES"].to_list(), molecules_df["inchi_key"].to_list(),molecules_df["Name"].to_list()

    # Generate a naming list
    lengeds_list = gen_func.combine_strings_in_lists([inchi_key_list,names_list])


    # Create the chunks of objects to process according to the max_mols_ppage parameter provided
    smiles_list_chunks = gen_func.split_list_in_chunks(smiles_list, max_mols_ppage)
    lengeds_list_chunks = gen_func.split_list_in_chunks(lengeds_list, max_mols_ppage)

    counter = 0 # used to number figures with chunks
    for index, list_item in enumerate(smiles_list_chunks):
        df_chunk  = pd.DataFrame({'SMILES':list_item})
        mols_item = [Chem.MolFromSmiles(smile) for smile in list_item]
        img = MolsToGridImage(mols=mols_item, legends=lengeds_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
        img.save(f"{output_path}/{table_name}_{counter}.png")
        counter+=1


def sanitize_smiles(smiles,retain_stereo):
    try: 
        repaired_smiles = repair_smiles_errors(smiles)
        largest_mol_smiles = get_largest_fragment(repaired_smiles)
        # Add here potential processing steps for cleaning/sanitizing the SMILES notation 
        #####
        #
        
        if retain_stereo == 0:
            clean_smiles_free_stereo = remove_smiles_stereo(largest_mol_smiles)
            return clean_smiles_free_stereo
        else:
            return largest_mol_smiles
        
    except:
        pass
    

def repair_smiles_errors(smiles):
    
    if "[nH]" in smiles:
        sanitized_smiles = smiles.replace("[nH]","[NH]")
        return repaired_smiles
    else:
        return smiles
    
def get_largest_fragment(smiles):   
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_frags = rdmolops.GetMolFrags(mol, asMols = True)
        largest_mol = max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())
        largest_mol_smiles = Chem.MolToSmiles(largest_mol)    
        return largest_mol_smiles
    except:
        return "Error in largest fragment computation"
    
def remove_smiles_stereo(smiles):
    smiles_free_stereo = smiles.replace("@","")
    if "[CH]" in smiles_free_stereo:
        smiles_free_stereo = smiles_free_stereo.replace("[CH]","C")
    return smiles_free_stereo