import shutil
import pandas as pd
from tidyscreen import database_interactions as db_ints
import sys
import ast
import json
import csv
import os
from tidyscreen import ligands_processing

def copy_file_to_dest(file,destination): 
    shutil.copy(file,destination)

def create_blob_object(filename):
    """
    This function will create a binary object from a filename provided that is aimed to be stored as a BLOB object within the database.
    ------
    Parameters:
    ------
    - filename: the filename in the corresponding format to be stored in the table containing BLOB objects.
    ------
    Returns:
    ------
    A blob object containing the ligand .pdbqt file in a compressed (.tar.gz) format.
    """
    with open(filename, 'rb') as file:
        blobData = file.read()
    file.close()
    
    return blobData

def input_chem_filters_to_db(conn,filters_file):
    # Input the specified file
    df_from_input_file = pd.read_csv(filters_file)

    # Check if the  filters database table is present
    chem_filters_present = db_ints.check_table_presence(conn,"chem_filters")

    # If the 'chem_filters' table is not present, create it prior to storing the corresponding filters
    if chem_filters_present == 0:
        print("Creating the 'chem_filters' table.")
        db_ints.create_chemfilters_table(conn)
    
    try:
        # Check if the columns are consistent with the required format
        if df_from_input_file.columns.to_list() == ['Filter_Name', 'SMARTS', 'Description']:
            print("Column name in the input file are ok")
            stored_filters_df = retrieve_stored_filters(conn)
            # Generate a list of stored filters
            stored_filters_list = stored_filters_df["Filter_Name"].to_list()
            
            # Evaluate every row of the input df to check if the filter is new (based on the name) or if it already exists
            for index, row in df_from_input_file.iterrows():
                # Check if the filter name is new
                current_filter_name = row["Filter_Name"]
                if current_filter_name not in stored_filters_list:
                    # Store the new row in the original df
                    stored_filters_df.loc[len(stored_filters_df)] = row
                else:
                    print(f"The filter '{current_filter_name}' is already in the database. Passing.")
                    
            # Store the updated df in the main database
            stored_filters_df["Filter_id"] = stored_filters_df.index + 1
            db_ints.store_df_as_table_in_db(conn,stored_filters_df,"chem_filters")
            
        else:
            print("columns are bad")
    
    except:
        print("Error storing the chem_filters")

def input_single_reactions_to_db(conn,reactions_file):
    # Input the specified file
    df_from_input_file = pd.read_csv(reactions_file)
    
    # Check if the single reactions table is present in TidyScreen main database
    single_reactions_table_present = db_ints.check_table_presence(conn,"single_reactions")
    
    # If the 'single_reactions' table is not present, create it prior to storing the corresponding reactions
    if single_reactions_table_present == 0:
        print("Creating the 'single_reactions' table.")
        db_ints.create_single_reactions_table(conn)
    
    try:
        # Check if the columns are consistent with the required format
        if df_from_input_file.columns.to_list() == ['Reaction_Name', 'Reaction_SMARTS', 'Reaction_creator']:
            print("Column name in the input file are ok")
            stored_single_reactions_df = retrieve_stored_reactions(conn)
            # Generate a list of stored filters
            single_reactions_names_list = stored_single_reactions_df["Reaction_Name"].to_list()
            
            # Evaluate every row of the input df to check if the filter is new (based on the name) or if it already exists
            for index, row in df_from_input_file.iterrows():
                # Check if the filter name is new
                current_reaction_name = row["Reaction_Name"]
                if current_reaction_name not in single_reactions_names_list:
                    # Store the new row in the original df
                    stored_single_reactions_df.loc[len(stored_single_reactions_df)] = row
                else:
                    print(f"The filter '{current_reaction_name}' is already in the database. Passing.")
                    
            # Store the updated df in the main database
            stored_single_reactions_df["Reaction_id"] = stored_single_reactions_df.index + 1
            db_ints.store_df_as_table_in_db(conn,stored_single_reactions_df,"single_reactions")
            
        else:
            print("columns are bad")
            
    except Exception as error:
        print(error)

def retrieve_stored_filters(conn):
    
    chem_filters_present = db_ints.check_table_presence(conn,"chem_filters")
    if chem_filters_present == 1:
        df = pd.read_sql_query(f"SELECT * FROM chem_filters", conn)
        return df
    else:
        print("No 'chem_filters' table was found in the main database.")
        sys.exit()
        
def retrieve_stored_reactions(conn):
    
    single_reactions_present = db_ints.check_table_presence(conn,"single_reactions")
    if single_reactions_present == 1:
        df = pd.read_sql_query(f"SELECT * FROM single_reactions", conn)
        return df
    else:
        print("No 'single_reactions' table was found in the main database.")
        sys.exit()

def store_filters_pipeline(project_db_conn,main_database_conn,filters_number_list,filters_instances_list,filters_criteria_list,filters_pipeline_name):
    
    existing_filters_df = pd.read_sql_query(f"SELECT * FROM filters_pipelines", project_db_conn)

    # get number of rows in the 'existing_filters_df'
    num_filters_pipeline = existing_filters_df.shape[0]
    
    if num_filters_pipeline == 0:
        last_filters_pipeline = 0
    else:
        # get the last filter id:
        last_filters_pipeline_row = existing_filters_df.iloc[-1]
        last_filters_pipeline = last_filters_pipeline_row["Filters_pipeline_id"]

    # Convert the corresponding lists to a json object to compatibilize it with storage in a SQL database

    filters_number_json = json.dumps(filters_number_list)
    filters_instances_json = json.dumps(filters_instances_list)
    filters_criteria_json = json.dumps(filters_criteria_list)

    # If the 'filters_pipeline_list' already exists in the database, inform and exit
    for index, row in existing_filters_df.iterrows():
        if filters_number_json == row["Filters_id"]:
            print(f"The filter pipeline: '{filters_number_list}' already exists in the filters pipeline table. Exit.")
            sys.exit()

    # If the above check does not exist as a consequence of the filter_list already existing, the storing of the filter pipeline will continue with the execution of the code below.
    
    try:
        current_filters_pipeline_id = last_filters_pipeline + 1
    except Exception as error:
        print(error)
    
    filters_smarts_json = create_filter_smarts_json(main_database_conn,filters_number_json)

    # Append the data the last row of 'existing_filters_df'
    existing_filters_df.loc[len(existing_filters_df)] = [current_filters_pipeline_id, filters_pipeline_name, filters_number_json, filters_instances_json, filters_criteria_json, filters_smarts_json]

    db_ints.store_df_as_table_in_db(project_db_conn,existing_filters_df,"filters_pipelines")

def store_reactions_scheme(project_db_conn,main_database_conn,reactions_id_list,reaction_scheme_name):
    
    existing_reaction_schemes_df = pd.read_sql_query(f"SELECT * FROM reaction_schemes", project_db_conn)
    
    # get number of rows in the 'existing_filters_df'
    num_reaction_schemes = existing_reaction_schemes_df.shape[0]
    
    if num_reaction_schemes == 0:
        last_reaction_scheme_nbr = 0
        
    else:
        # get the last filter id:
        last_reaction_scheme_row = existing_reaction_schemes_df.iloc[-1]
        last_reaction_scheme_nbr = last_reaction_scheme_row["Reaction_scheme_id"]
        
    # Convert the corresponding lists to a json object to compatibilize it with storage in a SQL database
    
    reactions_scheme_id_json = json.dumps(reactions_id_list)
    
    # If the 'reactions_id_list' already exists in the database, inform and exit
    
    for index, row in existing_reaction_schemes_df.iterrows():
        if reactions_scheme_id_json == row["Reactions_id"]:
            print(f"The reaction scheme: '{reactions_id_list}' already exists in the reactions schemes table. Exit.")
            sys.exit()
    
    # If the above check does not exist as a consequence of the reaction scheme already existing, the storing of the reaction scheme will continue with the execution of the code below.
    
    try:
        current_reaction_scheme_nbr = last_reaction_scheme_nbr + 1
    except Exception as error:
        print(error)
    
    reaction_scheme_smarts_json = create_reaction_scheme_smarts_json(main_database_conn,reactions_scheme_id_json)

    # Append the data the last row of 'existing_reaction_schemes_df'
    existing_reaction_schemes_df.loc[len(existing_reaction_schemes_df)] = [current_reaction_scheme_nbr, reaction_scheme_name, reactions_scheme_id_json, reaction_scheme_smarts_json]

    db_ints.store_df_as_table_in_db(project_db_conn,existing_reaction_schemes_df,"reaction_schemes")
    
def create_filter_smarts_json(main_database_conn,filters_number_json):

    cursor = main_database_conn.cursor()
    
    filters_id_list = json.loads(filters_number_json)
    
    filters_smarts_list = []

    for filter_nbr in filters_id_list:
        cursor.execute(f"""SELECT SMARTS FROM chem_filters WHERE Filter_id == {filter_nbr};""")
        current_smarts = cursor.fetchone()[0]
        filters_smarts_list.append(current_smarts)

    # This will dump the list as a json object that is storable within the SQL database. Requires a deserialization when retrieving to convert it back to a list.
    filters_smarts_json = json.dumps(filters_smarts_list)

    return filters_smarts_json

def create_reaction_scheme_smarts_json(main_database_conn,reactions_scheme_json):
    
    cursor = main_database_conn.cursor()
    
    reactions_id_list = json.loads(reactions_scheme_json)
    
    reaction_scheme_smarts_list = []
    
    for reaction_nbr in reactions_id_list:
        cursor.execute(f"""SELECT Reaction_SMARTS FROM single_reactions WHERE Reaction_id == {reaction_nbr};""")
        current_smarts = cursor.fetchone()[0]
        reaction_scheme_smarts_list.append(current_smarts)
        
    # This will dump the list as a json object that is storable within the SQL database. Requires a deserialization when retrieving to convert it back to a list.
    reaction_scheme_smarts_json = json.dumps(reaction_scheme_smarts_list)
    
    return reaction_scheme_smarts_json

def create_json_object_from_lists(keys, values):
    
    dict_for_json = {}
    for index,key in enumerate(keys):
        dict_for_json[key] = values[index]
        #print(key)
        #print(values[index])
        
    return dict_for_json

def check_id_smarts_pair_in_main_filters_db(main_existing_filters_df,filters_number_list,filters_SMARTS_list):
    
    # Generate the lists of stored filters numbers and SMARTS fitlers
    main_db_filters_number_list = main_existing_filters_df["Filter_id"].to_list()
    main_db_filters_SMARTS_list = main_existing_filters_df["SMARTS"].to_list()
    
    main_db_stored_tuples_list = generate_filter_nb_filter_smarts_list_of_tuples(main_db_filters_number_list,main_db_filters_SMARTS_list)
    provided_tuple_list = generate_filter_nb_filter_smarts_list_of_tuples(filters_number_list,filters_SMARTS_list)
    
    for item in provided_tuple_list:
        if item not in main_db_stored_tuples_list:
            print(f"The filter combination {item} is not present in the main database")
            sys.exit()
    
def generate_filter_nb_filter_smarts_list_of_tuples(filters_number_list,filters_SMARTS_list):
    
    # Generate tuples containing two elements: [0]: filter_number and [1]: filter SMARTS. The tuple will be checked for existence in the main database filters
    list_of_provided_tuples = []
    for index, filter_number in enumerate(filters_number_list):
        current_tuple = (filter_number,filters_SMARTS_list[index]) 
        list_of_provided_tuples.append(current_tuple)
    
    return list_of_provided_tuples

def evaluate_smarts_reaction(smarts):

    molecularity_nbr = len(smarts.split('>>')[0].split('.'))
    
    return molecularity_nbr

def check_tables_vs_smarts(table_names,smarts):
    number_of_tables = len(table_names)
    molecularity_nbr = len(smarts.split('>>')[0].split('.'))

    if number_of_tables != molecularity_nbr:
        print(f"The number of tables ({table_names} is not consistent with the reaction provided {smarts})")
        sys.exit()

def clean_reaction_scheme_product_df(df,reaction_scheme_id):
    # Get the column names corresponding to the column to mantain
    df = df[[df.columns[-1]]]
    df = df.rename(columns={df.columns[0]: 'SMILES'})
    df = ligands_processing.compute_inchi_key_for_whole_df(df)
    df["Name"] = f"Obtained_from_react_scheme_{reaction_scheme_id}"


    return df

def report_on_reaction_scheme(df):
    print("Count of resulting products in the synthetic scheme")
    non_null_unique_counts = df.nunique()
    print(non_null_unique_counts)

def check_folder_existence(folder):
    if os.path.exists(folder) and os.path.isdir(folder):
        return 1
    else:
        return 0

def create_folder(folder):
    os.mkdir(folder)
    print(f"The folder: \n {folder} \n was succesfully created.")

def split_list_in_chunks(list,nbr_items):
    """"
    This function will receive a list of items, and will return a list of lists each one containing the corresponding 'nbr_items'
    ------
    Parameters
    ------
    -list: input list to be processed in chunks
    -nbr_item: number of chunks included in each component list
    """
    list_of_chunks = [list[i:i+nbr_items] for i in range(0, len(list), nbr_items)]

    return list_of_chunks

def combine_strings_in_lists(list_of_lists):

    # Take the first list as reference
    reference_list = list_of_lists[0]

    # Check if some list in the list of list is not equal in length. In that case, inform and exit
    for item in list_of_lists:
        if len(reference_list) != len(item):
            print(f"The list: \n '{item}' \n is not equal in length to the rest of lists provided. Stopping.")
            sys.exit()

    # Append the string values to the resulting list
    appended_list = [' \n '.join(strings) for strings in zip(*list_of_lists)]

    return appended_list