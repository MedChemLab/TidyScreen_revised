import shutil
import pandas as pd
from tidyscreen import database_interactions as db_ints
import sys

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
            print("Columns are ok")
            stored_filters_df = retrieve_stored_filters(conn)
            # Generate a list of stored filters
            stored_filters_list = stored_filters_df["Filter_Name"].to_list()
            
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
            

def retrieve_stored_filters(conn):
    
    chem_filters_present = db_ints.check_table_presence(conn,"chem_filters")
    if chem_filters_present == 1:
        df = pd.read_sql_query(f"SELECT * FROM chem_filters", conn)
        return df
    else:
        print("No 'chem_filters' table was found in the main database.")
        sys.exit()