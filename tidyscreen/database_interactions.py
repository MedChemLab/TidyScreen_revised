import sqlite3
from distutils.sysconfig import get_python_lib
from termcolor import colored

def connect_to_database(database_full_path):
    conn = sqlite3.connect(f'{database_full_path}')
    return conn

def check_table_presence(conn,table_name):
    "Will return 1 if exists, otherwise returns 0"
    cursor = conn.cursor()

    cursor.execute("""SELECT name FROM sqlite_master WHERE type='table' AND name=?;""", (table_name,))

    table_exists = cursor.fetchone()

    if table_exists:
        return 1
    else:
        return 0

def create_ligand_pdbqt_table(conn,table_name):
    try:
        sql = f"""CREATE TABLE '{table_name}' 
                (SMILES TEXT NOT NULL, 
                Name TEXT NOT NULL,
                inchi_key TEXT NOT NULL, 
                stored_pdbqt_file TEXT NOT NULL, 
                pdbqt_blob_item BLOB NOT NULL,
                stored_pdb_file TEXT,
                pdb_blob_item BLOB);"""      
            
        conn.execute(sql)
    
    except sqlite3.Error as error:
        print(error)

def create_chemfilters_table(conn,table_name="chem_filters"):
    try:
        sql = f"""CREATE TABLE '{table_name}' 
                (Filter_id INTEGER, 
                Filter_Name TEXT NOT NULL, 
                SMARTS TEXT NOT NULL,
                Description TEXT NOT NULL);"""      
            
        conn.execute(sql)

        print(f"Table: '{table_name}' has been created.")
    except sqlite3.Error as error:
        print(error)
    
def store_df_as_table_in_db(conn,df,table_name,action="replace"):
    df.to_sql(con=conn, name=table_name,if_exists=action,index=None)
    
    
