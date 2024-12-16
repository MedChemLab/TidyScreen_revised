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
    sql = f"""CREATE TABLE '{table_name}' 
            (SMILES TEXT NOT NULL, 
            inchi_key TEXT NOT NULL, 
            stored_file TEXT NOT NULL, 
            blob_item BLOB NOT NULL,
            stored_pdb_file TEXT,
            pdb_blob_item BLOB);"""      
        
    conn.execute(sql)
    
    conn.close()

if __name__ == '__main__':
    pass



