import sqlite3
from distutils.sysconfig import get_python_lib
from termcolor import colored
import pandas as pd
import json
import sys

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

def create_single_reactions_table(conn,table_name="single_reactions"):
    try:
        sql = f"""CREATE TABLE '{table_name}' 
                (Reaction_id INTEGER, 
                Reaction_Name TEXT NOT NULL, 
                Reaction_SMARTS TEXT NOT NULL,
                Reaction_creator TEXT NOT NULL);"""      
            
        conn.execute(sql)
        print(f"Table: '{table_name}' has been created.")
        
    except sqlite3.Error as error:
        print(error)

def create_filters_pipelines_table(conn,table_name="filters_pipelines"):
    try:
        sql = f"""CREATE TABLE '{table_name}' 
                (Filters_pipeline_id INTEGER, 
                Filters_pipeline_name TEXT NOT NULL, 
                Filters_id TEXT NOT NULL,
                Filters_instances TEXT NOT NULL,
                Filters_criteria TEXT NOT NULL,
                Filters_SMARTS TEXT NOT NULL);"""      
            
        conn.execute(sql)

        print(f"Table: '{table_name}' has been created.")
    except sqlite3.Error as error:
        print(error)

def create_reaction_schemes_table(conn,table_name="reaction_schemes"):
    try:
        sql = f"""CREATE TABLE '{table_name}' 
                (Reaction_scheme_id INTEGER, 
                Reaction_scheme_name TEXT NOT NULL, 
                Reactions_id TEXT NOT NULL,
                Reactions_SMARTS TEXT NOT NULL);"""      
            
        conn.execute(sql)

        print(f"Table: '{table_name}' has been created.")
    except sqlite3.Error as error:
        print(error)

def store_df_as_table_in_db(conn,df,table_name,action="replace"):
    
    table_exists = check_table_presence(conn,table_name)

    if table_exists == 1:
        replace_table_action = input(f"The table named '{table_name}' already exists in the database. I will replace it, are you ok with that? (y/n): ")

        if replace_table_action == 'y':
            pass
        else:
            print(f"Quiting to safely retain table {table_name}")
    
    try:
        df.to_sql(con=conn, name=table_name,if_exists=action,index=None)
        
    except Exception as error:
        print(error)
    
def retrieve_filters_pipeline_smarts(conn,filters_pipeline_id):
    
    cursor = conn.cursor()
    
    sql_query = f"""SELECT Filters_pipeline_name, Filters_id,Filters_instances, Filters_criteria, Filters_SMARTS FROM filters_pipelines WHERE Filters_pipeline_id == {filters_pipeline_id};"""
    cursor.execute(sql_query)
    
    # Check if the cursor is empty and inform the error
    if cursor == None:
        print(f"An error occurred with the selection of filters pipeline with id {filters_pipeline_id}")
        sys.exit()

    # Extract the info from the cursor and split it into variables
    cursor_contained_data = cursor.fetchall()[0]
    # Separate retrieved items
    stored_filters_pipeline_name = cursor_contained_data[0]
    stored_filters_id_json = cursor_contained_data[1] 
    stored_instances_json = cursor_contained_data[2] 
    stored_criteria_json = cursor_contained_data[3]
    stored_smarts_json = cursor_contained_data[4]
    
    filters_id_list = json.loads(stored_filters_id_json)
    filters_instances_list = json.loads(stored_instances_json)
    filters_criteria_list = json.loads(stored_criteria_json)
    filters_SMARTS_list = json.loads(stored_smarts_json)
    
    return stored_filters_pipeline_name, filters_id_list, filters_instances_list, filters_criteria_list, filters_SMARTS_list

def retrieve_reaction_scheme_smarts(conn,reaction_scheme_id):
    
    cursor = conn.cursor()
    
    sql_query = f"""SELECT Reaction_scheme_name, Reactions_id,Reactions_SMARTS FROM reaction_schemes WHERE Reaction_scheme_id == {reaction_scheme_id};"""
    
    cursor.execute(sql_query)
    
    # Check if the cursor is empty and inform the error
    if cursor == None:
        print(f"An error occurred with the selection of filters pipeline with id {reaction_scheme_id}")
        sys.exit()

    # Extract the info from the cursor and split it into variables
    cursor_contained_data = cursor.fetchall()[0]
    # Separate retrieved items
    stored_reaction_scheme_name = cursor_contained_data[0]
    stored_reactions_id_json = cursor_contained_data[1] 
    stored_reactions_SMARTS_json = cursor_contained_data[2] 

    return stored_reaction_scheme_name, stored_reactions_id_json, stored_reactions_SMARTS_json

def delete_table_from_db(conn,table_name):
    cursor = conn.cursor()
    table_name = table_name
    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
    conn.commit()
    conn.close()

    print(f"The table: \n '{table_name} \n was deleted from the CS database.")