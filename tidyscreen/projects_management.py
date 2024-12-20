import os
import sqlite3
from termcolor import colored
#from distutils.sysconfig import get_python_lib
import shutil
from tidyscreen import database_interactions as db_ints
import pickle
from tidyscreen import datareader
from tidyscreen import general_functions
from tidyscreen import database_interactions as db_ints
from tidyscreen import ligands_processing
import sys
import tidyscreen_dbs
import pandas as pd
from pandarallel import pandarallel


def list_projects(list_rows=1,main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')):
        
        # Check if the database exists
        if os.path.exists(main_projects_db):
            try: 
                conn = db_ints.connect_to_database(main_projects_db)
                cursor = conn.cursor()
                cursor.execute("SELECT * FROM projects")
                
                rows = cursor.fetchall()
                
                # If the default value = 1, the rows with projects will be printed
                if list_rows == 1:
                    
                    if len(rows) == 0:
                         print("No projects available in database")
                    
                    else:
                        for row in rows:
                            print(f'Project: ', colored(f'{row[0]}','green'), 
                                '\n \t located at', colored(f'{row[1]}','green') )
                

                return rows

            except Exception as error:
                print(error)

        else:
            print("TidyScreen says: \n No project database exists. Create a new project.")
            return []
        
        try:
            conn = sqlite3.connect(main_projects_db)
            cur = conn.cursor()
            cur.execute('SELECT * FROM projects')
            projects =cur.fetchall()
            cur.close()
            for id, project in enumerate(projects):
                print(f"Project: {id}, \n \t Name: {project[0]},  \n \t Located at: {project[1]} \n \t Description: {project[2]}")
        except Exception as error:
            print(f'Error: {error}')

def delete_all_projects(main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')):
        
        try:
            confirm = input("Are you sure you want to delete ALL projects in database? (y/n): ")
            if confirm == "y":
                conn = sqlite3.connect(main_projects_db)
                cur = conn.cursor()
                cur.execute('SELECT name FROM projects')
                project_names = cur.fetchall()
                for name in project_names:
                    delete_single_project(name[0])
                
            else:
                print("Projects deletion canceled")
        
        except Exception as error:
            print(f'Error: {error}')

def delete_single_project(project_name,main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')):
        conn = sqlite3.connect(main_projects_db)
        # This will delete the project datafolder in the local machine
        project_path = get_project_path(conn,project_name)
        shutil.rmtree(project_path)
        # This will delete the project registry in the general database
        delete_project_row(conn,project_name)
        print(f"TidyScreen says: \n Project '{project_name}' has been deleted from TidyScreen")

def delete_project_row(conn,project_name):
        cur = conn.cursor()
        cur.execute(f"DELETE from projects WHERE name = '{project_name}'")
        conn.commit()

def get_project_path(conn,project_name):
        cur = conn.cursor()
        cur.execute(f"SELECT path FROM projects WHERE name = '{project_name}'")
        project_path = cur.fetchall()[0][0]
        return project_path

def check_main_database_presence(main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')):
    # Check if the database exists
    if os.path.exists(main_projects_db):
        return 1
    else:
        return 0 


class ActivateProject:
    
    def __init__(self, project_name):
        self.main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')
        self.project_name = project_name
        ActivateProject.retrieve_project_path(self)
        ActivateProject.retrieve_project_path_variables(self)
        
    def retrieve_project_path(self):
        conn = sqlite3.connect(self.main_projects_db)
        # This will delete the project datafolder in the local machine
        project_path = get_project_path(conn,self.project_name)
        self.project_path = project_path
        
    def retrieve_project_path_variables(self):
        with open(f"{self.project_path}/.project_vars/paths/paths.pkl", "rb") as file:
            paths_dict = pickle.load(file)
        file.close()
        
        self.paths_dict = paths_dict

class ChemSpaceActions:
    
    def __init__(self,project):
        
        self.project = project
        
    def process_raw_cpds_file(self,file):
        
        chemspace_raw_path = self.project.paths_dict['chemspace']['raw_data']
        chemspace_processed_path = self.project.paths_dict['chemspace']['processed_data']
        
        # If raw file not in project, store it
        filename = file.split('/')[-1]
        table_name = filename.split('.')[0]
        destination_file = chemspace_raw_path+'/'+filename
        
        self.project.file_to_process = destination_file
        self.project.table_name = table_name

        if os.path.isfile(destination_file):
            print(f"The raw: '{filename}' file has already been processed")
        else:
            print("The raw file was copied to the raw data directory of the project")
            general_functions.copy_file_to_dest(file,destination_file)
            project = self.project
            datareader.InputRawChemspaceFile(project)

    def list_ligands_tables(self):
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        print(f"Tables in database corresponding to project: '{self.project.project_name}':")
        for table in tables:
            print(table[0])

        conn.close()

    def delete_cpds_table(self,table_name):
        # The code below will delete the table from the project database
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        cursor = conn.cursor()
        
        # Check the presence of the table in the database
        presence = db_ints.check_table_presence(conn,table_name)
        
        if presence == 1:
        
            query = f"DROP TABLE IF EXISTS {table_name}"
            cursor.execute(query)
            conn.commit()
            #conn.close()
            print(f"The table named: '{table_name}' was deleted succesfully")
            
        else: 
            print(f"An error occurred deleting table: '{table_name}'. Check existence.")

        # The code below will remove the raw file from the project folder
        try:
            directory = self.project.paths_dict["chemspace"]["raw_data"]
            for filename in os.listdir(directory):
                if filename.startswith(table_name):
                    file_path = os.path.join(directory,filename)
                    if os.path.isfile(file_path): # This will ensure the object to be deleted is a file
                        os.remove(file_path)
                        print(f"Deleted raw file: '{filename}'")
    
        except:
            pass
    
    def list_cpds_in_table(self,table_name):
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        cursor = conn.cursor()
        cursor.execute(f"SELECT SMILES,Name FROM {table_name}")
        
        rows = cursor.fetchall()
        
        for row in rows:
            print(row)
        
        print(f"A total of {len(rows)} molecules are present in table: '{table_name}'")
        
    def generate_ligands_pdqbt_files(self,table_name):
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        # Check if the corresponding .pdbqt table already exists
        pdbqt_table = table_name + "_pdbqt_files"
        exists = db_ints.check_table_presence(conn,pdbqt_table)
        
        if exists == 1:
            print(f"Table '{table_name}' has already been processed to generate .pdbqt files")
            sys.exit()
        
        else:
            print("Computing .pdbqt files")
            db_ints.create_ligand_pdbqt_table(conn,pdbqt_table)
            ligands_processing.process_ligand_table_to_pdbqt(conn,table_name,pdbqt_table,chemspace_db)


    def enumerate_ligands_stereoisomers(self,table_name):
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        
        try:
            sql = f"""SELECT SMILES, Name FROM {table_name};"""
            df = pd.read_sql_query(sql,conn)
            # This will parallelize the computation of stereoisomer info, and return it into the original df as a column 'stereo_info' in the form of a dictionary
            pandarallel.initialize(progress_bar=True)
            df["stereo_info"] = df.parallel_apply(lambda row: ligands_processing.enumerate_stereoisomers_info(row["SMILES"],row["Name"]), axis=1)
            
            # This will explode the column containing the dictionary of stere info into a dataframe
            exploded_df = ligands_processing.explode_stereo_column(df)
            
            # Store the exploded df into the database adding the '_stereo_enum' suffix to the original table name
            db_ints.store_df_as_table_in_db(conn,exploded_df,table_name+'_stereo_enum')
            
        except:
            print("Error")
        

                       
if __name__ == '__main__':
     
     #list_projects()
     my_project = ActivateProject("test1")
     my_project_CS = ChemSpaceActions(my_project)
     #my_project_CS.process_raw_cpds_file("/home/fredy/Desktop/test_tidyscreen/files/test_smi.smi")
     #my_project_CS.list_ligands_tables()
     #my_project_CS.generate_ligands_pdqbt_files("test_smi")
     #my_project_CS.delete_cpds_table("test2_pdqbt_files")