import os
import sqlite3
from termcolor import colored
from distutils.sysconfig import get_python_lib
import shutil
from tidyscreen import database_interactions as db_ints
import pickle
import datareader
import general_functions
from tidyscreen import database_interactions as db_ints
import sys

def list_projects(list_rows=1):
        database_path = get_python_lib() + "/tidyscreen_dbs"
        database_file = f'{database_path}/projects_database.db'
        # Check if the database exists
        if os.path.exists(database_file):
            try: 
                conn = db_ints.connect_to_database(database_file)
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
        
        projects_database = get_python_lib() + "/src"
        try:
            conn = sqlite3.connect(f'{projects_database}/projects_database.db')
            cur = conn.cursor()
            cur.execute('SELECT * FROM projects')
            projects =cur.fetchall()
            cur.close()
            for id, project in enumerate(projects):
                print(f"Project: {id}, \n \t Name: {project[0]},  \n \t Located at: {project[1]} \n \t Description: {project[2]}")
        except Exception as error:
            print(f'Error: {error}')

def delete_all_projects():
        projects_database = get_python_lib() + "/tidyscreen_dbs"
        try:
            confirm = input("Are you sure you want to delete ALL projects in database? (y/n): ")
            if confirm == "y":
                conn = sqlite3.connect(f'{projects_database}/projects_database.db')
                cur = conn.cursor()
                cur.execute('SELECT name FROM projects')
                project_names = cur.fetchall()
                for name in project_names:
                    delete_single_project(name[0])
                
            else:
                print("Projects deletion canceled")
        
        except Exception as error:
            print(f'Error: {error}')

def delete_single_project(project_name):
        projects_database = get_python_lib() + "/tidyscreen_dbs"
        conn = sqlite3.connect(f'{projects_database}/projects_database.db')
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

def check_main_database_presence():
    database_path = get_python_lib() + "/tidyscreen_dbs"
    database_file = f'{database_path}/projects_database.db'
    # Check if the database exists
    if os.path.exists(database_file):
        return 1
    else:
        return 0 


class ActivateProject:
    
    def __init__(self, project_name):
        #print(" \n TidyScreen says: \n 'Avaliable projects in database:' \n")
        #list_projects()
        self.project_name = project_name
        ActivateProject.retrieve_project_path(self)
        ActivateProject.retrieve_project_path_variables(self)
        
    def retrieve_project_path(self):
        projects_database = get_python_lib() + "/tidyscreen_dbs"
        conn = sqlite3.connect(f'{projects_database}/projects_database.db')
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
        print("Tables in the database:")
        for table in tables:
            print(table[0])

        conn.close()

    def generate_ligands_pdqbt_files(self,table_name):
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        pdqbt_table = table_name + "_pdqbt_files"

        # Check if the corresponding .pdbqt table already exists
        exists = db_ints.check_table_presence(conn,pdqbt_table)
        
        if exists == 1:
            print(f"Table '{table_name}' has already been processed to generate .pdbqt files")
            sys.exit()
        
        else:
             print("Computing .pdbqt files")
             db_ints.create_ligand_pdbqt_table(conn,pdqbt_table)

            
if __name__ == '__main__':
     
     #list_projects()
     my_project = ActivateProject("test1")
     my_project_CS = ChemSpaceActions(my_project)
     #my_project_CS.process_raw_cpds_file("/home/fredy/Desktop/test_tidyscreen/files/test_smi.smi")
     #my_project_CS.list_ligands_tables()
     my_project_CS.generate_ligands_pdqbt_files("test_smi")