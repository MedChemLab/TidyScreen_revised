from pathlib import Path
from termcolor import colored
#from distutils.sysconfig import get_python_lib
import sqlite3
import pickle
import shutil
import os
import sys
#import tidyscreen functions
from tidyscreen import database_interactions as db_ints
from tidyscreen import projects_management as proj_mg
from tidyscreen import datareader
from tidyscreen import general_functions
import tidyscreen_dbs

class CreateProject:
    """
    This manages the creation and configuration of a new screening project
    """    
    
    def __init__(self):
        
        self.main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')
        self.project_config()
        self.create_proj_structure()
        self.store_project_in_database()
        # Inform that the project was created
        print(f"TidyScreen says: \n Project {self.project_name} created successfully")

    def project_config(self, name=None, path=None):
        if name == None:
            self.project_name = input("Provide the project name to be created: ")
            
            self.check_project_presence()

        else:
            self.project_name = name
        
        if path == None:
            base_dir = input("Provide the projects base directory: ")
            project_path = f'{base_dir}/{self.project_name}'
            self.project_path = project_path
        else:
            self.project_path = path

        self.project_description = input("Provide a brief project description: ")

    def check_project_presence(self):
        main_db_exists = proj_mg.check_main_database_presence()
        
        # If the main database exists, then check the presence
        
        if main_db_exists == 1:
            available_projects = proj_mg.list_projects(list_rows=0)
            for project in available_projects:
                existing_project_name = project[0]
                # Check if the new project name already exists, and in that case inform and exit
                if self.project_name == existing_project_name:
                    print(f"TidyScreen says: \n The project: '{self.project_name}' already exists. Stopping.")
                    sys.exit()
        else:
            print("TidyScreen says: \n the main database does not exist, I will create it.")

    def create_proj_structure(self):
        # Define the folder structure of a new project
        project_folders_structure = {'chemspace':["raw_data","processed_data","misc"],
                                         'docking':["docking_assays","docking_registers","params",'raw_data','receptors'],
                                         'dynamics':["md_assays","md_registers","md_params"],
                                         '.project_vars':["paths"]}

        # Initialize a dictionary to contain all project paths
        all_paths_dict = {}
        
        try: 
            for base_folder in project_folders_structure.keys():
                current_level_dict = {}
                for folder in project_folders_structure[base_folder]:
                    new_folder = f"{self.project_path}/{base_folder}/{folder}"
                    # Create the corresponding folder
                    Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
                    # Create a dictionary with folder for the current level of folder
                    current_level_dict[folder]=new_folder
                # Create a dictionary for the base level of folders
                all_paths_dict[base_folder]=current_level_dict
            
            # Associate the path general paths dictionaty to the object
            self.project_paths = all_paths_dict
            # Save the project paths to a pickle file
            with open(f"{self.project_path}/.project_vars/paths/paths.pkl","wb") as paths_file:
                pickle.dump(self.project_paths, paths_file)

        except Exception as error:
            print(f"TidyScreen says: \n Error creating project: '{self.project_name}' folder structure")

    def store_project_in_database(self):
        conn = db_ints.connect_to_database(self.main_projects_db)
        cur = conn.cursor()
        # Create the projects table only if it does not exists
        cur.execute('CREATE TABLE IF NOT EXISTS projects (name VARCHAR, path VARCHAR, description VARCHAR)')
        # Store the project register
        cur.execute(f'INSERT INTO projects (name, path, description) values (?,?,?)', (self.project_name, self.project_path,self.project_description))
        conn.commit()

class MainDbConfigs:

    def __init__(self):
        self.main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')

    def list_projects(self):
        
        # Check if the database exists
        if os.path.exists(self.main_projects_db):
            try: 
                conn = db_ints.connect_to_database(self.main_projects_db)
                cursor = conn.cursor()
                cursor.execute("SELECT * FROM projects")
                
                rows = cursor.fetchall()
                

                if len(rows) == 0:
                    print("No projects available in database")
                    
                else:
                    for row in rows:
                        print(f'Project: \n ', colored(f'{row[0]}','green'),  '\n \t located at', colored(f'{row[1]}','green'), '\n \t Description:', colored(f'{row[2]}','green'))


            except Exception as error:
                print(error)

        else:
            print("TidyScreen says: \n No project database exists. Create a new project.")
            return []
        
    def add_chem_filters_from_file(self,filters_file):
        print("Start adding chemical filters to main database \n")
        conn = sqlite3.connect(self.main_projects_db)
        general_functions.input_chem_filters_to_db(conn,filters_file)
        print("Finishing adding chemical filters to main database \n")
    
    def list_available_chem_filter(self):
        conn = sqlite3.connect(self.main_projects_db)
        filters_df = general_functions.retrieve_stored_filters(conn)
        
        for index, row in filters_df.iterrows():
            print(f"{row['Filter_id']}: {row['Filter_Name']} - {row['SMARTS']}")
        
    def add_single_reactions_from_file(self,reactions_file):
        print("Start adding single reactions to main database \n")
        conn = sqlite3.connect(self.main_projects_db)
        general_functions.input_single_reactions_to_db(conn,reactions_file)
        print("Finishing adding single reactions to main database \n")

    def list_available_single_reactions(self):
        conn = sqlite3.connect(self.main_projects_db)
        single_reactions_df = general_functions.retrieve_stored_reactions(conn)
        
        print("These are the reactions available in the main database: \n")
        
        for index, row in single_reactions_df.iterrows():
            print(f"{row['Reaction_id']}: {row['Reaction_Name']} - {row['Reaction_SMARTS']}")