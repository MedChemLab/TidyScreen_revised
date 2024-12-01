from pathlib import Path
from termcolor import colored
from distutils.sysconfig import get_python_lib
import sqlite3
import pickle
import shutil
import os
import sys
#import database_interactions as dbs_ints
from tidyscreen import database_interactions as db_ints

class CreateProject:
    """
    This manages the creation and configuration of a new screening project
    """    
    def __init__(self):
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
        available_projects = list_available_projects(list_rows=0) # Do not list rows, just obtain them
        for project in available_projects:
            existing_project_name = project[0]
            # Check if the new project name already exists, and in that case inform and exit
            if self.project_name == existing_project_name:
                print(f"TidyScreen says: \n The project: '{self.project_name}' already exists. Stopping.")
                sys.exit()

    def create_proj_structure(self):
        try: 
            # Define a path folder to store project relevant paths
            levels_paths_dict = {}

            # Create the project base level folder
            project_base_level_folders = ["chemspace","docking","dynamics"]
            project_base_level_folders_dict = {}
            for folder in project_base_level_folders:
                new_folder = f"{self.project_path}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)    
                project_base_level_folders_dict[folder] = new_folder
            
            # Store level_0 paths as a in the general path dict            
            levels_paths_dict["project_base_paths"] = project_base_level_folders_dict 
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

            # Create the project level 1 folder for 'chemspace'
            folder_level_chemspace = ["raw_data","processed_data","misc"]
            base_folder = self.project_paths['project_base_paths']['chemspace']
            level_chemspace_dict = {}
            for folder in folder_level_chemspace:
                new_folder = f"{base_folder}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
                level_chemspace_dict[folder] = new_folder
           
           # Store level_1 chemspace paths in the general path dict            
            levels_paths_dict["level_1_CS_paths"] = level_chemspace_dict
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

            # Create the project level 1 folder for 'docking'
            folder_level_docking = ["docking_assays","docking_registers","params","raw_data","receptors"]
            base_folder = self.project_paths['project_base_paths']['docking']
            level_1_docking_dict = {}
            for folder in folder_level_docking:
                new_folder = f"{base_folder}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
                level_1_docking_dict[folder] = new_folder
           
           # Store level_1 docking paths in the general path dict            
            levels_paths_dict["level_1_DOCK_paths"] = level_1_docking_dict
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

            # Create the project level 1 folder for 'dynamics'
            folder_level_dynamics = ["md_assays","md_registers","md_params"]
            base_folder = self.project_paths['project_base_paths']['dynamics']
            level_dynamics_dict = {}
            for folder in folder_level_dynamics:
                new_folder = f"{base_folder}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
                level_dynamics_dict[folder] = new_folder
           
           # Create a hidden 'project_vars' folder to store config files
            new_folder = f"{self.project_path}/project_vars"
            new_folder_hidden = f"{self.project_path}/.project_vars"
            Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
            os.system(f'mv {new_folder} {new_folder_hidden}') # Make the folder hidden
            

           # Store level_dynamics paths in the general path dict            
            levels_paths_dict["level_DYN_paths"] = level_dynamics_dict
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

        # Save the project paths to a pickle file
            with open(f"{self.project_path}/.project_vars/paths.pkl","wb") as paths_file:
                pickle.dump(self.project_paths, paths_file)

        except Exception as error:
            print(colored(f"## Error: \n '{error}' \n ocurred while creating project structure","red"))

    def store_project_in_database(self):
        projects_database = get_python_lib() + "/tidyscreen/projects_database.db"
        conn = db_ints.connect_to_database(projects_database)
        cur = conn.cursor()
        # Create the projects table only if it does not exists
        cur.execute('CREATE TABLE IF NOT EXISTS projects (name VARCHAR, path VARCHAR, description VARCHAR)')
        # Store the project register
        cur.execute(f'INSERT INTO projects (name, path, description) values (?,?,?)', (self.project_name, self.project_path,self.project_description))
        conn.commit()
   
class ProjectsManagement:
    """
    This will manage actions related to querying the projects database
    """
    def __init__(self):
        #self.list_projects()
        pass

    def list_projects():
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
        projects_database = get_python_lib() + "/src"
        try:
            confirm = input("Are you sure you want to delete ALL projects in database? (y/n): ")
            if confirm == "y":
                conn = sqlite3.connect(f'{projects_database}/projects_database.db')
                cur = conn.cursor()
                cur.execute('SELECT name FROM projects')
                project_names = cur.fetchall()
                for name in project_names:
                    ProjectsManagement.delete_single_project(name[0])
                
            else:
                print("Projects deletion canceled")
        
        except Exception as error:
            print(f'Error: {error}')

    def delete_single_project(project_name):
        projects_database = get_python_lib() + "/src"
        conn = sqlite3.connect(f'{projects_database}/projects_database.db')
        # This will delete the project datafolder in the local machine
        project_path = ProjectsManagement.get_project_path(conn,project_name)
        shutil.rmtree(project_path)
        # This will delete the project registry in the general database
        ProjectsManagement.delete_project_row(conn,project_name)

    def delete_project_row(conn,project_name):
        cur = conn.cursor()
        cur.execute(f"DELETE from projects WHERE name = '{project_name}'")
        conn.commit()

    def get_project_path(conn,project_name):
        cur = conn.cursor()
        cur.execute(f"SELECT path FROM projects WHERE name = '{project_name}'")
        project_path = cur.fetchall()[0][0]
        return project_path

def list_available_projects(list_rows=1):
    database_path = get_python_lib() + "/tidyscreen"
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
            
                for row in rows:
                    print(f'Project: ', colored(f'{row[0]}','green'), 
                        '\n \t located at', colored(f'{row[1]}','green') )
            
            return rows

        except Exception as error:
            print(error)

    else:
        print("TidyScreen says: \n No project database exists. Create a new project.")
        return []


# This section is to test the module locally

if __name__ == '__main__':

    #list_available_projects()
    project = CreateProject()
    #ProjectsManagement.list_projects()
    #ProjectsManagement.delete_all_projects()
    #ProjectsManagement.delete_single_project("test4")
    #ProjectsManagement.get_project_path("test4")


