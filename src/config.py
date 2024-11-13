from pathlib import Path
from termcolor import colored
from distutils.sysconfig import get_python_lib
import sqlite3
import pickle
import shutil

class CreateProject:
    """
    This manages the creation and configuration of a new screening project
    """    
    def __init__(self):
        self.project_config()
        self.store_project_in_database()

    def project_config(self, name=None, path=None):
        if name == None:
            self.project_name = input("Provide the project name to be created: ")
        else:
            self.project_name = name
        if path == None:
            base_dir = input("Provide the projects base directory: ")
            project_path = f'{base_dir}/{self.project_name}'
            self.project_path = project_path
        else:
            self.project_path = path

        self.project_description = input("Provide a brief project description: ")

        self.create_proj_structure()

    def create_proj_structure(self):
        try: 
            # Define a path folder to store usable paths
            levels_paths_dict = {}

            # Create the project level 0 folder
            folder_level_1 = ["chemspace","docking","dynamics"]
            level_0_folders_dict = {}
            for folder in folder_level_1:
                new_folder = f"{self.project_path}/{self.project_name}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)    
                level_0_folders_dict[folder] = new_folder
            
            # Store level_0 paths in the general path dict            
            levels_paths_dict["level_0_paths"] = level_0_folders_dict 
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

            # Create the project level 1 folder for 'chemspace'
            folder_level_chemspace = ["raw_data","processed_data","misc"]
            base_folder = self.project_paths['level_0_paths']['chemspace']
            level_1_chemspace_dict = {}
            for folder in folder_level_chemspace:
                new_folder = f"{base_folder}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
                level_1_chemspace_dict[folder] = new_folder
           
           # Store level_1 chemspace paths in the general path dict            
            levels_paths_dict["level_1_CS_paths"] = level_1_chemspace_dict
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

            # Create the project level 1 folder for 'docking'
            folder_level_docking = ["docking_assays","docking_registers","params","raw_data","receptors"]
            base_folder = self.project_paths['level_0_paths']['docking']
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
            base_folder = self.project_paths['level_0_paths']['dynamics']
            level_1_dynamics_dict = {}
            for folder in folder_level_dynamics:
                new_folder = f"{base_folder}/{folder}"
                Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
                level_1_dynamics_dict[folder] = new_folder
           
           # Create a hidden 'project_vars' folder to store config files
            new_folder = f"{self.project_path}/{self.project_name}/.project_vars"
            Path(f"{new_folder}").mkdir(parents=True, exist_ok=False)
            #os.system(f'mv {new_folder} .{new_folder}') # Make the folder hidden

           # Store level_1 docking paths in the general path dict            
            levels_paths_dict["level_1_DYN_paths"] = level_1_dynamics_dict
            # Associate the path general paths dictionaty to the object
            self.project_paths = levels_paths_dict

        # Save the project paths to a pickle file
            with open(f"{self.project_path}/{self.project_name}/.project_vars/paths.pkl","wb") as paths_file:
                pickle.dump(self.project_paths, paths_file)

        except Exception as error:
            print(colored(f"## Error: \n '{error}' \n ocurred while creating project structure","red"))

    def store_project_in_database(self):
        projects_database = get_python_lib() + "/src"
        conn = sqlite3.connect(f'{projects_database}/projects_database.db')
        self.store_project_in_table(conn)

    def store_project_in_table(self,conn):
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

# This section is to test the module locally

if __name__ == '__main__':

    #project = CreateProject()
    #ProjectsManagement.list_projects()
    ProjectsManagement.delete_all_projects()
    #ProjectsManagement.delete_single_project("test4")
    #ProjectsManagement.get_project_path("test4")


