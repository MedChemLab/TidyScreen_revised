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
import json


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
        self.main_projects_db = sys.modules['tidyscreen_dbs'].__file__.replace('__init__.py','projects_database.db')
        # Create a reference to the chemspace db
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        project_db_conn = db_ints.connect_to_database(chemspace_db)
        self.project_CS_conn = project_db_conn

    def process_raw_cpds_file(self,file):
        # This is the NEW CODE
        # Generate the table_name to be eventually stored in the database
        table_name = file.split('/')[-1].split('.')[0]
        # Check if the table to be generated already exists in the ChemSpace database
        table_exists = db_ints.check_table_presence(self.project_CS_conn,table_name)

        if table_exists == 1:
            action = input(f"The table named: \n '{table_name}' \n already exists in the ChemSpace database.Do you want to delete it in order to proceed analysis? (y/n): ")

            if action == 'y':
                db_ints.delete_table_from_db(self.project_CS_conn,table_name)
            else:
                print(f"Preventing overwritting of: '{table_name} \n table. Stopping.'")
                sys.exit()
        
        ## If the analysis is not stopped before, continue the data reading analysis
        # Copy the raw data file to the corresponding ChemSpace folder in order to back it up
        #print(destination_file)
        destination_file = f"{self.project.paths_dict['chemspace']['raw_data']}/{file.split('/')[-1]}"
        general_functions.copy_file_to_dest(file,destination_file)

        # Execute the data processing
        self.project.file_to_process = destination_file
        self.project.table_name = table_name
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

    def create_filters_pipeline(self):

        # This section will ask for filter and inclusion/exclusion criteria until the 'end' keyword is entered
        filters_number_list = []
        filters_instances_list = []
        filters_criteria_list = []
        
        while True:
            filter_number = input("Enter the filter number to add (integer; 'quit' to end): ")

            if filter_number == 'quit':
                print("Finishing filters addition to the pipeline")
                filters_pipeline_name = input("Input the filters pipeline name: ")
                break
            else:
                filter_instances = input("Enter the number of allowed instances (integer): ")
                filters_number_list.append(int(filter_number))
                filters_instances_list.append(int(filter_instances))
                while True:
                    filter_citeria = input("Indicate if selected filter should be inclusive ('I') or exclusive ('E'): ")    
                    if filter_citeria == 'I' or filter_citeria == 'E':
                        filters_criteria_list.append(filter_citeria)
                        break
                    else:
                        print("Only the keyword 'I' or 'E' are accepted. Try again.")
                        
        # Check the presence of filters pipelines table in the project database
        project_db_conn = self.project_CS_conn
        exists = db_ints.check_table_presence(project_db_conn,"filters_pipeline")

        # If the filters pipelines table does not exists, create the table before storing the pipeline
        if exists == 0:
            db_ints.create_filters_pipelines_table(project_db_conn)

        try:
            main_database_conn = sqlite3.connect(self.main_projects_db)
            general_functions.store_filters_pipeline(project_db_conn,main_database_conn,filters_number_list,filters_instances_list,filters_criteria_list,filters_pipeline_name)
            print(f"The filter: '{filters_pipeline_name}' was successfully stored.")

        except Exception as error:
            print(error)
            
    def list_project_filters_pipelines(self,table_name="filters_pipelines"):
        # Connect to the project CSA database
        conn = self.project_CS_conn
        
        sql = f"""SELECT * FROM {table_name};"""
        df = pd.read_sql_query(sql,conn)
        
        print(df)
    
    def apply_filters_pipeline(self,filters_pipeline_id,table_name):
        # Connect to the project CSA database
        conn = self.project_CS_conn
        
        # Construct the pipeline filtering lists: i.e. criteia and smarts
        stored_filters_pipeline_name, filters_id_list, filters_instances_list, filters_criteria_list, filters_SMARTS_list = db_ints.retrieve_filters_pipeline_smarts(conn, filters_pipeline_id)
        
        # Perform the filtering on the table
        ligands_df_filtered = ligands_processing.filter_ligands_table_with_smarts(conn,table_name,filters_instances_list, filters_criteria_list, filters_SMARTS_list)

        # Store the filtered df into the CS database
        filtered_table_name = table_name+f'_filter_id_{filters_pipeline_id}'
        db_ints.store_df_as_table_in_db(conn,ligands_df_filtered,filtered_table_name,action="replace")
        
        print(f" Table '{table_name}' filtered on pipeline id: '{filters_pipeline_id}' was stored as '{filtered_table_name}'")
        
    def export_filters_pipeline(self,filters_pipeline_id):
        # Connect to the project CSA database
        conn = self.project_CS_conn
        
        # Construct the pipeline filtering lists: i.e. criteia and smarts
        stored_filters_pipeline_name, filters_id_list, filters_instances_list, filters_criteria_list, filters_SMARTS_list = db_ints.retrieve_filters_pipeline_smarts(conn, filters_pipeline_id)
        
        filters_pipeline_keys = ["Filters_pipeline_name","Filters_id","Filters_instances","Filters_criteria","Filters_SMARTS"]
        filters_pipeline_data = [stored_filters_pipeline_name,filters_id_list,filters_instances_list,filters_criteria_list,filters_SMARTS_list]
        
        dict_for_json = general_functions.create_json_object_from_lists(filters_pipeline_keys,filters_pipeline_data)
        
        # Set the output path of the .json file to the /misc folder in Chemspace
        output_path = self.project.paths_dict["chemspace"]["misc"]
        json_filename  = f"filter_pipeline_{filters_pipeline_id}.json"
        with open(f"{output_path}/{json_filename}", "w") as json_file:
            json.dump(dict_for_json, json_file,indent=4)
        
        json_file.close()
        
        print(f"Filters_pipeline written to {output_path}/{json_filename}")
    
    def import_filters_pipeline_from_json_file(self,json_file):
        
        with open(json_file,"r") as file:
            filters_pipeline_dict = json.load(file)
            
        # This section will retrieve the info from the corresponding dict
        filters_pipeline_name = filters_pipeline_dict["Filters_pipeline_name"]
        filters_number_list = filters_pipeline_dict["Filters_id"]
        filters_instances_list = filters_pipeline_dict["Filters_instances"]
        filters_criteria_list = filters_pipeline_dict["Filters_criteria"]
        filters_SMARTS_list = filters_pipeline_dict["Filters_SMARTS"]
        
        # Check the presence of filters pipelines table in the project database
        project_db_conn = self.project_CS_conn
        exists = db_ints.check_table_presence(project_db_conn,"filters_pipeline")
        
        # If the filters pipelines does not exists, create the table before storing the pipeline
        if exists == 0:
            db_ints.create_filters_pipelines_table(project_db_conn)
        
        # Here we need to check that each Filters_SMARTS exists in the main TidyScreen database, and that the numbering of each SMARTS is consistant with that stored in the main DB of filters
        main_database_conn = sqlite3.connect(self.main_projects_db)
        main_existing_filters_df = general_functions.retrieve_stored_filters(main_database_conn)
        general_functions.check_id_smarts_pair_in_main_filters_db(main_existing_filters_df,filters_number_list,filters_SMARTS_list)
        
        try:
            general_functions.store_filters_pipeline(project_db_conn,main_database_conn,filters_number_list,filters_instances_list,filters_criteria_list,filters_pipeline_name)
            print(f"The filter: '{filters_pipeline_name}' was sucessfully stored.")
        
        except Exception as error:
            print(error)
    
    def create_reaction_scheme(self):
        
        # This section will ask for sequential reaction if lists 
        reactions_id_list = []
    
        while True:
            reaction_id = input("Enter the reaction id to add to the synthetic scheme (integer; 'quit' to end): ")

            if reaction_id == 'quit':
                print("Finishing reactions addition to the synthetic scheme")
                reaction_scheme_name = input("Input the reaction scheme name: ")
                break
            else:
                reactions_id_list.append(reaction_id)
            
        # Check the presence of reaction schemes table in the project database
        project_db_conn = self.project_CS_conn
        exists = db_ints.check_table_presence(project_db_conn,"reaction_schemes")
    
        # If the filters pipelines table does not exists, create the table before storing the pipeline
        if exists == 0:
            db_ints.create_reaction_schemes_table(project_db_conn)
    
        try:
            main_database_conn = sqlite3.connect(self.main_projects_db)
            general_functions.store_reactions_scheme(project_db_conn,main_database_conn,reactions_id_list,reaction_scheme_name)
            print(f"The reaction scheme: '{reaction_scheme_name}' was successfully stored.")

        except Exception as error:
            print(error)
    
    def list_available_reaction_schemes(self,table_name="reaction_schemes"):
        chemspace_db = f"{self.project.paths_dict['chemspace']['processed_data']}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        cursor = conn.cursor()
        cursor.execute(f"SELECT * FROM {table_name}")
        rows = cursor.fetchall()
        
        for row in rows:
            print(row)
        pass

    def perform_reaction_on_table(self,reaction_scheme_id,reactants_table_names_list,store_flag=0,products_table_name="none"):
        # Connect to the project CSA database
        conn = self.project_CS_conn

        # Retrieve the required info respect to the reaction to be applied
        stored_reaction_scheme_name, stored_reactions_id_json, stored_reactions_SMARTS_json = db_ints.retrieve_reaction_scheme_smarts(conn,reaction_scheme_id)
        
        
        stored_reactions_id_list = json.loads(stored_reactions_id_json)
        stored_reactions_SMARTS_list = json.loads(stored_reactions_SMARTS_json)

        # This function will construct the initial reactants dataframe, that will be afterwards incremented as a consequence of the 'for' loops associated to each reaction in the scheme
        #molecules_df = ligands_processing.retrieve_initial_reactants(conn,reactants_table_names_list)
        molecules_df = ligands_processing.retrieve_initial_reactants_updated(conn,reactants_table_names_list)

        for scheme_index, reaction_smarts in enumerate(stored_reactions_SMARTS_list):
            
            # Check if the number of reactants tables is consistant with the current SMARTS. If fails, will inform and exit.
            general_functions.check_tables_vs_smarts(reactants_table_names_list[scheme_index],reaction_smarts)
                
            try:
                current_reaction_step_table_names = reactants_table_names_list[scheme_index]
                
            except:
                print(f"Problem retrieving names of compounds for the given synthetic schema. Revise the following info: \n List reactants to be processed in the synthetic pipeline: \n{reactants_table_names_list} \n Reactions id's involved in the schema: \n {stored_reactions_id_list}")
                sys.exit()
            
            # # Evaluate the type of reaction involved by analyzing the given 'reaction_smarts'
            molecularity_nbr = general_functions.evaluate_smarts_reaction(reaction_smarts)

            if molecularity_nbr == 1:
                print(f"Applying the following unimolecular reaction \n {reaction_smarts} \n in step {scheme_index} of the synthetic scheme")
                
                # This is the original synthetic method based on the molecules df
                molecules_df = ligands_processing.apply_unimolecular_reaction(molecules_df,current_reaction_step_table_names[0],reaction_smarts,scheme_index)

            elif molecularity_nbr == 2:
                print(f"Applying a bimolecular reaction in step {scheme_index} of the synthetic scheme")
                molecules_df = ligands_processing.apply_bimolecular_reaction(conn,molecules_df,current_reaction_step_table_names,reaction_smarts,scheme_index)
            
            else:
                print(f"A reaction of molecularity {molecularity_nbr} is currently not supported")

        if store_flag == 0:
            general_functions.report_on_reaction_scheme(molecules_df)
            molecules_df.to_csv("/home/fredy/Desktop/test_tidyscreen/products.csv")

        elif store_flag == 1:
            # Perform a cleaning and data identification on the resulting df prior to storage
            #molecules_df = general_functions.clean_reaction_scheme_product_df(molecules_df,scheme_index)
            molecules_df = general_functions.clean_reaction_scheme_product_df(molecules_df,reaction_scheme_id)

            # Store the products dataframe as a table in the CS project database
            if products_table_name == "none":
                print("You need to provide a products table name expliciting the variable: 'products_table_name'. Stopping")
                sys.exit()

            db_ints.store_df_as_table_in_db(conn,molecules_df,products_table_name)
            print(f"Successfully applied and stored the reaction scheme {reaction_scheme_id} on table '{reactants_table_names_list[0]}'")
        
        else:
            print(f"The 'store_flag' option was not provided adequately: {store_flag}. Stopping.")
            sys.exit()

    def plot_table_of_molecules(self,table_name,max_mols_ppage=25):
        # Connect to the project CSA database
        conn = self.project_CS_conn
        
        # Check the existence of the 'table_name' provided
        table_exists = db_ints.check_table_presence(conn, table_name)
        if table_exists == 0:
            print(f"The table named: \n {table_name} \n does not exists. Stopping")
            sys.exit()

        # The .png files corresponding to the molecules will be written to the /misc folder in ChemSpace 
        output_folder = f'{self.project.paths_dict["chemspace"]["misc"]}/{table_name}'
        # Check the existence of the output folder
        folder_exists = general_functions.check_folder_existence(output_folder)
        
        if folder_exists == 1:
            action = input(f"The folder: \n {output_folder} \n already exists, do you want to delete it and proceed? (y/n): ")
            if action == 'y':
                shutil.rmtree(output_folder)
            else:
                print(f"Stopping analysis to prevent overwritting of: \n {output_folder}")
                sys.exit()
        
        # Create the folder for the ouput
        general_functions.create_folder(output_folder)

        # Generate the depictions for the given table
        ligands_processing.depict_ligands_table(conn,table_name,output_folder,max_mols_ppage)

        print(f"Files (.png) with plotted molecules were successfully written to: \n {output_folder}")