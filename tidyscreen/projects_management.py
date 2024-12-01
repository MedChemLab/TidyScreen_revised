import os
import sqlite3
from termcolor import colored
from distutils.sysconfig import get_python_lib
import shutil
from tidyscreen import database_interactions as db_ints

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

if __name__ == '__main__':
     list_projects()