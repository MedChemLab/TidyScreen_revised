import sqlite3
from distutils.sysconfig import get_python_lib
from termcolor import colored

def connect_to_database(path, database_name):
    conn = sqlite3.connect(f'{path}/{database_name}')
    return conn

def list_available_projects():
    database_path = get_python_lib() + "/src"
    conn = connect_to_database(database_path,"projects_database.db")
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM projects")
    
    rows = cursor.fetchall()
    for row in rows:
        print(f'Project: ', colored(f'{row[0]}','green'), 
              '\n \t located at', colored(f'{row[1]}','green') )

if __name__ == '__main__':
    list_available_projects()



