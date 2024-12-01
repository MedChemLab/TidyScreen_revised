import sqlite3
from distutils.sysconfig import get_python_lib
from termcolor import colored

def connect_to_database(database_full_path):
    conn = sqlite3.connect(f'{database_full_path}')
    return conn

if __name__ == '__main__':
    pass



