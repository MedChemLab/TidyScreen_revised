import pandas as pd
from termcolor import colored
from tidyscreen import database_interactions as db_ints

class InputRawChemspaceFile:
    """
    This will read a file and parse the content in order to store the SMILES and eventually their names.
    """
    def __init__(self, project):
        # This method will apply all the actions on the provided file
        self.project = project
        self.print_parsing_info()
        self.query_user_about_parsing()
        self.parse_data()
        self.store_df_in_database()

    def print_parsing_info(self):
        # Ask for user input on how to parse the file upon showing 2 lines:
        with open (self.project.file_to_process) as input_file:
            nbr_lines = len(input_file.readlines())
            input_file.seek(0) # Go to the start of the file
            if nbr_lines > 1:
                for i in range(2):
                    line = next(input_file).strip()
                    if i  == 0:
                        print(colored(f"Sample lines:\n","blue"), colored(f"{line}","green"))
                    else:
                        print(colored(f" {line}","green"))
            else:
                for i in range(1):
                    line = next(input_file).strip()
                    print(colored(f"Sample lines:\n","blue"), colored(f"{line}","green"))

    def query_user_about_parsing(self):
        
        # Indicate if header is present
        while True:
            header_info = input("Header presence (y/n): ")
            if header_info == "y" or header_info == "n":
                break
            else:
                print("Only 'Yes (y)' or 'No (n)' values are accepted. Input again")
        
        
        # Indicate the separator string
        separator_info = input("Indicate field separator string: ")
        
        # Indicate the filed number containing the SMILES
        while True:
            try:
                field_to_parse = int(input("Field number containing SMILES to parse (0-based index): "))
                break
            except:
                print("Only Integers describing fields are accepted. Input again")

        # Indicate if molecules names are present
        while True:
            names_info = input("Does the file contains molecule names? (y/n): ")
            if names_info == "y" or names_info == "n":
                # In case names are present, indicate the field to parse
                if names_info == "y":
                    try:
                        names_field_to_parse = int(input("Field number containing NAMES to parse (0-based index): "))
                    except:
                        print("Only Integers describing fields are accepted. Input again")
                else:
                    names_field_to_parse = -1 # This value indicated that no name is present
                
                break
            else:
                print("Only 'Yes (y)' or 'No (n)' values are accepted. Input again")


        parsing_options_dict = {"header": header_info,
                                "sep": separator_info,
                                "smi_field": field_to_parse,
                                "name_field": names_field_to_parse}
        
        self.parsing_options_dict = parsing_options_dict
        

    def parse_data(self):
        if self.parsing_options_dict["header"] == 'y':
            header_flag = 0
        else:
            header_flag = None
        
        df = pd.read_csv(self.project.file_to_process, header=header_flag, usecols=[self.parsing_options_dict["smi_field"],self.parsing_options_dict["name_field"]], names=['SMILES','Name'], sep=self.parsing_options_dict["sep"])
        
        self.raw_data_df = df

    def store_df_in_database(self):
        db_path = self.project.paths_dict["chemspace"]["processed_data"]
        chemspace_db = f"{db_path}/chemspace.db"
        conn = db_ints.connect_to_database(chemspace_db)
        
        # Store the dataframe in the corresponding database
        self.raw_data_df.to_sql(con=conn, name=self.project.table_name,if_exists="replace",index=None)
        
        

if __name__ == '__main__':
    
    #input_file = InputCsvFile("/home/fredy/Desktop/tidyscreen/files/test_smi.smi")
    input_file = InputRawChemspaceFile("/home/fredy/Desktop/tidyscreen/files/test_smi_no_name.smi")

