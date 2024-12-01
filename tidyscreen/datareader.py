import pandas as pd
from termcolor import colored

class input_csv_file:
    """
    This will read a .csv file and parse the SMILES field value.
    """
    def __init__(self, file):  
        # This method will apply all the actions on the provided file
        self.file = file
        self.user_parsing_info()
        self.query_user_about_parsing()
        self.parse_data()

    def user_parsing_info(self):
        # Ask for user input on how to parse the file upon showing 2 lines:
        with open (self.file) as input_file:
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
        while True:
            header_info = input("Header presence (y/n): ")
            if header_info == "y" or header_info == "n":
                break
            else:
                print("Only 'Yes (y)' or 'No (n)' values are accepted. Input again")
        
        separator_info = input("Indicate field separator string: ")
        
        while True:
            try:
                field_to_parse = int(input("Field number containing SMILES to parse (0-based index): "))
                break
            except:
                print("Only Integers describing fields are accepted. Input again")

        parsing_options_dict = {"header": header_info,
                                "sep": separator_info,
                                "smi_field": field_to_parse}
        
        self.parsing_options_dict = parsing_options_dict
        

    def parse_data(self):
        
        if self.parsing_options_dict["header"] == 'y':
            header_flag = 0
        else:
            header_flag = None
        
        df = pd.read_csv(self.file, header=header_flag, usecols=[self.parsing_options_dict["smi_field"]], names=['SMILES'], sep=self.parsing_options_dict["sep"])
        
        print(df)

    def store_datafile(self):
        pass
        


if __name__ == '__main__':
    
    input_file = input_csv_file("/home/fredy/Desktop/test_tidyscreen/test2.csv")

