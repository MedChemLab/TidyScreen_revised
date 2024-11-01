from pathlib import Path

class campaing_configurator():
    """
    This manages the creation and configuration of a new screening project
    """
    
    def __init__(self) -> None:
        pass

    def project_config(self, name=None, path=None):
        if name == None:
            project_name = input("Provide the project name to be created: ")
        else:
            self.project_name = project_name

        if path == None:
            project_path = input("Provide the projects base directory: ")
        else:
            self.project_path = project_path

    def create_structure(self):
        Path(f"{self.project_path}/{self.project_name}").mkdir(parents=True, exist_ok=False)


