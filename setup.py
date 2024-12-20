from setuptools import setup, find_packages

setup(name='tidyscreen',
      version='0.0.1',
      description='A package to manage data associated to drug screening campaigns',
      long_description=open('README.md').read().strip(),
      author="MedChemLab staff",
      author_email="aquevedo@unc.edu.ar",
      url='https://github.com/MedChemLab/TidyScreen_revised',
      license='GPLv3',
      python_requires='>=3.10',
      package_data = {'': ['tidyscreen','tidyscreen_dbs']},
      install_requires=[
        'pandas==2.2.1',
        'rdkit==2023.9.5',
        'openbabel-wheel',
        'pysqlite3',
        'termcolor',
	'pandarallel',
      ],
      packages=find_packages(),
      keywords='drug-discovery',
      classifiers=[
          "Programming Language :: Python :: 3.10",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
          "Topic :: virtual high-throughput screening",
      ],
      include_package_data=True
      )
