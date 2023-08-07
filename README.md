How to use:

1- create an environment with Prerequisites

You can use virtualenv package. An easy-to-use guide is available here:
https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/creating-a-virtual-environment

or you can use these commands:
Windows:
```
py -m pip install --user virtualenv
```

Unix:
```
python3 -m pip install --user virtualenv
```

Creating a virtual environment with:

Windows:
```
py -m venv env
```

Unix:
```
python3 -m venv env
```

env is the name of the folder/environment
activate the new environment with:

Windows
```
.\env\Scripts\activate
```

Unix
```
source env/bin/activate
```
extract the compressed file 

install required packages with:
```
py -m pip install -r requirements.txt
```

2-Using the Package:

2-1- Put your peptide sequences in the "input_peptide_sequences.txt" file inside the "classifier" folder.


2-2- Run this command in the classifier folder and check the results in the "results.txt" file

Windos CMD:
```
python classifier.py > results.txt
```
Windows Explorer:

Execute the "runThisOnWindows.bat" file.


Unix
```
python3 classifier.py > results.txt
```

