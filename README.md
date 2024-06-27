This code is directly taken from Jesper Kristensen. I have made some modifications such as converting the script from Python2 to Python3 and fixing some typos. 
SOURCE: https://github.com/PatrickJamesDev/Output-All-Tersoff-Parameters-of-a-System-in-LAMMPS-Format/edit/main/README.md
# Output-All-Tersoff-Parameters-of-a-System-in-LAMMPS-Format
Generating mixed Tersoff potentials automatically with a Python3 script. 

# USAGE
python: PRINT_MIXED_TERSOFF.py --input example1
Supply your own individual Tersoff potentials in your input file according to the example inputs, then the script will provide the mixed potential. I.e. supplying potentials for Si and Ge would output the Si-Ge potential.  
