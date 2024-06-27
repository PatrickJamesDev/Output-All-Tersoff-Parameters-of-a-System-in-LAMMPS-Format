This code is directly taken from Jesper Kristensen. I have made some modifications such as converting the script from Python2 to Python3 and fixing some typos. 
SOURCE: http://jespertoftkristensen.com/JTK/Blog/Entries/2013/12/14_Script_to_Output_any_Tersoff_Mixing_Parameters.html
# Output-All-Tersoff-Parameters-of-a-System-in-LAMMPS-Format
Generating mixed Tersoff potentials automatically with a Python3 script. 

# USAGE
python: PRINT_MIXED_TERSOFF.py --input example1
Supply your own individual Tersoff potentials in your input file according to the example inputs, then the script will provide the mixed potential. I.e. supplying potentials for Si and Ge would output the Si-Ge potential.  
