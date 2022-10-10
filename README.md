# About
* This code can be used to calculate survival probability of interfacial water(IW).
* IW are considered as those water molecule which are within 6 Angstorm from the protein backbone atom.


# Required input file
* A sample input file of a single frame is given as Input.pdb
* Arrange atom according to this pdb file using gmx trjconv or any other trajectory processing tools.

# Compilation
To run the code:
* gfortan residence-time.f95 -o a.out
* ./a.out
