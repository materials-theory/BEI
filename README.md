# Bond Elongation Index (BEI) 
Bond Elongation Index (BEI), a bond-length-based index, for characterizing polyhedral distortions of perovskite units with the specific purpose of quanterfying the anisotropic ferroelectric distortions using a Python script

# Examples/Tutorials – How to use
Simple example for new users.
Pre-generated VASP structure file: Tutorial_POSCAR.vasp, representing the tetragonal phase of BaTiO<sub>3</sub>.

Usage example : ```python3 bei_mtg.py -c Ti -v O -i BaTiO3_tetra.vasp```

- Octahedral center atom (Ti) and vertex atom (O) must be set by -c and -v option.

- Structure file must be set by -i option.

# Structure files
Crystal_structure_KNO_with_Strain directory contains all POSCAR converged structure files used in "Enhanced polarization in epitaxially strained monoclinic potassium niobate for lead-free electromechanical applications, _J. Mater. Chem. C_ **9**, 13420 (2021)" work

# Publication and Citing BEI
If you use BEI in your research work, please include a citation to this article:

W. Hwang, J.-H. Lee*, and A. Soon*, Enhanced polarization in epitaxially strained monoclinic potassium niobate for lead-free electromechanical applications, _J. Mater. Chem. C_ **9**, 13420 (2021).
