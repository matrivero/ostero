# Description #

Ostero is a finite element code that solves implicitly the equilibrium equation which governs the non-linear mechanics of a deformable body subjected to large deformations. In other words, Ostero allows to determine the response of a deformable solid body to an applied load. 

Ostero intends to be a didactic code, and its main objective is to allow the user to understand the very basic structure of a non-linear mechanics finite elements code and to provide a framework for the beta testing of different models (elastic material models, contact models, fracture, plasticity, etc.). 

Ostero is based on the solidz module of the [Alya code](http://www.bsc.es/alya) and so far it can solve 2D problems using triangles or quadrilateral elements and an isolinear or neo-hookean material model.

Ostero can read mesh files generated using [Gmsh](http://gmsh.info) and write ouputs in Vtk ASCII format, which can easily be postprocessed with [Paraview](http://www.paraview.org).

Ostero is a mixture between Python and Fortran, and tries to take advantage of the main properties of each code. The parsing of the inputs parameters as well as boundaries conditions and other options is done using a main Python code. The solver part is coded in Fortran and included in the Python code as an external library. 

# Usage #

In order to execute Ostero you only need a Python interpreter and the F2PY package, which is a Fortran to Python interface generator. Since 2007, F2PY is part of [Numpy](http://docs.scipy.org/doc/numpy-dev/f2py).

To start using Ostero, first you need to compile the external Fortran library by doing:

```
#!bash

f2py -c -m external external.f90
```

Once you have checked the external library was generated (external.so), you are ready to execute the example cases. First change the directory to the example folder i.e.:

```
#!bash

cd examples/square_with_hole
```

Then, to execute the example you must do:

```
#!bash

../../finite_strain.py input_file.dat boundary_file.dat
```

The main program of Ostero is *finite_strain.py*, while the first argument is the input file and the second is the boundary file. In the input file you must specify:

* the mesh path (keyword *$mesh_path*)

* the constitutive model (keyword *$constitutive_model*, options *ISOL*, *BELY*, *ZIEN* and *LAUR*)

* the time step (keyword *$time_step_size*)

* the total number of time steps (keyword *$total_steps*)

In the boundary file there are only three keywords:

* the volume definition and its mechanical properties: young modulus and poisson parameter (keyword *$VolumeDefinition*)

* displacement boundary conditions (Dirichlet) (keyword *$BoundaryConditionsDisplacement*)

* pressure boundary conditions (Neumann) (keyword *$BoundaryConditionsPressure*)

# Some results obtained with Ostero... #

Pressure (top) + Fixed displacement (bottom) - Triangles
![ostero_square_with_hole.png](https://bitbucket.org/repo/a69BrG/images/1213456489-ostero_square_with_hole.png)

Pressure (top) + Fixed displacement (bottom) - Quads
![ostero_square_quads.png](https://bitbucket.org/repo/a69BrG/images/857170256-ostero_square_quads.png)

Pressure + Finite displacement + Fixed displacement - Triangles
![ostero_complete.png](https://bitbucket.org/repo/a69BrG/images/1379224850-ostero_complete.png)

Two materials: Pressure (top) + Fixed displacement (bottom) - Triangles
![ostero_two_materials.png](https://bitbucket.org/repo/a69BrG/images/3766508879-ostero_two_materials.png)

Ostero is licensed under [GNU GPLv3](http://www.gnu.org/copyleft/gpl.html)

Ostero comes with ABSOLUTELY NO WARRANTY. 

Ostero is free software, and you are welcome to redistribute it under certain conditions.

**If you have any doubt, problem or suggestion please don't hesitate to contact me: matias.rivero(at)bsc.es **