Code for of the global optimization algorithm featured in the article submitted to Ind. Eng. Chem. Res. in June 2023.
"Global Optimization of QCPs using MIP relaxations with a base-2 logarithmic partitioning scheme".

You need to change the directories, specificed in the WUN.gms file. This is particularly relevant in line 266:
"file OptFile /C:\Users\Castro\Documents\GAMS files\GitHub\Pooling_tp\Cplex.op9/;"
It will decide the location of the solution pool solnpool.gdx file, which will be accessed in line 399.