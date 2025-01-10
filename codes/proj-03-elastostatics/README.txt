1. "driver.m"  is developed without gmsh. I used a square plate: ( 0 , 1 ) * ( 0 , 1 ).

2. "g.m" is the Dirichlet boundary condition function that would be used in the manufactured solution.

3."driver2.m" is developed with gmsh. Complement the small 1/4 circle, then the plate in: ( -1 , 1 ) * ( -1 , 1 ).

4. "driver.m" can solve problem with Dirichlet boundary condition but not Neumann. Because it is much easier to use gmsh with "LINES tag", so I did it in "driver2.m". "driver2.m" can solve problem with Dirichlet boundary condition and Neumann boundary condition.
