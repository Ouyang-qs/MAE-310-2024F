1. "driver.m"  is developed without gmsh. I used a square plate: ( 0 , 1 ) * ( 0 , 1 ).

2. "g.m" is the Dirichlet boundary condition function that would be used in the manufactured solution.

3."driver2.m" is developed with gmsh. Complement the small 1/4 circle, then the plate in: ( -1 , 1 ) * ( -1 , 1 ).

4. "driver.m" can solve problem with Dirichlet boundary condition but not Neumann. Because it is much easier to use gmsh with "LINES tag", so I did it in "driver2.m". "driver2.m" can solve problem with Dirichlet boundary condition and Neumann boundary condition.
the boundary conditions are "g_win.m" and "h_win.m"

5. function "stress_cartesian.m": input coordinates and R, output [stress_xx, stress_yy,  stress_xy], but Tx=1 in this function, so in later calculation in "driver2.m" need mutiply Tx

6. "driver3.m" is used for the last problem, actually it's almost the same with "driver2.m"

7. "new5.m" "new10.m" "new20.m" "new40.m" are produced by gmsh, the number are the one in the right hand side of the equals sign  "Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 21;"  // 点1到7的所有曲线段上，每个线段生成 21 个等距的网格点，把每条线段分成 20 个网格单元

8. "error1.m" calculates error in"driver.m";    "error2.m"  calculates error in"driver2.m"
