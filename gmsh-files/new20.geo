// Gmsh project created on Mon Jan 13 01:29:31 2025
R = 0.5;   				//定义变量值
L = 2;

Point(1) = {L, 0, 0};  		//指定点：右下角
Point(2) = {L, L, 0};   		//右上角
Point(3) = {0, L, 0}; 		//左上角
Point(4) = {0, 0, 0};		//左下角
Point(5) = {R, 0, 0};	//小正方形右下角
Point(6) = {0, R, 0};	//小正方形左上角
Point(7) = {Cos(Pi/4) * R, Sin(Pi/4) * R, 0};  //小正方形的右上角

Circle(1) = {5, 4, 7};		//指定圆弧：点4为圆心，点5画到点7
Circle(2) = {7, 4, 6};

Line(3) = {6, 3};			//指定线：点6-点3
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

Curve Loop(1) = {-3,-2,-7,-4}; 		//画封闭回路，line4-line7-circle2-line3
Plane Surface(1) = {1};  		//使用curve loop1 来生成一个面

Curve Loop(2) = {7, -1, -6, -5}; 	//画封闭回路，line7-circle1逆向-circle6逆向-line5逆向
Plane Surface(2) = {2};		   	//使用curve loop2 来生成一个面

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 21;  // 点1到7的所有曲线段上，每个线段生成21个等距的网格点，把每条线段分成 20 个网格单元

Transfinite Surface{1};  			//对surface1生成网格，以边上的节点数等距生成网格
Transfinite Surface{2};

Recombine Surface{1};			//*把三角形网格合成为的四边形网格*
Recombine Surface{2};

Mesh.ElementOrder = 1;		//generate linear mesh (1st-order elements)
Mesh.Algorithm = 8;			//use the Frontal-Delaunay algorithm to generate unstructured meshes, which is suitable for ...								//complex geometries, especially in 2D ( search through the Internet )

//explain the meaning of the geo file：这个文件先定义了图形的边界和形状（从点到面），然后等距生成网格，再把网格进行了优化

// EOF

//+
Physical Curve("up", 8) = {4};
//+
Physical Curve("down", 9) = {6};
//+
Physical Curve("left", 10) = {3};
//+
Physical Curve("right", 11) = {5};
//+
Physical Curve("arc", 12) = {2};
//+
Physical Curve("arc", 12) += {1};
//+
Physical Surface("blue", 13) = {1, 2};

