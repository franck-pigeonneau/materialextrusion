/* Déclaration des variables globales pour l'algorithme BDF */
const Float ABDF[7][7]={{     0.0,0., 0.0,   0.0,  0.0,  0.0,   0.0},
			{     1.0,1., 0.0,   0.0,  0.0,  0.0,   0.0},
			{     1.5,2.,-0.5,   0.0,  0.0,  0.0,   0.0},
			{  11./6.,3.,-1.5, 1./3.,  0.0,  0.0,   0.0},
			{ 25./12.,4.,-3.0, 4./3.,-0.25,  0.0,   0.0},
			{137./60.,5.,-5.0,10./3.,-1.25,1./5.,   0.0},
			{ 49./20.,6.,-7.5,20./3.,-3.75,6./5.,-1./6.}};

const Float BBDF[7][7]={{0.,0.,  0., 0.,  0.,0., 0.},
			{0.,1.,  0., 0.,  0.,0., 0.},
			{0.,2., -1., 0.,  0.,0., 0.},
			{0.,3., -3., 1.,  0.,0., 0.},
			{0.,4., -6., 4., -1.,0., 0.},
			{0.,5.,-10.,10., -5.,1., 0.},
			{0.,6.,-15.,20.,-15.,6.,-1.}};
