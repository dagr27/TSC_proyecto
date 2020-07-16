float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

float rIntegralMK(mesh m, int i){
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();
    
    float y1, y2, y3, y4;
    z1 = n1.getY();
    z2 = n2.getY();
    z3 = n3.getY();
    z4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();
    
    return ((5*y2*(3*pow(z2,2)+2*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2))+5*y3*(pow(z2,2)+z2*(2*z3+z4)+3*pow(z3,2)+2*z3*z4+pow(z4,2))+5*y4*(pow(z2,2)+z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2))+3*(pow(z2,2)+z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2)))/(180));
}


void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}

void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

/*Matriz A*/
void calculateGammaForA(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();
    
    float y1, y2, y3, y4;
    z1 = n1.getY();
    z2 = n2.getY();
    z3 = n3.getY();
    z4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    float a = ((8*x2*(3*pow(z2,2)+2*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2))+8*x3*(pow(z2,2)+z2*(2*z3+z4)+3*pow(z3,2)+2*z3*z4+pow(z4,2))+8*x4*(pow(z2,2)+z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2))+49*(pow(z2,2)+z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2)))/(2520));
    float b = ((16*x2*(6*pow(z2,2)+3*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2))+8*x3*(3*pow(z2,2) +2*z2*(2*z3+z4)+3*pow(z3,2) +2*z3*z4+pow(z4,2))+8*x4*(3*pow(z2,2)+2*z2*(z3+2*z4)+pow(z3,2) +2*z3*z4+3*pow(z4,2))+49*(3*pow(z2,2)+2*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,2)))/(2520));
    float c = ((8*x2*(3*pow(z2,2)+2*z2*(2*z3+z4)+3*pow(z3,2)+2*z3*z4+pow(z4,2))+16*x3*(pow(z2,2)+z2*(3*z3+z4)+6*pow(z3,2)+3*z3*z4+pow(z2,2))+8*x4*(pow(z2,2)+2*z2*(z3+z4)+3*pow(z4,2)+4*z3*z4+3*pow(z4,2))+49*(pow(z2,2)+z2*(2*z3+z4)+3*pow(z4,2)+2*z3*z4+pow(z4,2)))/(2520));
    float d = ((8*x2*(3*pow(z2,2)+2*z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2))+8*x3*(pow(z2,2)+2*z2*(z3+z4)+3*pow(z3,2)+4*z3*z4+3*pow(z4,2))+16*x4*(pow(z2,2)+z2*(z3+3*z4)+pow(z3,2)+3*z3*z4+6*pow(z4,2))+49*(pow(z4,2)+z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2)))/(2520));

	G1.at(0).at(0) = a;   	                                G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = b;                                     G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = c;                                     G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = d;                                     G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = a;   	                                G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = b;                                     G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = c; 	                                G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = d;                                     G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = a; 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = b; 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = c;
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = d; 
}

/*Matriz G*/
void calculateGammaForG(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();
    
    float y1, y2, y3, y4;
    z1 = n1.getY();
    z2 = n2.getY();
    z3 = n3.getY();
    z4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    float a = ((x2*(3*pow(z2,2)+2*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z2,2))+x3*(pow(z2,2)+z2*(2*z3+z4)+3*pow(z3,2)+2*z3*z4+pow(z4,2))+x4*(pow(z2,2)+z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2)))/(630));
    float b = ((2*x2*(6*pow(z2,2)+3*z2*(z3+z4)+pow(z3,2)+z3*z4+pow(z4,3))+x3*(3*pow(z2,2)+2*z2*(2*z3+z4)+3*pow(z3,2)+2*z3*z4+pow(z4,2))+x4*(3*pow(z2,2)+2*z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2)))/(630));
    float c = ((x2*(3*pow(z2,2)+2*z2*(2*z3+z4)+3*pow(z3,2)+2*z3*z4+pow(z4,2))+2*x3*(pow(z2,2)+z2*(3*z3+z4)+6*pow(z3,2)+3*z3*z4+pow(z4,2))+x4*(pow(z2,2)+2*z2*(z3+z4)+3*pow(z3,2)+4*z3*z4+3*pow(z4,2)))/(630));
    float d = ((x2*(3*pow(z2,2)+2*z2*(z3+2*z4)+pow(z3,2)+2*z3*z4+3*pow(z4,2))+x3*(pow(z2,2) +2*z2*(z3+z4)+3*pow(z3,2)+4*z3*z4+3*pow(z4,2))+2*x4*(pow(z3,2)+z2*(z3+3*z4)+pow(z3,2)+3*z3*z4+6*pow(z4,2)))/(630));

	G1.at(0).at(0) = a;   	                                G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = b;                                     G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = c;                                     G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = d;                                     G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = a;   	                                G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = b;                                     G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = c; 	                                G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = d;                                     G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = a; 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = b; 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = c;
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = d; 
}

/*Matriz D*/
void calculateGammaForD(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();
    
    float y1, y2, y3, y4;
    z1 = n1.getY();
    z2 = n2.getY();
    z3 = n3.getY();
    z4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    float a = ((y2*(4*pow(z2,3)+3*pow(z2,2)*(z3+z4)+2*z2*(pow(z3,2)+z3*z4+pow(z4,2))+pow(z3,3)+pow(z3,2)*z4+z3*pow(z4,2)+pow(z4,3))+y3*(pow(z2,3)+pow(z2,2)*(2*z3+z4)+z2*(3*pow(z3,2)+2*z3*z4+pow(z4,2))+4*pow(z3,2)+3*pow(z3,2)*z4+2*z3*pow(z4,2)+pow(z4,3))+y4*(pow(z2,3)+pow(z2,3)*(z3+2*z4)+z2*(pow(z3,2)+2*z3*z4+3*pow(z4,2))+pow(z3,2)+2*pow(z3,2)*z4+3*z3*pow(z4,2)+4*pow(z4,3)))/(6720));
    float b = ((2*y2*(10*pow(z2,3)+6* pow(z2,2)*(z3+z4)+3*z2*( pow(z3,2)+z3*z4+ pow(z4,2))+ pow(z3,3)+ pow(z3,2) *z4+z3*pow(z4,2)+pow(z2,3))+y3*(4*pow(z2,2)+3*pow(z2,2)*(2*z3+z4)+2*z2*(3*pow(z3,2)+2*z3*z4+pow(z4,2))+4*pow(z3,3)+3*pow(z3,2)*z4+2*z3*pow(z4,2)+pow(z4,3))+y4*(4*pow(z2,3)+3*pow(z2,2)*(z3+2*z4)+2*z2*(pow(z3,2)+2*z3*z4+3*pow(z4,2))+pow(z3,3)+2*pow(z3,2)*z4+3*z3*pow(z4,2)+4*pow(z4,3)))/(6720));
    float c = ((y2*(4*pow(z2,3)+3*pow(z2,2)*(2*z3+z4)+2*z2*(3*pow(z3,2)+2*z3*z4+pow(z4,2))+4*pow(z3,3)+3*pow(z3,2)*z4+2*z3*pow(z4,2)+pow(z4,3))+2*y3*(pow(z2,3)+pow(z2,2)*(3*z3+z4)+z2*(6*pow(z3,2)+3*z3*z4+pow(z4,2))+10*pow(z3,3)+6*pow(z3,2)*z4+3*z3*pow(z4,2)+pow(z4,3))+y4*(pow(z2,3)+2*pow(z2,2)*(z3+z4)+z2*(3*pow(z3,2)+4*z3*z4+3*pow(z4,2))+2*(2*pow(z3,3)+3*pow(z3,2)*z4+3*z3*pow(z4,2)+2*pow(z4,3))))/(6720));
    float d = ((y2*(4*pow(z2,3)+3*pow(z2,2)*(z3+2*z4)+2*z2*(pow(z3,2)+2*z3*z4+3*pow(z4,2))+pow(z3,3)+2*pow(z3,2)*z4+3*z3*pow(z4,2)+4*pow(z4,3))+y3*(pow(z2,3)+2*pow(z2,2)*(z3+z4)+z2*(3*pow(z3,2)+4*z3*z4+3*pow(z4,2))+2*(2*pow(z3,3)+3*pow(z3,2)*z4+3*z3*pow(z4,2)+2*pow(z4,3)))+2*y4*(pow(z2,3)+pow(z2,2)*(z3+3*z4)+z2*(pow(z3,2)+3*z3*z4+6*pow(z4,2))+pow(z3,3)+3*pow(z3,2)*z4+6*z3*pow(z4,2)+10*pow(z4,3)))/(6720));

	G1.at(0).at(0) = a;   	                                G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = b;                                     G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = c;                                     G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = d;                                     G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = a;   	                                G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = b;                                     G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = c; 	                                G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = d;                                     G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = a; 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = b; 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = c;
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = d; 
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixK,matrixG,matrixD,matrixC,matrixE,matrixW;
    float u_bar,nu,rho,Volumen,J,Determinant;
    /*
        | A+K  P |
        |  D   0 |
    */

    //Matrix A
    Matrix g_matrix, Alpha, Beta, gA_matrix,gG_matrix,gD_matrix;

    //u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    //Esto es para que me escupa la matriz A
    //Escalar que multiplica las componentes de A
    float real_A = (float) (J)/(Determinant);
    

    calculateGammaForA(gA_matrix,m,e);
    calculateGammaForG(gG_matrix,m,e);
    calculateGammaForD(gD_matrix,m,e);    

    /*Calculando y creando las matrices "Griegas"*/
    calculateGamma(g_matrix);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);

    productRealMatrix(real_A, productMatrixMatrix(gA_matrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);
    //Genera matrixA
    
    


    //Esto es para que me calcule la Matriz K
    Matrix Alpha_t,Beta_t;
    //EL volumen se vuelve un escalar, por lo tanto lo podemos poner en el numerador
    Volumen = rIntegralMK(m,e);
    float real_k = (float) (Volumen)/(Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);
   

	//Esta es la que me devuelve la matriz G
	Matrix Omega;
	calculateOmega(Omega);
	float real_G = (float) (J)/(Determinant);
	 
	productRealMatrix(real_G, productMatrixMatrix(gG_matrix,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixG);
    

    //Esta es la que nos da la matriz D
    Matrix g_matrix_t,Omega_t;
    float real_D = (float)(J/(Determinant));

    transpose(Omega, Omega_t);
    //Transponiendo mi matriz Resultante de componentes
    transpose(gD_matrix,g_matrix_t);
    //Efectuando producto de matrices
    productRealMatrix(real_D,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);
    

    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixG);
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

//Vamos a sacar el vector columna F
void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix;

    calculateF(f, m);

    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,24);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    
    return b;
}
