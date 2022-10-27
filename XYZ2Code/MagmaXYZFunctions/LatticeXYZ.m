// LatticeXYZ.m

/*
INPUT:
    S:= {p_1,...,p_s}, p_i primes in S
    p:= prime in S
    mu:= the precision on the p-adic field Q_p

OUTPUT:
    if p >= 5 and p0 is a primitive root (mod p):
        B:= [b_1,...,b_{s-2}, b_0], a matrix over ZZ with columns b_i,
            b_0,b_1,...,b_{s-2} form a basis for the reduced p-adic approximation lattice Gamma_mu* 

    if (p < 5) or (p >= 5 and p0 is not a primitive root (mod p)):
        B:= [b_1,...,b_{s-1}, b_0], a matrix over ZZ with columns b_i,
            b_0,b_1,...,b_{s-1} form a basis for the p-adic approximation lattice Gamma_mu
                            
    m0:= ord_p(log_p(p0))
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7];
    > p:= 2;
    > mu:= 20;
    > B,m0:= LatticeXYZ(S,p,mu);
    > B;
    [      1       0       0]
    [      0       1       0]
    [ 508661  260946 1048576] 
    > m0;
    2
    
    > S:= [2,7,11,13];
    > p:= 7;
    > mu:= 25;
    > B,m0:= LatticeXYZ(S,p,mu);
    > B;
    [                     1                      0                      0]
    [                     0                      1                      0]
    [ 654610875222193089255  353482029863608455592 1341068619663964900807]
    > m0;
    1
    
*/


function LatticeXYZ(S,p,mu)

    s:= #S;             // computes the length of S
    A,p0,m0,LatticeType := InitialSortXYZ(S,p,mu); 
    R:= pAdicField(p,mu+5);
    
    t:= [-pAdicLog(R!i,p)/pAdicLog(R!p0,p) : i in A];   // computes theta_i = -log_p(p_i)/log_p(p0) for i=1,...,s-2
    tm:= [ConvertpAdic(i,p,mu) : i in t];               // computes theta_i^{mu} = theta_i (mod p^mu) for i=1,...,s-2
    Append(~tm,p^mu);   // adjoins p^mu to tm
    
    row_tm := Matrix(IntegerRing(),1,s-1,tm);   // writes theta_i^{mu} values as a row matrix
    I:= ScalarMatrix(s-2,1);                    // (s-2)-identity matrix
    Z:= ZeroMatrix(IntegerRing(),s-2,1);        // creates (s-2)x1 zero matrix
    
    B := VerticalJoin(HorizontalJoin(I,Z),row_tm);      // creates matrix associated to basis of Gamma_{mu}, where B is an (s-1)x(s-1) matrix with b_0 = B[s-1]
                                                    
    // if (p < 5)  or (p >= 5 and p0 is not a primitive root (mod p)): Gamma_{mu}* = Gamma_{mu}, hence B does not change, else:
    if (p ge 5) and (LatticeType eq "reduced") then     // if p >= 5 and LatticeType == "reduced", compute matrix Gamma_{mu}* associated reduced basis
        alpha := [];
        for i in A do
            for j in [1..p-1] do
                if i mod p eq Modexp(p0,j,p) then   // for i = 1,...s-2, computes alpha_i, where p_i = p0^{alpha_i} mod p
                    Append(~alpha,j);   // stores alpha_i values
                    break;
                end if;                                
            end for;
        end for;    
        
        gamma0 := (p-1)/2;      // computes gamma_0* = (p-1)/2
        gamma0 := IntegerRing() ! gamma0;   // converts gamma0* into an integer
                                                         
        for i in [1..s-2] do
            gamma := (-alpha[i] - tm[i]) mod gamma0;    // computes -gamma_i* = -alpha_i - theta_i^{mu} mod (p-1)/2
            AddColumn(~B,gamma,s-1,i);  // computes b_i* = b_i + (-gamma_i*)b_0
        end for;
        
        MultiplyColumn(~B,gamma0,s-1);  // computes b_0* = gamma_0*b_0     
    end if;
   
    return B,m0;
end function;