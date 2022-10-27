// InitialSortXYZ.m

/*
INPUT:
    S:= {p_1,...,p_s}, p_i primes in S
    p:= prime in S
    mu:= the precision on the p-adic field Q_p

OUTPUT:
    if S\{p} contains a primitive root (mod p):
        p0:= prime of S such that ord_p(log_p(p0)) is minimal and a primitive root (mod p)
        LatticeType:= "reduced", indicates that a reduced lattice Gamma_mu* may be used in this case

    if S\{p} does not contain a primitive root (mod p):
        p0:= prime of S such that ord_p(log_p(p0)) is minimal
        LatticeType:= "-", indicates that the original lattice Gamma_mu must be used in place of Gamma_mu*
        
    m0:= ord_p(log_p(p0))
    A:= the primes S\{p,p0}, sorted
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7];
    > p:= 5;
    > mu:= 20;
    > A,p0,m0,LatticeType:= InitialSortXYZ(S,p,mu);
    > A;
    [ 3, 7 ]
    > p0;
    2
    > m0;
    1
    > LatticeType;
    reduced
    
    > S:= [2,7,11,13];
    > p:= 7;
    > mu:= 25;
    A,p0,m0,LatticeType:= InitialSortXYZ(S,p,mu);
    > A;
    [ 11, 13 ]
    > p0;
    2
    > m0;
    1
    > LatticeType;
    -
        
*/


function InitialSortXYZ(S,p,mu)

    s:= #S;             // computes the length of Sp

    R:= pAdicField(p,mu+5);
    A:= Exclude(S,p);           // removes p from S to rename the primes
    
    P0:= [Valuation(pAdicLog(R!i,p)) : i in A];         // computes ord_p(log_p(p_i)) for each prime p_i in S as a p-adic integer
    LatticeType:= "-";          // sets default LatticeType to "-", unless otherwise specified
    
    if (p ge 5) then
        if (PrimitiveRoot(p) in A) and P0[Index(A,PrimitiveRoot(p))] eq Min(P0) then    // if A contains a primitive root mod p and if index of primitive root mod p corresponds to minimal ord_p(log_p(p_i))
            p0index := Index(A,PrimitiveRoot(p));       // computes index of primiitive root mod p in S
            LatticeType:= "reduced";
            p0:= A[p0index];    // rename p0 such that ord_p(log_p(p0)) is minimal and p0 is a primitive root mod p
        end if;
    end if;
    
    if (p lt 5) or (LatticeType eq "-") then
        p0index0,p0index:= Min(P0);     // finds index of minimal ord_p(log_p(p_i)), p0index0; and it's location in index in P0
        p0:= A[p0index];
    end if;
    
    m0:= Min(P0);       // m0 = ord_p(log_p(p0))
    Exclude(~A,p0);
    Sort(~A);           // reorder prime set S by size of primes to include only p_1,...,p_{s-2}

    return A,p0,m0,LatticeType;
    
end function;   