// FPRestrictionsXYZ.m

/*
INPUT:
    S:= [p_1,...,p_s], p_i primes in S
    U:= [[f_1,p_1],...,[f_s,p_s]] such that 
        ord_{p_i}(xyz) <= f_i for each p_i in S, where x + y = z 
    Order:= [p_j,..,p_m], list of primes of S to be assessed
    C:= a constant 0 < C < 1 such that
        m:= f_i - m0 - ceil(C*f_i), where
            f_i:= element of [f_i,p_i] of U
            m:= the precision on the {p_i}-adic field Q_{p_i}), lattice invariant
            m0:= ord_{p_i}(log_{p_i}(p0)), lattice invariant  

OUTPUT:
    F:= [[f_1,p_1],...,[f_s,p_s]], updated such that 
        ord_{p_i}(xyz) <= f_i for each p_i in S 
    eqns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = z_i such that 
        ord_p(xyz) is in the range [m+m0, f_p], for a given prime p in Order
    
COMMENTS:
    Computes all S-unit equations x + y = z involving ord_p(xyz) in the range [m+m0, f_p], for given prime p in S:
        For a given prime p, to find all solutions (x,y,z) coprime with ord_p(xyz) <= f_p:
            Choose m < f_p - m0 and consider the lattice Gamma_{m}* = LatticeXYZ(S,p,m)
            If a solution (x,y,z) exists with ord_p(z) in [m+m0,f_p], then the vector (x_1,...,x_{s-2},x_0)^T with x_i = ord_{p_i}(x/y) for i = 0,...,s-2 is in the lattice
            The length of this vector is bounded by (f_{p_0}^2 + f_{p_1}^2 + ... + f_{p_(s-2)}^2)
    This algorithm uses the Fincke-Pohst algorithm to find all vectors of length at most (f_{p_0}^2 + f{p_1}^2 + ... + f{p_(s-2)}^2) and sorts them to give (x,y,z) for each p in Order

    p:= prime in Order, a subset of S
    m:= lattice coefficient Gamma_{m} = lattice(S,p,m), where m < f_p - m0
    m0:= invariant of the lattice Gamma_{m}
    f_{p_i}:= positive integer associated to p_i in S such that ord_{p_i}(xyz) <= f_{p_i} at start of iteration

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,5,11];
    > X0,U:= SmallestLatticeBoundXYZ(S);
    > X0;
    13
    > U;
    [
        [ 13, 2 ],
        [ 5, 5 ],
        [ 4, 11 ]
    ]
    > Order:= S;                                     
    > C:= 0.5;                      
    > F,eqns:= FPRestrictionsXYZ(S,U,Order,C);
    > F;
    [
        [ 13, 2 ],
        [ 5, 5 ],
        [ 2, 11 ]
    ]
    > eqns;
    []
    > U:=F;
    > U;
    [
        [ 13, 2 ],
        [ 5, 5 ],
        [ 2, 11 ]
    ]
    > F,eqns:= FPRestrictionsXYZ(S,U,Order,C);
    > F;
    [
        [ 13, 2 ],
        [ 3, 5 ],
        [ 2, 11 ]
    ]
    > eqns;
    []
    > U:=F;
    > U;
    [
        [ 13, 2 ],
        [ 3, 5 ],
        [ 2, 11 ]
    ]
    > F,eqns:= FPRestrictionsXYZ(S,U,Order,C);
    > F;
    [
        [ 5, 2 ],
        [ 3, 5 ],
        [ 2, 11 ]
    ]
    > eqns;
    []
    > U:=F;
    > U;
    [
        [ 5, 2 ],
        [ 3, 5 ],
        [ 2, 11 ]
    ]
    > F,eqns:= FPRestrictionsXYZ(S,U,Order,C);
    > F;
    [
        [ 5, 2 ],
        [ 1, 5 ],
        [ 2, 11 ]
    ]
    > eqns;
    [
        [ 121, 4, 125 ]
    ]

*/


function FPRestrictionsXYZ(S,U,Order,C)
    s:= #S;     // computes the size of S
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    F:= U;      // creates copy of U, to avoid changing U
    Sbound:= 10^(100);  // sets bound on the number of solutions in Short Vectors algorithm; this bound is intentionally chosen to be as large as possible
    eqns:= [];  // stores S-unit solutions x + y = z with ord_p(xyz) in [m+m0,f_p]
    LatticeDets:= [];   // stores the determinant of each matrix; this matrix is associated to a corresponding lattice
                        // the relevant short vectors are found within an ellipsoid within this lattice; the volume of this ellipsoid is 
                        // (( Pi(RealField())^(n/2) )/( Gamma(n/2+1) ))*(bound^(n/2)/Det(B)), where B is an n x n matrix
                        // hence (heuristically) the smaller the volume of the ellipsoid, the less short vectors to be found

    for p in Order do
        Indexp:= Index(S,p);    // computes the index of p in S
        fp:= F[Indexp][1];      // computes f_p value associated to p, where ord_p(xyz) <= fp
        if fp ne 0 then
        
            primes,p0,m0,LT:= InitialSortXYZ(S,p,100);      // computes primes:= {p_1,...,p_{s-2}}; p0 such that ord_p(log_p(p0)) is minimal; m0:= ord_p(log_p(p0)); LT:= lattice type
            
            Append(~primes,p0);     // writes primes:= {p_1,...,p_{s-2},p_0}
            Fp:= [F[Index(S,i)][1] : i in primes];  // for each p_i in primes, writes Fp[i]:= f_{p_i}, where ord_{p_i}(xyz) <= Fp[i]
            bound:= &+[i^2 : i in Fp];      // computes bound associated to p, f_{p_0}^2 + f_{p_1}^2 + ... + f_{p_(s-2)}^2
            
            if fp - m0 gt 0 then    // fail-safe ensuring m is chosen so that final fp value remains positive, despite C value
                if fp gt 8 then
                    m:= fp - m0 - Ceiling(C*fp);    // for "large" (i.e. > 8) values of fp, chooses m < fp - m0 according to input value, C
                else
                    m:= fp - m0 - 1;        // for "small" (i.e. <= 8) values of fp, decreases m only by 1
                end if;
            elif fp - m0 lt 0 then
                m:= 1;
            else
                m:= fp - m0;        // if fp - mu0 == 0, chooses m so that resulting fp = 0
            end if;
            
            B,m0:= LatticeXYZ(S,p,m);
            n:= NumberOfColumns(B);
            VolumeOfLattice:= (( Pi(RealField())^(n/2) )/( Gamma(n/2+1) ))*(bound^(n/2)/Determinant(B));        // computes the volume of the ellipsoid
            Append(~LatticeDets, < p, VolumeOfLattice > );       // stores the volume of the ellipsoid associated to the matrix B
        end if;
    end for;
    
    MinLatticeVol, IndexMinLatticeVol:= Min([d[2] : d in LatticeDets]);         // selects the prime p of Order such that the corresponding matrix B gives the smallest-volume ellipsoid
    
    p:= LatticeDets[IndexMinLatticeVol][1];
    Indexp:= Index(S,p);    // computes the index of p in S
    fp:= F[Indexp][1];      // computes f_p value associated to p, where ord_p(xyz) <= fp
    if (fp ne 0) then
    
        primes,p0,m0,LT:= InitialSortXYZ(S,p,100);      // computes primes:= {p_1,...,p_{s-2}}; p0 such that ord_p(log_p(p0)) is minimal; m0:= ord_p(log_p(p0)); LT:= lattice type
        
        Append(~primes,p0);     // writes primes:= {p_1,...,p_{s-2},p_0}
        Fp:= [F[Index(S,i)][1] : i in primes];  // for each p_i in primes, writes Fp[i]:= f_{p_i}, where ord_{p_i}(xyz) <= Fp[i]
        bound:= &+[i^2 : i in Fp];      // computes bound associated to p, f_{p_0}^2 + f_{p_1}^2 + ... + f_{p_(s-2)}^2
        
        if fp - m0 gt 0 then    // fail-safe ensuring m is chosen so that final fp value remains positive, despite C value
            if fp gt 8 then
                m:= fp - m0 - Ceiling(C*fp);    // for "large" (i.e. > 8) values of fp, chooses m < fp - m0 according to input value, C
            else
                m:= fp - m0 - 1;        // for "small" (i.e. <= 8) values of fp, decreases m only by 1
            end if;
        elif fp - m0 lt 0 then
            m:= 1;
        else
            m:= fp - m0;        // if fp - mu0 == 0, chooses m so that resulting fp = 0
        end if;
        
        cb:= [m + m0,fp];       // computes [m+m0,fp] for which ord_p(xyz) is required to be in
        B,m0:= LatticeXYZ(S,p,m);
        M,w,maxReached:= FinckePohst(B,bound,Sbound);   // computes the vectors of B whose norm is <= bound, enumerated using the Fincke-Pohst algorithm
                                                        // returns short lattice coefficient vectors, w, where Transpose(M*Transpose(w[i])) = ord_{p_i}(x/y)
                                                        // maxReached determines whether or not the number of solutions returned, #w, is less than SBound
        
        while maxReached eq true do
            Sbound:= Sbound*2;  // iterates FinckePosht code until the bound on the number of solutions in Short Vectors algorithm, Sbound, is larger than the number of solutions returned
                                // ie. until maxReached returns 'false'
            M,w,maxReached:= FinckePohst(B,bound,Sbound);       // computes the vectors of B whose norm is <= bound, enumerated using the Fincke-Pohst algorithm
                                                                // returns short lattice coefficient vectors, w, where Transpose(M*Transpose(w[i])) = ord_{p_i}(x/y)
                                                                // maxReached determines whether or not the number of solutions returned, #w, is less than SBound
        end while;
            
        if #w eq 0 then         
            soln:= [];  // if Fincke-Pohst returns zero vectors

        else
            soln:= [];  // if Fincke-Pohst returns => 1 vectors:

            j:= 1;
            while j le #w do    // checks each vector w[j] found in Fincke-Pohst, weeding them out accordingly
                w_j:= Transpose(M*Transpose(w[j]));
                if &and[Abs(w_j[1,i]) le Fp[i] : i in [1..s-1]] then   // checks if w_j[i]:= ord_{p_i}(w_j/y) <= Fp[i] for i = 1,...,s-1
                    
                    r:= &*[primes[i]^(w_j[1,i]) : i in [1..s-1]];      // computes x/y (or any other permutation of quotients of x,y,z)
                    rem1, s1:= SFactors(Numerator(r) + Denominator(r),S);       // computes rem1:= the remainder of (x/y) or (y/x), (ie. num = x,y, denom = x,y) after dividing by all primes of S
                                                                                   // computes s1:= [[ord_p(n),p] : p in S], where n = x/y or n = y/x
                    rem2, s2:= SFactors(Abs(Numerator(r) - Denominator(r)),S);  // may be > 0 or < 0; if negative: denom = z, num = x,y, representing (x/z) or (y/z)
                    
            // discards the vector w[j] if s1, s2 == 1 (since not in range [m+m0,f_p]), 
                // ie. z = 1, so ord_p(xyz) = 0
            // for each vector w[j], for s1 and s2 repectively, verifies that the following conditions hold:
                // 1. all factors of s1 (resp. s2) are in S
                    // ie. contains only primes of S
                    // verified by checking rem1 == 0 (resp. rem2 == 0)
                // 2. ord_p(s1) (resp. ord_p(s2)) is in the range [m+m0,f_p]
                    // hence also checks p|s1 (resp. p|s2)
                    // verified by checking s1[Indexp][1] in [m+m0..fp]
                // 3. ord_{p_i}(xyz) <= Fp[i] for each i = 0,...,s-2 where (p_i)|s1 (resp s2.)
                    // verified by checking &and[s1[i][1] le F[i][1] : i in [1..s]]
                    
                    if (Abs(rem1) eq 1) and (s1[Indexp][1] in [m+m0..fp]) then
                        if &and[s1[i][1] le F[i][1] : i in [1..s]] then
                            x:= Max(Numerator(r), Denominator(r));
                            y:= Min(Numerator(r), Denominator(r));
                            z:= Numerator(r) + Denominator(r);
                            Append(~soln,[x,y,z]);      //  stores solution num + denom = s1
                        end if;
                    end if;
                    
                    
                    if (Abs(rem2) eq 1) and (s2[Indexp][1] in [m+m0..fp]) then
                        if &and[s2[i][1] le F[i][1] : i in [1..s]] then
                            if (Numerator(r) - Denominator(r)) lt 0 then        // if s2 < 0
                                x:= Max(Numerator(r), Abs(Numerator(r) - Denominator(r)));
                                y:= Min(Numerator(r), Abs(Numerator(r) - Denominator(r)));
                                z:= Denominator(r);
                                Append(~soln,[x,y,z]);  // stores num + s2 = denom
                            else                                                // if s2 > 0
                                x:= Max(Denominator(r), Abs(Numerator(r) - Denominator(r))); 
                                y:= Min(Denominator(r), Abs(Numerator(r) - Denominator(r)));
                                z:= Numerator(r);
                                Append(~soln,[x,y,z]);  // stores s2 + denom = num
                            end if;
                        end if;
                    end if;
    
                end if;                                       
                j:= j + 1;
            end while;    
        end if;
        
        for sol in soln do
            if (sol in eqns) eq false then
                Append(~eqns, sol);     // stores each [x,y,z] in eqns vector
            end if;    
        end for;
        
        F[Indexp][1]:= m+m0-1;  // update f_p value, terminates algorithm
    end if;
    
    return F,eqns;
end function;
