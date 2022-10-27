// SmallestLatticeBoundXYZ.m

/*
INPUT:
    S:= [p_1,...,p_s], p_i primes in S

OUTPUT:
    X0:= mu + m0 - 1, smallest upper bound such that m(xyz) <= X0, where
        mu:= smallest lattice invariant such that m(xyz) <= X0
        m0:= ord_p(log_p(p0)), a fixed value
        m(xyz):= max(ord_{p_i}(xzy)) where max runs over all primes p_i in S, and x,y,z are relatively prime
            ie. bounds the largest possible exponent of the p_i appearing in x,y, or z, where x + y = z
    U:= [[u_1,p_1],...,[u_s,p_s]], where
        u_i:= mu + m0 - 1, and
        ord_{p_i}(xyz)<= u_i for p_i in S
    
COMMENTS:
    Lemma 3.14 (de Weger, Algorithms for diophantine equations) in this instance:
        If (c_1 = 0, c_2 = 1) l(Gamma_{mu}*) > sqrt(s-1)*X1, then x + y = z has no solutions with
            mu + m0 <= ord_p(z) <= m(xyz) <= X1
    Hence all solutions of x + y = z must have ord_p(xyz) <= m(xyz):= max(ord_{p_i}(xyz)) < mu + m0
    This result is applied successively for each p to reduce the bound to the smallest such bound X0 = mu + m0 - 1

    Lemma 3.4 (Lenstra, Lenstra and Lovasz):
        Let c_1,...,c_n be a reduced basis of the lattice Gamma. Then
            l(Gamma) >= 2^((1-n)/2)*|c_1|
    This result is used to compute l(Gamma_{mu}*)  

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
    
*/


function SmallestLatticeBoundXYZ(S)
    s:= #S;     // computes the size of S
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    X0:= LatticeBoundXYZ(S);    // computes initial upper bound
    const:= 2^(-(s-2)/2);       // computes constant such that l(Gamma_{mu}*) >= const*|c_1|      
    X:= [];     // stores X0 values such that l > sqrt(s-1)*X0

    fixX0:= true;
    while fixX0 eq true do
        U:= [];         // stores u_i = mu + m0 - 1 values

        for p in S do   // runs through each prime p in S, in the order provided
            mu:= Ceiling((s*Log(X0))/Log(p));   // estimates initial mu value for the given prime, p
            
            changeC:= true;
            while changeC eq true do
                B,m0 := LatticeXYZ(S,p,mu);     // computes the lattice Gamma_{mu}* associated to (S,p,mu) and the associated m0 value

                if s eq 3 then
                    l:= F3Approx(B);            // computes short lattice length based on F3Approx algorithm when S = {p_1,p_2,p_3} only
                else
                    C:= LLL(Transpose(B));      // for s > 3, computes lattice lenght l(Gamma_{mu}*) via LLL algorithm
                    x:= C[1];
                    l:= const*Sqrt(Norm(x));    // computes |c_1|, the norm of the first column of the reduced lattice
                end if;
                                                                                                                
                if l le Sqrt(s-1)*X0 then
                    mu:= mu + 5;        // changeC still True, returns to start of changeC with updated mu value, repeats until l = l(Gamma_{mu}*) > sqrt(s-1)*X0
                else
                    u:= mu + m0 - 1;    // application of Lemma 3.14: no solutions with mu + m0 <= ord_p(z), hence ord_p(z) <= mu + m0 + 1
                                        // invariant under permutation, so ord_p(x),ord_p(y) < mu + m0 also - ie. ord_p(xyz) <= mu + m0 - 1

                    Append(~U,[u,p]);   // stores u_i = mu + m0 - 1 values; p_i, for p_i, such that p_1 < p_2 <...< p_s
                    changeC:= false;    // exits changeC algorithm, continues to next value of p in S
                end if;
             end while;
        end for;        
        
        X0:= Max([u[1]: u in U]);       // selects largest u_i value, X0, so ord_p(xyz) <= X0 holds for all primes p in S
        Append(~X,X0);  // stores all X0 values
        if &+[1: x in X | x eq X0] eq 2 then    // computes number of times X0 appears in X; terminates algorithm if X0 twice - ie. lowest upper bound cannot be further improved
            fixX0 := false;
        end if;
    end while;
    
    return X0,U;
end function;                    