// LatticeBoundXYZ.m

/*
INPUT:
    S:= {p_1,...,p_s}, p_i primes in S

OUTPUT:
    The bound C8 such that m(xyz) <= C8, based on Theorem 6.1 of the Reference
        m(xyz):= max(ord_{p_i}(xzy)) where max runs over all primes p_i in S, and x,y,z are relatively prime
        ie. bounds the largest possible exponent of the p_i appearing in x,y, or z, where x + y = z

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7];
    > LatticeBoundXYZ(S); 
    1037919720087088125.94830455465
    
    > S:= [2,7,11];
    > LatticeBoundXYZ(S);
    1719613135665434941.04287726437
    
*/


function LatticeBoundXYZ(S)
    s:= #S;     // computes the size of S
    Sort(~S);   // orders primes by size p_1 < ... < p_s
    t:= Floor(2*s/3);
    P:= &*S;    // computes the product of the primes p_i in S
    
    q:= [];     // stores smallest primes Q such that Q does not divide divide p_i(p_i - 1), for each p_i in S
    for i in [1..s] do
        Q:= 3;
        while (S[i]*(S[i]-1)) mod Q eq 0 do
            Q:= NextPrime(Q);   // computes smallest prime Q such that Q does not divide p_i(p_i - 1)
        end while;
        Append(~q,Q);
    end for;
    q:= Max(q);
    
    C12n:= [768523, 476217, 373024, 318871, 284931, 261379, 2770008];   // C1(2,n) as in Lemma 2.6 CAN MAKE BETTER - THERE HAS SINCE BEEN IMPROVEMENTS ON THIS RESULT
    if t lt 8 then
        C12t:= C12n[t-1];       // C1(2,t)
        a1:= 56*Exp(1)/15;      // a1 as in Lemma 2.6, if t <= 7 (note 2 <= t always)
    else
        C12t:= C12n[7];
        a1:= 8*Exp(1)/3;        // a1 as in Lemma 2.6, if t>7
    end if;
    M:= Max([((p-1)*(2 + (1/(p-1)))^t)/((Log(p))^(t+2)) : p in S]);     // computes max_i as in the definition of U
    U:= (C12t)*(a1^t)*(t^(t+5/2))*(q^(2*t))*(q-1)*((Log(t*q))^2)*(M)*((Log(S[s]))^t)*(Log(4*Log(S[s])) + (Log(S[s]))/(8*t));
    C1:= U/(6*t);
    C2:= U*Log(4);
    
    V:= [];
    for i in [s-t+1..s] do
        vi:= Max([1,Log(S[i])]);        // computes V_i = max(1, log(p_i)) for i = s-t+1,...,s
        Append(~V,vi);
    end for;                                  
    Omega:= &*V;     // Omega
    
    C3:= (2^(9*t + 26))*(t^(t+4))*(Omega)*(Log(Exp(1)*V[#(V)-1]));            // computes C_i as in Theorem 6.1
    C4:= Max([7.4, (C1*Log(P/S[1])+C3)/Log(S[1])]);
    C5:= ((C2)*(Log(P/S[1]))+(C3)*(Log(Exp(1)*V[#V])) + 0.327)/Log(S[1]);
    C6:= Max([C5, ((C2)*Log(P/S[1])+Log(2))/Log(S[1])]);
    C7:= 2*(C6 + C4*Log(C4));
    C8:= Max([S[s], Log(2*(P/S[1])^(S[s]))/Log(S[1]), C2 + C1*Log(C7), C7]);
    return C8;
end function;    
