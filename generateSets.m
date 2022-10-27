/*
generateSets.m

<description>

Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    18 October 2022
*/

// RecursiveSets.m

/*
INPUT:
    S:= [2,p_1,...,p_{i-1}], a list of primes such that 2*p_1*...*p_{i-1} <= m
    nSets:= [S_1,...,S_{k-1}], where
        S_j:= [2,p_1,...,p_{n_j}] is a list of primes such that 2*p_1*...*p_{n_j} <= m
    i:= index of the set S of prime to be added
    n:= the number of elements of the largest possible set S_j, where
        S_j:= [2,p_1,...,p_{n-1}] and 2*p_1*...*p_{n-1} <= m
    m:= upper bound on the product of the primes of each set S_j of nSets

OUTPUT:
    S:= [2,p_1,...,p_{i-1}, p_i], a list of primes such that 2*p_1*...*p_i <= m
    nSets:= [S_1,...,S_{k-1}, S_k], where
        S_i:= [2,p_1,...,p_{n_i}] is a list of primes such that
            2*p_1*...*p_{n_i} <= m and S_i is not a subset of any other set S_j of nSets

COMMENTS:
    This is a recursive algorithm. The procedure re-enters itself until i > n, when the output is finally generated.

EXAMPLE:
     S:= [2,3,5,11];
     nSets:= [];
     i:= 5;
     n:= 7;
     m:= 10^6;
     RecursiveSets(~S,~nSets,i,n,m);
    > nSets;
    [
        [ 2, 3, 5, 11, 13, 17 ],

*/


/*
      Determines the prime divisors of a conductor N.

      Parameters
          N: <param>: <param type>
	      <param description>
      Returns
          <param>: <param type>
	      The set X as a string without whitespace.
*/
SetOutputFile(dir cat "/" cat N cat "XYZ2tmp.txt");
load "./Code/parseIO.m";

OutFile:=dir cat "/" cat N cat "XYZ2Forms.csv";
sN:=N;
N:=StringToInteger(N);
S:=PrimeDivisors(N);
if 2 notin S then
    Insert(~S,1,2);
end if;
fprintf OutFile,"%o,%o\n",sN,seqEnumToString(S);
print "Data for N = " cat sN cat " written to " cat OutFile;
exit;
