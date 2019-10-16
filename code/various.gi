invertMatrixEmbedding := function(A,m)

## Input: 2nx2n-matrix A over the natural integers,

## we assume m congruent 1 or 2 mod 4, so we have a ring of integers Z[w] with w = Sqrt(-m).
## Output: Its preimage A_original in GL_n(Z[w]).
local w, n, A_original, a_upper, a_lower, b_left, b_right, A11, A12, A21, A22, j, k;

if m=0 then return A;fi;
n := Size(A)/2;
A_original := NullMat(n,n);
for j in [1..n] do
  for k  in [1..n] do
    if m mod 4 = 3 then
    w := Sqrt(-m)/2 +1/2;
    A11 := A[2*j-1][2*k-1];
    A22 := A[2*j][2*k];
    A12 :=  A[2*j-1][2*k];
    A21 :=  A[2*j][2*k-1];
    if A12*(m+1)/(-4) = A21 then
      if A11 = A22-A12 then
        A_original[j][k] := A11 +A12*w;
      else
        Print("Error on identity block");
      fi;
    else
        Print("Error on imaginary block: A12 = ",A12,", A21 = ",A21);
    fi;
    else
    w := Sqrt(-m);
    ## Local inversion code in the case m congruent to 1 or 2 mod 4,
    ## w^2 = -m :
    ## Compute back a+b*w from a 2x2-matrix block,
    ## where 1 is represented in GL_2(Z) as the 2x2-identity matrix,
     ## and w is represented as [[0,1],[-m,0]].
    a_upper := A[2*j-1][2*k-1];
    a_lower := A[2*j][2*k];
    if a_upper = a_lower then
        A_original[j][k] := a_upper;
    #Print("A_",j,",",k," = ",a_upper,", ");
    else
        Print("Error on identity block");
    fi;
    b_right :=  A[2*j-1][2*k];
    b_left  :=  A[2*j][2*k-1]/(-m);
    if b_left = b_right then
        A_original[j][k] := a_upper+b_left*w;
    #Print("A_",j,",",k," = ",a_upper,"+",b_left,"*w, ");
    else
        Print("Error on imaginary block: b_right = ",b_right,", b_left = ",b_left);
    fi;
    fi;
  od;
od;
return A_original;
end;; 
##################################################
matrixEmbedding:=function(A,m)
## Input: nxn-matrix A in GL_n(Z[w])

## we assume m congruent 1 or 2 mod 4, so we have a ring of integers Z[w] with w = Sqrt(-m).
## Output: Its image A_emb in GL_2xn(Z).

local i,j, A_emb, w, re, im;
    w:=Sqrt(-m);
    n:=Size(A);
    A_emb:=NullMat(2*n,2*n);
    for i in [1..n] do
	for j in [1..n] do
	    re:=RealPart(A[i][j]);
	    im:=ImaginaryPart(A[i][j]);
	    A_emb[2*i-1][2*j-1]:=re;
            A_emb[2*i][2*j]:=re;
	    A_emb[2*i-1][2*j]:=im*m;
	    A_emb[2*i][2*j-1]:=-im;
	od;
    od;

    return A_emb;
end;
##################################################
Mult:=function(C,g,h)
local p;
p:=Position(C!.elts,C!.elts[g]*C!.elts[h]);
if p=fail then Add(C!.elts,C!.elts[g]*C!.elts[h]); return Length(C!.elts);fi;

return p;

end;
##############################################################
