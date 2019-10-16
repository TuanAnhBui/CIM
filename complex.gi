IsInGamma1:=function(x,prime)
;

    if not EuclideanRemainder(x[1][1]-1,prime)=0 then return false;fi;
    if not EuclideanRemainder(x[2][1],prime)=0 then return false;fi;
    if not EuclideanRemainder(x[3][1],prime)=0 then return false;fi;

    return true;
end;

IsInGamma2:=function(x,prime)

return IsInGamma1(x,prime) and (EuclideanRemainder(x[3][2],prime)=0);
end;
###################################################
reps1:=[
IdentityMat(3)^-1,
[[1,0,0],[1,1,0],[0,0,1]]^-1,
[[1,0,0],[1,1,0],[1,0,1]]^-1,
[[1,0,0],[0,1,0],[1,0,1]]^-1,
[[0,0,1],[1,0,0],[0,1,0]]^-1,
[[0,0,1],[1,0,0],[1,1,0]]^-1,
[[0,1,0],[0,0,1],[1,0,0]]^-1
];;

IsInG1:=function(x); return IsInGamma1(x,2);end;

reps2:=[
IdentityMat(3)^-1,
[[1,0,0],[0,0,-1],[0,1,0]]^-1,
[[1,0,0],[0,1,0],[0,1,1]]^-1
];;

IsInG2:=function(x); return IsInGamma2(x,2);end;
##################################################
C:=ContractibleGcomplex("SL3Zh");;
eq:=InducedMapFromComplexOfFiniteIndexSubgroup(C,reps1,IsInG1);
eq2:=InducedMapFromComplexOfFiniteIndexSubgroup(eq!.source,reps2,IsInG2);
M:=ModPCellularCohomologyWithSage(eq2,2,3);