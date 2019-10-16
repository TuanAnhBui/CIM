###################################################
OrbitTesselation:=function(S,Reps,CheckMembership)
local Stabs, Orbs,i,j, reps, g, orb, stab, x, w, len, reps_pos, elms, elm, id;
    reps:=StructuralCopy(Reps);
    Stabs:=[];
    Orbs:=[];
    len:=Length(Reps);
    reps_pos:=[];
    id:=One(S);
    elms:=[];
    while not len=0 do
	x:=reps[1];
        reps_pos[Position(Reps,x)]:=[Length(Orbs)+1,1];
	orb:=[x];
        stab:=[];
        elm:=[id];
	Remove(reps,1);
        for g in S do
	    j:=Length(reps);
            if CheckMembership(x*g^-1*x^-1)=true then
		Add(stab,x*g^-1*x^-1);
	    fi;
            while not j=0 do
                w:=reps[j]*g*x^-1;
                
                if CheckMembership(w)=true then
		    Add(orb,reps[j]);
                    Add(elm,w);
                    reps_pos[Position(Reps,reps[j])]:=[Length(Orbs)+1,Length(orb)];
		    Remove(reps,j);
		fi;
		j:=j-1; 
            od;
        od;
	Add(Stabs,Group(stab));
	Add(Orbs,orb);
        Add(elms,elm);
	len:=Length(reps);
    od;

    return [Stabs,Orbs,reps_pos,elms];
end;
#############

####################################################

InducedMapFromComplexOfFiniteIndexSubgroup:=function(C,reps,CheckMembership)
local Dimension, Boundary, n, i, N, IsInH, EltsG, Gword2Hword, x, y, fullOrbits,
      StabsH, OrbitsH, g, Action, L, EltsH, NOrbitsH, PositionCoset, mapgensRec,
      G2HRec, G2H, Pair2Triple, Triple2Pair, reps_pos,elms,posH,posG, IDEL,
      Mult, MultRec, Stabilizer, BoundaryRec, BoundaryR, D, DD, mapgens, map;

IsInH:=CheckMembership;
EltsH:=[One(C!.group)];
EltsG:=C!.elts;
BoundaryR:=C!.boundary;
##### Length of the complex#####
N:=0;
while C!.dimension(N)>0 do
    N:=N+1;
od;
N:=N-1;
################################
posH:=function(g)
local p;
    p:=Position(EltsH,g);
    if p=fail then Add(EltsH,g); p:=Length(EltsH);fi;
    return p;
end;
################################
posG:=function(g)
local p;
    p:=Position(EltsG,g);
    if p=fail then Add(EltsG,g); p:=Length(EltsG);fi;
    return p;
end;
################################
StabsH:=[];
fullOrbits:=[];
OrbitsH:=[];
NOrbitsH:=[];
reps_pos:=[];
elms:=[];
for n in [0..N] do
    StabsH[n+1]:=[];
    fullOrbits[n+1]:=[];
    OrbitsH[n+1]:=[];
    NOrbitsH[n+1]:=[];
    reps_pos[n+1]:=[];
    elms[n+1]:=[];
    for i in [1..C!.dimension(n)] do
	L:=OrbitTesselation(C!.stabilizer(n,i),reps,IsInH);
	StabsH[n+1][i]:=L[1];
        for x in L[1] do 
            for y in Elements(x) do 
                posH(y);
            od;
        od;
        fullOrbits[n+1][i]:=L[2];
	OrbitsH[n+1][i]:=List(L[2],w->w[1]);
	NOrbitsH[n+1][i]:=Length(OrbitsH[n+1][i]);
        reps_pos[n+1][i]:=L[3];
        elms[n+1][i]:=L[4];
    od;
od;


################################
Dimension:=function(n)
if n>N then return 0;fi;
return Sum(NOrbitsH[n+1]);
end;
################################
PositionCoset:=function(g)
local i;
    for i in [1..Length(reps)] do
	if IsInH(EltsG[g]*reps[i]^-1) then
	    return i;
        fi;
    od;
    Print("Representative list is not complete! \n");
    Print(1/0);
    return fail;
end;
################################
G2HRec:=[];

G2H:=function(g)
local p,t, h;
if not IsBound(G2HRec[g]) then
    p:=PositionCoset(g);
    h:=EltsG[g]*reps[p]^-1; # Check this one with the group action
    t:=Position(EltsH,h);
    if t=fail then Add(EltsH,h); t:=Length(EltsH);fi;
    G2HRec[g]:=[t,p]; 
fi;
return G2HRec[g];
end;
################################
Triple2Pair:=function(n,k,i)
local m;
    m:= Sum(NOrbitsH[n+1]{[1..AbsInt(k)-1]}); 
    return [n,SignInt(k)*(m+i)];
end;
###############################
Pair2Triple:=function(n,k)
local i, d, kk;
    kk:=AbsInt(k);
    i:=1;
    d:=0;
    while true do
        if d+NOrbitsH[n+1][i]<kk then
            d:=d+NOrbitsH[n+1][i];            
        else break;
        fi;
        i:=i+1;
    od;
    return [n,SignInt(k)*i,k-d];
end;
################################
Gword2Hword:=function(n,w)
local x,y,v,pos, g, elm, p;
    v:=[];
    for x in w do
        y:=G2H(x[2]);
        pos:=reps_pos[n+1][AbsInt(x[1])][y[2]];
        elm:=elms[n+1][AbsInt(x[1])][pos[1]][pos[2]];
        g:=EltsH[y[1]]*elm;
        p:=Position(EltsH,g);
        if p=fail then Add(EltsH,g); p:=Length(EltsH);fi;
        y:=[Triple2Pair(n,x[1],pos[1])[2],p];
        Add(v,y);
    od;
    return v;
end;
################################
MultRec:=List([1..Length(reps)],i->[]);

Mult:=function(i,j)
local x,r;
    if not IsBound(MultRec[i][j]) then
        x:=reps[i]*EltsG[j];
        r:=Position(EltsG,x);

        if r=fail then Add(EltsG,x); r:=Length(EltsG); fi;
        MultRec[i][j]:= r;
    fi;

    return MultRec[i][j];
end;
################################
Stabilizer:=function(n,k)
local p, x, S, B;
    p:=Pair2Triple(n,k);
    S:=ConjugateGroup(C!.stabilizer(p[1],p[2]),OrbitsH[p[1]+1][p[2]][p[3]]^-1);
    B:=Filtered(Elements(S),w->IsInH(w));
    return Group(B);
end;

################################
BoundaryRec:=[];
for i in [1..N] do
BoundaryRec[i]:=[];
od;

Boundary:=function(n,k)
local kk,w,x,f;
    kk:=AbsInt(k);
    if n<=0 then return [];fi;
    
    if not IsBound(BoundaryRec[n][kk]) then 
        x:=Pair2Triple(n,kk);
        f:=Position(reps,OrbitsH[n+1][x[2]][x[3]]);
        w:=StructuralCopy(BoundaryR(n,x[2]));
        w:=List(w, y->[y[1],Mult(f,y[2])]);
        BoundaryRec[n][kk]:= Gword2Hword(n-1,w);
    fi;
    if k>0 then return BoundaryRec[n][kk]; 
    else return NegateWord(BoundaryRec[n][kk]); 
    fi;

end;
################################
Action:=function(k,i,g) # We assume that this G-complex is rigid, 
                        # means the stabilizer fixes the cell pointwise
return 1;
end;
################################
D:= Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            positionCoset:=PositionCoset,
            G2H:=G2H,
            Pair2Triple:= Pair2Triple,
            Triple2Pair:= Triple2Pair,
            nOrbits:=NOrbitsH,
            Gword2Hword:= Gword2Hword,
            mult:=Mult,
            fullOrbits:=fullOrbits,
            stabH:=StabsH,
            orbits:=OrbitsH,
            eltsG:=EltsG,
            elts:= EltsH,
            group:= Group(EltsH),
            stabilizer:= Stabilizer,
            action:=Action,
            properties:=
            [["length",1000],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

###################################
mapgensRec:=[];
for n in [0..N] do
    mapgensRec[n+1]:=[];
    for i in [1..Dimension(n)] do
        mapgensRec[n+1][i]:=[];
    od;
od;

IDEL:=posH(One(C!.group));

mapgens:=function(x,m)
local u,a,y;

if not IsBound(mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]) then
    if not x[2]=IDEL then
        y:=StructuralCopy(mapgens([x[1],IDEL],m));
        Apply(y,b->[b[1],posG(EltsH[x[2]]*EltsG[b[2]])]);
        if x[1]>0 then 
            mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=y;
	else
	    mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=NegateWord(y);
	fi;
    return y;
    fi;

    u:=Pair2Triple(m,AbsoluteValue(x[1]));
    y:=[[SignInt(x[1])*u[2],posG(OrbitsH[m+1][u[2]][u[3]])]];
    if x[1]>0 then
      mapgensRec[m+1][x[1]][x[2]]:=y;
    else
      mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=NegateWord(y);
    fi;
    return y;
   

else if x[1]>0 then return mapgensRec[m+1][AbsoluteValue(x[1])][x[2]];
else return NegateWord(mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]); fi;
fi;
end;
###################################
map:=function(w,m)
local  u,v,x,y,z;

    v:=Collected(w);
#Print("v=",v,"\n");
    Apply(v,x->MultiplyWord(x[2],  mapgens(x[1],m)));
#Print("v=",v,"\n");
    v:= Concatenation(v);
#Print("v=",v,"\n");
    return AlgebraicReduction(v);

end;
###################################
return Objectify(HapEquivariantChainMap,
	   rec(
	    source:=D,
	    target:=C,
	    mapping:=map,
            mapgens:=mapgens,
#            originalHom:=fail,
	    properties:=
	    [["type","equivariantChainMap"],
	     ["characteristic",0]  ]));
end;

###############################################################
GetTorsionSubcomplex_alt:=function(C,p)
local N, n, i, torsionCells, reIndex,centerSize,trivialTorsion,
      Dimension, Stabilizer, Action, Boundary;

##### Length of the complex#####
N:=0;
while C!.dimension(N)>0 do
    N:=N+1;
od;
N:=N-1;
################################
centerSize:=Size(Intersection(List([1..C!.dimension(N)],i->C!.stabilizer(N,i))));
trivialTorsion:=Lcm(centerSize,p);
## Identify positions of the torsion cells #####

torsionCells:=[];
reIndex:=[];
for n in [0..N] do
    torsionCells[n+1]:=[];
    reIndex[n+1]:=[];
    for i in [1..C!.dimension(n)] do
        if (Size(C!.stabilizer(n,i))/trivialTorsion) mod p = 0 then
            Add(torsionCells[n+1],i);
            reIndex[n+1][i]:=Length(torsionCells[n+1]);
        fi;
    od;
od;
#################################

Dimension:=function(n)
    if n>N then return 0;fi;
    return Length(torsionCells[n+1]);
end;
#################################


Stabilizer:=function(n,i)   
    return C!.stabilizer(n,torsionCells[n+1][i]);
end;
#################################

Action:=function(n,i,g)
    return 1;
end;
#################################

Boundary:=function(n,i)
local t, bdr;
    t:=torsionCells[n+1][i];
    bdr:=StructuralCopy(C!.boundary(n,t));
    Apply(bdr,w->[SignInt(w[1])*reIndex[n][AbsInt(w[1])],w[2]]);

    return bdr;
end;
#################################
return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:= C!.elts,
            group:= C!.group,
            stabilizer:= Stabilizer,
            action:=Action,
            torsion:=p,
            torsionCells:=torsionCells,
            reIndex:=reIndex,
            originalCellcomplex:=C,
            properties:=
            [["length",1000],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));
end;


##############################################################
ChainMapOfTorsionSubcomplexes:=function(f,p)
local D,C,map, R, S, mapping;

    D:=f!.source;
    C:=f!.target;
    map:=f!.mapping;
###################################
#### Get the torsion subcomplexes ##
    R:=GetTorsionSubcomplex_alt(D,p);
    S:=GetTorsionSubcomplex_alt(C,p);

###################################
mapping:=function(w,m)
local v, x, y, fv, fw;
    v:=[];
    for x in w do
        y:=[SignInt(x[1])*R!.torsionCells[m+1][AbsInt(x[1])],x[2]];
        Add(v,y);
    od;

    fv:=map(v,m);

    fw:=[];
    for x in fv do
        y:=[SignInt(x[1])*S!.reIndex[m+1][AbsInt(x[1])],x[2]];
        Add(fw,y);
    od;

    return fw;
end;
###################################
return Objectify(HapEquivariantChainMap,
	   rec(
	    source:=R,
	    target:=S,
	    mapping:=mapping,
	    properties:=
	    [["type","equivariantChainMap"],
	     ["characteristic",0]  ]));
end;
#################################################################
QuotientByTorsionSubcomplex_alt:=function(C,p)
local T, Dimension, Boundary, Stabilizer, n, reIndex, quotientCells,
      N, a, b, i, j, k, Action;

    T:=GetTorsionSubcomplex_alt(C,p);

##################################
Dimension:=function(n)
    if n=0 then return C!.dimension(0)-T!.dimension(0)+1;
    else return C!.dimension(n)-T!.dimension(n);
    fi;
end;
##################################

##### Length of the complex#####
N:=0;
while C!.dimension(N)>0 do
    N:=N+1;
od;
N:=N-1;
################################

## Identify positions of the quotient cells #####

quotientCells:=[];
reIndex:=[];

### For 0-cells ####
quotientCells[1]:=[];
reIndex[1]:=[];
for i in [1..C!.dimension(0)] do
    if IsBound(T!.reIndex[1][i]) then
        reIndex[1][i]:=Dimension(0);
        if not IsBound(quotientCells[1][Dimension(0)]) then
            quotientCells[1][Dimension(0)]:=i;
        fi;
    else
        a:=0;
            for j in [1..i] do
                if not IsBound(T!.reIndex[1][j]) then
                    a:=a+1;
                fi;
            od; 
            
            quotientCells[1][a]:=i;
            reIndex[1][i]:=a;
    fi;
od;

####################
for n in [1..N] do
    quotientCells[n+1]:=[];
    reIndex[n+1]:=[];
    for i in [1..C!.dimension(n)] do
        if not IsBound(T!.reIndex[n+1][i]) then
            a:=0;
            for j in [1..i] do
                if not IsBound(T!.reIndex[n+1][j]) then
                    a:=a+1;
                fi;
            od; 
            
            quotientCells[n+1][a]:=i;
            reIndex[n+1][i]:=a;            
        fi;
    od;
od;    
###############################
Stabilizer:=function(n,k)

if n=0 and k=Dimension(0) then
    return C!.group;
fi;

return C!.stabilizer(n,quotientCells[n+1][k]);

end;
##############################
Boundary:=function(n,k)
local i, w, bdr, b;
    i:=quotientCells[n+1][k];
    bdr:=C!.boundary(n,i);
    b:=[];
    for w in bdr do
        if IsBound(reIndex[n][AbsInt(w[1])]) then
            Add(b,[SignInt(w[1])*reIndex[n][AbsInt(w[1])],w[2]]);
        fi;
    od;
    return b;
end;
##############################
Action:=function(n,i,g)
    return 1;
end;
######################################
return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:= C!.elts,
            group:= C!.group,
            stabilizer:= Stabilizer,
            action:=Action,
            torsion:=p,
            quotientCells:=quotientCells,
            reIndex:=reIndex,
            originalCellcomplex:=C,
            properties:=
            [["length",1000],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

    
end;
#################################################################
ChainMapOfQuotientsByTorsionSubcomplex:=function(f,p)
local D,C,map, R, S, mapping;

    D:=f!.source;
    C:=f!.target;
    map:=f!.mapping;
###################################
#### Get the torsion subcomplexes ##
    R:=QuotientByTorsionSubcomplex_alt(D,p);
    S:=QuotientByTorsionSubcomplex_alt(C,p);

###################################
mapping:=function(w,m)
local v, x, y, fv, fw;
    v:=[];
    for x in w do
        y:=[SignInt(x[1])*R!.quotientCells[m+1][AbsInt(x[1])],x[2]];
        Add(v,y);
    od;

    fv:=map(v,m);

    fw:=[];
    for x in fv do
        if IsBound(S!.reIndex[m+1][AbsInt(x[1])]) then
            y:=[SignInt(x[1])*S!.reIndex[m+1][AbsInt(x[1])],x[2]];
            Add(fw,y);
        fi;
    od;

    return fw;
end;
###################################
return Objectify(HapEquivariantChainMap,
	   rec(
	    source:=R,
	    target:=S,
	    mapping:=mapping,
	    properties:=
	    [["type","equivariantChainMap"],
	     ["characteristic",0]  ]));



end;
###############################################################
