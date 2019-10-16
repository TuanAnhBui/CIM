ModPCellularCohomology:=function(C,p)
local lnth, Cells, nrCells,x,w,j,s,id, Mat,t,M,b,h,r,rank,
      Mult, CLeftCosetElt, pos, hom, Elts, MatRank, y, BdMat, LocalSparseMat,diag,
      KillGAction, Reduction, Boundary, myOne;



myOne:=One(GF(p));
##################################################################
LocalSparseMat:=function(M,nrows,ncols)

return Objectify(HapSparseMat,
               rec( rows:=nrows,
                    cols:=ncols,
                    characteristic:=0,
                    mat:=M));
end;

MatRank:=function(g)
  local B;
if not (IsBound(g[1]) and IsBound(g[1][1])) then return 0;
else
  return Length(BaseIntMat(g));
fi;
end;


    lnth:=0;
    while not C!.dimension(lnth)=0 do
        lnth:=lnth+1;
    od;
    lnth:=lnth-1;

###########################
Boundary:=function(N)
local k,row, Mt, i, j, x, sum, bdr;
    
    k:=N+1;
    if k=0 then return List([1..C!.dimension(1)],i->[]);;fi;
 
    bdr:=List([1..C!.dimension(k-1)],i->[]);
    if C!.dimension(k)>0 then
        for i in [1..C!.dimension(k)] do
            for j in [1..C!.dimension(k-1)] do
                sum:=0;
                for x in C!.boundary(k,i) do
                    if AbsoluteValue(x[1])=j then
                        sum := sum + SignInt(x[1]);
                    fi;
                od;
                if not sum mod p =0 then
                    Add(bdr[j],[i,sum*myOne]);
                fi;
            od;
        od;
    fi;
    return LocalSparseMat(bdr,C!.dimension(k-1),C!.dimension(k));
end;
###################################################################

rank:=[];
for j in [2..lnth+1] do

  rank[j]:=Rank(Boundary(j-2));


od;

rank[1]:=0;
rank[lnth+2]:=0;
hom:=[];
for j in [1..lnth+1] do

    r:=C!.dimension(j-1) - rank[j]-rank[j+1];
    Add(hom,r);
od;


 return [hom,Boundary];

end;


############################################################################
##############################################################
ModPCellularCohomologyWithSage:=function(eq,p,degree)
local C,D,matC,matD, myOne, Boundary,map,x,temp,j,f,n,hom,
      i,s,b, name,exportFile,bdr,t, inducedMap, run_file;

myOne:=One(GF(p));
hom:=[];
#####################################
Boundary:=function(R,N)
local k,row, Mt, i, j, x, sum, bdr;
    
    k:=N+1;
    if k=0 then return List([1..R!.dimension(1)],i->[]);;fi;
 
    bdr:=List([1..R!.dimension(k-1)],i->[]);
    if R!.dimension(k)>0 then
        for i in [1..R!.dimension(k)] do
            for j in [1..R!.dimension(k-1)] do
                sum:=0;
                for x in R!.boundary(k,i) do
                    if AbsoluteValue(x[1])=j then
                        sum := sum + SignInt(x[1]);
                    fi;
                od;
                if not sum mod p =0 then
                    Add(bdr[j],[i,sum mod p]);
                fi;
            od;
        od;
    fi;
    return bdr;
end;
#####################################
C:=Source(eq);
D:=Target(eq);

##### Boundary matrix for C ########
for i in [0..degree-1] do
    name:=Concatenation("C",String(i),"_C",String(i+1));
    exportFile:=Concatenation("./tmp/",name,".sage");
    PrintTo(exportFile,"\"","C_",i,"--> C_",i+1,": Boundary Matrix of size ",String(C!.dimension(i)),"x",String(C!.dimension(i+1)),"\"\n");
    AppendTo(exportFile,"M",String(i)," = Matrix(GF(",p,"),",C!.dimension(i),",",C!.dimension(i+1),",sparse = True) \n");
    bdr:=Boundary(C,i);
    for s in [1..C!.dimension(i)] do
        b:=bdr[s];
        for t in b do
            AppendTo(exportFile,"M",String(i),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
        od;
        AppendTo(exportFile,"\n");
    od;
od;

for i in [0..degree-1] do
    name:=Concatenation("D",String(i),"_D",String(i+1));
    exportFile:=Concatenation("./tmp/",name,".sage");
    PrintTo(exportFile,"\"","D_",i,"--> D_",i+1,": Boundary Matrix of size ",String(D!.dimension(i)),"x",String(C!.dimension(i+1)),"\"\n");
    AppendTo(exportFile,"N",String(i)," = Matrix(GF(",p,"),",D!.dimension(i),",",D!.dimension(i+1),",sparse = True) \n");
    bdr:=Boundary(D,i);
    for s in [1..D!.dimension(i)] do
        b:=bdr[s];
        for t in b do
            AppendTo(exportFile,"N",String(i),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
        od;
        AppendTo(exportFile,"\n");
    od;
od;
####### Matrix of chain map #####################
#** for degree 0 **
for n in [0..0] do
f:=List([1..D!.dimension(n)],i->[]);

for i in [1..D!.dimension(n)] do
    for j in [1..C!.dimension(n)] do
        temp:=0;
        for x in List(eq!.mapping([[j,1]],n),y->y[1]) do
            if x=i then temp:=temp+1;fi;
            if x=-i then temp:=temp-1; fi;
        od;
        if not temp mod p = 0 then
            Add(f[i],[j,temp mod p]);
        fi;
    od;
od;

name:=Concatenation("f",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\"","D_",n,"--> C_",n,": Chain map of size ",D!.dimension(n),"x",C!.dimension(n),"\"\n");
AppendTo(exportFile,"f",String(n)," = Matrix(GF(",p,"),",D!.dimension(n),",",C!.dimension(n),",sparse = True) \n");

for s in [1..D!.dimension(n)] do
    b:=f[s];
    for t in b do
        AppendTo(exportFile,"f",String(n),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
    od;
    AppendTo(exportFile,"\n");
od;

name:=Concatenation("cmd",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\" Commands for computing cohomology of degree ",n,"\"\n");
AppendTo(exportFile,"from time import time, gmtime, strftime \n"); 
AppendTo(exportFile,"logfile = strftime(\"%d%b%Y%H%M\", gmtime()) \n");

AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('LOGFILE FOR THE COMPUTATION OF D'+str(n)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"total_run_time = time() \n");
AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"n = ",n,"\n");
AppendTo(exportFile,"load(\"./tmp/C",n,"_C",n+1,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the boundary map C'+str(n)+'_C'+str(n+1)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/D",n,"_D",n+1,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the boundary map D'+str(n)+'_D'+str(n+1)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/f",n,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the chain map f'+str(n)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"V2 = N",n,".left_kernel()\n");
AppendTo(exportFile,"B = V2.basis()\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the kernel of N has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"L = [x*f",n," for x in B]\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Mapping the basis of the kernel of N has been done!')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"V1 = M",n,".left_kernel()\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the kernel of M has been done!')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"LL = [list(V1.coordinate_vector(x)) for x in L]\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computing the coordinates of the images of f has been done!')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"with open(os.path.expanduser(\"./results/d",String(n),".g\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('return '+str(LL)+';')\n");

AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation the matrix d'+str(n)+' has been done!')\n");
AppendTo(exportFile,"    f.write('Totally, it took '+str(time()-total_run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

od;



#** for degree > 0 **
for n in [1..degree-1] do
f:=List([1..D!.dimension(n)],i->[]);

for i in [1..D!.dimension(n)] do
    for j in [1..C!.dimension(n)] do
        temp:=0;
        for x in List(eq!.mapping([[j,1]],n),y->y[1]) do
            if x=i then temp:=temp+1;fi;
            if x=-i then temp:=temp-1; fi;
        od;
        if not temp mod p = 0 then
            Add(f[i],[j,temp mod p]);
        fi;
    od;
od;

name:=Concatenation("f",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\"","D_",n,"--> C_",n,": Chain map of size ",D!.dimension(n),"x",C!.dimension(n),"\"\n");
AppendTo(exportFile,"f",String(n)," = Matrix(GF(",p,"),",D!.dimension(n),",",C!.dimension(n),",sparse = True) \n");

for s in [1..D!.dimension(n)] do
    b:=f[s];
    for t in b do
        AppendTo(exportFile,"f",String(n),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
    od;
    AppendTo(exportFile,"\n");
od;

name:=Concatenation("cmd",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");

PrintTo(exportFile,"\" Commands for computing cohomology of degree ",n,"\"\n");
AppendTo(exportFile,"from time import time, gmtime, strftime \n"); 
AppendTo(exportFile,"logfile = strftime(\"%d%b%Y%H%M\", gmtime()) \n");

AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('LOGFILE FOR THE COMPUTATION OF D'+str(n)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"n = ",n,"\n");
AppendTo(exportFile,"total_run_time = time() \n");
AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/C",n-1,"_C",n,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the boundary map C'+str(n-1)+'_C'+str(n)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/C",n,"_C",n+1,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the boundary map C'+str(n)+'_C'+str(n+1)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/D",n-1,"_D",n,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the boundary map D'+str(n-1)+'_D'+str(n)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/D",n,"_D",n+1,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the boundary map D'+str(n)+'_D'+str(n+1)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"load(\"./tmp/f",n,".sage\")\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Finished loading the chain map f'+str(n)+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"V2 = N",n,".left_kernel()\n");
AppendTo(exportFile,"B = V2.basis()\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the kernel of N'+str(n)+' has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"W2 = N",n-1,".image()\n");

AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the image of N'+str(n-1)+' has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"Q2 = V2/W2\n");
AppendTo(exportFile,"B = Q2.basis()\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the quotient has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"del V2\n");
AppendTo(exportFile,"del W2\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"L = [Q2.lift(x)*f",n," for x in B]\n");
AppendTo(exportFile,"del Q2\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Mapping the basis of the kernel of N has been done!')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"V1 = M",n,".left_kernel()\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the kernel of M'+str(n)+' has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"W1 = M",n-1,".image()\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the image of N'+str(n-1)+' has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");

AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"Q1 = V1/W1\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation of the quotient has been done! ')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"del V1\n");
AppendTo(exportFile,"del W1\n");
AppendTo(exportFile,"map = Q1.quotient_map()\n");


AppendTo(exportFile,"run_time = time() \n");
AppendTo(exportFile,"LL = [list(Q1.coordinate_vector(map(x))) for x in L]\n");
AppendTo(exportFile,"del Q1\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computing the coordinates of the images of f has been done!')\n");
AppendTo(exportFile,"    f.write('It took '+str(time()-run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


AppendTo(exportFile,"with open(os.path.expanduser(\"./results/d",String(n),".g\"),", "\"a+\") as f:");
AppendTo(exportFile,"    f.write('return '+str(LL)+';') \n");

AppendTo(exportFile,"with open(os.path.expanduser(\"./logs/log",String(n),".txt\"),", "\"a+\") as f:\n");
AppendTo(exportFile,"    f.write('Computation the matrix d'+str(n)+' has been done!')\n");
AppendTo(exportFile,"    f.write('Totally, it took '+str(time()-total_run_time)+' seconds!'+'");
AppendTo(exportFile,"""\n""");
AppendTo(exportFile,"')\n");


od;

for n in [0..degree-1] do
name:=Concatenation("cmd",String(n));
run_file:=Concatenation("sage ","./tmp/",name,".sage");
Exec(run_file);

inducedMap:=ReadAsFunction(Concatenation("./results/d",String(n),".g"));
Add(hom,TransposedMat(inducedMap()));
od;
return hom;
end;
######################################################################
##############################################################
IntegralCellularCohomologyWithSage:=function(eq,degree)
local C,D,matC,matD, myOne, Boundary,map,x,temp,j,f,n,hom,
      i,s,b, name,exportFile,bdr,t, inducedMap, run_file;

hom:=[];
#####################################
Boundary:=function(R,N)
local k,row, Mt, i, j, x, sum, bdr;
    
    k:=N+1;
    if k=0 then return List([1..R!.dimension(1)],i->[]);;fi;
 
    bdr:=List([1..R!.dimension(k-1)],i->[]);
    if R!.dimension(k)>0 then
        for i in [1..R!.dimension(k)] do
            for j in [1..R!.dimension(k-1)] do
                sum:=0;
                for x in R!.boundary(k,i) do
                    if AbsoluteValue(x[1])=j then
                        sum := sum + SignInt(x[1]);
                    fi;
                od;
                if not sum = 0 then
                    Add(bdr[j],[i,sum]);
                fi;
            od;
        od;
    fi;
    return bdr;
end;
#####################################
C:=Source(eq);
D:=Target(eq);

##### Boundary matrix for C ########
for i in [0..degree-1] do
    name:=Concatenation("C",String(i),"_C",String(i+1));
    exportFile:=Concatenation("./tmp/",name,".sage");
    PrintTo(exportFile,"\"","C_",i,"--> C_",i+1,": Boundary Matrix of size ",String(C!.dimension(i)),"x",String(C!.dimension(i+1)),"\"\n");
    AppendTo(exportFile,"M",String(i)," = Matrix(ZZ,",C!.dimension(i),",",C!.dimension(i+1),",sparse = True) \n");
    bdr:=Boundary(C,i);
    for s in [1..C!.dimension(i)] do
        b:=bdr[s];
        for t in b do
            AppendTo(exportFile,"M",String(i),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
        od;
        AppendTo(exportFile,"\n");
    od;
od;

for i in [0..degree-1] do
    name:=Concatenation("D",String(i),"_D",String(i+1));
    exportFile:=Concatenation("./tmp/",name,".sage");
    PrintTo(exportFile,"\"","D_",i,"--> D_",i+1,": Boundary Matrix of size ",String(D!.dimension(i)),"x",String(C!.dimension(i+1)),"\"\n");
    AppendTo(exportFile,"N",String(i)," = Matrix(ZZ,",D!.dimension(i),",",D!.dimension(i+1),",sparse = True) \n");
    bdr:=Boundary(D,i);
    for s in [1..D!.dimension(i)] do
        b:=bdr[s];
        for t in b do
            AppendTo(exportFile,"N",String(i),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
        od;
        AppendTo(exportFile,"\n");
    od;
od;
####### Matrix of chain map #####################
#** for degree 0 **
for n in [0..0] do
f:=List([1..D!.dimension(n)],i->[]);

for i in [1..D!.dimension(n)] do
    for j in [1..C!.dimension(n)] do
        temp:=0;
        for x in List(eq!.mapping([[j,1]],n),y->y[1]) do
            if x=i then temp:=temp+1;fi;
            if x=-i then temp:=temp-1; fi;
        od;
        if not temp = 0 then
            Add(f[i],[j,temp]);
        fi;
    od;
od;

name:=Concatenation("f",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\"","D_",n,"--> C_",n,": Chain map of size ",D!.dimension(n),"x",C!.dimension(n),"\"\n");
AppendTo(exportFile,"f",String(n)," = Matrix(ZZ,",D!.dimension(n),",",C!.dimension(n),",sparse = True) \n");

for s in [1..D!.dimension(n)] do
    b:=f[s];
    for t in b do
        AppendTo(exportFile,"f",String(n),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
    od;
    AppendTo(exportFile,"\n");
od;

name:=Concatenation("cmd",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\" Commands for computing cohomology of degree ",n,"\"\n");
AppendTo(exportFile,"load(\"./tmp/C",n,"_C",n+1,".sage\")\n");
AppendTo(exportFile,"load(\"./tmp/D",n,"_D",n+1,".sage\")\n");
AppendTo(exportFile,"load(\"./tmp/f",n,".sage\")\n");
AppendTo(exportFile,"V2 = N",n,".left_kernel()\n");
AppendTo(exportFile,"B = V2.basis()\n");
AppendTo(exportFile,"L = [x*f",n," for x in B]\n");
AppendTo(exportFile,"V1 = M",n,".left_kernel()\n");
AppendTo(exportFile,"LL = [list(V1.coordinate_vector(x)) for x in L]\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./results/d",String(n),".g\"),", "\"a+\") as f:");
AppendTo(exportFile,"    f.write('return '+str(LL)+';')");


od;



#** for degree > 0 **
for n in [1..degree-1] do
f:=List([1..D!.dimension(n)],i->[]);

for i in [1..D!.dimension(n)] do
    for j in [1..C!.dimension(n)] do
        temp:=0;
        for x in List(eq!.mapping([[j,1]],n),y->y[1]) do
            if x=i then temp:=temp+1;fi;
            if x=-i then temp:=temp-1; fi;
        od;
        if not temp = 0 then
            Add(f[i],[j,temp]);
        fi;
    od;
od;

name:=Concatenation("f",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\"","D_",n,"--> C_",n,": Chain map of size ",D!.dimension(n),"x",C!.dimension(n),"\"\n");
AppendTo(exportFile,"f",String(n)," = Matrix(ZZ,",D!.dimension(n),",",C!.dimension(n),",sparse = True) \n");

for s in [1..D!.dimension(n)] do
    b:=f[s];
    for t in b do
        AppendTo(exportFile,"f",String(n),"[",s-1,",",AbsInt(t[1])-1,"] = ",t[2],";");
    od;
    AppendTo(exportFile,"\n");
od;

name:=Concatenation("cmd",String(n));
exportFile:=Concatenation("./tmp/",name,".sage");
PrintTo(exportFile,"\" Commands for computing cohomology of degree ",n,"\"\n");
AppendTo(exportFile,"load(\"./tmp/C",n-1,"_C",n,".sage\")\n");
AppendTo(exportFile,"load(\"./tmp/C",n,"_C",n+1,".sage\")\n");
AppendTo(exportFile,"load(\"./tmp/D",n-1,"_D",n,".sage\")\n");
AppendTo(exportFile,"load(\"./tmp/D",n,"_D",n+1,".sage\")\n");
AppendTo(exportFile,"load(\"./tmp/f",n,".sage\")\n");
AppendTo(exportFile,"V2 = N",n,".left_kernel()\n");
AppendTo(exportFile,"W2 = N",n-1,".image()\n");
AppendTo(exportFile,"Q2 = V2/W2\n");
AppendTo(exportFile,"del V2\n");
AppendTo(exportFile,"del W2\n");
AppendTo(exportFile,"B = Q2.basis()\n");
AppendTo(exportFile,"L = [Q2.lift(x)*f",n," for x in B]\n");
AppendTo(exportFile,"del Q2\n");
AppendTo(exportFile,"V1 = M",n,".left_kernel()\n");
AppendTo(exportFile,"W1 = M",n-1,".image()\n");
AppendTo(exportFile,"Q1 = V1/W1\n");
AppendTo(exportFile,"del V1\n");
AppendTo(exportFile,"del W1\n");
AppendTo(exportFile,"map = Q1.quotient_map()\n");

AppendTo(exportFile,"LL = [list(Q1.coordinate_vector(map(x))) for x in L]\n");
AppendTo(exportFile,"del Q1\n");
AppendTo(exportFile,"with open(os.path.expanduser(\"./results/d",String(n),".g\"),", "\"a+\") as f:");
AppendTo(exportFile,"    f.write('return '+str(LL)+';')");

od;

for n in [0..degree-1] do
#name:=Concatenation("cmd",String(n));
#run_file:=Concatenation("sage ","~/tmp/",name,".sage");
#Exec(run_file);

#inducedMap:=ReadAsFunction(Concatenation("./results/d",String(n),".g"));
#Add(hom,TransposedMat(inducedMap()));
od;
return hom;
end;
######################################################################


