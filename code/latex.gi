################################################
Matrix2Latex:=function(M)
local n,m, i, j, s;
    n:=Length(M);
    m:=Length(M[1]);
    s:="";
    s:=Concatenation(s,"$\\begin{pmatrix}");
    for i in [1..n] do
        for j in [1..m] do
            s:=Concatenation(s,String(M[i][j]));
            if not j=m then
                s:=Concatenation(s,"&");
            fi;
        od;
        if not i=n then 
            s:=Concatenation(s,"\\\\");
        fi;
    od;
    s:=Concatenation(s,"\\end{pmatrix}$");
    return s;
end;
################################################
MatCol:=function(M,i) return List(M,w->w[i]);end;
  
OrbitStabilizer2Latex:=function(L)
local N,i,j, w, Stab, Orbits, s, a;
    Stab:=L[1];
    Orbits:=L[2];
    N:=Length(Stab);
    s:="";
    s:=Concatenation(s,"\\begin{itemize}\n");
    for i in [1..N] do
        s:=Concatenation(s,"\\item Isotropy ",StructureDescription(Stab[i]),",   Orbit ");
        a:="{";
        for j in [1..Length(Orbits[i])] do
            w:=Orbits[i][j];
            a:=Concatenation(a,Concatenation(List(TransposedMat(invertMatrixEmbedding(w^-1,0))[1],String)));
            if not j=Length(Orbits[i]) then
                a:=Concatenation(a,",");
            fi;
        od;
        a:=Concatenation(a,"}");
        s:=Concatenation(s,a);
        s:=Concatenation(s,"\n");
    od;
    s:=Concatenation(s,"\\end{itemize} \n");
    return s;
end;
################################################
PrintToLatexFile:=function(eq,filename)
local N, n, i, C, D, L;

C:=eq!.target;
D:=eq!.source;
##### Length of the complex#####
N:=0;
while C!.dimension(N)>0 do
    N:=N+1;
od;
N:=N-1;
################################
PrintTo(filename,"\n");

AppendTo(filename,"\\documentclass[a4paper]{extarticle} \n");
AppendTo(filename,"\\usepackage{amssymb,amsmath,amsthm,graphicx,color} \n");
AppendTo(filename,"\\begin{document} \n");
AppendTo(filename,"\\title{Induce cell complex for finite index subgroup} \n");
AppendTo(filename,"\\maketitle \n");

for n in [0..N] do
    AppendTo(filename,Concatenation("\\subsection{For ",String(n),"-cells}\n"));
    for i in [1..C!.dimension(n)] do
        AppendTo(filename,Concatenation("\\subsubsection{The ",String(i),"-th of ",String(n),"-cells}\n"));
        L:=[D!.stabH[n+1][i],D!.fullOrbits[n+1][i]];
        AppendTo(filename,OrbitStabilizer2Latex(L));
        AppendTo(filename,"\n");
    od;
od;
AppendTo(filename,"\\end{document} \n");
end;

##########################################################################################