" Commands for computing cohomology of degree 1"
from time import time, gmtime, strftime 
logfile = strftime("%d%b%Y%H%M", gmtime()) 
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('LOGFILE FOR THE COMPUTATION OF D'+String(n)\n')
total_run_time = time() 
run_time = time() 
load("./tmp/C0_C1.sage")
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Finished loading the boundary map C'+String(n-1)+'_C'+String(n)+'\n')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/C1_C2.sage")
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Finished loading the boundary map C'+String(n)+'_C'+String(n+1)+'\n')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/D0_D1.sage")
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Finished loading the boundary map D'+String(n-1)+'_D'+String(n)+'\n')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/D1_D2.sage")
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Finished loading the boundary map D'+String(n)+'_D'+String(n+1)+'\n')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/f1.sage")
    f.write('Finished loading the chain map f'+String(n)+'\n')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
V2 = N1.left_kernel()
B = V2.basis()
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation of the kernel of N'+String(n)+' has been done! ')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
W2 = N0.image()
B = Q2.basis()
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation of the image of N'+String(n-1)+' has been done! ')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
Q2 = V2/W2
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation of the quotient has been done! ')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
del V2
del W2
run_time = time() 
L = [Q2.lift(x)*f1 for x in B]
del Q2
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Mapping the basis of the kernel of N has been done!')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
V1 = M1.left_kernel()
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation of the kernel of M'+String(n)+' has been done! ')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
W1 = M0.image()
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation of the image of N'+String(n-1)+' has been done! ')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
run_time = time() 
Q1 = V1/W1
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation of the quotient has been done! ')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
del V1
del W1
map = Q1.quotient_map()
run_time = time() 
LL = [list(Q1.coordinate_vector(map(x))) for x in L]
del Q1
with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computing the coordinates of the images of f has been done!')
    f.write('It took '+String(time()-run_time)+' seconds!'+'\n')
with open(os.path.expanduser("./results/d1.g"),"w") as f:    f.write('return '+str(LL)+';')with open(os.path.expanduser("./logs/logfile.txt"),"w") as f:
    f.write('Computation the matrix d'+String(n)+' has been done!')
    f.write('Totally, it took '+String(time()-total_run_time)+' seconds!'+'\n')
