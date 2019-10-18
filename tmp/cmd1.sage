" Commands for computing cohomology of degree 1"
from time import time, gmtime, strftime 
logfile = strftime("%d%b%Y%H%M", gmtime()) 
n = 1
with open(os.path.expanduser("./logs/log1.txt"),"w+") as f:
    f.write('LOGFILE FOR THE COMPUTATION OF D'+str(n)+'\n')
total_run_time = time() 
run_time = time() 
load("./tmp/C0_C1.sage")
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Finished loading the boundary map C'+str(n-1)+'_C'+str(n)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/C1_C2.sage")
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Finished loading the boundary map C'+str(n)+'_C'+str(n+1)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/D0_D1.sage")
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Finished loading the boundary map D'+str(n-1)+'_D'+str(n)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/D1_D2.sage")
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Finished loading the boundary map D'+str(n)+'_D'+str(n+1)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/f1.sage")
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Finished loading the chain map f'+str(n)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
V2 = N1.left_kernel()
B = V2.basis()
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation of the kernel of N'+str(n)+' has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
W2 = N0.image()
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation of the image of N'+str(n-1)+' has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
Q2 = V2/W2
B = Q2.basis()
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation of the quotient has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
del V2
del W2
run_time = time() 
L = [Q2.lift(x)*f1 for x in B]
del Q2
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Mapping the basis of the kernel of N has been done!')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
V1 = M1.left_kernel()
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation of the kernel of M'+str(n)+' has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
W1 = M0.image()
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation of the image of N'+str(n-1)+' has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
Q1 = V1/W1
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation of the quotient has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
del V1
del W1
map = Q1.quotient_map()
run_time = time() 
LL = [list(Q1.coordinate_vector(map(x))) for x in L]
del Q1
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computing the coordinates of the images of f has been done!')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
with open(os.path.expanduser("./results/d1.g"),"w+") as f:    f.write('return '+str(LL)+';') 
with open(os.path.expanduser("./logs/log1.txt"),"a+") as f:
    f.write('Computation the matrix d'+str(n)+' has been done!')
    f.write('Totally, it took '+str(time()-total_run_time)+' seconds!'+'\n')
