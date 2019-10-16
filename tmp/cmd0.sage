" Commands for computing cohomology of degree 0"
from time import time, gmtime, strftime 
logfile = strftime("%d%b%Y%H%M", gmtime()) 
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('LOGFILE FOR THE COMPUTATION OF D'+str(n)+'\n')
total_run_time = time() 
run_time = time() 
n = 0
load("./tmp/C0_C1.sage")
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Finished loading the boundary map C'+str(n)+'_C'+str(n+1)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/D0_D1.sage")
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Finished loading the boundary map D'+str(n)+'_D'+str(n+1)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
load("./tmp/f0.sage")
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Finished loading the chain map f'+str(n)+'\n')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
V2 = N0.left_kernel()
B = V2.basis()
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Computation of the kernel of N has been done! ')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
L = [x*f0 for x in B]
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Mapping the basis of the kernel of N has been done!')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
V1 = M0.left_kernel()
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Computation of the kernel of M has been done!')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
run_time = time() 
LL = [list(V1.coordinate_vector(x)) for x in L]
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Computing the coordinates of the images of f has been done!')
    f.write('It took '+str(time()-run_time)+' seconds!'+'\n')
with open(os.path.expanduser("./results/d0.g"),"w") as f:
    f.write('return '+str(LL)+';')
with open(os.path.expanduser("./logs/log0.txt"),"w") as f:
    f.write('Computation the matrix d'+str(n)+' has been done!')
    f.write('Totally, it took '+str(time()-total_run_time)+' seconds!'+'\n')
