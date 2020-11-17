


import numpy as np



def readSCATTER(G): ## Read scattering data of discrete neutron energy groups = G
    
    file = open("input"+str(G)+"g.txt")
    numpy_array = np.loadtxt(file, delimiter = ",", skiprows = 19+G, max_rows = G)
    return numpy_array

def readOTHER(G): ## Read "other" data of discrete neutron energy groups = G
    
    file = open("input"+str(G)+"g.txt")
    numpy_array = np.loadtxt(file, delimiter = ",", skiprows = 11, max_rows = G)
    return numpy_array

def C(G): ## Create Sout,A,Sin,F matrices from nuclear data
    
    SCATTER = readSCATTER(G) ## Makes diagonal zero, meaning scatter from h to h cross section set to zero (for math purposes, not reality).
    for row in range(0,G):
        for col in range(0,G):
            if row == col:
                SCATTER[row][col] = 0
            else:
                1==1
    OTHER = readOTHER(G)
    
    diagonal = np.empty([G]) ## Sout block
    for col in range(0,G):
        diagonal[col] = np.sum(SCATTER[:,col]) 
    Sout = np.diag(diagonal)
    
    diagonal = OTHER[:,0] ## A block
    A = np.diag(diagonal)
    
    Sin = SCATTER ## Sin block
    for row in range(0,G): ## Makes diagonal zero, meaning scatter from h to h cross section set to zero (for math purposes, not reality).
        for col in range(0,G):
            if row == col:
                Sin[row][col] = 0
            else:
                1==1
    
    F = np.empty([G,G]) ## F block
    for row in range(0,G):
        F[row] = OTHER[:,2][row] * OTHER[:,1]
    
    storage = [Sout,A,Sin,F]
        
    # print(Sout)
    # print('\n')
    # print(A)
    # print('\n')
    # print(Sin)
    # print('\n')
    # print(F)
    # print('\n')
    return storage
    
def M(G): ## Mobility matrix
    Sout,A,Sin = C(G)[0],C(G)[1],C(G)[2]
    temp = np.add(A,Sout)
    M = np.subtract(temp,Sin)
    return M

def F(G): ## Fission matrix
    return C(G)[3]

def P(G): ## Problem matrix
    a = np.linalg.inv(M(G))
    b = F(G)
    return np.matmul(a,np.transpose(b))

def eigenvalues(G):
    values,vectors = np.linalg.eig(P(G))
    return values

def eigenvectors(G):
    values,vectors = np.linalg.eig(P(G))
    return vectors

def powerit(G,maxit):
    P_ = P(G)
    k = 1 ## initialize
    temp = np.ones(G)
    flux = np.transpose(temp) ## initialize
    for i in range(1,maxit):
        num = np.matmul(P_,flux)
        temp2 = np.matmul(P_,flux)
        denom = np.linalg.norm(temp2,2)
        flux = np.divide(num,denom)
        
        temp3 = np.matmul(P_,flux)
        temp4 = np.transpose(temp3)
        num = np.matmul(temp4,flux)
        temp5 = np.transpose(flux)
        denom = np.matmul(temp5,flux)
        k = np.divide(num,denom)
    storage = [k,flux]
    return storage

M2 = M(2)
F2 = F(2)
P2 = P(2)
a = (powerit(2,3)[0])
b = (powerit(2,3)[1])
c = (eigenvalues(2))
d = (eigenvectors(2))

M8 = M(8)
F8 = F(8)
P8 = P(8)
e = (powerit(8,3)[0])
f = (powerit(8,3)[1])
g = (eigenvalues(8))
h = (eigenvectors(8))

np.set_printoptions(linewidth = 150)

with open('input2g.txt', 'r') as file:
    text2g = file.read()

with open('output2g.txt', 'w') as file:
    file.write(str(text2g))
    file.write('\n' + '\n' + '\n' + '\n' + 'The following data corresponds to neutron balance output with 2 discrete neutron energy groups' + '\n' + '\n')
    file.write(str(M2) + '\n' + 'The above data depicts the Migration matrix' + '\n' + '\n')
    file.write(str(F2) + '\n' + 'The above data depicts the Fission matrix' + '\n' + '\n')
    file.write(str(P2) + '\n' + 'The above data depicts the Problem matrix' + '\n' + '\n')
    file.write(str(a) + '\n' + 'The above value is the dominant eigenvalue found by the Power Iteration method using three iterations' + '\n' + '\n')
    file.write(str(b) + '\n' + 'The above values make up the dominant eigenvector found by the Power Iteration method using three iterations' + '\n' + '\n')
    file.write(str(c) + '\n' + 'The above values make up the eigenvalues found by numpy.linalg.eig()' + '\n' + '\n')
    file.write(str(d) + '\n' + 'The above matrix has columns that make up the eigenvectors found by numpy.linalg.eig()' + '\n' + '\n')

with open('input8g.txt', 'r') as file:
    text8g = file.read()

with open('output8g.txt', 'w') as file:
    file.write(str(text8g))
    file.write('\n' + '\n' + '\n' + '\n' + 'The following data corresponds to neutron balance output with 8 discrete neutron energy groups' + '\n' + '\n')
    file.write(str(M8) + '\n' + 'The above data depicts the Migration matrix' + '\n' + '\n')
    file.write(str(F8) + '\n' + 'The above data depicts the Fission matrix' + '\n' + '\n')
    file.write(str(P8) + '\n' + 'The above data depicts the Problem matrix' + '\n' + '\n')
    file.write(str(e) + '\n' + 'The above value is the dominant eigenvalue found by the Power Iteration method using three iterations' + '\n' + '\n')
    file.write(str(f) + '\n' + 'The above values make up the dominant eigenvector found by the Power Iteration method using three iterations' + '\n' + '\n')
    file.write(str(g) + '\n' + 'The above values make up the eigenvalues found by numpy.linalg.eig()' + '\n' + '\n')
    file.write(str(h) + '\n' + 'The above matrix has columns that make up the eigenvectors found by numpy.linalg.eig()' + '\n' + '\n')

    
















