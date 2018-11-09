#!/usr/bin/env python
"""
Monte Carlo simulation of the 2D Ising model
"""

from scipy import *
from scipy import weave
from pylab import *

#These should go into __name part
#THEN YOU HAVE TO MODIFY THE INPUT OF SAMPLEPYTHON GIVING THESE TO IT.
Nitt = 1000000  # total number of Monte Carlo steps
N = 10          # linear dimension of the lattice, lattice-size= N x N
warm = 1000     # Number of warmup steps
measure=100     # How often to take a measurement


#DEFINE A FUNCTION TO COMPUTE THE ENERGY 
def CEnergy(latt):
    "Energy of a 2D Ising lattice at particular configuration"
    Ene = 0
    #GO THRU ALL THE SITES OF THE LATTICE:
    for i in range(len(latt)):
        for j in range(len(latt)):
            S = latt[i,j]
            #NEAREST NEIGHBOR %N=modulo N
            #(i-1)%N is ok in python but not in C++ (it gives back -1 and 9 as it should)
            WF = latt[(i+1)%N, j] + latt[i,(j+1)%N] + latt[(i-1)%N,j] + latt[i,(j-1)%N]
            #SUM UP TO ENE
            Ene += -0.5*WF*S # Each neighbor gives energy 1.0
    return Ene # Each par counted twice

def RandomL(N):
# you can replace by single line:
    #return sign(2*rand(N,N)-1).astype(int) #By chuck

    "Radom lattice, corresponding to infinite temerature"
    latt = zeros((N,N), dtype=int)
    for i in range(N):
        for j in range(N):
            latt[i,j] = sign(2*rand()-1)
    return latt

def SamplePython(Nitt, latt, PW,T):
    ###YOU CAN MODIFY THIS CODE TO WORK LIKE THE C++ CODE BELOW, USING ARRAY[5]!!
    ### EXERCISE:
    "Monte Carlo sampling for the Ising model in Pythons"
    ##COMPUTING THE ENERGY AND MAGNETIZATION USING FUNCTIONS:
    Ene = CEnergy(latt)  # Starting energy
    Mn=sum(latt)         # Starting magnetization
    
    #AVERAGED QUANTITY:
    Naver=0       # Measurements
    Eaver=0.0
    Maver=0.0
    E2aver=0.0
    M2aver=0.0
    
    N2 = N*N
    for itt in range(Nitt):
        #SELECT ONE OF THE RANDOM SITE IN THE LATTICE (N2 sites)
        t = int(rand()*N2)
        #TRIAL SITE
        (i,j) = (t % N, t/N)
        #COMPUTE ENERGY AND MAGNETIZATION
        S = latt[i,j]
        #WEISS FIELD:
        WF = latt[(i+1)%N, j] + latt[i,(j+1)%N] + latt[(i-1)%N,j] + latt[i,(j-1)%N]
        #PROBABILITY!
        P = PW[4+S*WF] #S*WF = DeltaE/2
        #METROPOLIS ALGORITHM!!
        ## if P>=1 or P>rand()
        if P>rand(): # flip the spin ## ACCEPT
            #FLIP THE SPIN:
            latt[i,j] = -S
            #UPDATE THE ENERGY (see notes why 2*)
            Ene += 2*S*WF
            #MAGNETIZATINO newM - oldM = S- -S=2S
            Mn -= 2*S
            
        #MAKE THE MEASUREMENTS:
        #THE MEASUREMENTS HAVE TO BE DONE INSIDE THE LOOP. EACH 
        #FOR EVERY LOOP > WARMUP AND EVERY "measure"
        if itt>warm and itt%measure==0:
            Naver += 1
            Eaver += Ene
            Maver += Mn
            E2aver += Ene*Ene;
            M2aver += Mn*Mn
    aE=Eaver/Naver
    aM=Maver/Naver
    cv = (E2aver/Naver -(Eaver/Naver)**2)/T**2
    chi = (M2aver/Naver-(Maver/Naver)**2)/T

    return (aM, aE, cv, chi)




def SampleCPP(Nitt, latt, PW, T):
    "The same Monte Carlo sampling in C++"
    Ene = float(CEnergy(latt))  # Starting energy
    Mn = float(sum(latt))       # Starting magnetization

    #DEFINE AN ARRAY OF 5 STORING ALL THE QTIES, THE ARRAY IS DFLOAT
    # Measurements
    aver = zeros(5,dtype=float) # contains: [Naver, Eaver, Maver]
    
    code="""
    using namespace std;
    int N2 = N*N;
    for (int itt=0; itt<Nitt; itt++){
        int t = static_cast<int>(drand48()*N2);
        int i = t % N;
        int j = t / N;
        int S = latt(i,j);
        int WF = latt((i+1)%N, j) + latt(i,(j+1)%N) + latt((i-1+N)%N,j) + latt(i,(j-1+N)%N);
        double P = PW(4+S*WF);
        if (P > drand48()){ // flip the spin
            latt(i,j) = -S;
            Ene += 2*S*WF;
            Mn -= 2*S;
        }
        if (itt>warm && itt%measure==0){
            aver(0) += 1;
            aver(1) += Ene;
            aver(2) += Mn;
            aver(3) += Ene*Ene;
            aver(4) += Mn*Mn;
        }
    }
    """
    #USE WEAVE TO RUN C++ CODE
    weave.inline(code, ['Nitt','latt','N','PW','Ene','Mn','warm', 'measure', 'aver'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    #GET ALL THE AVERAGED QTIES USIGN THE AVER ARRAY:
    aE = aver[1]/aver[0]
    aM = aver[2]/aver[0]
    cv = (aver[3]/aver[0]-(aver[1]/aver[0])**2)/T**2
    chi = (aver[4]/aver[0]-(aver[2]/aver[0])**2)/T
    return (aM, aE, cv, chi)





if __name__ == '__main__':
    #RANDOM LATTICE:
    latt = RandomL(N)
    #TO TEST YOU CAN USE ONES ARRAY
    #latt=ones(N,N),dtype=int)etc...
    #DECLARE THE ARRAY:
    PW = zeros(9, dtype=float)

    #LOOP IN TEMPERATURE:
    wT = linspace(4,0.5,10)
    WMEAS=[]
    #EMPTY ARRAYS
    wMag=[]
    wEne=[]
    wCv=[]
    wChi=[]
    for T in wT:
        #THESE ARE THE EXPONENTS STORED IN AN ARRAY: 0,2,4,6,8
        # Precomputed exponents
        PW[4+4] = exp(-4.*2/T)
        PW[4+2] = exp(-2.*2/T)
        PW[4+0] = exp(0.*2/T)
        PW[4-2] = exp( 2.*2/T)
        PW[4-4] = exp( 4.*2/T)

        #YOU CAN DO THIS WAY:
        #for i in range(-4,6,2):
        #PW(4+i) = exp(-i*2/T)
        #print PW
        
        #SAMPLE WITH PYTHON!!!
        (aM, aE, cv, chi)=SamplePython(Nitt, latt, PW,T)
        #SAMPLE WITH C++
        # (aM, aE, cv, chi) = SampleCPP(Nitt, latt, PW, T)        
        wMag.append( aM/(N*N) )
        wEne.append( aE/(N*N) )
        wCv.append( cv/(N*N) )
        wChi.append( chi/(N*N) )
        
        print T, aM/(N*N), aE/(N*N), cv/(N*N), chi/(N*N)
        
    plot(wT, wEne, label='E(T)')
    plot(wT, wCv, label='cv(T)')
    plot(wT, wMag, label='M(T)')
    xlabel('T')
    legend(loc='best')
    show()
    plot(wT, wChi, label='chi(T)')
    xlabel('T')
    legend(loc='best')
    show()
    
