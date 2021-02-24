from libc.math cimport exp,pi
#from libcpp.vector cimport vector

cdef exponentFunc(pointRelate , pointMain , double a ,int n , double k , double t):
    sum=0
    cdef int N=n
    for i in range(N):
        sum=sum+(-(pointRelate[i]-pointMain[i])**(n-1)/(t*(2*a)**(n-1))-k**((n-1))*t**1)
    return exp(abs(sum)*(-1))

cdef coefficientFunc(double a,int n,double t):
    #resultTemp=8*np.pi*t*(a**(3*(n-1)/2))*((np.pi*t)**(1/2))
    cdef double resultTemp=0
    cdef double resultTemp2=0
    resultTemp=(2*a*(pi*t)**(1/2))**n
    resultTemp2=1/resultTemp
    return resultTemp2

def GaussianFunc( pointRelate,pointMain,double a,int n,double k,double t):
    return coefficientFunc(a,n,t)*exponentFunc(pointRelate,pointMain,a,n,k,t)