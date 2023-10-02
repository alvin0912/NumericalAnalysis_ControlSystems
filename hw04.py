from PIL import Image
import numpy as np
from matplotlib import pyplot as plt

def printImg(img):
    plt.figure()
    plt.imshow(img)
    plt.axis('off')
    plt.show()

class TransferFunction:
    def show(self):
        for i in range(len(self.Numerator)):
            print(self.Numerator[i],' s^',len(self.Numerator)-i-1, end=' ',sep='')
        print('')
        print('----------------------------')
        for i in range(len(self.Denominator)):
            print(self.Denominator[i],'s^',len(self.Denominator)-i-1, end=' ',sep='')
        print('')

    def poly_mul(self,p1,p2):
        self.tmplen=(len(p1)-1)+(len(p2)-1)+1
        self.tmpdata=np.zeros(self.tmplen)
        self.pos=0
        for i in range(len(p1)):
                for j in range(len(p2)):
                    self.pos=self.tmplen-(len(p1)-i+len(p2)-j-2)-1
                    self.tmpdata[self.pos]=self.tmpdata[self.pos]+p1[i]*p2[j]
        return self.tmpdata
    
    def poly_plus(self,p1,p2):
        self.tmplen=len(p1)if len(p1)>len(p2) else len(p2)
        p2=np.concatenate([np.zeros(len(p1)-len(p2)),p2])
        self.tmpdata=np.zeros(self.tmplen)
        for i in range(self.tmplen):
            self.tmpdata[i]=p1[i]+p2[i]
        return self.tmpdata

    def Series(self,G2):
        self.Numerator=self.poly_mul(self.Numerator,G2.Numerator)
        self.Denominator=self.poly_mul(self.Denominator,G2.Denominator)

    def Parallel(self,G2):
        self.Numerator=self.poly_plus(self.poly_mul(self.Numerator,G2.Denominator),\
            self.poly_mul(self.Denominator,G2.Numerator))
        self.Denominator=self.poly_mul(self.Denominator,G2.Denominator)

    def Nfeedback(self,H):
        self.tmpPoly=self.poly_mul(self.Numerator,H.Denominator)
        self.Denominator=self.poly_plus(self.poly_mul(self.Denominator,H.Denominator),\
            self.poly_mul(self.Numerator,H.Numerator))
        self.Numerator=self.tmpPoly
            
    def __init__(self,*args):
        if len(args)==2:
            self.Numerator=np.ones(1)*float(args[0])
            self.Denominator=np.ones(1)*float(args[1])
        elif len(args)==1:
            self.Numerator=np.ones(1)*float(args[0])
            self.Denominator=np.ones(1)
        else:
            for j in range(2):
                print('input Numerator') if j==0 else print('input Denominator')
                self.input=input()
                self.start=0
                self.end=0
                self.num_data=0
                self.data=np.zeros(20)
                for i in range(20):
                    self.end=self.input[self.start:].find(' ')+self.start
                    if self.end==self.start-1:
                        self.data[i]=self.input[self.start:]
                        self.num_data=self.num_data+1
                        break
                    self.data[i]=float(self.input[self.start:self.end])
                    self.start=self.end+1
                    self.num_data=self.num_data+1
                if j==0:
                    self.Numerator=self.data[:self.num_data]
                else:
                    self.Denominator=self.data[:self.num_data]

class System:
    def __init__(self):
        print('new System')
        self.T=TransferFunction(1,1)
    
    def TF_SS(self):
        self.num_x=len(self.T.Denominator)-1
        self.A=np.zeros([self.num_x,self.num_x])
        for i in range(self.num_x-1):
            self.A[i][i+1]=1
        for i in range(self.num_x):
            self.A[self.num_x-1][i]=self.T.Denominator[self.num_x-i]*-1

        self.B=np.zeros(self.num_x)
        self.B[self.num_x-1]=self.T.Numerator[0]

        self.C=np.zeros(self.num_x)
        self.C[0]=1
    
    def setInput(self,t):
        print('Select Input')
        print('(1)for unit step')
        print('(2)for unit ramp')
        arg=input()
        if arg=='1':
            if len(self.T.Numerator)==2:
                self.U=np.zeros(len(t))
            else:
                self.U=np.ones(len(t))
        else :
            if len(self.T.Numerator)==2:
                self.U=np.ones(len(t))
            else:
                self.U=t
    
    def f(self,X,t):
        return self.A.dot(X)+self.B*t

    def EulerMethod(self):
        X0=np.zeros(self.num_x)
        start=0
        end=self.end
        n=100
        h=(end-start)/n
        t=np.linspace(start,end,n)
        y=np.zeros(n)
        X=np.zeros([n,self.num_x])
        X[0,:] = X0
        self.setInput(t)

        for i in range(1,n):
            X[i,:] = X[i-1,:]+h*self.f(X[i-1],self.U[i-1])
            y[i]=self.C.dot(X[i,:])

        plt.plot(t,y,'o')
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        plt.title("Approximation Solution with Euler's Method")
        plt.grid()
        plt.show()

    def run(self):
        self.setTransferFunction()
        self.TF_SS()
        self.EulerMethod()
        
    
class HWsystem(System):
    def setTransferFunction(self):
        self.end=0.05
        print('Input G:')
        G=TransferFunction()
        H=TransferFunction(1,1)
        print('Input k')
        k=float(input())
        K=TransferFunction(k)
        G.Series(K)
        G.Nfeedback(H)
        self.T=G

class SIMsystem(System):
    def setTransferFunction(self):
        self.end=6
        print('Input G:')
        G=TransferFunction()
        self.T=G

    def f(self,X,t):
        return self.A.dot(X)+self.B*t

class SIMsystem_withC(System):
    def setTransferFunction(self):
        self.end=6
        print('Input G:')
        G=TransferFunction()
        G.show()
        print('Input C:')
        C=TransferFunction()
        G.Series(C)
        self.T=G

    def f(self,X,t):
        return self.A.dot(X)+self.B*t

print('Please select a system:')
img=Image.open('System.jpg')
printImg(img)
print('(1)for System 1')
print('(2)for System 2')
print('(3)for System 3')
arg=input()

if arg=='1':
    s1=HWsystem()
    s1.run()
elif arg=='2':
    s1=SIMsystem()
    s1.run()
else:
    s3=SIMsystem_withC()
    s3.run()