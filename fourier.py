"""
Author: A.Alexeevski
"""

import numpy as np
import math
import random
import matplotlib.pyplot as plt
import sys

class Fourier(object):
    """
    
    """
    def __init__(self):
        return(None)

    def compile_func(self, parameters, interval):
        """
        Compile a function as a sum of gausses
        parameters = [ (lambda, beta, gamma), ...]
        gauss = lambda*exp(-(beta^2)*(X-gamma)^2) here x^2 = x*x
        interval = [from,to,number of poins]        
        """
        
        X = np.array( np.arange(interval[0],interval[1],(interval[1]-interval[0])/float(interval[2]) ) )
        Y = np.array(len(X)*[0.0])
        for z in parameters:            
            Y += float(z[0])*np.exp(-(float(z[1])**2)*(X - len(X)*[float(z[2])])** 2) 
        return (X,Y)

    def write_func(self,X,Y,file):
        file.write("#X\tY\n")
        for i in range(len(X)):
            file.write("%s\t%s\n" % (X[i], Y[i] ) )

    def read_func(self, file):
            """
            """
            X=[]
            Y=[]
            for line in file:                      
                if len(line.strip() ) == 0:
                    continue
                if line.strip()[0] == "#":
                    continue
                line_split = line.strip().split()
                try:
                    X.append( float(line_split[0]) )
                    Y.append( float(line_split[1]) )
                except:
                    raise Exception("Wrong line in file")
                    sys.exit()
                
            return(np.array(X),np.array(Y) )

    def plot_func(self, X, Y, picture):
        plt.plot(X,Y, 'k')
        plt.savefig(picture)
        
    def fourier(self, Y):
        """        
        f = [ (k, F_k, phi_k), ... ] is a full ordered list of  Fourier coefficients 
        """
        n = len(Y)
        F=np.fft.fft(Y)

        f=[]
        
        f.append( (0,F[0].real/n,0)  )
        
        for i in range(1,int(n/2)-1):
            A =  F[i].real/n 
            B =  F[i].imag/n 
            f.append( (i, math.sqrt(A**2+B**2), np.arctan2(-B,A) )   )
            
        return(f)

    def handmade_f(self, Y): 
        """
        len(Y) = 2n+1 (required)
        f(x) is defined at [0, 2pi] in (2n+1) points with step 2pi/(2n+1)
        f(x) = A_0 + 2*sum k_0^N (A_k*cos(kx) + B_k*sin(kx))
        A_k = 1/(2n+1) sum s_0^2n (f(2pi*s/(2n+1) ) cos(  (2pi*k*/(2n+1))*s ) k = 0, ..., n
        B_k = 1/(2n+1) sum s_0^2n (f(2pi*s/(2n+1) ) sin(  (2pi*k*/(2n+1))*s ) k = 0, ..., n
        """
        f = []
        n = int((len(Y)-1)/2)
        for k in range(n+1):
            A=0
            B=0
            for s in range(2*n+1):
                A += (1/float(2*n+1))*Y[s]*np.cos( 2*np.pi*k*s/float(2*n+1) )
                B += (1/float(2*n+1))*Y[s]*np.sin( 2*np.pi*k*s/float(2*n+1) )
            if k==0:
                f.append( (0,A,0) )
                continue
            f.append(   ( k, math.sqrt(A**2+B**2), np.arctan2(B,A) )   )         
        return(f)
                      
    def write_f(self, f, file):
        file.write("#no.\tF\tphi\n")
        for x in f:           
            file.write("%s\t%s\t%s\n" % (x[0], x[1], x[2]) )
    def randomize_F(self, ft, sigma_percent):
        """
        Add gaussian noice with sigma = sigma_percet of ampl. F 
        """
        for i in range(len(ft)):            
            new_F= random.gauss(ft[i][1],(sigma_percent/float(100))*ft[i][1] )
            ft[i] = (ft[i][0],new_F,ft[i][2])
        return(ft)    

    def randomize_phi(self, ft, sigma_percent):
        """
        Add gaussian noice with sigma = sigma_percet of phase phi 
        """
        for i in range(len(ft)):            
            new_phi = random.gauss(ft[i][2],(sigma_percent/float(100))*ft[i][2] )
            ft[i] = (ft[i][0],ft[i][1],new_phi)
        return(ft)    
            
    def read_f(self, file):
        """
        """
        f=[]
        for line in file:                      
            if len(line.strip() ) == 0:
                continue
            if line.strip()[0] == "#":
                continue
            line_split = line.strip().split()
            try:
                f.append( (int(line_split[0]), float(line_split[1]), float(line_split[2]) ) )
            except:
                raise Exception("Wrong line in file")
        return(f)

    def reduce(self, f, garmonics):
        """
        garmonics = [n_0,n_2,...]
        """
        garmonics.sort()
        if garmonics[-1] > len(f)-1:
            raise Exception("n of garmonics > len(fourier)")
            sys.exit()
        f_reduced=[]
        for k in garmonics:
            f_reduced.append(f[k])          
        return(f_reduced)
                
    def recover(self, X, ft_reduced):
        """
        """
        Yrec =[]
        n = len(X)
        
        for s in range(n):
            y=0
            for i in range(len(ft_reduced)):
                z=ft_reduced[i]
                k = z[0]
                
                if k == 0:
                    y += z[1]                    
                else:                    
                    y +=  2*z[1]*np.cos(2*np.pi*(s/float(n))*k-z[2])
                #print m,n,k, y
            #print y
            Yrec.append(y)
        return( np.array(Yrec) )


    def write_two(self,X,Y1,Y2,file):
        """
        Write two functions Y1, Y2 in the same points x from X
        """
        file.write("#X\t  F1\t  F2\n")
        for i in range(len(X)):           
            file.write("%s\t%s\t%s\n" % (X[i], Y1[i], Y2[i]) )
            
    def plot_func_1(self, X, Y):
        plt.plot(X,Y, 'k')
        plt.show()            
        
    
    def plot_2(self, X, Y, Yrec, picture1):
        plt.plot(X,Y, 'k', X, Yrec, 'k:')
        plt.savefig(picture1)

