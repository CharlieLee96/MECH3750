# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 13:54:56 2018

@author: s4395102
"""
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
#data=[]
#with open('data2.txt','r') as dat:
#    for line in dat:
#        data.append(float(line))
#    print(len(data))
def plot_fft(text):
    """applies dft to a set of data points and shows the result on polar plane"""
    data = np.loadtxt(text)
    plt.figure(1)
    plt.plot(data,'-')
    plt.show()
    ft = np.fft.fft(data)
    rl = ft.real

    img = ft.imag
    print(img)
    plt.scatter(np.linspace(1,160,160,endpoint=True),list(img))
    plt.ylabel('imaginary part of fourier series')
#    plt.xlabel('real part of fourier series')
    plt.title('imaginary solution of fourier series')
    plt.show()
    
    plt.scatter(np.linspace(1,160,160,endpoint=True),list(rl))
    plt.ylabel('real part of fourier series')
#    plt.xlabel('real part of fourier series')
    plt.title('real solution of fourier series')
    plt.show()
    
    plt.scatter(list(rl),list(img))
    plt.ylabel('imaginary part of fourier series')
    plt.xlabel('real part of fourier series')
    plt.title('real vs imaginary parts of series')
    plt.show()
    mag=[]
    
    for i in ft:
        mag.append(np.absolute(i)) #magnitude of each signal
#    print(mag)
    mag2=mag[0:80] #chopped up version of whole because it just mirrors itself after 80th point.

#    ano = np.ma.anomalies(np.asarray(mag))
#    print ("ano ",ano)
#    print(len(ano))
    std=np.std(mag)
#    print("std ", std)
    outliers={}
    """figuring out outliers"""
    for i in mag2:
        if i > (stat.median(mag2))+1.5*std:
            outliers[str(mag2.index(i))]=i
        elif i < (stat.median(mag2))-1.5*std:
            outliers[str(mag2.index(i))]=i
    
#    plt.xlim((min(rl)*1.5,max(rl)*1.5))
#    plt.ylim((min(img)*1.5,max(img)*1.5))
#    plt.ylabel('Imaginary')
#    plt.xlabel('Real')
    plt.polar(rl,img,'go')
    plt.title('polar imagery of each signal')
    plt.show()
    print(len(mag2))
    plt.ylabel('Sound Magnitude')
    plt.xlabel('Frequency (Hz)')
    plt.title('frequency vs magnitude of each signal')
    plt.plot(np.linspace(1,len(mag2)/40,len(mag2)),mag2)
    plt.show()
    print(outliers)
 
    """compression of signals"""
    mag3=[] #compressed list, clipping out all weak signals
    mag3_index = []
    for i in mag2:
        if i >= 0.1*max(mag2):
            mag3.append(i)
            mag3_index.append(mag2.index(i))
    plt.ylabel('Sound Magnitude')
    plt.xlabel('Frequency (Hz)')
    print(mag3_index)
    plt.title('compressed frequency vs magnitude')
    mag3space=[]
    for i in mag3_index:
        mag3space.append(i/40)
    print(mag3space)
    plt.plot(mag3space,mag3)
    plt.show()
    return ft,rl,img,mag2,outliers
a=plot_fft('data2.txt')

a
#print(a[0])