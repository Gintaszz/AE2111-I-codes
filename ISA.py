from math import exp
#CONSTS

# [[(hmin,hmax),a],....]
layer_data = [[(0,11000),-0.0065],
        [(11000,20000),0],
        [(20000,32000),0.001],
        [(32000,47000),0.0028],
        [(47000,51000),0],
        [(51000,71000),-0.0028],
        [(71000,86000),-0.002]]

sea_level=[288.15,101325]
g=9.80665
R=287

#FUNCS
def ISA(h):
    f=1
    i=0
    pT=[0,sea_level[0]]
    while layer_data[i][0][0]<h: #while in relevant layers
        pT=frac(pT[1],h,layer_data[i]) #[p1/p0,T]
        f=f*pT[0]
        i=i+1
    return [f*sea_level[1],pT[1],(f*sea_level[1]/(R*pT[1]))]    #[p,t,rho]            
                       
def frac(T,h,l_data): #check if in bounds of layer
    #l_data [(hmin,hmax),a]
    if h>l_data[0][1]:
        u=l_data[0][1]-l_data[0][0]
    else:
        u=h-l_data[0][0]
    if l_data[1] == 0:
        return [exp((-(g)/(R*T))*(u)),T] #[p1/p0,T]
    else:
        return [((T+l_data[1]*(u))/T)**((-g)/(l_data[1]*R)),T+l_data[1]*(u)] #[p1/p0,T]
    