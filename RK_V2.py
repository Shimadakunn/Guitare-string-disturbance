#Résolution de l'équation d'onde par Runge Kutta d'ordre 4

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from pylab import*
import pdb
import time as tm

####################################################################################

"""Parametres de discrétisation/Résolution"""
#caractéristique physique 

with open("variables.txt", "r+") as file: #introduction des variables
  dx = float(file.readline()) #pas entre deux points de l'espace
  Duree = float(file.readline()) #Durée de la mesure de l'onde
  c = float(file.readline()) #m.s^-1   célérité de l'onde
  L = float(file.readline()) #Longueur de la corde en m
  file.close()

#Paramètres discrétisation de l'espace - maillage spatiale - indice i
Nx=int(L/dx) #Nombre de points de l'espace
X=linspace(0,L,Nx)

#Paramètre de discrétisation du temps -maillage temporel - indice n
dt=0.000001  #pas dans le temps
Nt=int(Duree/dt)  #Nombre de points dans l'espace
T=linspace(0,Duree,Nt)

#Caractéristique spec
alpha=(c/dx)**2

#Conditions aux limites
u0l=0     
unl=0

MODE=3 #Nombre de modes spatiaux

A=2
x0=L/2

#########################################################################################

#Fonctions des conditions initiales
def y0x(i):
  """int -> float
  Fonction de condition initiale y(0,x)
  onde sinusoïdale"""
  return sin((MODE*pi/L)*(i*dx))

def y0x2(i):
  """int -> float
  Fonction de condition initiale
  onde pincée triangulaire"""

  if((i*dx)<=x0):
    return A*(i*dx)/x0
  
  else:
    return A*(L-i*dx)/(L-x0)


def y0x3(i):
    """int -> float
  Fonction de condition initiale
    Corde pincée triangulaire"""
    
    return A*(i*dx)*(L-i*dx)/(x0*L)


#Initialise la matrice colonne avec les valeurs à l'état initial t=0
#Ici,on considère que U(x,0)=sin(x)
def init_U0(n):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=0 """

  U0=zeros((n,1),float)

  for i in range(0,n):
    U0[i,0]=y0x(i) #Modifier pour changer le profil initial
  return U0
  

#Initialise la matrice colonne avec les conditions initiales de dU/dt en t=0
def init_dU0(n):
    """int -> array(n,1)
    Retourne une matrice colonne composée des valeurs de la dérivée de la fonction d'onde à l'état initial t=0"""
    
    dU0=zeros((n,1),float)
    
    for i in range(0,n):
        dU0[i,0]=0
    
    #U0[n-1,0]=0
    return dU0


"""Définition des fonctions f et g"""
def f(n,i,U):
    """int + int + array(n,1) -> float
    Retourne la valeur de f(n,i,U)"""
    
    if(i==0):
        return alpha* (u0l - 2*U[i,0] + U[i+1,0])
    
    elif(i==n-1):
        return alpha* (U[i-1,0] - 2*U[i,0] + unl)
        
    else:
        return alpha* (U[i-1,0] - 2*U[i,0] + U[i+1,0] )
    
def g(z):
    """float -> float
    Retourne la valeur de g(n,i,U)=g(z)=z"""
    
    return z
 
 
#FONCTIONS F ET G
def K0i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*Z[i,0]
    
def L0i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*f(Nx,i,U)
    
def K1i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*(Z[i,0] + 1/2*L0i(i,U,Z))
    
def L1i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    if(i==0):
        T=alpha* (u0l - 2*(U[i,0] +1/2*K0i(i,U,Z))  + (U[i+1,0]+1/2*K0i(i+1,U,Z)))
        
    elif(i==Nx-1):
        T= alpha* ((U[i-1,0]+1/2*K0i(i-1,U,Z)) - 2*(U[i,0]+1/2*K0i(i,U,Z)) + unl)
        
    else:
        T= alpha* ((U[i-1,0]+1/2*K0i(i-1,U,Z)) - 2*(U[i,0]+1/2*K0i(i,U,Z)) + (U[i+1,0]+1/2*K0i(i+1,U,Z)) )
    
    return dt*T
   
def K2i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*(Z[i,0] + 1/2*L1i(i,U,Z))

def L2i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    if(i==0):
        T=alpha* (u0l - 2*(U[i,0] +1/2*K1i(i,U,Z))  + (U[i+1,0]+1/2*K1i(i+1,U,Z)) )
        
    elif(i==Nx-1):
        T= alpha* ((U[i-1,0]+1/2*K1i(i-1,U,Z)) - 2*(U[i,0]+1/2*K1i(i,U,Z)) + unl)
        
    else:
        T= alpha* ((U[i-1,0]+1/2*K1i(i-1,U,Z)) - 2*(U[i,0]+1/2*K1i(i,U,Z)) + (U[i+1,0]+1/2*K1i(i+1,U,Z)) )
    
    return dt*T
    
def K3i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*(Z[i,0] + L2i(i,U,Z) )   
 
def L3i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    if(i==0):
        T=alpha* (u0l+ - 2*(U[i,0] +K1i(i,U,Z))  + (U[i+1,0]+K1i(i+1,U,Z)) )
        
    elif(i==Nx-1):
        T= alpha* ((U[i-1,0]+K2i(i-1,U,Z)) - 2*(U[i,0]+K2i(i,U,Z)) + unl)
        
    else:
        T= alpha* ((U[i-1,0]+K2i(i-1,U,Z)) - 2*(U[i,0]+K2i(i,U,Z)) + (U[i+1,0]+K2i(i+1,U,Z)) )
    
    return dt*T

#######################################################################################################

start=tm.process_time()
 
#Création de la matrice résultat
U=zeros((Nx,Nt),float)
U[:,[0]]=init_U0(Nx)

z=zeros((Nx,Nt),float)
z[:,[0]]=init_dU0(Nx)

#Paramètres d'itération
t=1
i=0

print(Nx)
print(Nt)
tm.sleep(1)

while(t<Nt):
    i=0
    
    while(i<Nx):
        #calculs de k0 et l0
        k0,l0=K0i(i,U[:,[t-1]],z[:,[t-1]]) , L0i(i,U[:,[t-1]],z[:,[t-1]])
        #print("k0=",k0,"l0=",l0)
        #calculs de k1 et l1
        k1,l1=K1i(i,U[:,[t-1]],z[:,[t-1]]), L1i(i,U[:,[t-1]],z[:,[t-1]])

        #calculs de k2 et l2
        k2,l2=K2i(i,U[:,[t-1]],z[:,[t-1]]), L2i(i,U[:,[t-1]],z[:,[t-1]])
        
        #calculs de k3 et l3
        k3,l3=K3i(i,U[:,[t-1]],z[:,[t-1]]), L3i(i,U[:,[t-1]],z[:,[t-1]])
        #print("k3=",k3,"l3=",l3)
        #Calculs de u et de z pour i et t
        U[i,t]=U[i,t-1]+ 1/6*(k0 + 2*k1 + 2*k2 + k3)
        z[i]=z[i,t-1]+ 1/6*(l0 + 2*l1 + 2*l2 + l3)
        
        i+=1
        #pdb.set_trace()
    print("t=",t,"Nt=",Nt)
    t+=1
     

end=tm.process_time()
print("Le temps d'exécution est de ",end-start,"secondes")

#############################################################################################

# Affichage de la solution Sans animation
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position (x)', fontsize = 16)
ax.set_ylabel('temps (t)', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Amplitude en fonction de x et t")
ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,U,cmap = 'plasma')       
plt.show()

#Animation graphique 2D en fonction de la position
fig = plt.figure()
line, = plot([],[])
plt.xlim(0, L)
plt.ylim(-1, 1)
plt.xlabel("position")
plt.ylabel("amplitude")
plt.title("Graphique 2D en fonction de la positon et à t fixé")
plt.grid("True")

def animate(i):
  line.set_data(X,U[:,[i]])
  return line,

ani = animation.FuncAnimation(fig, animate, frames=100, blit=True, interval=20, repeat=False)
plt.show()

##Animation graphique 2D en fonction du temps
fig = plt.figure()
line, = plot([],[])
plt.xlim(0, 0.001)
plt.ylim(-1, 1)
plt.xlabel("Temps")
plt.ylabel("amplitude")
plt.title("Graphique 2D en fonction du temp et à x fixé ")
plt.grid("True")

def animate(i):
  line.set_data(T,U[i,:])
  return line,

ani = animation.FuncAnimation(fig, animate, frames=100, blit=True, interval=20, repeat=False)
plt.show()




