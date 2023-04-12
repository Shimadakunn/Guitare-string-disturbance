#Résolution de l'équation de d'Alembert par la méthode des différences finies

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from pylab import*
import time as tm


###################################################################################################################
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
dt=dx/c  #pas dans le temps
Nt=int(Duree//dt)  #Nombre de points dans l'espace
T=linspace(0,Duree,Nt)

#Conditions aux limites
u0l=0     
unl=0

#Parametre de calcul
alpha =(c*dt/dx)    
MODE=3 #Nombre de modes spatiaux

#x0
x0=L/2
A=2
###################################################################################################################
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

def dy0x(i):
  """int ->float
  Fonction de condition initiale d(0,x)"""
  return 0

#Création de la matrice tridiagonale  
def M(n):
  """int -> array(n,n)
 Renvoie matrice diagonale de dimension n*n contenant les coefficients du systeme"""

  M = zeros((n,n),float)

  for i in range(0,n):
    if(i==0):                                   #pour la première ligne de la matrice
      M[i,i], M[i,i+1]=(2-2*alpha**2),alpha**2
    elif(i==n-1):
      M[i,i-1], M[i,i]=alpha**2,(2-2*alpha**2)
    else:                                       #Pour la dernière ligne de la matrice
      M[i,i-1], M[i,i], M[i,i+1] = alpha**2,(2-2*alpha**2),alpha**2

  return M

###Fonctions des conditions initiales de la matrices

#Initialise la matrice colonne avec les valeurs à l'état initial t=0
#Ici,on considère que U(x,0)=sin(x)
def init_U0(n,y):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=0 """

  U0=zeros((n,1),float)

  for i in range(0,n):
    U0[i,0]=y(i)
  #U0[n-1,0]=0
  return U0


#Initialise la matrice colonne avec les valeurs à l'état t=1 
def init_U1(n,dy,y):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=1 """
  #Initialisation des tableaux
  U0=init_U0(n,y)
  U1=zeros((n,1),float)

  for i in range (0,n):
    U1[i]= dt*dy(i) + U0[i]
    #U1[i]= alpha**2*(U0[i-1] - 2*U0[i]+U0[i+1])+2*U0[i]  #Approximation de U1
  #U1[0]=U1[n-1]=0.0 #Conditions initiales à modifier 
  
  return U1

##############################################################################################

start=tm.process_time()

#Création de la matrice triagonale D
D=M(Nx)

#Création de la matrice U des résultats
U=zeros((Nx,Nt),float)

#Conditions initiales
U0=init_U0(Nx,y0x)        #A modifier pour changer le profil initial
U1=init_U1(Nx,dy0x,y0x)   #A modifier pour changer le profil initial
U[:,[0]]=U0
U[:,[1]]=U1


#Matrices des conditions limites
Cl=zeros((Nx,1),float)
Cl[0]=u0l
Cl[Nx-1]=unl

#Calcul du reste des valeurs
for j in range(1,Nt-1):
  U[:,[j+1]]=dot(D,U[:,[j]])-U[:,[j-1]]+Cl
  print("j=",j , "Nt=",Nt-1)


end=tm.process_time()
print("Le temps d'exécution est de ",end-start,"secondes")

#######################################################################################

#Affichage de la solution Sans animation
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position (x)', fontsize = 16)
ax.set_ylabel('temps (t)', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Amplitude en fonction de x et t")
ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,U,cmap = 'viridis')       
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
