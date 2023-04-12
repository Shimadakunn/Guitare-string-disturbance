from tkinter import *
import matplotlib.pyplot as plt

def explicite():
    dx = dx_entry.get()
    dt = dt_entry.get()
    c = c_entry.get()
    L = L_entry.get()
    with open("variables.txt", "w+") as file:
        file.write(dx + "\n")
        file.write(dt + "\n")
        file.write(c + "\n")
        file.write(L + "\n")
    exec(open(".\Euler_explicite.py").read(), globals())

def implicite():
    dx = dx_entry.get()
    dt = dt_entry.get()
    c = c_entry.get()
    L = L_entry.get()
    with open("variables.txt", "w+") as file:
        file.write(dx + "\n")
        file.write(dt + "\n")
        file.write(c + "\n")
        file.write(L + "\n")
    exec(open(".\Euler_implicite.py").read(), globals())

def rk():
    dx = dx_entry.get()
    dt = dt_entry.get()
    c = c_entry.get()
    L = L_entry.get()
    with open("variables.txt", "w+") as file:
        file.write(dx + "\n")
        file.write(dt + "\n")
        file.write(c + "\n")
        file.write(L + "\n")
    exec(open(".\RK_V2.py").read(), globals())

window = Tk()
#Setting up window
window.title("Euler explicite")
window.geometry("480x360")
window.config(background='black')

#creer la frame
frame = Frame(window, bg ='black')

#Titre
label_title = Label(frame, text = "Variations des conditions", font=("Arial", 18), bg='black', fg='white')
label_title.pack()

#dx à faire varier
dx_label = Label(frame, text = "dx", font=("Arial", 12), bg='black', fg='white')
dx_label.pack()
dx_entry = Entry(frame, text = "dx", font=("Arial", 12), bg='white', fg='black')
dx_entry.pack()

#dt à faire varier
dt_label = Label(frame, text = "dt", font=("Arial", 12), bg='black', fg='white')
dt_label.pack()
dt_entry = Entry(frame, text = "dt", font=("Arial", 12), bg='white', fg='black')
dt_entry.pack()

#c à faire varier
c_label = Label(frame, text = "c", font=("Arial", 12), bg='black', fg='white')
c_label.pack()
c_entry = Entry(frame, text = "c", font=("Arial", 12), bg='white', fg='black')
c_entry.pack()

#L à faire varier
L_label = Label(frame, text = "L", font=("Arial", 12), bg='black', fg='white')
L_label.pack()
L_entry = Entry(frame, text = "L", font=("Arial", 12), bg='white', fg='black')
L_entry.pack()

#button
bouton1 = Button(frame, text='explicite', font=("Arial",12), command = explicite)
bouton2 = Button(frame, text='implicite', font=("Arial",12), command = implicite)
bouton3 = Button(frame, text='RK', font=("Arial",12), command = rk)
bouton1.pack()
bouton2.pack()
bouton3.pack()

frame.pack(expand=YES)

window.mainloop()
