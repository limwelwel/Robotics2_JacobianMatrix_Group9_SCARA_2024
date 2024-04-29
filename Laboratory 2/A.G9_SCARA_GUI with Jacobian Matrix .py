from tkinter import *
from tkinter import messagebox
from tkinter import PhotoImage
import numpy as np
import roboticstoolbox as rtb
from roboticstoolbox import DHRobot, RevoluteDH, PrismaticDH
import spatialmath
from spatialmath import SE3
import matplotlib
matplotlib.use('TkAgg')

# Window and Title
Robotics2_G9_SCARA = Tk()
Robotics2_G9_SCARA.title("SCARA Calculator")
Robotics2_G9_SCARA.resizable(False, False)
Robotics2_G9_SCARA.configure(bg="tan")

def reset():
    a1_E.delete(0, END)
    a2_E.delete(0, END)
    a3_E.delete(0, END)
    a4_E.delete(0, END)
    a5_E.delete(0, END)

    T1_E.delete(0, END)
    T2_E.delete(0, END)
    d3_E.delete(0, END)

    X_E.delete(0, END)
    Y_E.delete(0, END)
    Z_E.delete(0, END)

def f_k():
    # link lengths in cm
    a1 = float(a1_E.get())
    a2 = float(a2_E.get())
    a3 = float(a3_E.get())
    a4 = float(a4_E.get())
    a5 = float(a5_E.get())

    T1 = float(T1_E.get())
    T2 = float(T2_E.get())
    d3 = float(d3_E.get())

    # degree to radian
    T1 = (T1/180.0)*np.pi
    T2 = (T2/180.0)*np.pi

    # Parametric Table (theta, alpha, r, d)
    PT = [[T1,(0.0/180.0)*np.pi,a2,a1],
          [T2,(180.0/180.0)*np.pi,a4,a3],
          [(0.0/180.0)*np.pi,(0.0/180.0)*np.pi,0,a5+d3]]
    
    # HTM formulae
    i = 0
    H0_1 = [[np.cos(PT[i][0]),-np.sin(PT[i][0])*np.cos(PT[i][1]),np.sin(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.cos(PT[i][0])],
            [np.sin(PT[i][0]),np.cos(PT[i][0])*np.cos(PT[i][1]),-np.cos(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.sin(PT[i][0])],
            [0,np.sin(PT[i][1]),np.cos(PT[i][1]),PT[i][3]],
            [0,0,0,1]]

    i = 1
    H1_2 = [[np.cos(PT[i][0]),-np.sin(PT[i][0])*np.cos(PT[i][1]),np.sin(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.cos(PT[i][0])],
            [np.sin(PT[i][0]),np.cos(PT[i][0])*np.cos(PT[i][1]),-np.cos(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.sin(PT[i][0])],
            [0,np.sin(PT[i][1]),np.cos(PT[i][1]),PT[i][3]],
            [0,0,0,1]]

    i = 2
    H2_3 = [[np.cos(PT[i][0]),-np.sin(PT[i][0])*np.cos(PT[i][1]),np.sin(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.cos(PT[i][0])],
            [np.sin(PT[i][0]),np.cos(PT[i][0])*np.cos(PT[i][1]),-np.cos(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.sin(PT[i][0])],
            [0,np.sin(PT[i][1]),np.cos(PT[i][1]),PT[i][3]],
            [0,0,0,1]]
    
    H0_1 = np.matrix(H0_1)
    H1_2 = np.matrix(H1_2)
    H2_3 = np.matrix(H2_3)

    H0_2 = np.dot(H0_1,H1_2)
    H0_3 = np.dot(H0_2,H2_3)

    X0_3 = H0_3[0,3]
    X_E.delete(0,END)
    X_E.insert(0,np.around(X0_3,3))

    Y0_3 = H0_3[1,3]
    Y_E.delete(0,END)
    Y_E.insert(0,np.around(Y0_3,3))

    Z0_3 = H0_3[2,3]
    Z_E.delete(0,END)
    Z_E.insert(0,np.around(Z0_3,3))

    ## Jacobian Matrix Program

    J_sw = Toplevel()
    J_sw.title("Velocity Calculator")
    J_sw.resizable(False,False)

    #1. Linear / Translational Vectors
    Z_1 = [[0],[0],[1]] # The [0,0,1] vector

    #Row 1 to 3, Column 1
    J1a = [[1,0,0],
           [0,1,0],
           [0,0,1]] #R0_0
    J1a = np.dot(J1a,Z_1)

    J1b_1 =H0_3[0:3,3:] 
    J1b_1 = np.matrix(J1b_1)

    J1b_2 = [[0],[0],[0]] # position vector at the base
    J1b_2 = np.matrix(J1b_2)

    J1b = J1b_1 - J1b_2

    J1 = [[(J1a[1,0]*J1b[2,0])-(J1a[2,0]*J1b[1,0])],
          [(J1a[2,0]*J1b[0,0])-(J1a[0,0]*J1b[2,0])],
          [(J1a[0,0]*J1b[1,0])-(J1a[1,0]*J1b[0,0])]]
    J1 = np.matrix(J1)

    #Row 1 to 3, Column 2
    J2a = H0_1[0:3,0:3]
    J2a = np.dot(J2a,Z_1)

    J2b_1 = H0_3[0:3,3:]
    J2b_1 = np.matrix(J2b_1)

    J2b_2 = H0_1[0:3,3:]
    J2b_2 = np.matrix(J2b_2)

    J2b = J2b_1 - J2b_2

    J2 = [[(J2a[1,0]*J2b[2,0])-(J2a[2,0]*J2b[1,0])],
          [(J2a[2,0]*J2b[0,0])-(J2a[0,0]*J2b[2,0])],
          [(J2a[0,0]*J2b[1,0])-(J2a[1,0]*J2b[0,0])]]
    J2 = np.matrix(J2)

    #Row 1 to 3, Column 3
    J3 = H0_2[0:3,0:3]
    J3 = np.dot(J3,Z_1)
    J3 = np.matrix(J3)

    #2. Rotation / Orientation Vectors
 
    #Row 4 to 6, Column 1
    J4 = J1a
    J4 = np.matrix(J4)

    #Row 4 to 6, Column 2
    J5 = J2a
    J5 = np.matrix(J5)

    #Row 4 to 6, Column 3
    J6 = [[0],[0],[0]]
    J6 = np.matrix(J6)

    #3. Concatenated Jacobian Matrix
    JM1 = np.concatenate((J1,J2,J3),1)
    JM2 = np.concatenate((J4,J5,J6),1)

    J = np.concatenate((JM1,JM2),0)
    J = np.matrix(J)

    #Velocity Slide Update
    def update_velocity():
        T1p = T1_slider.get()
        T2p = T2_slider.get()
        d3p = d3_slider.get()

        q = np.array([[T1p],[T2p],[d3p]])
        E = np.dot(J,q)

        xp_e = E[0,0]
        x_entry.delete(0,END)
        x_entry.insert(0,str(xp_e))

        yp_e = E[1,0]
        y_entry.delete(0,END)
        y_entry.insert(0,str(yp_e))

        zp_e = E[2,0]
        z_entry.delete(0,END)
        z_entry.insert(0,str(zp_e))

        ωx_e = E[3,0]
        ωx_entry.delete(0,END)
        ωx_entry.insert(0,str(ωx_e))

        ωy_e = E[4,0]
        ωy_entry.delete(0,END)
        ωy_entry.insert(0,str(ωy_e))

        ωz_e = E[5,0]
        ωz_entry.delete(0,END)
        ωz_entry.insert(0,str(ωz_e))

    #Jacobian Sliders
    T1_velo = Label(J_sw,text=("θ1* ="),font=(5))
    T1_slider = Scale(J_sw,from_=0,to_=3.1416,orient=HORIZONTAL,length=100,sliderlength=10)
    T1_unit = Label(J_sw,text=("rad/s"),font=(5))

    T2_velo = Label(J_sw,text=("θ2* ="),font=(5))
    T2_slider = Scale(J_sw,from_=0,to_=3.1416,orient=HORIZONTAL,length=100,sliderlength=10)
    T2_unit = Label(J_sw,text=("rad/s"),font=(5))

    d3_velo = Label(J_sw,text=("d3* ="),font=(5))
    d3_slider = Scale(J_sw,from_=0,to_=30,orient=HORIZONTAL,length=100,)
    d3_unit = Label(J_sw,text=("cm/s"),font=(5))

    T1_velo.grid(row=0,column=0)
    T1_slider.grid(row=0,column=1)
    T1_unit.grid(row=0,column=2)

    T2_velo.grid(row=1,column=0)
    T2_slider.grid(row=1,column=1)
    T2_unit.grid(row=1,column=2)

    d3_velo.grid(row=2,column=0)
    d3_slider.grid(row=2,column=1)
    d3_unit.grid(row=2,column=2)

    #Jacobian Entries and Labels
    x_velo = Label(J_sw,text=("x* = "),font=(5))
    x_entry = Entry(J_sw,width=10,font=(10))
    x_unit = Label(J_sw,text=("cm/s"),font=(5))
    x_velo.grid(row=3,column=0)
    x_entry.grid(row=3,column=1)
    x_unit.grid(row=3,column=2)

    y_velo = Label(J_sw,text=("y* = "),font=(5))
    y_entry = Entry(J_sw,width=10,font=(10))
    y_unit = Label(J_sw,text=("cm/s"),font=(5))
    y_velo.grid(row=4,column=0)
    y_entry.grid(row=4,column=1)
    y_unit.grid(row=4,column=2)

    z_velo = Label(J_sw,text=("z* = "),font=(5))
    z_entry = Entry(J_sw,width=10,font=(10))
    z_unit = Label(J_sw,text=("cm/s"),font=(5))
    z_velo.grid(row=5,column=0)
    z_entry.grid(row=5,column=1)
    z_unit.grid(row=5,column=2)

    ωx_velo = Label(J_sw,text=("ωx* = "),font=(5))
    ωx_entry = Entry(J_sw,width=10,font=(10))
    ωx_unit = Label(J_sw,text=("rad/s"),font=(5))
    ωx_velo.grid(row=6,column=0)
    ωx_entry.grid(row=6,column=1)
    ωx_unit.grid(row=6,column=2)

    ωy_velo = Label(J_sw,text=("ωy* = "),font=(5))
    ωy_entry = Entry(J_sw,width=10,font=(10))
    ωy_unit = Label(J_sw,text=("rad/s"),font=(5))
    ωy_velo.grid(row=7,column=0)
    ωy_entry.grid(row=7,column=1)
    ωy_unit.grid(row=7,column=2)

    ωz_velo = Label(J_sw,text=("ωz* = "),font=(5))
    ωz_entry = Entry(J_sw,width=10,font=(10))
    ωz_unit = Label(J_sw,text=("rad/s"),font=(5))
    ωz_velo.grid(row=8,column=0)
    ωz_entry.grid(row=8,column=1)
    ωz_unit.grid(row=8,column=2)

    #Update Button
    update_but = Button(J_sw,text="Update",bg="green",fg="white",command=update_velocity)
    update_but.grid(row=9,column=0)

    # Create links
    # [robot_variable]=DHRobot([RevoluteDH(d,r,alpha,offset)])
    SCARA = DHRobot ([
        RevoluteDH(a1/100,0,(0.0/180.0)*np.pi,(0.0/180.0)*np.pi,qlim=[-np.pi/2,np.pi/2]),
        PrismaticDH(0,a2/100,0,0,qlim=[0,0]),
        RevoluteDH(a3/100,0,(180.0/180.0)*np.pi,(0.0/180.0)*np.pi,qlim=[-np.pi/2,np.pi/2]),
        PrismaticDH(0,a4/100,0,0,qlim=[0,0]),
        PrismaticDH(0,0,(0.0/180.0)*np.pi,a5/100,qlim=[0,d3]),
        ], name="SCARA")
    
    #plot joints
    q1 = np.array([T1,0,T2,0,d3/100])

    #plot scale
    x1 = -0.5
    x2 = 0.5
    y1 = -0.5
    y2 = 0.5
    z1 = 0.0
    z2 = 0.5

    #Plot command
    SCARA.plot(q1,limits=[x1,x2,y1,y2,z1,z2],block=True)

def i_k():
# Inverse Kinematics using Graphical Method
    # link lengths in cm
    a1 = float(a1_E.get())
    a2 = float(a2_E.get())
    a3 = float(a3_E.get())
    a4 = float(a4_E.get())
    a5 = float(a5_E.get())

    # Position Vector in cm
    xe = float(X_E.get())
    ye = float(Y_E.get())
    ze = float(Z_E.get())

    #try & except
    try:
        phi2 = np.arctan(ye/xe)
    except:
        phi2 = -1 #NAN Error
        messagebox.showerror(title="DivideZero Error", message="Undefine error if X=0")
    
    # To solve for Theta 1 or Th1
    phi2 = np.arctan(ye/xe) #1
    r1 = np.sqrt(ye**2 + xe**2) #2
    phi1 = np.arccos((a4**2 - a2**2 - r1**2) / (-2 * a2 * r1)) #3
    T1 = phi2 - phi1 #4

    # To solve for Theta 2 or Th2
    phi3 = np.arccos((r1**2 - a2**2 - a4**2) / (-2 * a2 * a4)) #5
    T2 = np.pi - phi3 #6

    # To solve for d3
    d3 = a1 + a3 - a5 - ze #7

    T1_E.delete(0,END)
    T1_E.insert(0,np.around(T1*180/np.pi,3))

    T2_E.delete(0,END)
    T2_E.insert(0,np.around(T2*180/np.pi,3))

    d3_E.delete(0,END)
    d3_E.insert(0,np.around(d3,3))

    # Create links
    # [robot_variable]=DHRobot([RevoluteDH(d,r,alpha,offset)])
    SCARA = DHRobot ([
        RevoluteDH(a1/100,0,(0.0/180.0)*np.pi,(0.0/180.0)*np.pi,qlim=[-np.pi/2,np.pi/2]),
        PrismaticDH(0,a2/100,0,0,qlim=[0,0]),
        RevoluteDH(a3/100,0,(180.0/180.0)*np.pi,(0.0/180.0)*np.pi,qlim=[-np.pi/2,np.pi/2]),
        PrismaticDH(0,a4/100,0,0,qlim=[0,0]),
        PrismaticDH(0,0,(0.0/180.0)*np.pi,a5/100,qlim=[0,d3]),
        ], name="SCARA")
    
    #plot joints
    q1 = np.array([T1,0,T2,0,d3/100])

    #plot scale
    x1 = -0.5
    x2 = 0.5
    y1 = -0.5
    y2 = 0.5
    z1 = 0.0
    z2 = 0.5

    #Plot command
    SCARA.plot(q1,limits=[x1,x2,y1,y2,z1,z2],block=True)

# Link Lengths and Joint Variables Frame
FI = LabelFrame(Robotics2_G9_SCARA,text="Link Lengths and Joint Variables",font=(S),bg="tan",fg="black")
FI.grid(row=0,column=0)

# link lengths label
a1 = Label(FI,text="a1 = ",font=(10),bg="antiquewhite",fg="black")
a1_E = Entry(FI,width=5,font=(10),bg="white",fg="black")
cm1 = Label(FI,text="cm",font=(10),bg="wheat",fg="black")

a2 = Label(FI,text="a2 = ",font=(10),bg="antiquewhite",fg="black")
a2_E = Entry(FI,width=5,font=10,bg="white",fg="black")
cm2 = Label(FI,text="cm",font=(10),bg="wheat",fg="black")

a3 = Label(FI,text="a3 = ",font=(10),bg="antiquewhite",fg="black")
a3_E = Entry(FI,width=5,font=10,bg="white",fg="black")
cm3 = Label(FI,text="cm",font=(10),bg="wheat",fg="black")

a4 = Label(FI,text="a4 = ",font=(10),bg="antiquewhite",fg="black")
a4_E = Entry(FI,width=5,font=10,bg="white",fg="black")
cm4 = Label(FI,text="cm",font=(10),bg="wheat",fg="black")

a5 = Label(FI,text="a5 = ",font=(10),bg="antiquewhite",fg="black")
a5_E = Entry(FI,width=5,font=10,bg="white",fg="black")
cm5 = Label(FI,text="cm",font=(10),bg="wheat",fg="black")

a1.grid(row=0,column=0)
a1_E.grid(row=0,column=1)
cm1.grid(row=0,column=2)

a2.grid(row=1,column=0)
a2_E.grid(row=1,column=1)
cm2.grid(row=1,column=2)

a3.grid(row=2,column=0)
a3_E.grid(row=2,column=1)
cm3.grid(row=2,column=2)

a4.grid(row=3,column=0)
a4_E.grid(row=3,column=1)
cm4.grid(row=3,column=2)

a5.grid(row=4,column=0)
a5_E.grid(row=4,column=1)
cm5.grid(row=4,column=2)

# Joint Variables label
T1 = Label(FI,text=" T1 = ",font=(10),bg="antiquewhite",fg="black")
T1_E = Entry(FI,width=5,font=10,bg="white",fg="black")
deg1 = Label(FI,text="deg",font=(10),bg="wheat",fg="black")

T2 = Label(FI,text=" T2 = ",font=(10),bg="antiquewhite",fg="black")
T2_E = Entry(FI,width=5,font=10,bg="white",fg="black")
deg2 = Label(FI,text="deg",font=(10),bg="wheat",fg="black")

d3 = Label(FI,text=" d3 = ",font=(10),bg="antiquewhite",fg="black")
d3_E = Entry(FI,width=5,font=10,bg="white",fg="black")
cm6 = Label(FI,text="cm",font=(10),bg="wheat",fg="black")

T1.grid(row=0,column=3)
T1_E.grid(row=0,column=4)
deg1.grid(row=0,column=5)

T2.grid(row=1,column=3)
T2_E.grid(row=1,column=4)
deg2.grid(row=1,column=5)

d3.grid(row=2,column=3)
d3_E.grid(row=2,column=4)
cm6.grid(row=2,column=5)

# Button Frame
BF = LabelFrame(Robotics2_G9_SCARA,text="Forward Kinematics",font=(5),bg="tan",fg="black")
BF.grid(row=1,column=0)

# Buttons
FK = Button(BF,text="Forward",font=(10),bg="cornflowerblue",fg="black",command=f_k)
rst = Button(BF,text="RESET",font=(10),bg="darkkhaki",fg="black",command=reset)
IK = Button(BF,text="↑ Inverse",font=(10),bg="plum",fg="black",command=i_k)

FK.grid(row=0,column=0)
rst.grid(row=0,column=1)
IK.grid(row=0,column=2)

# Position Vectors Frame
PV = LabelFrame(Robotics2_G9_SCARA,text="Position Vectors",font=(5),bg="tan",fg="black")
PV.grid(row=2,column=0)

# Position Vectors Label
X = Label(PV,text="X = ",font=(10),bg="antiquewhite",fg="black")
X_E = Entry(PV,width=5,font=10,bg="white",fg="black")
cm7 = Label(PV,text="cm",font=(10),bg="wheat",fg="black")

Y = Label(PV,text="Y = ",font=(10),bg="antiquewhite",fg="black")
Y_E = Entry(PV,width=5,font=10,bg="white",fg="black")
cm8 = Label(PV,text="cm",font=(10),bg="wheat",fg="black")

Z = Label(PV,text="Z = ",font=(10),bg="antiquewhite",fg="black")
Z_E = Entry(PV,width=5,font=10,bg="white",fg="black")
cm9 = Label(PV,text="cm",font=(10),bg="wheat",fg="black")

X.grid(row=0,column=0)
X_E.grid(row=0,column=1)
cm7.grid(row=0,column=2)

Y.grid(row=1,column=0)
Y_E.grid(row=1,column=1)
cm8.grid(row=1,column=2)

Z.grid(row=2,column=0)
Z_E.grid(row=2,column=1)
cm9.grid(row=2,column=2)

#display image
img = PhotoImage(file="SCARA.png")
img = img.subsample(3,3)
PI = Label(Robotics2_G9_SCARA,image=img)
PI.grid(row=3,column=0)

Robotics2_G9_SCARA.mainloop()
