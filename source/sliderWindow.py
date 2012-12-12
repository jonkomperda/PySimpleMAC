from mpi4py import *
from Tkinter import *

comm = MPI.COMM_WORLD
myid = comm.Get_rank()

class App(Frame):
    """docstring for App"""
    def __init__(self, master):
        self.master = master
        Frame.__init__(self, master)
        
        self.slider = Scale(master, from_=100, to=10000, resolution=100,\
                            orient=HORIZONTAL, command= self.updateRe)
        self.slider.pack()
    
    def updateRe(self,master):
        re = self.slider.get()
        comm.isend(obj=re, dest=0, tag=11)
        


if __name__ == '__main__':
    root = Tk()
    app = App(root)
    app.mainloop()
        