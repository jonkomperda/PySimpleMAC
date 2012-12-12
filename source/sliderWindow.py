class App(Frame):
    """
    The main GUI window which controls the plot window. Should be called from a separate 
    MPI process as the plot window.
    """
    def __init__(self, master):
        """
        Initializes the frame with all buttons and parameters. Called upon creation of an App.
        """
        self.master = master
        Frame.__init__(self, master)
        master.title("SimpleMAC GUI")
        
        self.topFrame = LabelFrame(master,text='Reynolds Number',padx=5,pady=5)
        self.topFrame.pack(fill=BOTH,padx=10,pady=10)
        self.botFrame = Frame(master)
        self.botFrame.pack()
        
        self.slider = Scale(self.topFrame, from_=100, to=10000, resolution=100,\
                            orient=HORIZONTAL, command= self.updateRe, width=25)
        self.slider.pack(fill=X)
        
        self.start = Button(self.botFrame,text="Start",command=self.startRun,width=10)
        self.start.pack(side=LEFT)
        self.stop  = Button(self.botFrame,text="Stop", command=self.stopRun,width=10)
        self.stop.pack(side=LEFT)
        self.reset = Button(self.botFrame,text="Reset", command=self.reset,width=10)
        self.reset.pack(side=LEFT)
        
        menubar = Menu(master)
        filemenu = Menu(menubar)
        menubar.add_cascade(label='File', menu=filemenu)
        filemenu.add_command(label='Start simulation', command=self.startRun)
        filemenu.add_command(label='Stop simulation', command=self.stopRun)
        filemenu.add_separator()
        filemenu.add_command(label='Save window...', command=self.saveit)
        filemenu.add_separator()
        filemenu.add_command(label='Reset simulation', command=self.reset)
        filemenu.add_separator()
        filemenu.add_command(label='Exit', command=self.quitIt)
        viewmenu = Menu(menubar)
        menubar.add_cascade(label='View', menu=viewmenu)
        viewmenu.add_command(label='U-Velocity',command=self.viewU)
        viewmenu.add_command(label='V-Velocity',command=self.viewV)
        viewmenu.add_command(label='Pressure',command=self.viewP)
        viewmenu.add_command(label='Vorticity',command=self.viewW)
        master.config(menu=menubar)
    
    def updateRe(self,master):
        """
        Function which communicates with the computational thread to change the Reynolds number
        in the solver. Recieves value from slider and uses non-blocking send to main thread.
        """
        re = self.slider.get()
        comm.isend(obj=re, dest=0, tag=11)
        
    def startRun(self):
        """
        Tells the compute thread to begin computation. Should be called by a button.
        """
        start = True
        comm.isend(obj=start, dest=0, tag=12)
    
    def stopRun(self):
        """
        Tells the compute thread to end computation. Should be called by a button.
        """
        start = False
        comm.isend(obj=start,dest=0, tag=12)
    
    def quitIt(self):
        """
        Sends a quit message to both threads. Proper way to exit the program and ensure both
        threads recieve an exit message.
        """
        comm.isend(obj=False,dest=0, tag=13)
        self.master.quit()
    
    def reset(self):
        """
        Sends reset message to computational thread. Restores initial conditions to the solver.
        """
        comm.isend(obj=True,dest=0,tag=14)
    
    def viewU(self):
        """
        Sets plot to display U-Velocity
        """
        comm.isend(obj='u',dest=0,tag=25)
    
    def viewV(self):
        """
        Sets plot to display V-Velocity
        """
        comm.isend(obj='v',dest=0,tag=25)
    
    def viewP(self):
        """
        Sets plot to display Pressure
        """
        comm.isend(obj='p',dest=0,tag=25)
    
    def viewW(self):
        """
        Sets plot to display Vorticity
        """
        comm.isend(obj='w',dest=0,tag=25)
    
    def saveit(self):
        """
        Uses a communicator to have computational thread save a jpeg
        """
        comm.isend(obj='1',dest=0,tag=44)
