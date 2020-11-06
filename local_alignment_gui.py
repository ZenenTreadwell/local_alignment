#!/usr/bin/python

import os
import numpy as np
from tkinter import filedialog
from tkinter import *

# Algorithm Execution Variables
S1=""
S2=""
match=0
mismatch=0
gapOpen=0
gapExtend=0
scoreMatrix=None
gapMatrix=None
opened=False
highestPos=0

#Calculates the score of an entry in the matrix
def score (matrix,gapMatrix, i, j):
	global opened
	#Flag for if it was a match or mismatch
	x=0 
	if S1[i-1]==S2[j-1]:
		x=match
	else:
		x=mismatch
	extendedL=False
	extendedU=False
	upLeft = matrix[i - 1][j - 1] + x
	#Were just opening a gap now
	if gapMatrix[i-1][j]==0:
		up   = matrix[i - 1][j] + gapOpen
	else:
		up   = matrix[i - 1][j] + gapExtend
		extendedU=True
	if gapMatrix[i][j-1]==0:
		left = matrix[i][j-1] + gapOpen
	else:
		#A gap is already open
		left = matrix[i][j-1] + gapExtend
		extendedL=True
	#pull from second matrix at smae coordinates to determine if its an extension or new gap
	if max(0,upLeft,up,left)==up  and extendedU==True:
		#Extended a gap
		return max(0,upLeft,up,left), 1
	elif  max(0,upLeft,up,left)==left and extendedL==True:
		return max(0,upLeft,up,left), 2
	elif (max(0,upLeft,up,left)==up or max(0,upLeft,up,left)==left) and extendedL==False and extendedU==False:
		#New gap added
		return max(0,upLeft,up,left), 3
	return max(0,upLeft,up,left),0

# Creates the matrix of distance scores
def makeMatrix(row,col):
	matrix = [[0 for i in range(col)] for j in range(row)]
	gapMatrix = [[0 for i in range(col)] for j in range(row)]
	maxScore=0
	k=(0,0)
	for i in range (1,row):
		for j in range (1, col):
			nextt=score(matrix,gapMatrix,i,j)
			if nextt[0]> maxScore:
				maxScore=nextt[0]
				k=(i,j)
			matrix[i][j] = nextt[0]
			gapMatrix[i][j]=nextt[1]
	#print (np.matrix(gapMatrix))
	return matrix, k, gapMatrix

def localAlignment(sequence1,sequence2,matchValue,mismatchValue,gapOvalue,gapEvalue):
	global S1
	global S2
	global match
	global mismatch
	global gapOpen
	global gapExtend
	global scoreMatrix
	global gapMatrix
	global highestPos
	S1=sequence1
	S2=sequence2
	match=matchValue
	mismatch=mismatchValue
	gapOpen=gapOvalue
	gapExtend=gapEvalue
	scoreMatrix,highestPos, gapMatrix = makeMatrix(len(S1)+1,len(S2)+1)

def traceback():
    row = highestPos[0]
    col=highestPos[1]
    ans1=""
    ans2=""
    ans3=""
    nextt=0
    # print (highestPos)
    # print (np.matrix(scoreMatrix))
    while scoreMatrix[row][col]!=0:
    	flag1=0
    	flag2=0
    	flag3=0
    	flag4=0
    	flag5=0
    	flag6=0
    	if scoreMatrix[row][col]==scoreMatrix[row-1][col-1]+match:
    		flag1=scoreMatrix[row-1][col-1]
    	if scoreMatrix[row][col]== scoreMatrix[row-1][col-1]+mismatch:
    		flag2=scoreMatrix[row-1][col-1]
    	if scoreMatrix[row][col]==scoreMatrix[row-1][col]+gapOpen:
    		flag3=scoreMatrix[row-1][col]
    	if scoreMatrix[row][col]==scoreMatrix[row][col-1]+gapOpen:
    		flag4=scoreMatrix[row][col-1]
    	if scoreMatrix[row][col]==scoreMatrix[row-1][col]+gapExtend and gapMatrix[row-1][col]==1:
    		flag5=scoreMatrix[row-1][col]
    	if scoreMatrix[row][col]==scoreMatrix[row][col-1]+gapExtend and gapMatrix[row][col-1]==2:
    		flag6=scoreMatrix[row][col-1]
    	maxx=max(flag1,flag2,flag3,flag4,flag5,flag6)
    	if maxx==flag1:
    		ans1+=S1[row-1]
    		ans2+=S2[col-1]
    		row-=1
    		col-=1
    		ans3+="|"
    	elif maxx==flag2:		
    		ans1+=S1[row-1]
    		ans2+=S2[col-1]
    		row-=1
    		col-=1
    		ans3+=" "
    	elif maxx==flag3:
    		#If we backtrace up
    		ans1+=S1[row-1]
    		ans2+="-"
    		ans3+=" "
    		row-=1
    	elif maxx==flag4:
    		#If we backtrace up
    		ans2+=S2[col-1]
    		ans1+="-"
    		ans3+=" "
    		col-=1
    	elif maxx==flag5:
    		#Loop and add until we hit the gap open
    		while gapMatrix[row][col]==1:
    			#If we backtrace up
    			ans1+=S1[row-1]
    			ans2+="-"
    			ans3+=" "
    			row-=1
    	elif maxx ==flag6:
    		while gapMatrix[row][col]==2:
    			#If we backtrace up
    			ans2+=S2[col-1]
    			ans1+="-"
    			ans3+=" "
    			col-=1

    output = {
        's1_opt' : ans1[::-1],
        'mid' : ans3[::-1],
        's2_opt' : ans2[::-1]
    }

    return output

# Master GUI Class
class GUI:
    def __init__(self, master):
        self.master = master
        master.title("Sequence Local Alignment with Affine Gap")
        #master.geometry("550x250")

        self.upperFrame = Frame(master,bd = 10)
        self.input1 = Input(self.upperFrame,1)
        self.input2 = Input(self.upperFrame,2)
        self.upperFrame.pack()

        self.lowerFrame = Frame(master,bd=10)

        self.affineGapBox = AffineGetBox(self.lowerFrame)
        self.affineGapBox.frame.pack(side=LEFT,anchor=E)

        self.buttonFrame = Frame(master, bd = 10)
        self.show_matrix_button = Button(self.buttonFrame, text = "Show Matrix", command = self.matrixWindow)
        self.show_matrix_button.pack()
        self.show_gap_matrix_button = Button(self.buttonFrame, text = "Show Gap Matrix", command = self.gapMatrixWindow)
        self.show_gap_matrix_button.pack()
        self.show_alignment_button = Button(self.buttonFrame, text = "Show Optimal Alignment", command = self.resultWindow)
        self.show_alignment_button.pack()

        self.buttonFrame.pack(side=RIGHT)

        self.compute_button = Button(self.lowerFrame, text="COMPUTE", command=self.compute,bg="red",height=5,width=8)
        self.compute_button.pack(side=RIGHT)
        self.lowerFrame.pack(fill=BOTH)

        self.output1 = None
        self.align = None
        self.output1 = None

    def compute(self):
        input1 = self.input1.entry.get()
        input2 = self.input2.entry.get()

        if (len(input1) < len(input2)): # The larger string should be input1
            temp = input1
            input1 = input2
            input2 = temp

        specs = self.affineGapBox.getVals()

        localAlignment(input1, input2, specs['match'],-specs['mismatch'],-specs['gap_open'],-specs['gap_extend'])
        output_strings = traceback()
        self.output1 = output_strings['s1_opt']
        self.align = output_strings['mid']
        self.output2 = output_strings['s2_opt']

        self.resultWindow()

    def matrixWindow(self):
        subroot = Toplevel(self.master)
        if (scoreMatrix != None):
            text = np.matrix(scoreMatrix)
            window = TextWindow(subroot,"Score Matrix",text)

    def gapMatrixWindow(self):
        subroot = Toplevel(self.master)
        if (gapMatrix != None):
            text = np.matrix(gapMatrix)
            window = TextWindow(subroot,"Gap Matrix",text)
        
    def resultWindow(self):
        subroot = Toplevel(self.master)
        if (self.output1 != None):
            window = TextWindow(subroot,"Optimal Local Alignment",self.output1+'\n'+self.align+'\n'+self.output2)

class TextWindow:
    def __init__(self,master,title,txt):
        self.master = master
        master.title(title)
        self.label = Label(self.master,text=txt, font=("Courier",12))
        self.label.pack()

class Input:
    def __init__(self,master,num):
        self.master = master
        self.frame = Frame(master)
        
        self.label = Label(self.frame,text="Input {}: ".format(num))
        self.label.pack(side=LEFT)

        self.entry = Entry(self.frame,width=30)
        self.entry.pack(side=LEFT)

        self.file = Button(self.frame,text="Import from file",command=self.getFile)
        self.file.pack(side=RIGHT)

        self.frame.pack()


    def getFile(self):
        filename = filedialog.askopenfilename(title = "Select Gene Sequence",filetypes = (("Text Files","*.txt"),("All Files","*.*")))
        with open(filename,"r") as f:
            self.entry.delete(0,END)
            string = f.read()
            self.entry.insert(0,string[:-1])

class IntSelect:
    def __init__(self,master,txt,intvar,low,high):
        self.master = master
        self.frame = Frame(master)

        self.label = Label(self.frame,text=txt)
        self.label.pack(side=LEFT)

        nums = []
        for i in range(low,high):
            nums.append(i)

        self.intMenu = OptionMenu(self.frame,intvar, *nums)
        self.intMenu.pack(side=RIGHT)
    
    def change_dropdown(*args):
        print(intvar.get())

class AffineGetBox:
    def __init__(self,master):
        self.master = master
        self.frame = Frame(master,bd=5,bg="grey")
        self.label = Label(self.frame,text="Algorithm Specifications",bg="grey")
        self.label.pack(anchor=N)

        self.msInt = IntVar(self.frame)
        self.msInt.set(3)
        self.matchselect = IntSelect(self.frame,"Match Score:",self.msInt,0,10)
        self.matchselect.frame.pack(anchor=E,fill=X)

        self.nmsInt = IntVar(self.frame)
        self.nmsInt.set(1)
        self.nonmatchselect = IntSelect(self.frame,"Non-Match Penalty:",self.nmsInt,0,10)
        self.nonmatchselect.frame.pack(anchor=E,fill=X)

        self.gosInt = IntVar(self.frame)
        self.gosInt.set(1)
        self.gapselect = IntSelect(self.frame,"Gap Opening Penalty:",self.gosInt,0,10)
        self.gapselect.frame.pack(anchor=W,fill=X)

        self.gesInt = IntVar(self.frame)
        self.gesInt.set(1)
        self.gapselect = IntSelect(self.frame,"Gap Extend Penalty:",self.gesInt,0,10)
        self.gapselect.frame.pack(anchor=W,fill=X)

    def getVals(self):
        output = {
            "match" : self.msInt.get(),
            "mismatch" : self.nmsInt.get(),
            "gap_open" : self.gosInt.get(),
            "gap_extend" : self.gesInt.get(),
        }
        return output

root = Tk()
my_gui = GUI(root)
root.mainloop()

