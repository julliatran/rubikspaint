# Name : Jullia Tran 
# Course: 15-112 Carnegie Mellon University
# Date: 11/27/18 
# Program : RubiksPaint

from math import pi, sin, cos
import sys, os, copy

import tkFileDialog as filedialog
import numpy as np

import direct.directbase.DirectStart
from direct.gui.OnscreenImage   import OnscreenImage
from direct.gui.OnscreenText    import OnscreenText
from direct.gui.DirectGui       import *
from direct.task                import Task 
from direct.showbase.ShowBase import ShowBase
from direct.showbase.DirectObject import DirectObject
from direct.interval.IntervalGlobal import LerpHprInterval, Func, Sequence
from panda3d.core import *
from panda3d.core import TransparencyAttrib

from imageToRubiks import imageToRubiks


base.setBackgroundColor(0.12157,0.47059,0.47451,1)
base.disableMouse()
base.camera.setPos(0, -10, 0)

# You can't normalize inline so this is a helper function
# Took this function from panda3d sample code: procedural-cube 
# by Kwasi Mensah 
def normalized(*args):
    myVec = LVector3(*args)
    myVec.normalize()
    return myVec

# Helper function to make a square given the Lower-Left-Hand and
# Upper-Right-Hand corners
# Use similar approach from panda3d sample code: procedural-cube,
# by Kwasi Mensah to make a square but I added ox,oy,oz
def makeSquare(x1, y1, z1, x2, y2, z2, ox = 0, oy = 0, oz = 0, col=(1,0,0)):
    format = GeomVertexFormat.getV3n3cpt2()
    vdata = GeomVertexData('square', format, Geom.UHDynamic)

    vertex = GeomVertexWriter(vdata, 'vertex')
    normal = GeomVertexWriter(vdata, 'normal')
    color = GeomVertexWriter(vdata, 'color')
    texcoord = GeomVertexWriter(vdata, 'texcoord')

    # reposition vertices based on offset
    x1 += ox
    y1 += oy
    z1 += oz
    x2 += ox
    y2 += oy
    z2 += oz

    # make sure we draw the sqaure in the right plane
    if x1 != x2:
        vertex.addData3(x1, y1, z1)
        vertex.addData3(x2, y1, z1)
        vertex.addData3(x2, y2, z2)
        vertex.addData3(x1, y2, z2)

        normal.addData3(normalized(2 * x1 - 1, 2 * y1 - 1, 2 * z1 - 1))
        normal.addData3(normalized(2 * x2 - 1, 2 * y1 - 1, 2 * z1 - 1))
        normal.addData3(normalized(2 * x2 - 1, 2 * y2 - 1, 2 * z2 - 1))
        normal.addData3(normalized(2 * x1 - 1, 2 * y2 - 1, 2 * z2 - 1))

    else:
        vertex.addData3(x1, y1, z1)
        vertex.addData3(x2, y2, z1)
        vertex.addData3(x2, y2, z2)
        vertex.addData3(x1, y1, z2)

        normal.addData3(normalized(2 * x1 - 1, 2 * y1 - 1, 2 * z1 - 1))
        normal.addData3(normalized(2 * x2 - 1, 2 * y2 - 1, 2 * z1 - 1))
        normal.addData3(normalized(2 * x2 - 1, 2 * y2 - 1, 2 * z2 - 1))
        normal.addData3(normalized(2 * x1 - 1, 2 * y1 - 1, 2 * z2 - 1))

    color.addData4f(col[0], col[1], col[2], 1.0)
    color.addData4f(col[0], col[1], col[2], 1.0)
    color.addData4f(col[0], col[1], col[2], 1.0)
    color.addData4f(col[0], col[1], col[2], 1.0)

    texcoord.addData2f(0.0, 1.0)
    texcoord.addData2f(0.0, 0.0)
    texcoord.addData2f(1.0, 0.0)
    texcoord.addData2f(1.0, 1.0)

    tris = GeomTriangles(Geom.UHDynamic)
    tris.addVertices(0, 1, 3)
    tris.addVertices(1, 2, 3)

    square = Geom(vdata)
    square.addPrimitive(tris)
    return square

def offsetToFaces(ox, oy, oz):
    result = []
    if ox == -1: result.append("left")
    elif ox == 1: result.append("right")
    if oy == -1: result.append("front")
    elif oy == 1: result.append("back")
    if oz == -1: result.append("bottom")
    elif oz == 1: result.append("top")
    return result

#return the coordinate and orientation for 27 cubes
def makeCubes(colors, faceToColor, faces, cubeToFaces):
    for ox in [-1, 0, 1]:
        for oy in [-1, 0, 1]:
            for oz in [-1, 0, 1]:
                square0 = makeSquare(-0.5, -0.5, -0.5,  0.5, -0.5,  0.5,
                                     ox, oy, oz, colors[faceToColor["front"]])
                square1 = makeSquare(-0.5,  0.5, -0.5,  0.5,  0.5,  0.5,
                                     ox, oy, oz, colors[faceToColor["back"]])
                square2 = makeSquare(-0.5,  0.5,  0.5,  0.5, -0.5,  0.5,
                                     ox, oy, oz, colors[faceToColor["top"]])
                square3 = makeSquare(-0.5,  0.5, -0.5,  0.5, -0.5, -0.5,
                                     ox, oy, oz, colors[faceToColor["bottom"]])
                square4 = makeSquare(-0.5, -0.5, -0.5, -0.5,  0.5,  0.5,
                                     ox, oy, oz, colors[faceToColor["left"]])
                square5 = makeSquare(0.5,  -0.5, -0.5,  0.5,  0.5,  0.5, 
                                     ox, oy, oz, colors[faceToColor["right"]])
                snode = GeomNode('square%d%d%d' % (ox%3,oy%3,oz%3))
                snode.addGeom(square0)
                snode.addGeom(square1)
                snode.addGeom(square2)
                snode.addGeom(square3)
                snode.addGeom(square4)
                snode.addGeom(square5)

                cube = render.attachNewNode(snode)
                cube.setTwoSided(True)

                cubeFaces = offsetToFaces(ox, oy, oz)
                cubeToFaces[cube] = set(cubeFaces)
                for face in cubeFaces: faces[face].append(cube)

class FaceSolver(DirectObject):

    def __init__(self):
        self.setupCameraControls()

        ####################
        ## OBJECT GLOBALS ##
        ####################
        self.theta       = 90
        self.phi         = 0
        self.rotateTime  = 0.7
        self.finalFace   = ["white","white","white","white","white",
                            "white","white","white","white",]
        self.cubeToFaces = None

        #####################
        ## CLASS CONSTANTS ##
        #####################

        # RGB values for Rubik's cube colors
        #https://www.schemecolor.com/rubik-cube-colors.php
        self.colors = {
            "blue"   : (0,0.31765,0.72941),
            "white"  : (1,1,1),
            "yellow" : (1,0.83529,0),
            "green"  : (0,0.61961,0.37647),
            "orange" : (1,0.34510,0),
            "red"    : (0.76863,0.11765,0.22745)
        }
        #for the top face
        self.colorTofaceToColor = {
            "white":  {"front": "blue", "back": "green", "top": "white",
                       "bottom": "yellow", "left": "red", "right": "orange"},
            "blue":   {"front": "yellow", "back": "white", "top": "blue",
                       "bottom": "green", "left": "red", "right": "orange"}, 
            "red":    {"front": "blue", "back": "green", "top": "red",
                       "bottom": "orange", "left": "yellow", "right": "white"},
            "green":  {"front": "orange", "back": "red", "top": "green",
                       "bottom": "blue", "left": "white", "right": "yellow"},
            "yellow": {"front": "blue", "back": "green", "top": "yellow",
                       "bottom": "white", "left": "orange", "right": "red"},
            "orange": {"front": "yellow", "back": "white", "top": "orange",
                       "bottom": "red", "left": "blue", "right": "green"},
        }

        # Mapping from faces to direction (0 is x, 1 is y, 2 is z)
        self.startingFaceToDir = {"left": 1, "right": 1, "top": 2,
                                  "bottom": 2, "front": 0, "back": 0}

        # h,p,r rotate around the z,x,y axes respectively      
        self.faceHpr = {
            "left"   : LVecBase3f(0.,90.,0.),
            "right"  : LVecBase3f(0.,-90.,0.),
            "top"    : LVecBase3f(-90.,0.,0.),
            "bottom" : LVecBase3f(90.,0.,0.),
            "front"  : LVecBase3f(0.,0.,90.),
            "back"   : LVecBase3f(0.,0.,-90.),
        }

        # Describes the faces that are going to be affected if a rotation in 
        # this orientation is called, in clock-wise order
        self.faceOrder = {
            "left"   : ["top","front","bottom","back"],
            "right"  : ["back","bottom","front","top"],
            "top"    : ["back","right","front","left"],
            "bottom" : ["left","front","right","back"],
            "front"  : ["top","right","bottom","left"],
            "back"   : ["left","bottom","right","top"]
        }

    ################
    ## USER INPUT ##
    ################

    def setupCameraControls(self):
        self.accept("arrow_left", self.cameraRotate,[-6,0])
        self.accept("arrow_right", self.cameraRotate,[6,0])
        self.accept("arrow_down", self.cameraRotate,[0,-6])
        self.accept("arrow_up", self.cameraRotate,[0,6])

    #for open cube
    def _setupKeyControls(self):
        # Set up rotation keys
        self.accept("l", self.addRotation, ["left"])
        self.accept("r", self.addRotation, ["right"])
        self.accept("t", self.addRotation, ["top"])
        self.accept("b", self.addRotation, ["bottom"])
        self.accept("f", self.addRotation, ["front"])
        self.accept("k", self.addRotation, ["back"])
        
        self.accept("shift-l", self.addRotation, ["left",-1])
        self.accept("shift-r", self.addRotation, ["right",-1])
        self.accept("shift-t", self.addRotation, ["top",-1])
        self.accept("shift-b", self.addRotation, ["bottom",-1])
        self.accept("shift-f", self.addRotation, ["front",-1])
        self.accept("shift-k", self.addRotation, ["back",-1])
        self.accept("enter", self._startSequence)


    #for open cube
    def _disableKeyControls(self):
        #Disable rotation keys
        self.ignore("l")
        self.ignore("r")
        self.ignore("t")
        self.ignore("b")
        self.ignore("f")
        self.ignore("k")
        self.ignore("enter")

        self.ignore("shift-l")
        self.ignore("shift-r")
        self.ignore("shift-t")
        self.ignore("shift-b")
        self.ignore("shift-f")
        self.ignore("shift-k")

    def setupKeyControls(self):
        self.accept("space", self.startSequence)

    def disableKeyControls(self):
        self.ignore("space")

    ####################
    ## CAMERA CONTROL ##
    ####################


    def cameraRotate(self,dphi,dtheta):
        self.phi += dphi
        self.theta += dtheta
        if self.theta < 0: self.theta = 0
        elif self.theta > 180: self.theta = 180

        phiR = self.phi * (pi / 180.0)
        thetaR = self.theta * (pi / 180.0)

        camera.setPos(10*sin(thetaR)*sin(phiR),
                     -10*sin(thetaR)*cos(phiR),
                     -10*cos(thetaR))
        
        camera.setHpr(self.phi, -self.theta + 90, 0)

    ############################
    ## SETUP AND REMOVE CUBES ##
    ############################

    # mapping each cube to the color of the sides
    def setupCubesToColor(self):
        for cube in self.cubeToFaces:
            self.cubeToColor[cube] = ["","",""]
            cubeFaces = list(self.cubeToFaces[cube])
            for face in cubeFaces:
                self.cubeToColor[cube][self.startingFaceToDir[face]] = \
                  self.faceToColor[face]

    def setupCubes(self):
        # mapping from faces to cubes within the face
        self.faces = {
            "front"  : [],
            "left"   : [],
            "right"  : [],
            "top"    : [],
            "bottom" : [],
            "back"   : []
        }
        self.cubeToFaces = {}
        self.cubeToColor = {}
        self.pivot = render.attachNewNode("pivot")
        
        # Generate the 27 cubes
        makeCubes(self.colors, self.faceToColor, self.faces, self.cubeToFaces)
        self.setupCubesToColor()

    def removeCubes(self):
        if self.cubeToFaces != None:
            for cube in self.cubeToFaces:
                cube.removeNode()

    ##############################################
    ## INITIATE OR MODIFY SEQUENCE OF ROTATIONS ##
    ##############################################

    # Adds a face rotation to the sequence.
    # This adds both the rotation function, which updates
    # data structures to track information regarding the cube,
    # as well as the rotation animation function LerpHprInterval
    # to the sequence.
    def addRotation(self, rotatingFace, dir = 1.0):
        self.sequence.append(Func(self.rotate, rotatingFace, dir))
        self.sequence.append(LerpHprInterval(self.pivot, self.rotateTime,
                             self.faceHpr[rotatingFace] * dir))

    # Begins executing a sequence of rotations.
    # Same approach as how Epihaius from this page uses sequence
    # from https://discourse.panda3d.org/t/rubiks-cube-in-panda/15586/3
    def startSequence(self):
        self.disableKeyControls()
        self.removeCubes()
        self._setSolution()
        self.sequence.append(Func(self.setupKeyControls))
        self.sequence.start()
        self.sequence = Sequence()

    # start sequence for open cube mode
    def _startSequence(self):
        self._disableKeyControls()
        self.sequence.append(Func(self._setupKeyControls))
        self.sequence.start()
        self.sequence = Sequence()
    # Compute and execute solution for a given final face color configuration
    def _setSolution(self):
        # Mapping from faces to colors
        self.faceToColor = self.colorTofaceToColor[self.finalFace[4]]

        # Setup virtual cubes for computing solution
        self.setupCubes()
        
        # Compute solution by simulating moves on virtual cubes
        self.solution = []
        self.solveCross()
        self.solveCorners()
        self.cleanSolution()

        # Remove virtual cubes and setup actual cubes to render
        self.removeCubes()
        self.setupCubes()

        # Initialize sequence
        self.sequence = Sequence()

        # Add solution to sequence
        for (rotatingFace, dir) in self.solution:
            self.addRotation(rotatingFace, dir)
        self.solution = []

    def setSolution(self, finalFace):
        self.finalFace = finalFace
        self.removeCubes()
        self._setSolution()
        self.setupKeyControls()

    ########################################################
    ## PERFORM A FACE ROTATION AND UPDATE DATA STRUCTURES ##
    ########################################################

    def rotate(self, rotatingFace, dir = 1):
        # update the set of faces each cube on the rotating face is in
        faceOrder = self.faceOrder[rotatingFace]
        if dir == -1:
            faceOrder = faceOrder[::-1]
        oldCubeFaces = {}
        for cube in self.faces[rotatingFace]:
            oldCubeFaces[cube] = list(self.cubeToFaces[cube])
            
            newCubeFaces = []
            for face in oldCubeFaces[cube]:
                if face in faceOrder:
                    idx = faceOrder.index(face)
                    newFace = faceOrder[(idx+1)% len(faceOrder)]
                    newCubeFaces.append(newFace)
                else: newCubeFaces.append(face)

            self.cubeToFaces[cube] = set(newCubeFaces)

        # update what cubes are in each face
        for cube in self.faces[rotatingFace]:
            for face in oldCubeFaces[cube]:
                if face != rotatingFace:
                    self.faces[face].remove(cube)
            for face in self.cubeToFaces[cube]:
                if face != rotatingFace:
                    self.faces[face].append(cube) 

        # Update color orientation of the cube
        for cube in self.faces[rotatingFace]:
            if rotatingFace == "left" or rotatingFace == "right":
                self.cubeToColor[cube][0], self.cubeToColor[cube][2] = \
                      self.cubeToColor[cube][2], self.cubeToColor[cube][0]
            elif rotatingFace == "top" or rotatingFace == "bottom":
                self.cubeToColor[cube][0], self.cubeToColor[cube][1] = \
                      self.cubeToColor[cube][1], self.cubeToColor[cube][0]
            elif rotatingFace == "front" or rotatingFace == "back":
                self.cubeToColor[cube][1], self.cubeToColor[cube][2] = \
                      self.cubeToColor[cube][2], self.cubeToColor[cube][1]
                
        # Updates the pivot to have rotatingFace cubes as children.
        # Used the same approach as how Epihaius from this page uses pivots
        # from https://discourse.panda3d.org/t/rubiks-cube-in-panda/15586/3
        # he/she uses 6 pivots, I use 1
        children = self.pivot.getChildren()
        children.wrtReparentTo(render)
        self.pivot.clearTransform()
        for cube in self.faces[rotatingFace]:
            cube.wrtReparentTo(self.pivot)

    #######################################
    ## COMPUTE SOLUTION HELPER FUNCTIONS ##
    #######################################

    ############################
    ## COMPUTE CROSS SOLUTION ##
    ############################
    
    #I used the beginner's method to do the cube, first the cross, then the
    #corners or rubiks cube. 

    #Check if the edge is already the correct color
    def checkEdge(self, face, targetColor):
        for cube in self.faces["top"]:
            if face in self.cubeToFaces[cube] and \
                    len(self.cubeToFaces[cube]) == 2 and \
                    self.cubeToColor[cube][2] == targetColor:
               return True
        return False

    # Find an edge with the given color
    def findEdge(self, color):
        for cube in self.cubeToColor:
            if ("" in self.cubeToColor[cube]) and \
               (color in self.cubeToColor[cube]) and \
               (self.cubeToColor[cube].count("") != 2) and \
               (not "top" in self.cubeToFaces[cube]): 
                return cube

    # Find the edge on the top face with the given face
    def findTargetEdge(self, face):
        for cube in self.faces["top"]:
            if face in self.cubeToFaces[cube] and \
              len(self.cubeToFaces[cube]) == 2:
                return cube;

    def solveCross(self):
        targetEdgeColors = [(self.finalFace[1], "back"),
                            (self.finalFace[3], "left"),
                            (self.finalFace[5], "right"),
                            (self.finalFace[7], "front")]

        for targetColor, face in targetEdgeColors:
            # Check if the piece at that point is already the correct color,
            # if so, moves on
            if self.checkEdge(face, targetColor):
                continue
            
            originalEdge = self.findEdge(targetColor)

            if originalEdge == None:
                self.solution.append((face,1))
                self.rotate(face,1)
                originalEdge = self.findEdge(targetColor)

            originalFaces = self.cubeToFaces[originalEdge]
            
        
            # Edge is on bottom face
            if "bottom" in originalFaces:
                originalEdgeColors = self.cubeToColor[originalEdge]

                while not face in self.cubeToFaces[originalEdge]:
                    self.solution.append(("bottom",1))
                    self.rotate("bottom",1)
                
                # Case where color is on bottom
                if targetColor == originalEdgeColors[2]:
                    self.solution.append((face,1))
                    self.rotate(face,1)
                    self.solution.append((face,1))
                    self.rotate(face,1)

                # Case where color is on side of bottom
                else:
                    self.solution.append((face,1))
                    self.rotate(face,1)

                    self.solution.append(("top",1))
                    self.rotate("top",1)

                    faceOrder = self.faceOrder["top"]
                    idx = faceOrder.index(face)
                    nextFace = faceOrder[(idx+1)% len(faceOrder)]

                    self.solution.append((nextFace,-1))
                    self.rotate(nextFace,-1)

                    self.solution.append(("top",-1))
                    self.rotate("top",-1)

            # Edge is not on bottom face
            else:
                targetEdge = self.findTargetEdge(face)

                # color is on the left side of face
                if ("left" in originalFaces and "back" in originalFaces and \
                      self.cubeToColor[originalEdge][1] == targetColor) or \
                   ("back" in originalFaces and "right" in originalFaces and \
                      self.cubeToColor[originalEdge][0] == targetColor) or \
                   ("right" in originalFaces and "front" in originalFaces and \
                      self.cubeToColor[originalEdge][1] == targetColor) or \
                   ("front" in originalFaces and "left" in originalFaces and \
                      self.cubeToColor[originalEdge][0] == targetColor):
                   
                    originalFace = "left"
                    if "back" in originalFaces and "right" in originalFaces:
                        originalFace = "back"
                    elif "right" in originalFaces and "front" in originalFaces:
                        originalFace = "right"
                    elif "front" in originalFaces and "left" in originalFaces:
                        originalFace = "front"

                    counter = 0
                    while not originalFace in self.cubeToFaces[targetEdge]:
                        self.solution.append(("top",1))
                        self.rotate("top",1)
                        counter += 1

                    self.solution.append(("top",1))
                    self.rotate("top",1)

                    faceOrder = self.faceOrder["top"]
                    idx = faceOrder.index(originalFace)
                    nextFace = faceOrder[(idx+1)% len(faceOrder)]

                    self.solution.append((nextFace,-1))
                    self.rotate(nextFace,-1)

                    self.solution.append(("top",-1))
                    self.rotate("top",-1)

                    for i in range(counter):
                        self.solution.append(("top",-1))
                        self.rotate("top",-1)

                # Color is on the right side of face
                else:

                    originalFace = "left"
                    if "front" in originalFaces and "right" in originalFaces:
                        originalFace = "front"
                    elif "right" in originalFaces and "back" in originalFaces:
                        originalFace = "right"
                    elif "back" in originalFaces and "left" in originalFaces:
                        originalFace = "back"

                    counter = 0
                    while not originalFace in self.cubeToFaces[targetEdge]:
                        self.solution.append(("top",1))
                        self.rotate("top",1)
                        counter += 1

                    self.solution.append(("top",-1))
                    self.rotate("top",-1)

                    faceOrder = self.faceOrder["bottom"]
                    idx = faceOrder.index(originalFace)
                    nextFace = faceOrder[(idx+1)% len(faceOrder)]

                    self.solution.append((nextFace,1))
                    self.rotate(nextFace,1)

                    self.solution.append(("top",1))
                    self.rotate("top",1)

                    for i in range(counter):
                        self.solution.append(("top",-1))
                        self.rotate("top",-1)  

    ##############################
    ## COMPUTE CORNERS SOLUTION ##
    ##############################

    # Check if piece is already completed
    def checkCorners(self, face, targetColor):
        for cube in self.faces["top"]:
            if face[0] in self.cubeToFaces[cube] and\
               face[1] in self.cubeToFaces[cube] and \
               self.cubeToColor[cube][2] == targetColor:
               return True
        return False

    def findCorners(self,color):
        for cube in self.cubeToColor:
            if (color in self.cubeToColor[cube]) and \
               (len(self.cubeToFaces[cube]) == 3) and \
               (not "top" in self.cubeToFaces[cube]): 
                return cube


    def solveCorners(self):
        targetCornerColors = [ (self.finalFace[0], ("back","left")),
                               (self.finalFace[2], ("back","right")),
                               (self.finalFace[6], ("front","left")),
                               (self.finalFace[8], ("front","right"))] 
        

        for targetColor, face in targetCornerColors:
            #check if the corner is already done
            if self.checkCorners(face,targetColor):
                continue
            
            originalCorner = self.findCorners(targetColor)
            
            # Brings the piece up top down to bottom.
            # Used the algorithms Ri,Di,R from this website: 
            # https://www.rubiks.com/blog/how-to-solve-the-rubiks-cube-stage-3
            # but I case it differently and came up with these rotations by
            # looking at a rubiks cube.
            if originalCorner == None:
                
                if "left" in face:
                    if "back" in face:
                        self.solution.append((face[1],-1))
                        self.rotate(face[1],-1)
                        self.solution.append(("bottom",-1))
                        self.rotate("bottom",-1)
                        self.solution.append((face[1],1))
                        self.rotate(face[1],1)
                    elif "front" in face:
                        self.solution.append((face[1],1))
                        self.rotate(face[1],1)
                        self.solution.append(("bottom",-1))
                        self.rotate("bottom",-1)
                        self.solution.append((face[1],-1))
                        self.rotate(face[1],-1)
                    originalCorner = self.findCorners(targetColor)
                elif "right" in face:
                    if "back" in face:
                        self.solution.append((face[1],1))
                        self.rotate(face[1],1)
                        self.solution.append(("bottom",1))
                        self.rotate("bottom",1)
                        self.solution.append((face[1],-1))
                        self.rotate(face[1],-1)
                    elif "front" in face:
                        self.solution.append((face[1],-1))
                        self.rotate(face[1],-1)
                        self.solution.append(("bottom",1))
                        self.rotate("bottom",1)
                        self.solution.append((face[1],1))
                        self.rotate(face[1],1)
                    originalCorner = self.findCorners(targetColor)
            originalFaces = self.cubeToFaces[originalCorner]
            
            #if the corner is at the bottom
            if "bottom" in originalFaces:
                originalCornerColors = self.cubeToColor[originalCorner]
                
                while (not face[0] in self.cubeToFaces[originalCorner] or not\
                                 face[1] in self.cubeToFaces[originalCorner]):
                    self.solution.append(("bottom",-1))
                    self.rotate("bottom",-1)
                        
                if "back" in face:
                    # Case when the color is facing along the y-axis
                    if targetColor == originalCornerColors[1]:
                        if "left" in face:
                            self.solution.append((face[1],-1))
                            self.rotate(face[1],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[1],1))
                            self.rotate(face[1],1)
                        elif "right" in face:
                            self.solution.append((face[1],1))
                            self.rotate(face[1],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[1],-1))
                            self.rotate(face[1],-1)
                    # Case when the color is facing along the x-axis
                    elif targetColor == originalCornerColors[0]:
                        if "left" in face:
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                        elif "right" in face:
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                    # Case when the color is facing along the z-axis
                    elif targetColor == originalCornerColors[2]:
                        if "left" in face:
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                        elif "right" in face:
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)

                #front corners
                elif "front" in face:
                    # Case when the color is facing along the y-axis
                    if targetColor == originalCornerColors[1]:
                        if "left" in face:
                            self.solution.append((face[1],1))
                            self.rotate(face[1],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[1],-1))
                            self.rotate(face[1],-1)
                        elif "right" in face:
                            self.solution.append((face[1],-1))
                            self.rotate(face[1],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[1],1))
                            self.rotate(face[1],1)
                    # Case when the color is facing along the y-axis
                    elif targetColor == originalCornerColors[0]:
                        if "left" in face:
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                        elif "right" in face:
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                    # Case when the color is facing along the z-axis
                    elif targetColor == originalCornerColors[2]:
                        if "left" in face:
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                        elif "right" in face:
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)
                            self.solution.append(("bottom",-1))
                            self.rotate("bottom",-1)
                            self.solution.append((face[0],1))
                            self.rotate(face[0],1)
                            self.solution.append(("bottom",1))
                            self.rotate("bottom",1)
                            self.solution.append((face[0],-1))
                            self.rotate(face[0],-1)

    ###################################################
    ## CLEAN UP FINAL SOLUTION (REMOVE EXCESS MOVES) ##
    ###################################################

    def removeBackAndForth(self,solution):
        result = []

        # remove back and forth rotations
        i = 0
        while i < len(solution):
            if i < len(solution)-1:
                face1, dir1 = solution[i]
                face2, dir2 = solution[i+1]
                if face1 == face2 and dir1 != dir2:
                    i += 2
                else:
                    result.append(solution[i])
                    i += 1
            else:
                result.append(solution[i])
                i += 1

        return result

    def reduceTriples(self,solution):
        result = []

        # reduce triple rotation
        i = 0
        while i < len(solution):
            if i < len(solution)-2:
                face1, dir1 = solution[i]
                face2, dir2 = solution[i+1]
                face3, dir3 = solution[i+2]
                if face1 == face2 and dir1 == dir2 and face2 == face3 \
                  and dir2 == dir3:
                    result.append((face1, -dir1))
                    i += 3
                else:
                    result.append(solution[i])
                    i += 1
            else:
                result.append(solution[i])
                i += 1
        return result

    # Clean solution (remove excess moves)
    def cleanSolution(self):
        solution = self.solution
        solution = self.reduceTriples(solution)
        solution = self.removeBackAndForth(solution)
        solution = self.reduceTriples(solution)
        solution = self.removeBackAndForth(solution)
        solution = self.removeBackAndForth(solution)
        self.solution = solution

class UserInterface(DirectObject):

    def __init__(self, faceSolver):
        self.rgbToColor = {
            (196,30,58)  : "red",
            (0,158,96)   : "green",
            (0,81,186)   : "blue",
            (255,88,0)   : "orange",
            (255,213,0)  : "yellow",
            (255,255,255): "white"
        }

        ## OBJECT GLOBAL VARIABLES ##
        self.faceSolver = faceSolver
        self.maxNumCubes = 100

        self.previewMinX = -0.45
        self.previewMaxX = 0.45
        self.previewMinY = -0.3
        self.previewMaxY = 0.6

        self.backdropMinX = -0.9
        self.backdropMaxX = 0.9
        self.backdropMinY = -0.9
        self.backdropMaxY = 0.9

        self.squares = []
        self.squareBBoxs = []
        self.highlightSquares = []
        self.squareStep = 0.0

        self.mouseClickX = 0.0
        self.mouseClickY = 0.0

        self.filepath = ""
        self.processedFilepath = ""
        self.colorMatrix = None
        self.myImage = None

        ##### Main Frame #####

        #mainFrame globals
        self.mainFrame = None
        self.titleName = None
        self.importButton = None
        self.retryDialog = None
        self.startButton = None
        self.cubeAmount = None
        

        #used panda3d manual for frames, images and texts

        self.mainFrame = DirectFrame(frameColor = (0.12157,0.47059,0.47451,1),
                                     frameSize = (-2,2,-2,2),
                                     pos=(0,0,0))
        
        self.logo = OnscreenImage("images/rubiks.png", scale = 0.70, pos = (0, 0, 0.11))
        self.logo.setTransparency(TransparencyAttrib.MAlpha)
        self.logo.reparentTo(self.mainFrame)
        self.titleRubiks = OnscreenImage("images/RubiksTitle.png", scale = 0.5, pos = (-0.44, 0, 0.7))
        self.titleRubiks.setTransparency(TransparencyAttrib.MAlpha)
        self.titleRubiks.reparentTo(self.mainFrame)
        self.titlePaint = OnscreenImage("images/PaintTitle.png", scale = 0.5, pos = (0.5, 0, 0.7))
        self.titlePaint.setTransparency(TransparencyAttrib.MAlpha)
        self.titlePaint.reparentTo(self.mainFrame)


        self.importButton = DirectButton(text = "Import Image",
                                         pos = (0,0,-0.5),
                                         text_scale = (0.1,0.1),
                                         text_fg = (0,0.31765,0.37647,1),
                                         scale = 0.8, command=self.importImage)                 
        self.importButton.reparentTo(self.mainFrame)

        self.startButton = DirectButton(text = "Start!", pos = (0,0,-0.8),\
                             text_scale = (0.1,0.1),text_fg =(0,0.31765,0.37647,1),\
                             scale = 0.8, command = self.toPictureFrame)
        self.startButton.reparentTo(self.mainFrame)
        self.startButton.hide()
        
        self.openCubeButton = DirectButton(text = "Open Cube", pos = (0,0,-0.8),\
                             text_scale = (0.1,0.1),text_fg =(0,0.31765,0.37647,1),\
                             scale = 0.8, command = self.toOpenCubeFrame)
        self.openCubeButton.reparentTo(self.mainFrame)

        self.cubeAmount = DirectEntry(text = "", scale=0.05,
                command=self.setMaxNumCubes, initialText = "Enter amount of cubes (1-2000):",
                text_fg =(0,0.31765,0.37647,1),numLines = 2,focus=0,
                 focusInCommand=self.clearText, pos=(-0.25,0,-0.8))
        self.cubeAmount.reparentTo(self.mainFrame)

        self.cubeSpeed = DirectOptionMenu(text= "Speed", scale=0.07,
            items=["Slow","Normal","Fast","Ludicrous"],initialitem=1,
            command=self.setCubeSpeed, textMayChange=1,
            pos=(-1,0,-0.9), text_fg = (0,0.31765,0.37647,1))
        self.cubeSpeed.reparentTo(self.mainFrame)

        self.speedText = OnscreenText(text = 'Speed:',
                                    pos=(-1.0,-0.8), scale=0.07,fg = (1.0,1.0,1.0,1),
                                    align=TextNode.ALeft)
        self.speedText.reparentTo(self.mainFrame)
        ##### Picture Frame #####

        # pictureFrame globals 
        self.pictureFrame = DirectFrame(frameColor = (0.12157,0.47059,0.47451,1),
                                     frameSize = (-2,2,-2,2),
                                     pos=(0,0,0))
        
        self.mainMenuButton = DirectButton(text = "Main Menu", pos = (-1.1,0,0.85),\
                             text_scale = (0.1,0.1),text_fg =(0,0.31765,0.37647,1),\
                             scale = 0.5, command = self.toMainMenuFrame)
        self.mainMenuButton.reparentTo(self.pictureFrame)
        
        

        ##### Cube Frame #####

        self.cubeFrame = DirectFrame()
        self.toPictureFrameButton = DirectButton(text = "Back to Picture",
                                         pos = (-1.0,0,0.85),
                                         text_scale = (0.1,0.1),
                                         text_fg = (0,0.31765,0.37647,1),
                                         scale = 0.5, command=self.toPictureFrame)
        self.toPictureFrameButton.reparentTo(self.cubeFrame)
        
        self.spaceText = OnscreenText(text = 'Press spacebar to solve\nUse arrow keys to rotate camera',
                                    pos=(-1.3,-0.9), scale=0.05,fg = (1.0,1.0,1.0,1),
                                    align=TextNode.ALeft)
        self.spaceText.reparentTo(self.cubeFrame)


        #### Open Cube Frame #####
        self.openCubeFrame = DirectFrame()
        self.toMainMenuButtonFromOpenCube = DirectButton(text = "Back to Menu",
                                         pos = (-1.1,0,0.85),
                                         text_scale = (0.1,0.1),
                                         text_fg = (0,0.31765,0.37647,1),
                                         scale = 0.5, command=self._toMainMenuFrame)
        self.toMainMenuButtonFromOpenCube.reparentTo(self.openCubeFrame)
        self.cubeAmount.hide()

        self.arrowsText = OnscreenText(text = 
        'Use keys l, r, t, b, f, k to add rotations\nUse shift + key to negate rotation\nHit enter to start rotation\nUse arrow keys to rotate camera',
                        pos=(-1.3,-0.8), scale=0.05,fg=(1.0,1.0,1.0,1),
                        align=TextNode.ALeft, wordwrap= 20.0)
        self.arrowsText.reparentTo(self.openCubeFrame)      

        self.pictureFrame.hide()
        self.cubeFrame.hide()
        self.openCubeFrame.hide()
        
        
        self.ignore("mouse1")
        self.faceSolver._disableKeyControls()
    
    def setCubeSpeed(self,arg):
        if arg == "Slow":
            self.faceSolver.rotateTime = 2.0
        elif arg == "Normal":
            self.faceSolver.rotateTime = 1.0
        elif arg == "Fast":
            self.faceSolver.rotateTime = 0.5
        elif arg == "Ludicrous":
            self.faceSolver.rotateTime = 0.07

    # Setting the face color of the rubiks according to the selected face.
    # Pass in the input on this face to compute the solution
    def selectFace(self):
        self.ignore("mouse1")

        w = len(self.colorMatrix[0])//3 if self.colorMatrix != None else 0

        x, y = None, None
        if base.mouseWatcherNode.hasMouse(): 
            x, y = base.mouseWatcherNode.getMouseX(), base.mouseWatcherNode.getMouseY()

        if x != None and y != None:
            for i in range(len(ui.squares)):
                bbox = ui.squareBBoxs[i]
                squareNode = ui.squares[i]
                ll = render2d.getRelativePoint(squareNode, Point3(bbox[0], 0, bbox[1]))
                ur = render2d.getRelativePoint(squareNode, Point3(bbox[2], 0, bbox[3]))
                l, b, r, t = ll.getX(), ll.getZ(), ur.getX(), ur.getZ()

                #within the box
                if x >= l and y >= b and x <= r and y <= t:
                    brow = 3 * (i // w)
                    bcol = 3 * (i % w)

                    #one rubiks face is 9 squares --> object row and ocject col for top cube
                    self.finalFace = []
                    for orow in range(3):
                        for ocol in range(3):
                            color = self.rgbToColor[tuple(self.colorMatrix[brow+orow,bcol+ocol])]
                            self.finalFace.append(color)

                    self.cubeFrame.show()
                    for highlightSquare in ui.highlightSquares:
                                highlightSquare.removeNode()
                    taskMgr.remove('CheckMousePosition')
                    self.pictureFrame.hide()

                    self.faceSolver.setSolution(self.finalFace)

    # Draws the image on screen
    def drawImage(self, minX, minY, maxX, maxY, colorMatrix, frame):
        for squareNode in self.squares:
            squareNode.removeNode()
        self.squares = []
        self.squareBBoxs = []

        h, w = tuple(colorMatrix.shape[0:2])

        previewW = maxX - minX
        previewH = maxY - minY

        step = 1.0
        startX = minX
        startY = minY
        if w > h:
            step = previewW / float(w)
            startY = 0.5 * previewH + minY - 0.5 * step * h
        else:
            step = previewH / float(h)
            startX = 0.5 * previewW + minX - 0.5 * step * w
        self.squareStep = step

        for brow in range(0,h,3):
            for bcol in range(0,w,3):

                snode = GeomNode('target_face_%d%d' % (brow, bcol))

                bbox = [startX + step * bcol,
                        startY + step * (h-brow-1),
                        startX + step * (bcol+1),
                        startY + step * (h-brow)]

                for orow in range(3):
                    for ocol in range(3):
                        row = brow + orow
                        col = bcol + ocol

                        x0 = startX + step * col
                        y0 = startY + step * (h-row-1)
                        x1 = startX + step * (col+1)
                        y1 = startY + step * (h-row)

                        if x0 < bbox[0]: bbox[0] = x0
                        if y0 < bbox[1]: bbox[1] = y0
                        if x1 > bbox[2]: bbox[2] = x1
                        if y1 > bbox[3]: bbox[3] = y1

                        color = tuple([i/255.0 for i in colorMatrix[row,col]])
        
                        square = makeSquare(x0, 0, y0, x1, 0, y1, col=color)                
                        snode.addGeom(square)
                        
                squareNode = aspect2d.attachNewNode(snode)
                squareNode.reparentTo(frame)
        
                self.squares.append(squareNode)
                self.squareBBoxs.append(bbox)

    # Command for button to show picture frame
    def toPictureFrame(self): 
        self.mainFrame.hide()
        self.cubeFrame.hide()
        self.openCubeFrame.hide()
        self.pictureFrame.hide()

        self.colorMatrix, self.processedFilepath = imageToRubiks(self.filepath,\
                                                             self.maxNumCubes)

        self.drawImage(self.backdropMinX, self.backdropMinY,
                       self.backdropMaxX, self.backdropMaxY,
                       self.colorMatrix, self.pictureFrame)
        taskMgr.add(checkMousePosition, 'CheckMousePosition')
        
       
        self.accept("mouse1",self.selectFace)
        
        self.pictureFrame.show()

    # Command for button to show menu frame
    def toMainMenuFrame(self):
        self.ignore("mouse1")

        for highlightSquare in ui.highlightSquares:
                    highlightSquare.removeNode()
        taskMgr.remove('CheckMousePosition')
        self.maxNumCubes = 100
        
        self.startButton.hide()

        self.pictureFrame.hide()
        self.cubeFrame.hide()
        self.openCubeFrame.hide()

        self.logo.show()
        self.openCubeButton.show()
        self.mainFrame.show()
    
    # Back to main menu from open cube frame
    def _toMainMenuFrame(self):
        self.ignore("mouse1")
        self.faceSolver.removeCubes()
        
        self.openCubeFrame.hide()
        self.cubeFrame.hide()
        self.pictureFrame.hide()

        self.faceSolver._disableKeyControls()

        self.mainFrame.show()

    # Command for button to show open cube frame
    def toOpenCubeFrame(self):
        self.faceSolver.disableKeyControls()
        self.faceSolver.removeCubes()

        self.ignore("mouse1")
        self.mainFrame.hide()
        self.cubeFrame.hide()
        self.pictureFrame.hide()

        self.faceSolver.faceToColor = self.faceSolver.colorTofaceToColor["white"]
        self.faceSolver.setupCubes()
        self.faceSolver._setupKeyControls()
        self.faceSolver.sequence = Sequence()
        
        self.openCubeFrame.show()
        
    
    #command for retryDialog
    def itemSel(self, arg):
        if arg == 1:
            self.importImage()
        else:
            self.retryDialog.hide()
            self.mainFrame.show()
    
    # Direct entry commands
    def clearText(self):
        self.cubeAmount.enterText('')

    def setMaxNumCubes(self, textEntered):
        if 0 < int(textEntered) <= 2000:
            self.maxNumCubes = int(textEntered)
            self.cubeAmount.hide()
            self.startButton.show()


    def importImage(self):
        self.startButton.hide()
        # show an "Open" dialog box and return the path to the selected file
        #https://pythonspot.com/tk-file-dialogs/
        self.filepath = filedialog.askopenfilename() 
        fileExt = os.path.splitext(self.filepath)[1]
        isJpg = fileExt == ".jpg" or fileExt == ".JPG"

        if self.retryDialog != None: self.retryDialog.hide()

        if isJpg:
            self.logo.hide()
            self.openCubeButton.hide()
            #set amount of cubes according to user input (between 1-1000)
            self.cubeAmount.show()

            self.colorMatrix, self.processedFilepath = imageToRubiks(self.filepath,\
                                                             self.maxNumCubes)
            self.drawImage(self.previewMinX, self.previewMinY,
                           self.previewMaxX, self.previewMaxY,
                           self.colorMatrix, self.mainFrame)

        elif self.filepath != "":
            self.retryDialog = RetryCancelDialog(dialogName="RetryNoCancelDialog",\
                     text="Please retry with a .jpg image", text_fg =(0,0.31765,0.37647,1),\
                     command = self.itemSel)
            self.retryDialog.reparentTo(self.mainFrame)
            self.retryDialog.show()
   
# Create frames
faceSolver = FaceSolver()
ui = UserInterface(faceSolver)

# Similar approach to https://discourse.panda3d.org/t/mouse-over-click-on-text/3320
def getMousePosition():
    if base.mouseWatcherNode.hasMouse(): 
        return base.mouseWatcherNode.getMouseX(), base.mouseWatcherNode.getMouseY()
    return None, None

# Similar approach to https://discourse.panda3d.org/t/mouse-over-click-on-text/3320
# Draws the black square box used to select the cube on the image.
def checkMousePosition(task): 
    # Get mouse position
    x, y = getMousePosition()

    inside = False
    if x != None and y != None:
        for i in range(len(ui.squares)):
            bbox = ui.squareBBoxs[i]
            squareNode = ui.squares[i]
            ll = render2d.getRelativePoint(squareNode, Point3(bbox[0], 0, bbox[1]))
            ur = render2d.getRelativePoint(squareNode, Point3(bbox[2], 0, bbox[3]))
            l, b, r, t = ll.getX(), ll.getZ(), ur.getX(), ur.getZ()

            if x >= l and y >= b and x <= r and y <= t:
                inside = True
                for highlightSquare in ui.highlightSquares:
                    highlightSquare.removeNode()

                step = ui.squareStep / 4

                squareLeft   = makeSquare(bbox[0]-step, 0, bbox[1]-step, bbox[0], 0, bbox[3]+step, col=(0,0,0))
                squareBottom = makeSquare(bbox[0]-step, 0, bbox[1]-step, bbox[2]+step, 0, bbox[1], col=(0,0,0))
                squareRight  = makeSquare(bbox[2], 0, bbox[1]-step, bbox[2]+step, 0, bbox[3]+step, col=(0,0,0))
                squareTop    = makeSquare(bbox[0]-step, 0, bbox[3], bbox[2]+step, 0, bbox[3]+step, col=(0,0,0))

                snodeLeft    = GeomNode('highlight_square_left')
                snodeBottom  = GeomNode('highlight_square_bottom')
                snodeRight   = GeomNode('highlight_square_right')
                snodeTop     = GeomNode('highlight_square_top')

                snodeLeft.addGeom(squareLeft)
                snodeBottom.addGeom(squareBottom)
                snodeRight.addGeom(squareRight)
                snodeTop.addGeom(squareTop)

                squareLeftNode  = aspect2d.attachNewNode(snodeLeft)
                squareBottomNode = aspect2d.attachNewNode(snodeBottom)
                squareRightNode  = aspect2d.attachNewNode(snodeRight)
                squareTopNode    = aspect2d.attachNewNode(snodeTop)

                ui.highlightSquares = [squareLeftNode, squareBottomNode,
                                       squareRightNode, squareTopNode]
    
    if not inside:
        for highlightSquare in ui.highlightSquares:
                    highlightSquare.removeNode()

    return Task.cont

run()
