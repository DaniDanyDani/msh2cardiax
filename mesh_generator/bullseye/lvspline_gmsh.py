#!python
import sys
import os
import time
from math import sin, cos, pi, tan

# Variaveis globais
malhaTet = True

# GMSH STUFF
gmsh_line_id = 0
gmsh_vertex_id = 0
gmsh_surface_id = 0
gmsh_volume_id = 0

geoFile = None

# ------------------------------------------------------------------------------

class Block:
    def __init__(self, id, a1_1, a1_2, a1_3, a1_4, a2_1, a2_2, a2_3, a2_4, t0, t1, p0, p1):
        self.num = id
        self.phi0 = p0
        self.phi1 = p1
        self.theta0 = t0
        self.theta1 = t1
        self.apex = False

        if (a1_3 == None and a1_4 == None and a2_3 == None and a2_4 == None):
            a1_3, a1_4, a2_3, a2_4 = -1, -1, -1, -1
            self.apex = True

        self.a_endo1, self.a_endo2 = a1_1, a1_2
        self.a_endo3, self.a_endo4 = a1_3, a1_4
        self.a_epi1, self.a_epi2 = a2_1, a2_2
        self.a_epi3, self.a_epi4 = a2_3, a2_4

    def __str__(self):
        a, b, c, d = self.a_endo1, self.a_endo2, self.a_endo3, self.a_endo4
        e, f, g, h = self.a_epi1, self.a_epi2, self.a_epi3, self.a_epi4
        s = 'Bloco %d ;\n' % self.num
        s += '  a_endo = (%f, %f, %f, %f);\n' % (a, b, c, d)
        s += '  a_epi  = (%f, %f, %f, %f);\n' % (e, f, g, h)
        s += '  phi    = (%f, %f);\n' % (self.phi0, self.phi1)
        s += '  theta  = (%f, %f);\n' % (self.theta0, self.theta1)
        s += '  apex   = %s' % self.apex
        return s

    def getNum(self):
        return int(self.num)

    def getEndoData(self):
        a11, a12, a13, a14 = self.a_endo1, self.a_endo2, self.a_endo3, self.a_endo4
        return a11, a12, a13, a14

    def getEpiData(self):
        a21, a22, a23, a24 = self.a_epi1, self.a_epi2, self.a_epi3, self.a_epi4
        return a21, a22, a23, a24

    def getTheta(self):
        return self.theta0, self.theta1

    def getPhi(self):
        return self.phi0, self.phi1

    def inside(self, teta, phi):
        if (phi > 0):
            if (teta >= self.theta0 and teta <= self.theta1 and
                        phi >= self.phi0 and phi <= self.phi1):
                return True
        elif (phi < 0):
            if (teta >= self.theta0 and teta <= self.theta1 and
                        phi >= self.phi0 and phi <= self.phi1):
                return True
        else:
            return False

    def isApex(self):
        return self.apex


# ------------------------------------------------------------------------------

def findBlock(blocklist, teta, phi):
    myblock, blockid = None, None
    for i in range(len(blocklist)):
        res = blocklist[i].inside(teta, phi)
        if (res):
            myblock = blocklist[i]
            blockid = i + 1
    if (blockid is None):
        print("ERROR BLOCK NOT FOUND")
        return None
    return myblock

# ------------------------------------------------------------------------------

def createVertex(x, y, z):
    global gmsh_vertex_id, geoFile
    gmsh_vertex_id += 1
    geoFile.write("Point(%d) = {%e, %e, %e, 1.0};\n" % (gmsh_vertex_id,x,y,z))
    return gmsh_vertex_id

# ------------------------------------------------------------------------------

def createLine(nodes):
    global gmsh_line_id, geoFile
    gmsh_line_id += 1
    geoFile.write("Line(%d) = {%d, %d};\n" % (gmsh_line_id,nodes[0],nodes[1]))
    return gmsh_line_id

# ------------------------------------------------------------------------------

def createSpline(vertexList):
    global gmsh_line_id, geoFile
    gmsh_line_id += 1
    cmdstr = 'Spline(%d) = {' % gmsh_line_id
    for i in range(len(vertexList)-1):
        cmdstr += ' %d,' % vertexList[i]
    cmdstr += ' %d' % vertexList[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)
    return gmsh_line_id

# ------------------------------------------------------------------------------

def createSurface(curvesList):
    cmdstr = 'create surface curve'
    for i in range(len(curvesList)):
        cmdstr += ' %d' % curvesList[i]
    cubit.cmd(cmdstr)
    surfID = cubit.get_last_id('surface')
    return surfID

# ------------------------------------------------------------------------------

def createCompoundSurface(surfID):

    global gmsh_surface_id, geoFile

    gmsh_surface_id +=1

    cmdstr = 'Compound Surface(%d) = {' % gmsh_surface_id
    for i in range(len(surfID)-1):
        cmdstr += ' %d,' % surfID[i]
    cmdstr += ' %d' % surfID[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)

    return gmsh_surface_id


def createPhysicalSurface(surfID, marker):

    global geoFile

    if isinstance(surfID,int):
        geoFile.write('Physical Surface(%d) = {%d};\n' % (marker, surfID) )

    else:
        cmdstr = 'Physical Surface(%d) = {' % marker
        for i in range(len(surfID)-1):
            cmdstr += ' %d,' % surfID[i]
        cmdstr += ' %d' % surfID[-1]
        cmdstr += ' };\n'
        geoFile.write("%s" % cmdstr)


def createPhysicalVolume(volID, marker):
    global geoFile
    geoFile.write('Physical Volume(%d) = {%d};\n' % (marker, volID) )


# ------------------------------------------------------------------------------


def createSurfaceBase(endoBase, epiBase):

    global gmsh_line_id, gmsh_surface_id, geoFile

    
    gmsh_surface_id += 1

    cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
    for i in range(len(endoBase)-1):
        cmdstr += ' %d,' % endoBase[i]

    cmdstr += ' %d' % endoBase[-1] 
    cmdstr += ' };\n'

    geoFile.write("%s" % cmdstr)

    cmdstr2 = "Plane Surface(%d) = {%d" % (gmsh_surface_id,gmsh_surface_id)
    bID = gmsh_surface_id
    
    gmsh_surface_id += 1
    cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
    for i in range(len(epiBase)-1):
        cmdstr += ' %d,' % epiBase[i]
    
    cmdstr += ' %d' % epiBase[-1]  
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)

    cmdstr2+=", %d};\n" % (gmsh_surface_id)
    geoFile.write("%s" % cmdstr2)

    gmsh_surface_id += 1
    geoFile.write('Surface Loop(%d) = {%d};\n' % (gmsh_surface_id, bID) )


    return gmsh_surface_id, bID


def createSurfaceSkin(splineList, splineVertex):
    # sp1 = splineList[0]
    # sp2 = splineList[-1]
    # cubit.cmd('create surface skin curve %d to %d' % (sp1, sp2))
    # surfID = cubit.get_last_id('surface')
    # return surfID
    global gmsh_line_id, gmsh_surface_id, geoFile

    baseSplines1 = []
    baseSplines2 = []
    baseSplines3 = []
    baseSplines4 = []
    surfaceLoop = []

    for i in range(len(splineList[0])-1):

        sid1 = createSpline([splineVertex[0][i], splineVertex[0][i+1]])
        baseSplines1.append(sid1)
        sid2 = createSpline([splineVertex[1][i], splineVertex[1][i+1]])
        baseSplines2.append(sid2)
        sid3 = createSpline([splineVertex[2][i], splineVertex[2][i+1]])
        baseSplines3.append(sid3)
        sid4 = createSpline([splineVertex[3][i], splineVertex[3][i+1]])
        baseSplines4.append(sid4)
        sid = [sid1, sid2, sid3, sid4]

        for j in range(3):
            gmsh_surface_id += 1
            cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
            cmdstr += ' %d,' % sid[j]
            cmdstr += ' %d,' % splineList[j][i+1]
            cmdstr += ' %d,' % -sid[j+1]
            cmdstr += ' %d' % -splineList[j][i]
            cmdstr += ' };\n'
            geoFile.write("%s" % cmdstr)
            geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
            surfaceLoop.append(gmsh_surface_id)


        gmsh_surface_id += 1
        cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
        cmdstr += ' %d,' % sid[3]
        cmdstr += ' %d,' % splineList[3][i+1]
        cmdstr += ' %d' % -splineList[3][i]
        cmdstr += ' };\n'
        geoFile.write("%s" % cmdstr)
        geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
        surfaceLoop.append(gmsh_surface_id)



    sid1 = createSpline([splineVertex[0][-1], splineVertex[0][0]])
    baseSplines1.append(sid1)
    sid2 = createSpline([splineVertex[1][-1], splineVertex[1][0]])
    baseSplines2.append(sid2)
    sid3 = createSpline([splineVertex[2][-1], splineVertex[2][0]])
    baseSplines3.append(sid3)
    sid4 = createSpline([splineVertex[3][-1], splineVertex[3][0]])
    baseSplines4.append(sid4)
    sid = [sid1, sid2, sid3, sid4]

    for j in range(3):
        gmsh_surface_id += 1
        cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
        cmdstr += ' %d,' % sid[j]
        cmdstr += ' %d,' % splineList[j][0]
        cmdstr += ' %d,' % -sid[j+1]
        cmdstr += ' %d' % -splineList[j][-1]
        cmdstr += ' };\n'
        geoFile.write("%s" % cmdstr)
        geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
        surfaceLoop.append(gmsh_surface_id)


    gmsh_surface_id += 1
    cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
    cmdstr += ' %d,' % sid[3]
    cmdstr += ' %d,' % splineList[3][0]
    cmdstr += ' %d' % -splineList[3][-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)
    geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
    surfaceLoop.append(gmsh_surface_id)

    surfID = gmsh_surface_id


    gmsh_surface_id +=1
    cmdstr = 'Surface Loop(%d) = {' % gmsh_surface_id
    for i in range(len(surfaceLoop)-1):
        cmdstr += ' %d,' % surfaceLoop[i]
    cmdstr += ' %d' % surfaceLoop[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)



    return gmsh_surface_id, surfaceLoop, baseSplines1


# ------------------------------------------------------------------------------

def createSurfaceSkin2(splineList):
    msg = 'create surface skin curve '
    for i in len(splineList):
        msg += '%d ' % splineList[i]
    cubit.cmd(msg)
    surfID = cubit.get_last_id('surface')
    return surfID

# ------------------------------------------------------------------------------

def createVolume(surfList):
    endoSID, epiSID, baseSID = surfList
    global gmsh_volume_id, gmsh_surface_id, geoFile
    gmsh_volume_id +=1

    cmdstr = 'Volume(%d) = {%d, %d, %d};\n' % (gmsh_volume_id, endoSID, baseSID, epiSID)

    geoFile.write("%s" % cmdstr)

    return gmsh_volume_id

# ------------------------------------------------------------------------------

def createVolume2(surfList):
    global geoFile
    cmdstr = 'create volume surface'
    for i in range(len(surfList)):
        cmdstr += ' %d' % surfList[i]
    cmdstr += " heal"
    geoFile.write(cmdstr)
    cubit.cmd(cmdstr)
    volID = cubit.get_last_id('volume')
    return volID


# ------------------------------------------------------------------------------

def volumeList():
    cubit.cmd('group "allVolumes" add vol all')
    groupID = cubit.get_id_from_name("allVolumes")
    volList = cubit.get_group_volumes(groupID)
    return volList


# ------------------------------------------------------------------------------

def uniteVolumes(vlist):
    print("Unindo volumes")
    cmdstr = 'unite volume'
    for v in vlist:
        cmdstr += ' %d' % int(v)
    cubit.cmd(cmdstr)


# ------------------------------------------------------------------------------

def ellipsoid(a, c, phi, theta):
    """
    Parametric equations for the ellipsoid
    """
    PHI0 = pi / 6

    #translate base to 0
    offset = c * sin(PHI0)

    x = a * cos(phi) * cos(theta)
    y = a * cos(phi) * sin(theta)
    z = c * sin(phi) + offset
    return (x, y, z)


# ------------------------------------------------------------------------------

def printBullseye(b):
    for i in range(len(b)):
        print("Segmento %d: %.3f cm" % (i + 1, b[i]))


# ------------------------------------------------------------------------------

def readBullseye(inputfile):
    f = open(inputfile)
    f.readline()
    f.readline()
    a1 = float(f.readline()) / 2

    f.readline()
    c1 = float(f.readline())

    # antigo
    # c1 = float(f.readline())/(1.565)  #C1 Diastole
    # c1 = float(f.readline())/(1.75)  #C1 Sistole

    bull=[]

    for i in range(19):
        f.readline()
        bull.append( float( f.readline() ) ) 

    printBullseye(bull)

    return a1, c1, bull


# ------------------------------------------------------------------------------

def lerp(r0, r1, theta0, theta1, theta):
    """
    Linear interpolation
    """
    r = r0 + ((r1 - r0) / (theta1 - theta0)) * (theta - theta0)
    return r


# ------------------------------------------------------------------------------

def createSurfaceListVertex(lst):
    msg = 'create surface vertex '
    for i in range(len(lst)):
        msg += '%d ' % lst[i]
    cubit.cmd(msg)
    surfID = cubit.get_last_id('surface')
    return surfID


# ------------------------------------------------------------------------------

def createTriMesh(lst):
    msg = 'mesh surface '
    for i in range(len(lst)):
        msg += '%d ' % lst[i]
    cubit.cmd(msg)
    MeshID = cubit.get_last_id('mesh')
    return MeshID


# -------------------------------------------------------------------------------

def createBlocks(a1, bull):
    nblocks = 0
    lblocks = []

    inc = 7 * pi / 30

    # AQUI
    phi0 = -pi / 6
    phi1 = -pi / 6 + inc
    phi2 = -pi / 6 + 2 * inc
    phi3 = -pi / 6 + 3 * inc

    print("Modificacoes Base phi2->phi3")
    print("Gerando blocos da regiao basal")
    for i in range(6):
        nblocks += 1
        print("Bloco %d" % (nblocks))
        theta0 = (0 + i) * pi / 3.0
        theta1 = (1 + i) * pi / 3.0
        if (i < 5):
            a2_0 = a1 + (bull[0 + i] / 2);
            a2_1 = a1 + (bull[1 + i] / 2)
            a2_2 = a1 + (bull[6 + i] / 2);
            a2_3 = a1 + (bull[7 + i] / 2)
            a1_0 = a1 - (bull[0 + i] / 2);
            a1_1 = a1 - (bull[1 + i] / 2)
            a1_2 = a1 - (bull[6 + i] / 2);
            a1_3 = a1 - (bull[7 + i] / 2)

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi0, phi1)
            lblocks.append(b)
        else:
            a2_0 = a1 + (bull[5 + (i - 5)] / 2);
            a2_1 = a1 + (bull[0 + (i - 5)] / 2)
            a2_2 = a1 + (bull[11 + (i - 5)] / 2);
            a2_3 = a1 + (bull[6 + (i - 5)] / 2)
            a1_0 = a1 - (bull[5 + (i - 5)] / 2);
            a1_1 = a1 - (bull[0 + (i - 5)] / 2)
            a1_2 = a1 - (bull[11 + (i - 5)] / 2);
            a1_3 = a1 - (bull[6 + (i - 5)] / 2)

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi0, phi1)
            lblocks.append(b)

    print("Gerando blocos da regiao medial")
    for i in range(6):
        nblocks += 1
        print("Bloco %d" % (nblocks))
        theta0 = (0 + i) * pi / 3
        theta1 = (1 + i) * pi / 3
        if (i < 5):
            a2_0 = a1 + (bull[i + 6] / 2);
            a2_1 = a1 + (bull[i + 7] / 2)
            a2_2 = a1 + (bull[i + 12] / 2);
            a2_3 = a1 + (bull[i + 13] / 2)
            a1_0 = a1 - (bull[i + 6] / 2);
            a1_1 = a1 - (bull[i + 7] / 2)
            a1_2 = a1 - (bull[i + 12] / 2);
            a1_3 = a1 - (bull[i + 13] / 2)

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi1, phi2)
            lblocks.append(b)
        else:
            a2_0 = a1 + (bull[5 + (i + 6 - 5)] / 2)
            a2_1 = a1 + (bull[0 + (i + 6 - 5)] / 2)
            a2_2 = a1 + (bull[11 + (i + 6 - 5)] / 2)
            a2_3 = a1 + (bull[6 + (i + 6 - 5)] / 2)
            a1_0 = a1 - (bull[5 + (i + 6 - 5)] / 2)
            a1_1 = a1 - (bull[0 + (i + 6 - 5)] / 2)
            a1_2 = a1 - (bull[11 + (i + 6 - 5)] / 2)
            a1_3 = a1 - (bull[6 + (i + 6 - 5)] / 2)

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi1, phi2)
            lblocks.append(b)

    print("Gerando blocos da regiao apical")
    for i in range(6):
        nblocks += 1
        print("Bloco %d" % (nblocks))
        block = i
        theta0 = (0 + i) * pi / 3
        theta1 = (1 + i) * pi / 3
        if (i < 5):
            a2_0 = a1 + (bull[i + 12] / 2);
            a2_1 = a1 + (bull[i + 13] / 2)
            a1_0 = a1 - (bull[i + 12] / 2);
            a1_1 = a1 - (bull[i + 13] / 2)

            b = Block(nblocks, a1_0, a1_1, None, None, a2_0, a2_1, None, None,
                      theta0, theta1, phi2, phi3)
            lblocks.append(b)
        else:
            a2_0 = a1 + (bull[17] / 2);
            a2_1 = a1 + (bull[12] / 2)
            a1_0 = a1 - (bull[17] / 2);
            a1_1 = a1 - (bull[12] / 2)

            b = Block(nblocks, a1_0, a1_1, None, None, a2_0, a2_1, None, None,
                      theta0, theta1, phi2, phi3)
            lblocks.append(b)

    # imprime lista de blocos
    for b in lblocks:
        print(b)

    print("Criacao de blocos OK")

    return lblocks


# ------------------------------------------------------------------------------

def geraLVSplines(blocks, a1, c1, c2, phi0):
    # inicializacao

    ncirc = 18
    nlong = 18

    phi0, phi1 = -pi / 6, pi / 2
    theta0, theta1 = 0, 2 * pi

    dt = (theta1) / (ncirc - 1)
    dp = (phi1 - phi0) / (nlong - 1)

    print("Delta theta = %f" % dt)
    print("Delta phi   = %f" % dp)

    print("Gmsh .geo file generation;\n;\n;\n")

    # cria nos do apex - endo e epi
    x, y, z = ellipsoid(0, c1, phi1, 0)
    apexEndo = createVertex(x, y, z)
    x, y, z = ellipsoid(0, c2, phi1, 0)
    apexEpi = createVertex(x, y, z)

    LEpiBase = []
    LEndoBase = []
    ZEndoBase = []


    # epicardium
    theta = theta0
    
    splst1 = []
    splst2= []
    splst3 = []
    splst4 = []
    EpisplstVertex1 = []
    EpisplstVertex2 = []
    EpisplstVertex3 = []
    EpisplstVertex4 = []

    for j in range(ncirc-1):
        phi = phi0
        vlist = []
        for i in range(nlong - 1):
            myblock = findBlock(blocks, theta, phi)
            a1_0, a1_1, a1_2, a1_3 = myblock.getEpiData()

            bphi0, bphi1 = myblock.getPhi()
            btet0, btet1 = myblock.getTheta()

            if (not myblock.isApex()):
                ai = lerp(a1_0, a1_2, bphi0, bphi1, phi)
                af = lerp(a1_1, a1_3, bphi0, bphi1, phi)
                a = lerp(ai, af, btet0, btet1, theta)
            else:
                a = lerp(a1_0, a1_1, btet0, btet1, theta)

            x, y, z = ellipsoid(a, c2, phi, theta)
            
            
            
            vid = createVertex(x, y, z)

            

            vlist.append(vid)

            if phi == phi0:
                LEpiBase.append(vid)
                ZEndoBase.append(z)

            phi = phi + dp

        vlist.append(apexEpi)


        sid = createSpline(vlist[0:7])
        splst1.append(sid)
        EpisplstVertex1.append(vlist[0])
        sid = createSpline(vlist[6:11])
        splst2.append(sid)
        EpisplstVertex2.append(vlist[6])
        sid = createSpline(vlist[10:15])
        splst3.append(sid)
        EpisplstVertex3.append(vlist[10])
        sid = createSpline(vlist[14:])
        splst4.append(sid)
        EpisplstVertex4.append(vlist[14])


        theta = theta + dt

    epiSID, surfID, epiBaseSplines = createSurfaceSkin([splst1, splst2, splst3, splst4], [EpisplstVertex1, EpisplstVertex2, EpisplstVertex3, EpisplstVertex4])

    ID=createCompoundSurface(surfID)
    createPhysicalSurface(ID, 40)
    

    # endocardium
    theta = theta0
    splst1 = []
    splst2= []
    splst3 = []
    splst4 = []
    EndosplstVertex1 = []
    EndosplstVertex2 = []
    EndosplstVertex3 = []
    EndosplstVertex4 = []

    for j in range(ncirc - 1):
        phi = phi0
        vlist = []
        for i in range(nlong - 1):
            # encontra o bloco atual
            myblock = findBlock(blocks, theta, phi)

            # recupera informacoes do bloco
            a1_0, a1_1, a1_2, a1_3 = myblock.getEndoData()

            bphi0, bphi1 = myblock.getPhi()
            btet0, btet1 = myblock.getTheta()

            if (not myblock.isApex()):
                ai = lerp(a1_0, a1_2, bphi0, bphi1, phi)
                af = lerp(a1_1, a1_3, bphi0, bphi1, phi)
                a = lerp(ai, af, btet0, btet1, theta)
            else:
                a = lerp(a1_0, a1_1, btet0, btet1, theta)

            # cria ponto
            x, y, z = ellipsoid(a, c1, phi, theta)
            
            if phi == phi0: z=ZEndoBase[j]

            #print 'z: ', z, ' phi: ', phi, 'theta: ', theta
            
            vid = createVertex(x, y, z)
            # adiciona vertice na lista
            vlist.append(vid)

            if phi == phi0: 
                LEndoBase.append(vid)

            # incrementa angulo phi
            phi = phi + dp

        vlist.append(apexEndo)



        #print("Criando spline %d" % sid)
        # adiciona curva (spline) na lista
        sid = createSpline(vlist[0:7])
        splst1.append(sid)
        EndosplstVertex1.append(vlist[0])
        sid = createSpline(vlist[6:11])
        splst2.append(sid)
        EndosplstVertex2.append(vlist[6])
        sid = createSpline(vlist[10:15])
        splst3.append(sid)
        EndosplstVertex3.append(vlist[10])
        sid = createSpline(vlist[14:])
        splst4.append(sid)
        EndosplstVertex4.append(vlist[14])


        
        # incrementa angulo theta
        theta = theta + dt


    # TODO: aprender a criar a surface no GMSH, problemas com Line Loop
    #print("Criando superficie do endocardio")
    endoSID, surfID, endoBaseSplines = createSurfaceSkin([splst1, splst2, splst3, splst4], [EndosplstVertex1, EndosplstVertex2, EndosplstVertex3, EndosplstVertex4])

    
    ID = createCompoundSurface(surfID)
    createPhysicalSurface(ID, 30)

    # pega id da curva da base no endocardio
    #endoCID = cubit.get_last_id('curve')
    #
    # print("Curva da base endocardio Id = %d" % endoCID)
    # print("Criando superficie da cavidade")
    # endoBaseSID = createSurface([endoCID])
    #
    # print("Criando volume da cavidade")
    # cstr = "create volume surface %d %d heal keep" % (int(endoBaseSID), int(endoSID))
    # print(cstr)
    # cubit.cmd(cstr)
    # CavidadeVID = cubit.get_last_id('volume')
    # print(CavidadeVID, type(CavidadeVID))
    #
    # # calcula e mostra o volume da cavidade
    # cubit.cmd('set echo on')
    # cubit.cmd('list volume %d geometry' % CavidadeVID)
    # cubit.cmd('set echo off')
    #
    # # limpa cavidade
    # cubit.cmd('delete surface %d' % endoBaseSID)
    # cubit.cmd('delete volume %d' % CavidadeVID)
    #

    baseSID, surfID = createSurfaceBase(endoBaseSplines, epiBaseSplines)

    print endoBaseSplines
    print epiBaseSplines

    createPhysicalSurface(surfID, 10)
    
    surflist = [endoSID, epiSID, baseSID]
    veID = createVolume(surflist)

    createPhysicalVolume(veID, 0)




    # print("Criando superficie do epicardio")
    # epiSID = createSurfaceSkin(splst)
    # epiLID = cubit.get_last_id('curve')
    #
    # print("Criando superficie da base")
    # lBaseLines = []
    # if len(LEpiBase) == len(LEndoBase):
    #     for i in range(len(LEpiBase)):
    #         a = LEndoBase[i]
    #         b = LEpiBase[i]
    #         linha = [a, b]
    #         c = createLine(linha)
    #         lBaseLines.append(c)
    # else:
    #     print("Erro")
    #     sys.exit(1)
    # baseSID = createSurfaceSkin(lBaseLines)
    #
    # print("Criando volume do VE")
    # surflist = [endoSID, epiSID, baseSID]
    # veID = createVolume(surflist)

    # print("Limpando")
    # cubit.cmd("delete vertex all")
    # cubit.cmd("delete curve all")
    # cubit.cmd("delete surface all")
    #
    # cubit.cmd('set echo on')
    #
    # tam1 = c2 / 3
    # tam2 = 2 * c2 / 3
    #
    # print("Ajeitando")
    # cubit.cmd("imprint all")
    # cubit.cmd("merge all")
    #
    # if (malhaTet):
    #     criaMalhaTetraedro([veID])
    # else:
    #     criaMalhaHexaedro()
    #
    # criaMalhaSurperficieEndocardio()
    # criaMalhaSurperficieEpicardio()
    # salvaMalha()


# -----------------------------------------------------------------------------

def salvaMalha():
    cubit.cmd('set abaqus precision 6')
    cubit.cmd('export abaqus "saida.inp" dimension 3 everything overwrite')
    cubit.cmd('export STL binary "saida.stl" mesh overwrite')


# -----------------------------------------------------------------------------

def criaMalhaSurperficieEndocardio():
    cubit.cmd('set duplicate block elements off')
    cubit.cmd('block 3333 surface 1')


# -----------------------------------------------------------------------------

def criaMalhaSurperficieEpicardio():
    cubit.cmd('set duplicate block elements off')
    cubit.cmd('block 4444 surface 5')


# -----------------------------------------------------------------------------

def criaMalhaHexaedro():
    pass


# -----------------------------------------------------------------------------

def criaMalhaTetraedro(baseName):

    cmdstr = '/home/joventino/gmsh/bin/gmsh ' + baseName + '.geo -3 -o ' + baseName + '.msh'

    os.system(cmdstr)
    #os.system('gmsh ./lv_teste.geo -3 -o ./teste.vtk')

# -----------------------------------------------------------------------------

def criaFibra(baseName):

	from dolfin import *
	from fiberrules import *
	from dolfin_utils.meshconvert import meshconvert

	parameters["form_compiler"]["quadrature_degree"] = 2
	parameters.allow_extrapolation = True

	#convert mesh to dolfin format in order to use in Fiberrules
	ifilename = baseName + '.msh'
	ofilename = baseName + '_fenics.xml'
	iformat = 'gmsh'
	meshconvert.convert2xml(ifilename, ofilename, iformat=iformat)

	#Load mesh for fiberrules
	mesh = Mesh(ofilename)
	domains = MeshFunction("size_t", mesh, baseName + '_fenics_facet_region.xml')

	for facet in facets(mesh):
		mesh.domains().set_marker((facet.index(), domains[facet]), 2)

	#Fibers varying transmurally 60 to -60
	fiber_angle_epi = 150 #150 - 90 = 60
	fiber_angle_endo = 30 #30 - 90 = -60

	fiber_space = FunctionSpace(mesh, 'DG',0)
	#fiber_space = FunctionSpace(mesh, "Lagrange", 1)

	fibers, sheet_normals, cross_sheet = dolfin_fiberrules(mesh, fiber_space, fiber_angle_epi, fiber_angle_endo)

	dolfin_to_vtk(fibers, baseName)
	dolfin_to_carpfile(mesh, baseName)
	fibers_to_carpfile(mesh, fibers, baseName)


# ------------------------------------------------------------------------------

def convert2cardiax(material_params, pvloop_params, baseName):
	import cardiax_gmsh2xml as c

	c.gmsh2xml(baseName + '.msh', baseName + '_cardiax.xml', material_params, pvloop_params, 'bullseye/pvloop_data.txt')



# ------------------------------------------------------------------------------



def run(a1, c1, bull, material_params, pvloop_params, baseName):


    print baseName, bull

    global geoFile

    fileName = baseName + '.geo'
    geoFile = open(fileName,'w')


    c2 = c1 + bull[18]

    reduction = c1*sin(pi/6)/(1.0 + sin(pi/6))
    c1 = c1 - reduction
    a1 = a1/2.0 # diameter/2
    print 'C1: ',c1

    # teste
    reduction = c2*sin(pi/6)/(1.0 + sin(pi/6))
    c2 = c2 - reduction

    print 'bull[18]: ', bull[18]
    print 'C2: ', c2

    # cria lista com os blocos e suas informacoes
    blocklist = createBlocks(a1, bull)

    # define geometria e gera malha
    geraLVSplines(blocklist, a1, c1, c2, 0)

    geoFile.close()

    criaMalhaTetraedro(baseName)

    criaFibra(baseName)

    convert2cardiax(material_params, pvloop_params, baseName)




# MAIN -------------------------------------------------------------------------

# Diastole1
#arq = 'parametrosNovos2Diastole.txt'
#arq = 'parametrosBaiDiastole.txt'
#arq = 'parametrosRegular.txt'
#val = 1.59
#val = 1.0

## Sistole
##arq = 'parametrosNovos1Sistole_metro.txt'

## le arquivos de parametros
##a1, c1, bull = readBullseye(arq)

##print a1, c1, bull

##run(a1, c1, bull, 'lvTeste')
