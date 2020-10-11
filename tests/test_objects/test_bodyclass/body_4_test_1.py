
import bodyClass_V1
import numpy as np
from scipy import interpolate
import scipy.io as sio

def readData(file):
    matFile = sio.loadmat(file) 
    keys = list(matFile.keys())[-1]
    datas = matFile[keys]
    return datas

body_4_1 = bodyClass_V1.BodyClass('rm3.h5') #regular B2B_Case1:b2b = 0

w = np.conj(np.transpose(readData("./testData/body_4_test/w.mat"))) 
waveDir = [0]
CIkt = 601
CTTime = readData("./testData/body_4_test/CTTime.mat")[0]
dt = 0.1
rho = 1000
g = 9.81
waveType = 'regular'
waveAmpTime = np.conj(np.transpose(readData("./testData/body_4_test/waveAmpTime.mat"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 0
nlHydro = 0
B2B = 0
numFreq = []
body_4_1.bodyNumber = 1
body_4_1.bodyTotal = [2]
body_4_1.readH5file()
body_4_1.hydroStiffness = np.zeros((6, 6))
body_4_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_4_1.linearDamping = np.zeros((6, 6))
body_4_1.mass = 'equilibrium'
body_4_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
adjMassWeightFun = 5
body_4_1.momOfInertia = [20907301, 21306090.66, 37085481.11];  
body_4_1.cg = np.array([0 , 0, -10])
body_4_1.checkinputs()

# result1 = readData("./testData/body_4_test/body1_fAddedMass.mat")
# np.testing.assert_allclose(body_4_1.hydroForce['fAddedMass'], result1)
