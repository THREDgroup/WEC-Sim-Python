
import bodyClass_V1
import numpy as np
from scipy import interpolate
import scipy.io as sio

def readData(file):
    matFile = sio.loadmat(file) 
    keys = list(matFile.keys())[-1]
    datas = matFile[keys]
    return datas

body_5_1 = bodyClass_V1.BodyClass('rm3.h5') #regular B2B_Case1:b2b = 0
body_5_2 = bodyClass_V1.BodyClass('rm3.h5') #regular B2B_Case1:b2b = 0

w = np.conj(np.transpose(readData("./testData/body_5_test/w.mat"))) 
waveDir = [0]
CIkt = 601
CTTime = readData("./testData/body_5_test/CTTime.mat")[0]
dt = 0.1
rho = 1000
g = 9.81
waveType = 'regular'
waveAmpTime = np.conj(np.transpose(readData("./testData/body_5_test/waveAmpTime.mat"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 0
nlHydro = 0
B2B = 1
numFreq = []
body_5_1.bodyNumber = 1
body_5_1.bodyTotal = [2]
body_5_1.readH5file()
body_5_1.hydroStiffness = np.zeros((6, 6))
body_5_1.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_5_1.linearDamping = np.zeros((6, 6))
body_5_1.mass = 'equilibrium'
body_5_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
adjMassWeightFun = 5
body_5_1.momOfInertia = [20907301, 21306090.66, 37085481.11]

iBod = 2 #body number
body_5_2.bodyNumber = 2
body_5_2.bodyTotal = [2]
body_5_2.readH5file()
body_5_2.hydroStiffness = np.zeros((6, 6))
body_5_2.viscDrag = {'Drag':np.zeros((6, 6)),'cd':np.zeros(6),'characteristicArea':np.zeros(6)}
body_5_2.linearDamping = np.zeros((6, 6))
body_5_2.mass = 'equilibrium'
body_5_2.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
adjMassWeightFun = 5
body_5_2.momOfInertia = [94419614.57, 94407091.24, 28542224.82]
AAAA=body_5_1.hydroForce['fAddedMass']
acc =np.array([[80, 33, 90, 62, 34, 89],[72, 14, 56, 93, 22, 19]])
iBod =2
fam = np.zeros(np.shape(acc))
for i in range(6):
    tmp = np.zeros(np.size(acc[:,i]))
    for j in range(6):
        if B2B == 1:
            jj = (iBod-1)*6+j
        else:
            jj = j
        iam = AAAA[i,jj]
        tmp = tmp + acc[:,j]* iam
    fam[:,i] = tmp

# result1 = readData("./testData/body_5_test/body1_fAddedMass_adjusted.mat")
# np.testing.assert_allclose(body_5_1.hydroForce['fAddedMass'], result1)

# result2 = readData("./testData/body_5_test/body2_fAddedMass_adjusted.mat")
# np.testing.assert_allclose(body_5_2.hydroForce['fAddedMass'], result2)