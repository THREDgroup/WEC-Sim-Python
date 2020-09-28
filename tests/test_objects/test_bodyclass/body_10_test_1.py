
import bodyClass_V1
import numpy as np
from scipy import interpolate
import scipy.io as sio

def readData(file):
    matFile = sio.loadmat(file) 
    keys = list(matFile.keys())[-1]
    datas = matFile[keys]
    return datas

body_10_1 = bodyClass_V1.BodyClass('ellipsoid.h5')

w = np.conj(np.transpose(readData("./testData/body_10_test/w.mat"))) 
waveDir = [0]
CIkt = 1201
CTTime = readData("./testData/body_10_test/CTTime.mat")[0]
dt = 0.05
rho = 1025
g = 9.81
waveType = 'regular'
waveAmpTime = np.conj(np.transpose(readData("./testData/body_10_test/waveAmpTime.mat"))) 
iBod = 1 #body number
numBod = [] #later change it to 2 to check
ssCalc = 0
nlHydro = 2
B2B = 0
numFreq = []
body_10_1.bodyNumber = 1
body_10_1.bodyTotal = [1]
body_10_1.readH5file()
body_10_1.hydroStiffness = np.zeros((6, 6))
body_10_1.viscDrag = {'Drag':np.zeros((6, 6)),
                      'cd':np.array([1,0,1,0,1,0]),
                      'characteristicArea':np.array([25,0,np.pi*5**2,0,np.pi*5**5,0])}
body_10_1.linearDamping = np.zeros((6, 6))
body_10_1.mass = 'equilibrium'
body_10_1.bodyGeo("./testData/body_10_test/elipsoid.stl")

# if body_10_1.mass == 'equilibrium':
#     body_10_1.massCalcMethod = body_10_1.mass
#     if nlHydro == 0:
#         body_10_1.mass = body_10_1.hydroData['properties']['disp_vol'] * rho
#     else:
#         cg_tmp = body_10_1.hydroData['properties']['cg'][0]
#         z = np.conj(np.transpose(np.array(body_10_1.bodyGeometry['center'])))[2] + cg_tmp[2]
#         zr = [0 if x > 0 else x for x in z]
#         area = np.conj(np.transpose(body_10_1.bodyGeometry['area']))[0]
#         av = np.array([area,area,area])*-1*np.conj(np.transpose(body_10_1.bodyGeometry['norm']))
#         tmp = rho*np.array([zr, zr, zr])*-1*av
#         body_10_1.mass = sum(tmp[2])
#         AAA = body_10_1.mass
# elif body_10_1.mass == 'fixed':
#     body_10_1.massCalcMethod = body_10_1.mass
#     body_10_1.mass = 999
#     body_10_1.momOfInertia = [999, 999, 999]
# else:
#     body_10_1.massCalcMethod = 'user'

body_10_1.hydroForcePre(w,waveDir,CIkt,CTTime,numFreq,dt,rho,g,waveType,waveAmpTime,iBod,numBod,ssCalc,nlHydro,B2B)
result1 = readData("./testData/body_10_test/body10_linearHydroRestCoef.mat")
np.testing.assert_allclose(body_10_1.hydroForce['linearHydroRestCoef'], result1)
result2 = readData("./testData/body_10_test/body10_visDrag.mat")
np.testing.assert_allclose(body_10_1.hydroForce['visDrag'], result2)
result3 = readData("./testData/body_10_test/body10_linearDamping.mat")
np.testing.assert_allclose(body_10_1.hydroForce['linearDamping'], result3)
result4 = readData("./testData/body_10_test/body10_userDefinedFe.mat")
np.testing.assert_allclose(body_10_1.hydroForce['userDefinedFe'], result4)
result5 = readData("./testData/body_10_test/body10_re.mat")[0]
np.testing.assert_allclose(body_10_1.hydroForce['fExt']['re'], result5)
result6 = readData("./testData/body_10_test/body10_im.mat")[0]
np.testing.assert_allclose(body_10_1.hydroForce['fExt']['im'], result6)
result7 = readData("./testData/body_10_test/body10_md.mat")[0]
np.testing.assert_allclose(body_10_1.hydroForce['fExt']['md'], result7)
result8 = readData("./testData/body_10_test/body10_fAddedMass.mat")
np.testing.assert_allclose(body_10_1.hydroForce['fAddedMass'], result8)
result9 = readData("./testData/body_10_test/body10_fDamping.mat")
np.testing.assert_allclose(body_10_1.hydroForce['fDamping'], result9)
result10 = readData("./testData/body_10_test/body10_vertex.mat")
np.testing.assert_allclose(body_10_1.bodyGeometry['vertex'], result10)
result11 = readData("./testData/body_10_test/body10_face.mat")
np.testing.assert_allclose(body_10_1.bodyGeometry['face'], result11)
result12 = readData("./testData/body_10_test/body10_norm.mat")
np.testing.assert_allclose(body_10_1.bodyGeometry['norm'], result12)
result13 = readData("./testData/body_10_test/body10_area.mat")
np.testing.assert_allclose(body_10_1.bodyGeometry['area'], result13)
result14 = readData("./testData/body_10_test/body10_center.mat")
np.testing.assert_allclose(body_10_1.bodyGeometry['center'], result14)
#body_10_1.bodyGeometry['numVertex'], 1442)
#body_10_1.bodyGeometry['numFace'], 2880)
 