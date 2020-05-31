#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:28:49 2020

@author: logical
"""
import os
import waveClass
class ParaviewClass:
        def waveElevationGrid(self, t, X, Y):
        """
        

        Parameters
        ----------
        t : int
            the current time
        X : list
            (m x n) matrix of X coordinates at which to calculate
            the wave elevation
        Y : list
            (m x n) matrix of Y coordinates at which to calculate
            the wave elevation.

        Returns
        -------
        Z - (m x n) matrix of Z coordinates of the wave elevation

        """
        if (self.wType == 'noWave') and (self.wType == 'noWaveCIC') and (self.wType == 'etaImport'):
            Z = np.zeros((np.size(X),np.size(X)))
        elif (self.wType == 'regular') and (self.wType == 'regularCIC'):
            Xt = X*np.cos(self.waveDir*np.pi/180) + Y*np.sin(self.waveDir*np.pi/180)
            Z = self.A*np.cos(-1*self.k*Xt + self.wF*t)
        elif (self.wType == 'irregular') and (self.wType == 'spectrumImport'):
            Z = np.zeros((np.size(X),np.size(X)))
            Xt = X*np.cos(self.waveDir*np.pi/180) + Y*np.sin(self.waveDir*np.pi/180)
            for iw in range(1,len(self.wF)+1):
                Z = Z + np.sqrt(self.wA(iw)*self.dw(iw))*np.cos(-1*self.k(iw)*Xt + self.wF(iw)*t + self.phase(iw))
        return(Z)
    
    def write_paraview_vtp_wave(self, t, numPointsX, numPointsY, domainSize, model, simdate, mooring):
        """
        Write vtp files for visualization using Paraview        

        Parameters
        ----------
        t : TYPE
            DESCRIPTION.
        numPointsX : TYPE
            DESCRIPTION.
        numPointsY : TYPE
            DESCRIPTION.
        domainSize : TYPE
            DESCRIPTION.
        model : TYPE
            DESCRIPTION.
        simdate : TYPE
            DESCRIPTION.
        mooring : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        #ground plane
        hDir = os.getcwd()
        os.mkdir('vtk')
        vDir = os.path.join(hDir, 'vtk')
        os.chdir(vDir)
        filename = 'ground.txt'
        fid = open(filename,"w+")
        fid.write(str(domainSize) + "\n")
        fid.write(str(self.waterDepth) + "\n")
        fid.write(str(mooring) + "\n")
        fid.close()
        #wave
        x = np.linspace(-domainSize, domainSize, numPointsX)
        y = np.linspace(-domainSize, domainSize, numPointsY)
        [X,Y] = np.meshgrid(x,y)
        lx = np.size(x,1)
        ly = np.size(y,1)
        numVertex = lx * ly
        numFace = (lx-1) * (ly-1)
        for it in range(0,len(t)):
            #open file
            os.mkdir('waves')
            wDir = os.path.join(hDir, 'waves')
            os.chdir(wDir)
            filenames = 'waves_'+str(it)+'.vtp'
            fids = open(filenames,"w+")
            # calculate wave elevation
            Z = self.waveElevationGrid(t[it], X, Y)
            # write header
            fids.write('<?xml version="1.0"?>\n')
            fids.write('<!-- WEC-Sim Visualization using ParaView -->\n')
            fids.write('<!--   model: ' , model , ' - ran on ' , simdate , ' -->\n')
            fids.write('<!--   wave:  ' , self.wType ,' -->\n')
            fids.write('<!--   time:  ' , str(t(it)) , ' -->\n')
            fids.write('<VTKFile type="PolyData" version="0.1">\n')
            fids.write('  <PolyData>\n')
            # write wave info
            fids.write('    <Piece NumberOfPoints="' + str(numVertex) + '" NumberOfPolys="' + str(numFace) + '">\n')
            # write points
            fids.write('      <Points>\n')
            fids.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
            for jj in range(1,len(y)+1):
                for ii in range(1,len(x)+1):
                    pt = [X[jj][ii], Y[jj][ii], Z[jj][ii]]
                    """
                    if not work try np.array on X,Y,Z
                    """
                    fids.write('          {:5.5f} {:5.5f} {:5.5f}\n'.format(pt[0],pt[1],pt[2]))
            fids.write('        </DataArray>\n')
            fids.write('      </Points>\n')
            # write squares connectivity
            fids.write('      <Polys>\n')
            fids.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
            for jj in range(1,ly):
                for ii in range(1,lx):
                    p1 = (jj-1)*lx + (ii-1)
                    p2 = p1+1
                    p3 = p2 + lx
                    p4 = p1 + lx
                    fids.write('          {} {} {} {}\n'.format(p1,p2,p3,p4))
            fids.write('        </DataArray>\n')
            fids.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
            fids.write('         ')
            for ii in range(1,numFace+1):
                n = ii * 4
                fids.write(' {}'.format(n))
            fids.write('\n')
            fids.write('        </DataArray>\n')
            fids.write('      </Polys>\n')
            # end file
            fids.write('    </Piece>\n')
            fids.write('  </PolyData>\n')
            fids.write('</VTKFile>')
            #close file
            fids.close()
     