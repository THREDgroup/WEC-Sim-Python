# BEMIO Class
BEMIO class reads data generated from BEM solver and creates HDF5 file.
## Procedure:
1. Use CAD software (e.g. Rhinoceros 3D, SOLIDWORKS, and etc.) to generate .gdf or .stl file. Rhinoceros 3D is more preferable than SOLIDWORKS.
2. Run meshConverter.py to change .gdf(Rhinoceros 3D) or .stl(SOLIDWORKS) file to .gdf(WAMIT), .dat(NEMOH), .out(AQWA) type file. 
3. Use BEM solver.
4. Place necessary files into BEMIO/data directory.
5. Run BEMIO.py to generate HDF5 file.

# Simulation Class
To be updated.

# Wave Class
Wave Class creates wave and its related varaibles.
## Procedure:
1. Import necessay properties from SimulationClass.
2. Generate a WaveClass (wave type).
3. Run waveSetup method to produce wave.
4. Check produced wave by running a plotETA method.
5. Check wave spectrum for irregular wave by running plotSpectrum method.

# Body Class
To be updated.

# Note

# BEM Solver Guide
## Procedure (WAMIT):
To be updated.

## Procedure (NEMOH-Windows):
1. Download NEMOH file from 
2. Create run.bat file on working directory
3. To be updated.

## Procedure (NEMOH-Linux, MAC):
1. Download and extract NEMOH file
2. Go to the NEMOH directory and run make file
3. Copy all files located in the tutorial directory
4. Go to the bin directory and paste files
5. Run preProc
6. Run solver
7. Run postProc
8. Move FKForce.tac and ..... to BEMIO/data directory
9. To be updated

## Procedure (AQWA):
To be updated 

## API
@software{trimesh,
	author = {{Dawson-Haggerty et al.}},
	title = {trimesh},
	url = {https://trimsh.org/},
	version = {3.2.0},
	date = {2019-12-8},
}

