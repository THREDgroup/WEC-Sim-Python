import numpy as np
import warnings
from stl import mesh
import trimesh
import time
def import_stl_fast(filename,mode,*eol):
    """
    # function to import ASCII STL files into matlab. 
    """ 
    # method using numpy-stl
    # your_mesh = mesh.Mesh.from_file(filename)
    # v = your_mesh.vectors.reshape(int(np.size(your_mesh.vectors)/3),3)
    # n = your_mesh.normals
    # a = np.transpose(your_mesh.areas)
    
    #method using trimesh
    your_mesh = trimesh.load_mesh(filename)
    n = your_mesh.face_normals
    c = your_mesh.triangles_center
    a = your_mesh.area_faces
    v = your_mesh.triangles.reshape(int(np.size(your_mesh.triangles)/3),3)
    [vertex,faces] = np.unique(v,return_inverse=True,axis=0)
    face = faces.reshape(int(np.size(faces)/3),3)+1
    # ver = your_mesh.vertices
    # vert = ver[np.lexsort((ver[:,2], ver[:,1],ver[:,0]))]
    
    if mode == 1:
        vout = [a,c,n,vertex,face]
    elif mode == 2:
        vout = [a,c,n,v]
    else:
        warnings.warn("invalid mode",DeprecationWarning)
        vout = None
    return(vout)
tick = time.clock()
[area,center,normal,vertex,face]= import_stl_fast("elipsoid.stl",1,1)


tock = time.clock()
results = tock-tick

your_mesh_v1 = trimesh.load_mesh("elipsoid.stl")

# face = your_mesh_v1.faces_unique_edges
# faces = face[np.lexsort((face[:,0], face[:,1],face[:,2]))]


tick2 = time.clock()
your_mesh = mesh.Mesh.from_file("elipsoid.stl")
v2 = your_mesh.vectors.reshape(int(np.size(your_mesh.vectors)/3),3)
n = your_mesh.normals
a = np.transpose(your_mesh.areas)
tock2 = time.clock()
results2 = tock2-tick2
# your_meshs = mesh.Mesh.from_file("elipsoid.stl")
# volume, cog, inertia = your_meshs.get_mass_properties()
#mesh2 = trimesh.load_mesh("elipsoid.stl")
# %open file
# fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
# if fid == -1
#     error('File could not be opened, check name or path.')
# end

# %=====================================================

# fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.

# if eol == 1
#     fmt = '%*s %*s %f32 %f32 %f32 \n %*s %*s \n %*s %f32 %f32 %f32 \n %*s %f32 %f32 %f32 \n %*s %f32 %f32 %f32 \n %*s \n %*s \n';
# elseif eol == 2
#     fmt = '%*s %*s %f32 %f32 %f32 \r\n %*s %*s \r\n %*s %f32 %f32 %f32 \r\n %*s %f32 %f32 %f32 \r\n %*s %f32 %f32 %f32 \r\n %*s \r\n %*s \r\n';
# end
# C=textscan(fid, fmt, 'HeaderLines', 1);
# fclose(fid);

# %extract normal vectors and vertices
# tnorm = cell2mat(C(1:3));
# tnorm = double(tnorm);

# v1 = cell2mat(C(4:6));
# v2 = cell2mat(C(7:9));
# v3 = cell2mat(C(10:12));

# if isnan(C{1}(end))
#     tnorm = tnorm(1:end-1,:); %strip off junk from last line
# end

# if isnan(C{4}(end))
#     v1 = v1(1:end-1,:); %strip off junk from last line
#     v2 = v2(1:end-1,:); %strip off junk from last line
#     v3 = v3(1:end-1,:); %strip off junk from last line
# end

# v_temp = [v1 v2 v3]';
# v = zeros(3,numel(v_temp)/3);

# v(:) = v_temp(:);
# v = v';

#     varargout = cell(1,nargout);
#     switch mode
#         case 1
#             [p,t]=fv2pt(v,length(v)/3);%gets points and triangles

#             varargout{1} = p;
#             varargout{2} = t;
#             varargout{3} = tnorm;
#         case 2
#             varargout{1} = v;
#             varargout{2} = tnorm;
#     end
# end

# %%
# function [p,t]=fv2pt(v,fnum)

# %gets points and triangle indexes given vertex and facet number
# c=size(v,1);

# %triangles with vertex id data
# t=zeros(3,fnum);
# t(:)=1:c;

# %now we have to keep unique points fro vertex
# [p,~,j]=unique(v,'rows'); %now v=p(j) v(i)=p;
# t(:)=j(t(:));
# t=t';

# end