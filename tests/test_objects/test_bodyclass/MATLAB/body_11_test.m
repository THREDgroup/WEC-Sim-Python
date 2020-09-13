clear all, close all, clc
% %mode1 = import_stl_fast('ellipsoid/elipsoid.stl',1)
% %error
% mode2 = import_stl_fast('ellipsoid/elipsoid.stl',2)
% %mode1_eol1 = import_stl_fast('ellipsoid/elipsoid.stl',1,1)
% %error
% mode2_eol1 = import_stl_fast('ellipsoid/elipsoid.stl',2,1)
% %mode1_eol2 = import_stl_fast('ellipsoid/elipsoid.stl',1,2)
% %error on linux
% %mode2_eol2 = import_stl_fast('ellipsoid/elipsoid.stl',2,2)
% %error on linux
filename = 'ellipsoid/elipsoid.stl';
%open file
fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
if fid == -1
    error('File could not be opened, check name or path.')
end

%=====================================================

fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.

fmt = '%*s %*s %f32 %f32 %f32 \n %*s %*s \n %*s %f32 %f32 %f32 \n %*s %f32 %f32 %f32 \n %*s %f32 %f32 %f32 \n %*s \n %*s \n';

C=textscan(fid, fmt, 'HeaderLines', 1);
fclose(fid);

%extract normal vectors and vertices
tnorm = cell2mat(C(1:3));
tnorm = double(tnorm);

v1 = cell2mat(C(4:6));
v2 = cell2mat(C(7:9));
v3 = cell2mat(C(10:12));

if isnan(C{1}(end))
    tnorm = tnorm(1:end-1,:); %strip off junk from last line
end

if isnan(C{4}(end))
    v1 = v1(1:end-1,:); %strip off junk from last line
    v2 = v2(1:end-1,:); %strip off junk from last line
    v3 = v3(1:end-1,:); %strip off junk from last line
end

v_temp = [v1 v2 v3]';
v = zeros(3,numel(v_temp)/3);

v(:) = v_temp(:);
v = v';
varargout = cell(1,2);
varargout{1} = v;
varargout{2} = tnorm;
%     varargout = cell(1,nargout);
%     switch mode
%         case 1
%             [p,t]=fv2pt(v,length(v)/3);%gets points and triangles
% 
%             varargout{1} = p;
%             varargout{2} = t;
%             varargout{3} = tnorm;
%         case 2
%             varargout{1} = v;
%             varargout{2} = tnorm;
%     end
    
 function [p,t]=fv2pt(v,fnum)

%gets points and triangle indexes given vertex and facet number
c=size(v,1);

%triangles with vertex id data
t=zeros(3,fnum);
t(:)=1:c;

%now we have to keep unique points fro vertex
[p,~,j]=unique(v,'rows'); %now v=p(j) v(i)=p;
t(:)=j(t(:));
t=t';

end


