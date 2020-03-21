clc;clear all;close all;
a = [1 3 4 7];
b = [4 5 6 8];
c = [4 5 5 8];
[mesh_x,mesh_y] = meshgrid(linspace(min(a),max(a),30), ...
    linspace(min(b),max(b),30));
zz = griddata(a,b,c,mesh_x(:),mesh_y(:),'v4');
zz = reshape(zz,size(mesh_x));
mesh(mesh_x,mesh_y,zz);
plot3(a,b,c,'ro',mesh_x,mesh_y,zz,'b','LineWidth',2);
%plot3(a,b,c,'Marker','o');
