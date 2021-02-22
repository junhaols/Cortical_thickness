function [Cortical_Thickness,Pial_area,White_area] = LJH_Cortical_Thickness_area(surf_pial,surf_white)
%% This function is used to calculate the thickness and the areas.
% The surface format should be freesurfer
% surface such as '?h.pial','?h.white' or gifti format as 'xxxx.gii'.
% For freesurfer surface,the function "read_surf.m" should be set path,this
% function is in the path of you freesurfer:~/freesurfer/matlab/read_surf.m
% For gifti surface,the function "gifti.m" should be set path,this function
% is in the gifti package.

%% input
%    surf_pial:the pial surface
%    surf_white:the whote surface
%% output
%    Cortical_Thickness:The Cortical Thickness.
%    Pial_area:The area of pial surface.
%    White_area:The area of white surface.

tic
%surf_pial='/Users/Junhao/data/ToolBox/freesurfer_matlab/lh.pial';
%surf_white='/Users/Junhao/data/ToolBox/freesurfer_matlab/lh.white';

%% check the files format
[~,~,suffix1]=fileparts(surf_pial);
[~,~,suffix2]=fileparts(surf_white);

if strcmp(suffix1,'.pial')&&strcmp(suffix2,'.white')
    [vertex1, face1,~] = read_surf(surf_pial);
    [vertex2, face2,~] = read_surf(surf_white);
    min_index=min(face1);
    if ismember(0,min_index) %% freesurfer surface:The vertex index in the faces start from 0.We should add 1 to all the vertex index.
        face1=face1+1;
        face2=face2+1;
    end
elseif strcmp(suffix1,'.gii')&&strcmp(suffix2,'.gii')
    surf1=gifti(surf_pial);
    vertex1=surf1.vertices;%% gifti surface:The vertex index start from 1.
    face1=surf1.faces;
    
    surf2=gifti(surf_white);
    vertex2=surf2.vertices;
    face2=surf2.faces;
else
    disp('The format of surface files should be "xxx.pial,xxx.white or xxx.gii"');
end
    
     
%% for thickness
nVert=length(vertex1);
thickness=zeros(nVert,1);
%%
%nParts=100;
%flag=mod(nVert,nParts);

for i=1:nVert
    
    % the dist1
    data1=vertex1(i,:);
    data1=repmat(data1,[nVert,1]);
    det_data1=data1-vertex2;  
    Dist1=sqrt(sum(det_data1.^2,2));
    minDist1=min(Dist1);
    
    thickness1=minDist1;
    
    % the dist2
 
    data2=vertex2(i,:);
    data2=repmat(data2,[nVert,1]);
    det_data2=data2-vertex1;  
    Dist2=sqrt(sum(det_data2.^2,2));
    minDist2=min(Dist2);   
    thickness2=minDist2;  
    
    % thickness
    thickness(i)=0.5*(thickness1+thickness2);
end



%% for areas.

% A1=f1(:,1);
% B1=f1(:,2);
% C1=f1(:,3);
%% pial-surf
A1=vertex1(face1(:,1),:);
B1=vertex1(face1(:,2),:);
C1=vertex1(face1(:,3),:);
a1 = sqrt((A1(:,1)-B1(:,1)).^2+(A1(:,2)-B1(:,2)).^2+(A1(:,3)-B1(:,3)).^2);
b1 = sqrt((C1(:,1)-B1(:,1)).^2+(C1(:,2)-B1(:,2)).^2+(C1(:,3)-B1(:,3)).^2);
c1 = sqrt((A1(:,1)-C1(:,1)).^2+(A1(:,2)-C1(:,2)).^2+(A1(:,3)-C1(:,3)).^2);
p1 = (a1+b1+c1)/2;
%% S1 is the area of all the triangles (pial-surf).
S1= sqrt(p1.*(p1-a1).*(p1-b1).*(p1-c1));

%% white-surf
A2=vertex2(face2(:,1),:);
B2=vertex2(face2(:,2),:);
C2=vertex2(face2(:,3),:);
a2 = sqrt((A2(:,1)-B2(:,1)).^2+(A2(:,2)-B2(:,2)).^2+(A2(:,3)-B2(:,3)).^2);
b2 = sqrt((C2(:,1)-B2(:,1)).^2+(C2(:,2)-B2(:,2)).^2+(C2(:,3)-B2(:,3)).^2);
c2 = sqrt((A2(:,1)-C2(:,1)).^2+(A2(:,2)-C2(:,2)).^2+(A2(:,3)-C2(:,3)).^2);
p2 = (a2+b2+c2)/2;
%% S2 is the area of all the triangles (white-surf).
S2= sqrt(p2.*(p2-a2).*(p2-b2).*(p2-c2));

% Each vertex gets one third of the area of each triangle it is a part of.
S1_vert=zeros(length(vertex1),1);
S2_vert=zeros(length(vertex2),1);
for i=1:nVert
    
    % for pial surface
      % find the tringles that include the point.
      index1=find(face1==i);
      tri1=zeros(length(index1),1);
     
      for j=1:length(index1)
          det1=index1(j)-length(face1);
          if det1<=0
              tri1(j)=index1(j);
          elseif det1>0 && det1<=length(face1)
              tri1(j)=index1(j)-length(face1);
          else
              tri1(j)=index1(j)-2*length(face1);
          end
      end
    
      S1_vert(i)=1/3*sum(S1(tri1));% the vertice-wise area is 1/3 of each triangle it is a part of.
      
     % for white surface
     
      index2=find(face2==i);
      tri2=zeros(length(index2),1);
     
      for j=1:length(index2)
          det2=index2(j)-length(face2);
          if det2<=0
              tri2(j)=index2(j);
          elseif det2>0 && det2<=length(face2)
              tri2(j)=index2(j)-length(face2);
          else
              tri2(j)=index2(j)-2*length(face2);
          end
      end
      S2_vert(i)=1/3*sum(S2(tri2));   
      
end
      
 Cortical_Thickness=thickness;
 Pial_area=S1_vert;
 White_area=S2_vert;

      
toc              
          
    
    
    
    
    









