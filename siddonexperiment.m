
% the points
p1x =0 ;
p1y =0 ;
p1z = 10;

p2x = 0;
p2y = 0;
p2z = -3;

P1 = [p1x p1y p1z ];
P2 = [p2x p2y p2z];

% the number of planes
Nx =4 ;
Ny = 4;
Nz = 4;

% separation between planes
dx = 1;
dy = 1;
dz = 1;
d=1;

% position of first plane (in x or y or z)
X1 = -2;
Y1 = -2;
Z1 = -2;

% define i j k
i= 1:1:Nx;
j= 1:1:Ny;
k= 1:1:Nz;
 
 plot3([p1x p2x], [p1y p2y], [ p1z p2z],'r');
 grid on 
% draw the grid

    X_start = X1;
    X_end = X1 + (Nx-1) * dx;
    
    Y_start = Y1;
    Y_end = Y1 + (Ny-1) * dy;
    
    Z_start = Z1;
    Z_end  = Z1 + (Nz-1) * dz;
    
    

    %figure;
    %hold on;
    for i = 1:Nx
        for j = 1:Ny
            for k=1:Nz
         
                X(i) = X1 + (i-1)*dx; 
                Y(j) = Y1 + (j-1)*dy;
                Z(k) = Z1 + (k-1)*dz;
        
         
            end
        end
    end
 
  
%  find alpha_min and alpha_max
 if (p2x-p1x ~= 0)
 
alpha_x1 = (X1 - p1x)/(p2x - p1x); % hit first x plane
alpha_xend = (X1 + (Nx-1)*dx- p1x) / (p2x - p1x); % hit last x plane
 else
     disp('alpha_x is undefined')
 end
 
  if (p2y-p1y ~= 0)
 
alpha_y1 = (Y1 - p1y)/(p2y - p1y); % hit first y plane
alpha_yend = (Y1 + (Ny-1)*dy - p1y) / (p2y - p1y); % hit last y plane
  else
      disp('alpha_y is undefined')
  end
  
  if (p2z-p1z ~= 0)
 
alpha_z1 = (Z1 - p1z)/(p2z - p1z); % hit first z plane
alpha_zend = (Z1 + (Nz-1)*dz - p1z) / (p2z - p1z); % hit last z plane
  else
      disp('alpha_z is undefined')
  end
  

% calculate range of parametric values 

alpha_xmin = min([alpha_x1, alpha_xend]); % the first x plane we hit
alpha_xmax = max([alpha_x1, alpha_xend]); % the last x plane we hit

alpha_ymin = min([alpha_y1, alpha_yend]); % the first y plane we hit
alpha_ymax = max([alpha_y1, alpha_yend]); % the last y plane we hit

alpha_zmin = min([alpha_z1, alpha_zend]); % the first z plane we hit
alpha_zmax = max([alpha_z1, alpha_zend]); % the last z plane we hit

alpha_min = max([0,alpha_xmin, alpha_ymin, alpha_zmin]); % the first plane we hit (x or y or z)
alpha_max = min([1,alpha_xmax, alpha_ymax, alpha_zmax]); % the lane plane we hit (x or y or z)

if (alpha_min < alpha_max) 
    disp('Intersecting voxel space...');
else 
    disp('The ray does not intersecting voxel space...');
end
% 
% calculate range of indicates

if ( p1x <= p2x ) 
    i_min=Nx -(X(Nx)-alpha_min*(p2x-p1x)-p1x)/dx;
else 
    i_min = Nx-(X(Nx)-alpha_max*(p2x-p1x)-p1x)/dx;
end

if ( p1x <= p2x) 
    i_max =1+(p1x+alpha_max*(p2x-p1x)-X(1))/dx;
else
    i_max = 1+(p1x+alpha_min*(p2x-p1x)-X(1))/dx;
end

if (p1y <= p2y) 
    j_min = Ny-(X(Nx)-alpha_min*(p2y-p1y)-p1y)/dy;
else 
    j_min = Ny-(X(Nx)-alpha_max*(p2y-p1y)-p1y)/dy;
end

if (p1y <= p2y) 
    j_max = 1+(p1y+alpha_max*(p2y-p1y)-Y(1))/dy;
else
    j_max = 1+(p1y+alpha_min*(p2y-p1y)-Y(1))/dy;
end

if (p1z <= p2z)
    k_min = Nz-(Z(Nz)-alpha_min*(p2z-p1z)-p1z)/dz;
else 
    k_min = Nz-(Z(Nz)-alpha_max*(p2z-p1z)-p1z)/dz;
end

if (p1z <= p2z) 
    k_max = 1+(p1z+alpha_max*(p2z-p1z)-Z(1))/dz;
else
    k_max = 1+(p1z+alpha_min*(p2z-p1z)-Z(1))/dz;
end


for i = 1:Nx
        for j = 1:Ny
            for k=1:Nz
     
                alphax(i) = (X(i)-p1x)/(p2x-p1x);
                alphay(j) = (Y(j)-p1y)/(p2y-p1y);
                alphaz(k)= (Z(k)-p1z)/(p2z-p1z);  
            end
        end
end
 
alpha_x=sort(alphax);
alpha_y=sort(alphay);
alpha_z=sort(alphaz);
  
  Alphaset=[alpha_min,alpha_x ,alpha_y,alpha_z,alpha_max];
  alphaSET= sort(Alphaset);
  

 
  
  n = floor((i_max - i_min + 1) + (j_max - j_min + 1)+(k_max -k_min+1)+1);
%   m = 1:n+1;
  
 for m = 1:n
     
     Alphaset=[alpha_min,alpha_x ,alpha_y,alpha_z,alpha_max];
     alphaSET= sort(unique(Alphaset));
    
 end 

d12=sqrt((p2x-p1x)^2+(p2y-p1y)^2+(p2z-p1z)^2);

for m=1:1:n
    
l(m)=d12*(alphaSET(m+1)-alphaSET(m));
alpha_mid(m)= (alphaSET(m+1)+alphaSET(m))/2;

i(m)= floor(1+ (p1x +alpha_mid(m)*(p2x-p1x)-X1)/dx);
j(m)=floor(1+ (p1y +alpha_mid(m)*(p2y-p1y)-Y1)/dy);
k(m)=floor(1+ (p1z +alpha_mid(m)*(p2z-p1z)-Z1)/dz);


% i(m)= (1+ (p1x +alpha_mid(m)*(p2x-p1x)-p1x)/dx);
% j(m)=(1+ (p1y +alpha_mid(m)*(p2y-p1y)-p1y)/dy);
% k(m)=(1+ (p1z +alpha_mid(m)*(p2z-p1z)-p1z)/dz);
% 

disp(['(i): ', num2str(i(m)),'  ','(j): ', num2str(j(m)) , '  ', '(k): ', num2str(k(m)), '  ', 'VoxelLenght: ', num2str(l(m))]);
 
end


  
  
  