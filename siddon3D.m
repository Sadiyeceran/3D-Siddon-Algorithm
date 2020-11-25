 
% the points
p1x =0 ;
p1y =0 ;
p1z =10;

p2x = 0;
p2y = 0;
p2z = -3;

Dx=0;
Dy=0;
Dz=-3;
Ddktr=[Dx Dy Dz];

P1 = [p1x p1y p1z ];
P2 = [p2x p2y p2z];

% the number of planes
Nx =4 ; %boy
Ny = 4;%geniþlik
Nz = 4  ;%yükseklik

% separation between planes
%1 cm =1 br
dx = 1;
dy = 1;
dz = 1;
d=1;

% position of first plane (in x or y or z)
X1 = Nx;
Y1 = Ny;
Z1 = Nz;
%KÜP

 plot3([X1 X1],[Y1 Ny/2],[Z1 Z1],'m'); hold on; grid on;
 plot3([X1 Nx/2],[Y1 Y1],[Z1 Z1],'m'); hold on; grid on;
 plot3([X1 Nx/2],[Ny/2 Ny/2],[Z1 Z1],'m'); hold on; grid on;
 plot3([Nx/2 Nx/2],[Y1 Ny/2],[Z1 Z1],'m'); hold on; grid on;
 %küpün yan yüzeyleri 
 plot3([X1 X1],[Y1 Y1],[Z1 Nz/2],'m'); hold on; grid on;
 plot3([X1 X1],[Ny/2 Ny/2],[Z1 Nz/2],'m'); hold on; grid on;
 plot3([Nx/2 Nx/2],[Ny/2 Ny/2],[Z1 Nz/2],'m'); hold on; grid on;
 plot3([Nx/2 Nx/2],[Y1 Y1],[Z1 Nz/2],'m'); hold on; grid on; 
 %küpün üst kýsmý
 plot3([X1 X1],[Y1 Ny/2],[Nz/2 Nz/2],'m'); hold on; grid on; 
 plot3([X1 Nx/2],[Y1 Y1],[Nz/2 Nz/2],'m'); hold on; grid on;
 plot3([X1 Nx/2],[Ny/2 Ny/2],[Nz/2 Nz/2],'m'); hold on; grid on;
 plot3([Nx/2 Nx/2],[Y1 Ny/2],[Nz/2 Nz/2],'m'); hold on; grid on;

% define i j k
i= 1:1:Nx;
j= 1:1:Ny;
k= 1:1:Nz;
 
 plot3([p1x p2x], [p1y p2y], [ p1z p2z],'r');
 grid on 
 
rx=0;
ry=0;
rz=1;

ref_point=[rx ry rz];
plot3(rx,ry,rz,'*');
hold on
grid on

for theta=0:30:30
   for a =(-Nx/2-1 ):1:(Nx/2)
    for b =(-Ny/2-1 ):1:(Ny/2)
        
        
        p2x=a+0.5;
        p2y=b+0.5;
        p2z=Z1-1;
        
        P2 =[p2x  p2y p2z]; 
        
        start_p=[p1x-rx p1y-ry p1z-rz 1] ;
        end_p=[p2x-rx p2y-ry p2z-rz 1];
         
         start_pX=[p1x-rx p1y-ry p1z-rz 1] ;%for rotation  centre red line
         end_pX=[Dx-rx Dy-ry Dz-rz 1];%for rotation  centre red line
        

    RotX=[1 0 0 0; 0 cosd(theta) sind(theta) 0; 0 -sind(theta) cosd(theta) 0; 0 0 0 1];
    RotY=[cosd(theta) 0 sind(theta) 0; 0 1 0 0; -sind(theta) 0 cosd(theta) 0; 0 0 0 1];
    RotZ=[cosd(theta) -sind(theta) 0 0; sind(theta) cosd(theta) 0 0; 0 0 1 0; 0  0 0 1];
    
%          %X EKSENI
% %         
          start_p2= start_p*RotX;
          end_p2= end_p*RotX;
          
          start_p2x= start_pX*RotX;%for rotation  centre red line
          end_p2x= end_pX*RotX;%for rotation  centre red line
            
            
              %Y EKSENI
        
%           start_p2= start_p*RotY;
%           end_p2= end_p*RotZ;
% 
%               Z EKSENI
%         
%           start_p2= start_p*RotZ;
%           end_p2= end_p*RotZ;
%                 
%           start_p2x= start_pX*RotZ;%for rotation  centre red line
%           end_p2x= end_pX*RotZ;%for rotation  centre red lin


  plot3([start_p2(1)+rx end_p2(1)+rx],[start_p2(2)+ry end_p2(2)+ry],[start_p2(3)+rz end_p2(3)+rz],'b');
%   x=sqrt((end_p2(1)-start_p2(1))^2+(end_p2(2)-start_p2(2))^2+(end_p2(3)-start_p2(3))^2) 
  hold on
  grid on
        
  plot3([start_p2x(1)+rx end_p2x(1)+rx],[start_p2x(2)+ry end_p2x(2)+ry],[start_p2x(3)+rz end_p2x(3)+rz],'r');
  hold on; 
  grid on;
  disp(P2)     

    X_start = X1;
    X_end = X1 + (Nx/2-1) * dx;
    
    Y_start = Y1;
    Y_end = Y1 + (Ny/2-1) * dy;
    
    Z_start = Z1;
    Z_end  = Z1 + (Nz/2-1) * dz;
    
    
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
alpha_xend = (X1 + (Nx/2-1)*dx- p1x) / (p2x - p1x); % hit last x plane
 else
     disp('alpha_x is undefined')
 end
 
  if (p2y-p1y ~= 0)
 
alpha_y1 = (Y1 - p1y)/(p2y - p1y); % hit first y plane
alpha_yend = (Y1 + (Ny/2-1)*dy - p1y) / (p2y - p1y); % hit last y plane
  else
      disp('alpha_y is undefined')
  end
  
  if (p2z-p1z ~= 0)
 
alpha_z1 = (Z1 - p1z)/(p2z - p1z); % hit first z plane
alpha_zend = (Z1 + (Nz/2-1)*dz - p1z) / (p2z - p1z); % hit last z plane
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
  
if ( p1x <= p2x ) 
    i_min=Nx/2 -(X(Nx/2)-alpha_min*(p2x-p1x)-p1x)/dx;
else 
    i_min = Nx/2-(X(Nx/2)-alpha_max*(p2x-p1x)-p1x)/dx;
end

if ( p1x <= p2x) 
    i_max =1+(p1x+alpha_max*(p2x-p1x)-X(1))/dx;
else
    i_max = 1+(p1x+alpha_min*(p2x-p1x)-X(1))/dx;
end

if (p1y <= p2y) 
    j_min = Ny/2-(Y(Ny/2)-alpha_min*(p2y-p1y)-p1y)/dy;
else 
    j_min = Ny/2-(Y(Ny/2)-alpha_max*(p2y-p1y)-p1y)/dy;
end

if (p1y <= p2y) 
    j_max = 1+(p1y+alpha_max*(p2y-p1y)-Y(1))/dy;
else
    j_max = 1+(p1y+alpha_min*(p2y-p1y)-Y(1))/dy;
end

if (p1z <= p2z)
    k_min = Nz/2-(Z(Nz/2)-alpha_min*(p2z-p1z)-p1z)/dz;
else 
    k_min = Nz/2-(Z(Nz/2)-alpha_max*(p2z-p1z)-p1z)/dz;
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

  n = floor((i_max - i_min + 1) + (j_max - j_min + 1)+(k_max -k_min+1)+1);
%   m = 1:n+1;
 if (alpha_max > alpha_min) 
      disp('Intersecting voxel space...');
    
 for m = 1:n
         
     Alphaset=[alpha_min,alpha_x ,alpha_y,alpha_z,alpha_max];
     alphaSET= sort((Alphaset));
     %alphaSET = sort(unique(Alphaset));
  
 end 
 
 
d12=sqrt((p2x-p1x)^2+(p2y-p1y)^2+(p2z-p1z)^2);
    
for m=1:n
    
l(m)=d12*(alphaSET(m+1)-alphaSET(m));
alpha_mid(m)= (alphaSET(m+1)+alphaSET(m))/2;

i(m)= floor(1+ (p1x +alpha_mid(m)*(p2x-p1x)-X(1)/dx));
j(m)=floor(1+ (p1y +alpha_mid(m)*(p2y-p1y)-Y(1)/dy));
k(m)=floor(1+ (p1z +alpha_mid(m)*(p2z-p1z)-Z(1)/dz));

% 
% i(m)= (1+ (p1x +alpha_mid(m)*(p2x-p1x)-p1x)/dx);
% j(m)=(1+ (p1y +alpha_mid(m)*(p2y-p1y)-p1y)/dy);
% k(m)=(1+ (p1z +alpha_mid(m)*(p2z-p1z)-p1z)/dz);



disp(['(i): ', num2str(i(m)),'  ','(j): ', num2str(j(m)) , '  ', '(k): ', num2str(k(m)), '  ', 'VoxelLenght: ', num2str(l(m))]);


end

 else
     disp(' The ray does not intersecting voxel space...');
 end
    end
   end  
end 


