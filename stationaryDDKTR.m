rx=0;
ry=0;
rz=-4;

ref_point=[rx ry rz];
plot3(rx,ry,rz,'*');
hold on
grid on

Nx=4;
Ny=4;
Nz=6;

X_start=-Nx/2;
Y_start=-Ny/2;
Z_start=-Nz/2;

Xend=Nx/2;
Yend=Ny/2;
Zend=Nz/2;


P1x=0;
P1y=0;
P1z=100;
P1=[0 0 100];

Dx=0;
Dy=0;
Dz=-4;
Ddktr=[Dx Dy Dz];

plot3([P1x Dx],[P1y Dy],[P1z Dz],'r'); %centre line
hold on; 
grid on;



for theta=-20:5:20
   for a =(-Nx/2-1 ):1:(Nx/2)
    for b =(-Ny/2-1 ):1:(Ny/2)
        
       
        P2x=a+0.5;
        P2y=b+0.5;
        P2z=Z_start-1;
        
        P2 =[P2x  P2y P2z]; 
       
        start_p=[P1x-rx P1y-ry P1z-rz 1] ;
        end_p=[P2x-rx P2y-ry P2z-rz 1];
        
        start_pX=[P1x-rx P1y-ry P1z-rz 1] ;%for rotation  centre red line
        end_pX=[Dx-rx Dy-ry Dz-rz 1];%for rotation  centre red line
        


    RotX=[1 0 0 0; 0 cosd(theta) sind(theta) 0; 0 -sind(theta) cosd(theta) 0; 0 0 0 1];
    RotY=[cosd(theta) 0 sind(theta) 0; 0 1 0 0; -sind(theta) 0 cosd(theta) 0; 0 0 0 1];
    RotZ=[cosd(theta) -sind(theta) 0 0; sind(theta) cosd(theta) 0 0; 0 0 1 0; 0  0 0 1];
    
         %X EKSENI
%         
          start_p2= start_p*RotX;
           end_p2= end_p*RotX;
%          end_p2= end_p;
          
          start_p2x= start_pX*RotX;%for rotation  centre red line 
          end_p2x= end_pX*RotX;%for rotation centre red line
%               end_p2x= end_pX;
       
            
              %Y EKSENI
        
%           start_p2= start_p*RotY;
%           end_p2= end_p*RotZ;

              %Z EKSENI
        
%           start_p2= start_p*RotZ;
%           end_p2= end_p*RotZ;


   plot3([start_p2(1)+rx end_p2(1)+rx],[start_p2(2)+ry end_p2(2)+ry],[start_p2(3)+rz end_p2(3)+rz],'b');
   x=sqrt((end_p2(1)-start_p2(1))^2+(end_p2(2)-start_p2(2))^2+(end_p2(3)-start_p2(3))^2);  
   hold on
   grid on
        
   plot3([start_p2x(1)+rx end_p2x(1)+rx],[start_p2x(2)+ry end_p2x(2)+ry],[start_p2x(3)+rz end_p2x(3)+rz],'r');     
   hold on; 
   grid on;
   disp(P2)     
   
        
        
    end 
  end
end



