
%{
    Derivadas Parciales X,Y y Laplaciano
    Ricardo Lopez R. A01066515
    09.22.2019

%}

%main()
function partialDerivates()
    close all
    h = .05;
    x = -2:h:2;
    
    [X,Y] = meshgrid(x);
    F = exp(-X.^2-Y.^2);
    dimension = size(F);
    dfdx = zeros(dimension(1),dimension(2));
    dfdy = zeros(dimension(1),dimension(2));
    d2fdx = zeros(dimension(1),dimension(2));
    d2fdy = zeros(dimension(1),dimension(2));
    
    figure(1);
    surf(X,Y,F);
    zlabel('F(x,y)');
    xlabel('x') 
    ylabel('y')
       
    partialDerivateX(F,dimension,dfdx,h,X,Y);
    partialDerivateY(F,dimension,dfdy,h,X,Y);
    
    laplacian = partialSeconDerivateX(F,dimension,d2fdx,h,X,Y) + ...
        partialSecondDerivateY(F,dimension,d2fdy,h,X,Y);
    
    figure(6);
    surf(X,Y,laplacian);
    zlabel('Laplacian');
    xlabel('x') 
    ylabel('y')
    
end

function partialDerivateX (F,dimension,dfdx,h,X,Y)

    %Partial Derivate with respect to X
    for i = 1:dimension(1)
        for j = 2:dimension(2)-1
            dfdx(i,j) = (F(i,j+1) - F(i,j-1))/(2*h);            
        end 
        
    end
    
    figure(2)
    surf(X,Y,dfdx);
    zlabel('dFdx');
    xlabel('x');
    ylabel('y');
end

function partialDerivateY (F,dimension,dfdy,h,X,Y)

    %Partial Derivate with respect to y
    for i = 1:dimension(1)
        for j = 2:dimension(2)-1
            dfdy(j,i) = (F(j+1,i) - F(j-1,i))/(2*h);            
        end         
    end
    
    figure(3)
    surf(X,Y,dfdy);
    zlabel('dFdy');
    xlabel('x');
    ylabel('y');    
    
end

function [secondDerivateX] = partialSeconDerivateX(F,dimension,d2fdx,h,X,Y)
    %%Second derivate with respect to X
    for i = 1:dimension(1)
        for j = 2:dimension(2)-1
            d2fdx(i,j) = (F(i,j+1) -2*F(i,j) + F(i,j-1))/(h^2);            
        end 
        
    end
    
    secondDerivateX = d2fdx;
    figure(4)
    surf(X,Y,d2fdx);
    zlabel('d2Fdx2');
    xlabel('x');
    ylabel('y');
    
end

function [secondDerivateY] = partialSecondDerivateY(F,dimension,d2fdy,h,X,Y)
    %Second derivate with respect to y
    for i = 1:dimension(1)
        for j = 2:dimension(2)-1
            d2fdy(j,i) = (F(j+1,i)-2*F(j,i) + F(j-1,i))/(h^2);            
        end         
    end 
    
    secondDerivateY = d2fdy;
    figure(5)
    surf(X,Y,d2fdy);
    zlabel('d2Fdy2');
    xlabel('x');
    ylabel('y');   
    
end
