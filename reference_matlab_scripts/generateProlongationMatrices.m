


% -- negative z
% --      x
% --      +  -
% -- y +  1  0
% --   -  2  3
% 
% -- positive z
% --      x
% --      +  -
% -- y +  5  4
% --   -  6  7
% -  [[1, 0.5, 0.25, 0.5, 0.5, 0.25, 0.125, 0.25],
% --   [0.5, 1, 0.5, 0.25, 0.25, 0.5, 0.25, 0.125],
% --   [0.25, 0.5, 1, 0.5, 0.125, 0.25, 0.5, 0.25],
% --   [0.5, 0.25, 0.5, 1, 0.25, 0.125, 0.25, 0.5],
% --   [0.5, 0.25, 0.125, 0.25, 1, 0.5, 0.25, 0.5],
% --   [0.25, 0.5, 0.25, 0.125, 0.5, 1, 0.5, 0.25],
% --   [0.125, 0.25, 0.5, 0.25, 0.25, 0.5, 1, 0.5],
% --   [0.25, 0.125, 0.25, 0.5, 0.5, 0.25, 0.5, 1]]

w = zeros(8,8,8);
sz = 1;

fprintf('let prolongationWeightsSmall :[8][8][8]f64 =[\n');

for anchorNode = 0:7 % which cell are we looking at?
    c_anchor = getCoords(anchorNode);
    
    for sendingNode = 0:7 %which node is sending
        
        csend = 2.* getCoords(sendingNode);
        
        
        for recvNode = 0:7 %which node is recieving
            co = c_anchor + getCoords(recvNode);
            d=getDist(csend,co);
            
            startR = (sz*recvNode)+1;
            startS = (sz*sendingNode)+1;
            
            endR = sz*(recvNode+1);
            endS = sz*(sendingNode+1);
        
            
            w(startR:endR,startS:endS,anchorNode+1) = diag(repmat(getWeigth(d),sz,1));
            
        end
    end
    
    
    printFutMat(w(:,:,anchorNode+1));
    if (anchorNode~=7)
        fprintf(',\n');
    end
end
fprintf(']\n');



function [coord] = getCoords(li)
switch li
    case 0
        coord = [0,1,0];
    case 1
        coord = [1,1,0];
    case 2
        coord = [1,0,0];
    case 3
        coord = [0,0,0];
    case 4
        coord = [0,1,1];
    case 5
        coord = [1,1,1];
    case 6
        coord = [1,0,1];
    case 7
        coord = [0,0,1];
    otherwise
        coord = [-3,-3,-3]
end
end

function printFutMat(m)
fprintf('[');
for i=1:size(m,1)
    fprintf('[');
    for j=1:size(m,2)
        
        fprintf('%g',m(i,j));
        
        if (j~=size(m,2))
            fprintf(',');
        end
    end
    if (i~=size(m,1))
        fprintf('],\n');
    else
        fprintf(']');
    end
    
end
fprintf(']');
end

function [d] = getDist(c1,c2)
c = c1 - c2;
d = (c(1)).^2 + (c(2)).^2 + (c(3)).^2;
end

function [w] = getWeigth(d)
switch d
    case 0
        w=1;
    case 1
        w=0.5;
    case 2
        w=0.25;
    case 3
        w=0.125;
    otherwise
        w = 0;  
end
end
