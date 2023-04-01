
fid = fopen('keConstants_new.fut','w');

l = 1;
%kepre = Ke3D();
kepre = Ke3DNum();

k_old = getKEPreIntegration(2);
k_test = getKEPreIntegration_shift(2);

fprintf(fid,'type localMatrix = {xx: f64, xy: f64, xz: f64, yx: f64, yy: f64, yz: f64, zx: f64, zy: f64, zz: f64}\n\n');
fprintf(fid,'let getke_l0(recv :i32, send :i32)  :localMatrix = \n');
fprintf(fid,'match (recv,send)\n');
for i=1:8
    for j=1:8
        fprintf(fid,'case (%d,%d) -> ', i-1,j-1);
        kel = kepre(3*i-2:3*i,3*j-2:3*j);
        

        fprintf(fid,'{xx=(%25.25f),xy=(%25.25f),xz=(%25.25f),yx=(%25.25f),yy=(%25.25f),yz=(%25.25f),zx=(%25.25f),zy=(%25.25f),zz=(%25.25f)}\n',...
            kel(1,1), kel(1,2), kel(1,3), ...
            kel(2,1), kel(2,2), kel(2,3), ...
            kel(3,1), kel(3,2), kel(3,3));
    end
end
fprintf(fid,'case _ -> {xx=0,xy=0,xz=0,yx=0,yy=0,yz=0,zx=0,zy=0,zz=0}\n');



fprintf(fid,'\n\nlet keconst :[24][24]f64 = \n');
fprintf(fid,'[');
for i=1:24
    fprintf(fid,'[');
    for j=1:24
        fprintf(fid,'%25.25f',kepre(i,j));
        if (j~=24)
            fprintf(fid,',');
        end
    end
    fprintf(fid,']');
    if (i~=24)
        fprintf(fid,',');
    end
    fprintf(fid,'\n');
    
end
fprintf(fid,']\n');

fclose(fid);
fid = fopen('keLevelConstants.fut','w');

for l = 1:4
    kepre = getKEPreIntegration_shift(l);
    nl = pow2(l-1);
    
    fprintf(fid,'\n\n\n');
    fprintf(fid,'let keslices_l%d :[%d][%d][%d][8][3][24]f64 = [\n',l-1,nl,nl,nl);
    for ic=1:nl
        fprintf(fid,'[');
        for jc=1:nl
            fprintf(fid,'[');
            for kc=1:nl
                fprintf(fid,'[');
                for n=1:8
                    fprintf(fid,'[');
                    for s=1:3
                        fprintf(fid,'[');
                        i = 3*(n-1) + s;
                        for j=1:24
                            fprintf(fid,'%25.25f',kepre{ic,jc,kc}(i,j));
                            if (j~=24)
                                fprintf(fid,',');
                            end
                        end
                        fprintf(fid,']');
                        if (s~=3)
                            fprintf(fid,',');
                        end
                    end
                    fprintf(fid,']');
                    if (n~=8)
                        fprintf(fid,',');
                    end
                end
                if (kc~=nl)
                    fprintf(fid,'],\n');
                else
                    fprintf(fid,']');
                end
            end
            if (jc~=nl)
                fprintf(fid,'],\n');
            else
                fprintf(fid,']');
            end
        end
        if (ic~=nl)
            fprintf(fid,'],\n');
        else
            fprintf(fid,']');
        end
    end
    
    fprintf(fid,']\n');
end


fclose(fid);






function printFutVec(a)
fprintf('[');
for i=1:size(a)


        fprintf('%f15f32',a(i));
       
        if (i~=size(a))
        fprintf(',');
        end
end
fprintf(']');
end


function KE = Ke3D()
nu = 0.3;
C = [2/9 1/18 1/24 1/36 1/48 5/72 1/3 1/6 1/12];
A11 = [-C(1) -C(3) -C(3) C(2) C(3) C(3); -C(3) -C(1) -C(3) -C(3) -C(4) -C(5);...
    -C(3) -C(3) -C(1) -C(3) -C(5) -C(4); C(2) -C(3) -C(3) -C(1) C(3) C(3);...
    C(3) -C(4) -C(5) C(3) -C(1) -C(3); C(3) -C(5) -C(4) C(3) -C(3) -C(1)];
B11 = [C(7) 0 0 0 -C(8) -C(8); 0 C(7) 0 C(8) 0 0; 0 0 C(7) C(8) 0 0;...
    0 C(8) C(8) C(7) 0 0; -C(8) 0 0 0 C(7) 0; -C(8) 0 0 0 0 C(7)];
A22 = [-C(1) -C(3) C(3) C(2) C(3) -C(3); -C(3) -C(1) C(3) -C(3) -C(4) C(5);...
    C(3) C(3) -C(1) C(3) C(5) -C(4); C(2) -C(3) C(3) -C(1) C(3) -C(3);...
    C(3) -C(4) C(5) C(3) -C(1) C(3); -C(3) C(5) -C(4) -C(3) C(3) -C(1)];
B22 = [C(7) 0 0 0 -C(8) C(8); 0 C(7) 0 C(8) 0 0; 0 0 C(7) -C(8) 0 0;...
    0 C(8) -C(8) C(7) 0 0; -C(8) 0 0 0 C(7) 0; C(8) 0 0 0 0 C(7)];
A12 = [C(6) C(3) C(5) -C(4) -C(3) -C(5); C(3) C(6) C(5) C(3) C(2) C(3);...
    -C(5) -C(5) C(4) -C(5) -C(3) -C(4); -C(4) C(3) C(5) C(6) -C(3) -C(5);...
    -C(3) C(2) C(3) -C(3) C(6) C(5); C(5) -C(3) -C(4) C(5) -C(5) C(4)];
B12 = [-C(9) 0 -C(9) 0 C(8) 0; 0 -C(9) -C(9) -C(8) 0 -C(8); C(9) C(9) -C(9) 0 C(8) 0;...
    0 -C(8) 0 -C(9) 0 C(9); C(8) 0 -C(8) 0 -C(9) -C(9); 0 C(8) 0 -C(9) C(9) -C(9)];
A13 = [-C(4) -C(5) -C(3) C(6) C(5) C(3); -C(5) -C(4) -C(3) -C(5) C(4) -C(5);...
    C(3) C(3) C(2) C(3) C(5) C(6); C(6) -C(5) -C(3) -C(4) C(5) C(3);...
    C(5) C(4) -C(5) C(5) -C(4) -C(3); -C(3) C(5) C(6) -C(3) C(3) C(2)];
B13 = [0 0 C(8) -C(9) -C(9) 0; 0 0 C(8) C(9) -C(9) C(9); -C(8) -C(8) 0 0 -C(9) -C(9);...
    -C(9) C(9) 0 0 0 -C(8); -C(9) -C(9) C(9) 0 0 C(8); 0 -C(9) -C(9) C(8) -C(8) 0];
A14 = [C(2) C(5) C(5) C(4) -C(5) -C(5); C(5) C(2) C(5) C(5) C(6) C(3);...
    C(5) C(5) C(2) C(5) C(3) C(6); C(4) C(5) C(5) C(2) -C(5) -C(5);...
    -C(5) C(6) C(3) -C(5) C(2) C(5); -C(5) C(3) C(6) -C(5) C(5) C(2)];
B14 = [-C(9) 0 0 -C(9) C(9) C(9); 0 -C(9) 0 -C(9) -C(9) 0; 0 0 -C(9) -C(9) 0 -C(9);...
    -C(9) -C(9) -C(9) -C(9) 0 0; C(9) -C(9) 0 0 -C(9) 0; C(9) 0 -C(9) 0 0 -C(9)];
A23 = [C(2) C(5) -C(5) C(4) -C(5) C(5); C(5) C(2) -C(5) C(5) C(6) -C(3);...
    -C(5) -C(5) C(2) -C(5) -C(3) C(6); C(4) C(5) -C(5) C(2) -C(5) C(5);...
    -C(5) C(6) -C(3) -C(5) C(2) -C(5); C(5) -C(3) C(6) C(5) -C(5) C(2)];
B23 = [-C(9) 0 0 -C(9) C(9) -C(9); 0 -C(9) 0 -C(9) -C(9) 0; 0 0 -C(9) C(9) 0 -C(9);...
    -C(9) -C(9) C(9) -C(9) 0 0; C(9) -C(9) 0 0 -C(9) 0; -C(9) 0 -C(9) 0 0 -C(9)];
KE = 1/(1+nu)/(2*nu-1)*([A11 A12 A13 A14; A12' A22 A23 A13'; A13' A23' A22 A12'; A14' A13 A12 A11] +...
    nu*[B11 B12 B13 B14; B12' B22 B23 B13'; B13' B23' B22 B12'; B14' B13 B12 B11]);
end

%% FUNCTION Ke3D - ELEMENT STIFFNESS MATRIX
function KE = Ke3DNum()

nu = 0.3;
a = 0.5;
b = a;
c = a;

xx = [-a  b -c ...
       a  b -c ...
       a -b -c ...
      -a -b -c ...
      -a  b  c ...
       a  b  c ...
       a -b  c ...
      -a -b  c];


xpts = [-1/sqrt(3), 1/sqrt(3)];
ypts = [-1/sqrt(3), 1/sqrt(3)];
zpts = [-1/sqrt(3), 1/sqrt(3)];

C = getC(nu);
KE = zeros(24);

for xi = xpts
    for eta = ypts
        for zeta = zpts
         
            [B,jdet] = getB([xi, eta, zeta], xx);

            KE = KE +jdet*(B'*C*B);
        end
    end
end
end


function [C] = getC(nu)
temp1 = (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
temp2 = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
temp3 = 1.0 / (2.0 * (1.0 + nu));

C = zeros(6);
C(1, 1) = temp1;
C(2, 2) = temp1;
C(3, 3) = temp1;
C(4, 4) = temp3;
C(5, 5) = temp3;
C(6, 6) = temp3;
C(1, 2) = temp2;
C(2, 1) = temp2;
C(1, 3) = temp2;
C(3, 1) = temp2;
C(2, 3) = temp2;
C(3, 2) = temp2;
end

function [B,jdet] = getB(iso,xe)

    xi = iso(1);
    eta = iso(2);
    zeta = iso(3);
    
    n1xi   = -0.125 * (1 + eta) * (1 - zeta);
    n1eta  =  0.125 * (1 - xi) * (1 - zeta);
    n1zeta = -0.125 * (1 - xi) * (1 + eta);
    n2xi   =  0.125 * (1 + eta) * (1 - zeta);
    n2eta  =  0.125 * (1 + xi) * (1 - zeta);
    n2zeta = -0.125 * (1 + xi) * (1 + eta);

    n3xi   =  0.125 * (1 - eta) * (1 - zeta);
    n3eta  = -0.125 * (1 + xi) * (1 - zeta);
    n3zeta = -0.125 * (1 + xi) * (1 - eta);
    n4xi   = -0.125 * (1 - eta) * (1 - zeta);
    n4eta  = -0.125 * (1 - xi) * (1 - zeta);
    n4zeta = -0.125 * (1 - xi) * (1 - eta);

    n5xi   = -0.125 * (1 + eta) * (1 + zeta);
    n5eta  =  0.125 * (1 - xi) * (1 + zeta);
    n5zeta =  0.125 * (1 - xi) * (1 + eta);
    n6xi   =  0.125 * (1 + eta) * (1 + zeta);
    n6eta  =  0.125 * (1 + xi) * (1 + zeta);
    n6zeta =  0.125 * (1 + xi) * (1 + eta);

    n7xi   =  0.125 * (1 - eta) * (1 + zeta);
    n7eta  = -0.125 * (1 + xi) * (1 + zeta);
    n7zeta =  0.125 * (1 + xi) * (1 - eta);
    n8xi   = -0.125 * (1 - eta) * (1 + zeta);
    n8eta  = -0.125 * (1 - xi) * (1 + zeta);
    n8zeta =  0.125 * (1 - xi) * (1 - eta);

    L = zeros(6,9);
    jac = zeros(3);
    jacinvt = zeros(9);
    Nt = zeros(9,24);
    
    L(1, 1) = 1.0;
    L(2, 5) = 1.0;
    L(3, 9) = 1.0;
    L(4, 2) = 1.0;
    L(4, 4) = 1.0;
    L(5, 6) = 1.0;
    L(5, 8) = 1.0;
    L(6, 3) = 1.0;
    L(6, 7) = 1.0;

    Nt(1, 1)  = n1xi;
    Nt(2, 1)  = n1eta;
    Nt(3, 1)  = n1zeta;
    Nt(1, 4)  = n2xi;
    Nt(2, 4)  = n2eta;
    Nt(3, 4)  = n2zeta;
    Nt(1, 7)  = n3xi;
    Nt(2, 7)  = n3eta;
    Nt(3, 7)  = n3zeta;
    Nt(1, 10)  = n4xi;
    Nt(2, 10)  = n4eta;
    Nt(3, 10)  = n4zeta;
    Nt(1, 13) = n5xi;
    Nt(2, 13) = n5eta;
    Nt(3, 13) = n5zeta;
    Nt(1, 16) = n6xi;
    Nt(2, 16) = n6eta;
    Nt(3, 16) = n6zeta;
    Nt(1, 19) = n7xi;
    Nt(2, 19) = n7eta;
    Nt(3, 19) = n7zeta;
    Nt(1, 22) = n8xi;
    Nt(2, 22) = n8eta;
    Nt(3, 22) = n8zeta;

    Nt(4, 2)  = n1xi;
    Nt(5, 2)  = n1eta;
    Nt(6, 2)  = n1zeta;
    Nt(4, 5)  = n2xi;
    Nt(5, 5)  = n2eta;
    Nt(6, 5)  = n2zeta;
    Nt(4, 8)  = n3xi;
    Nt(5, 8)  = n3eta;
    Nt(6, 8)  = n3zeta;
    Nt(4, 11) = n4xi;
    Nt(5, 11) = n4eta;
    Nt(6, 11) = n4zeta;
    Nt(4, 14) = n5xi;
    Nt(5, 14) = n5eta;
    Nt(6, 14) = n5zeta;
    Nt(4, 17) = n6xi;
    Nt(5, 17) = n6eta;
    Nt(6, 17) = n6zeta;
    Nt(4, 20) = n7xi;
    Nt(5, 20) = n7eta;
    Nt(6, 20) = n7zeta;
    Nt(4, 23) = n8xi;
    Nt(5, 23) = n8eta;
    Nt(6, 23) = n8zeta;

    Nt(7, 3)  = n1xi;
    Nt(8, 3)  = n1eta;
    Nt(9, 3)  = n1zeta;
    Nt(7, 6)  = n2xi;
    Nt(8, 6)  = n2eta;
    Nt(9, 6)  = n2zeta;
    Nt(7, 9)  = n3xi;
    Nt(8, 9)  = n3eta;
    Nt(9, 9)  = n3zeta;
    Nt(7, 12) = n4xi;
    Nt(8, 12) = n4eta;
    Nt(9, 12) = n4zeta;
    Nt(7, 15) = n5xi;
    Nt(8, 15) = n5eta;
    Nt(9, 15) = n5zeta;
    Nt(7, 18) = n6xi;
    Nt(8, 18) = n6eta;
    Nt(9, 18) = n6zeta;
    Nt(7, 21) = n7xi;
    Nt(8, 21) = n7eta;
    Nt(9, 21) = n7zeta;
    Nt(7, 24) = n8xi;
    Nt(8, 24) = n8eta;
    Nt(9, 24) = n8zeta;


    jac(1, 1) = n1xi * xe(1) + n2xi * xe(4) + n3xi * xe(7) + n4xi * xe(10) +...
        n5xi * xe(13) + n6xi * xe(16) + n7xi * xe(19) + n8xi * xe(22);
    jac(2, 1) = n1eta * xe(1) + n2eta * xe(4) + n3eta * xe(7) + n4eta * xe(10) +...
        n5eta * xe(13) + n6eta * xe(16) + n7eta * xe(19) + n8eta * xe(22);
    jac(3, 1) = n1zeta * xe(1) + n2zeta * xe(4) + n3zeta * xe(7) + n4zeta * xe(10) + n5zeta * xe(13) + n6zeta * xe(16) +n7zeta * xe(19) + n8zeta * xe(22);

    jac(1, 2) = n1xi * xe(2) + n2xi * xe(5) + n3xi * xe(8) + n4xi * xe(11) + n5xi * xe(14) + n6xi * xe(17) +n7xi * xe(20) + n8xi * xe(23);
    jac(2, 2) = n1eta * xe(2) + n2eta * xe(5) + n3eta * xe(8) + n4eta * xe(11) + n5eta * xe(14) + n6eta * xe(17) +n7eta * xe(20) + n8eta * xe(23);
    jac(3, 2) = n1zeta * xe(2) + n2zeta * xe(5) + n3zeta * xe(8) + n4zeta * xe(11) + n5zeta * xe(14) + n6zeta * xe(17) +n7zeta * xe(20) + n8zeta * xe(23);

    jac(1, 3) = n1xi * xe(3) + n2xi * xe(6) + n3xi * xe(9) + n4xi * xe(12) + n5xi * xe(15) + n6xi * xe(18) +n7xi * xe(21) + n8xi * xe(24);
    jac(2, 3) = n1eta * xe(3) + n2eta * xe(6) + n3eta * xe(9) + n4eta * xe(12) + n5eta * xe(15) + n6eta * xe(18) +n7eta * xe(21) + n8eta * xe(24);
    jac(3, 3) = n1zeta * xe(3) + n2zeta * xe(6) + n3zeta * xe(9) + n4zeta * xe(12) + n5zeta * xe(15) + n6zeta * xe(18) +n7zeta * xe(21) + n8zeta * xe(24);

    jdet = det(jac);
    ijac = inv(jac);
    
    jacinvt(1:3,1:3) = ijac;
    jacinvt(4:6,4:6) = ijac;
    jacinvt(7:9,7:9) = ijac;

    B = (L * jacinvt * Nt);
end

%% FUNCTION getKEPreIntegration - preintegrate KE for cell structure
function [KEpre] = getKEPreIntegration(l)

nu = 0.3;
ncell = 2^(l-1);
int_points = 25;
C = getC(nu);

a = 0.5;
b = a;
c = a;

xx = [-a  b -c ...
       a  b -c ...
       a -b -c ...
      -a -b -c ...
      -a  b  c ...
       a  b  c ...
       a -b  c ...
      -a -b  c];

spacing = 2/ncell/int_points;
subCellVolume = spacing*spacing*spacing;

%pre- integrate matrices, to speed up product.
KEpre = cell(ncell,ncell,ncell);
for ii = 1:ncell
    for kk = 1:ncell
        for jj = 1:ncell
            KEpre{ii,jj,kk} = zeros(24);
            
            starti = -1 + spacing/2 + 2/ncell*(ii-1);
            endi   = 1 - spacing/2 - 2/ncell*(ncell-ii);
            ipts = starti:spacing:endi;
            
            startj = -1 + spacing/2 + 2/ncell*(jj-1);
            endj   = 1 - spacing/2 - 2/ncell*(ncell-jj);
            jpts = startj:spacing:endj;
            
            startk = -1 + spacing/2 + 2/ncell*(kk-1);
            endk   = 1 - spacing/2 - 2/ncell*(ncell-kk);
            kpts = startk:spacing:endk;
            
            for xi = ipts
                for eta = jpts
                    for zeta = kpts
                        [B,jdet] = getB([xi,eta,zeta],xx);
                        KEpre{ii,jj,kk} = KEpre{ii,jj,kk} + (jdet * subCellVolume) * (B'*C*B);
                    end
                end
            end
        end
    end
end
end

%% FUNCTION getKEPreIntegration - preintegrate KE for cell structure
function [KEpre] = getKEPreIntegration_shift(l)

nu = 0.3;
ncell = 2^(l-1);
C = getC(nu);

a = 0.5 * ncell;
b = a;
c = a;

xx = [-a  b -c ...
       a  b -c ...
       a -b -c ...
      -a -b -c ...
      -a  b  c ...
       a  b  c ...
       a -b  c ...
      -a -b  c];

%pre- integrate matrices, to speed up product.
KEpre = cell(ncell,ncell,ncell);
for ii = 1:ncell
    for kk = 1:ncell
        for jj = 1:ncell
            KEpre{ii,jj,kk} = zeros(24);
            
            ax = -1 + 2/ncell*(ii-1);
            bx = -1 + 2/ncell*(ii);
            weighti = (bx-ax) / 2;
            
            ay = -1 + 2/ncell*(jj-1);
            by = 1 - 2/ncell*(ncell-jj);
            weightj = (by-ay) / 2;
            
            az = -1 + 2/ncell*(kk-1);
            bz   = 1 - 2/ncell*(ncell-kk);
            weightk = (bz-az) / 2;
            
            xpts = [-1/sqrt(3), 1/sqrt(3)];
            ypts = [-1/sqrt(3), 1/sqrt(3)];
            zpts = [-1/sqrt(3), 1/sqrt(3)];

            for xi = xpts
                for eta = ypts
                    for zeta = zpts
                        xi_new = weighti * xi + ((ax+bx)/2);
                        eta_new = weighti * eta + ((ay+by)/2);
                        zeta_new = weighti * zeta + ((az+bz)/2);
                        
                        [B,jdet] = getB([xi_new, eta_new, zeta_new], xx);
                        KEpre{ii,jj,kk} = KEpre{ii,jj,kk} +(jdet*weighti*weightj*weightk)*(B'*C*B);
                    end
                end
            end
        end
    end
end
end