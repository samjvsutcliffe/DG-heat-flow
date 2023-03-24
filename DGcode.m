function DGcode
% Create the mesh ---------------------------------------------------------
dx = 0.1;
[etpl,ftlp,nnls,nls,nfcs,coord] = makemesh(dx);                             % Make the intial mesh
nen=4;nov = (nen);neDoF = (nov)^2;                                          % number of nodes and degrees of freedom for an element

% Shape functions ---------------------------------------------------------
fNR = @(X)[(1-X(1)).*(1-X(2)),(1-X(1)).*(1+X(2)),(1+X(1)).*(1+X(2)),(1+X(1)).*(1-X(2))]./4;
fdNr = @(X)[X(2) - 1, - X(2) - 1, X(2) + 1, 1 - X(2); X(1) - 1, 1 - X(1), X(1) + 1, - X(1) - 1]./4;


% volume integral ---------------------------------------------------------
ed = zeros(nnls,size(etpl,2));F = zeros(nnls*nov,1);
krow = zeros(nnls*neDoF,1);kcol = krow; kval= krow;                         % dof matrix, and column storage variables
a  = 1/sqrt(3); gpv = [0]'; wp = 1; ngp=1;                % Gauss point locations
kloc = 0;
for nel = 1:nnls                                                            % Looping over all the elements
    ke = zeros(nov);                                                        % Local stiffness matrix
    fe = zeros(nov,1);
    for gp = 1:ngp                                                          % Gauss point loop
        dNr = fdNr(gpv(gp,:));
        JT = dNr*coord(etpl(nel,:),:);  detJ = det(JT);                     % Jocabian
        dNx = JT\dNr;                                                       % Shape function derivative
        ke = ke+dNx.'*dNx*detJ*wp;                                          % Local stiffness matrix summation
        x = fNR(gpv(gp,:))*coord(etpl(nel,:),:);
        fe = fe + 2*(pi^2)*fNR(gpv(gp,:))'*sin(x(1)*pi)*sin(x(2)*pi)*detJ*wp;
    end
    ed(nel,1:nov) = (1:nov) + max(ed(:));                                   % Creating degrees of freedom for current element
    F(ed(nel,1:nov))=fe;
    kloc = (1:neDoF) + max(kloc);                                           % Steering vector index
    krow(kloc) = reshape(ed(nel,1:nov).'*ones(1,nov),neDoF,1);              % Steering vector row values
    kcol(kloc) = reshape(ones(nov,1)*ed(nel,1:nov)  ,neDoF,1);              % Steering vector column values
    kval(kloc) = reshape(ke,neDoF,1);                                       % Vector of local stiffness values
end


% Surface integral --------------------------------------------------------
a  = sqrt(3/5); wp = [5/9 8/9 5/9]; ngp=3;                                  % Gauss weights
gpf(1).gp = [-1 -1 -1;-a 0  a]'; gpf(2).gp = [-a 0 a  ; 1 1 1]';            % Gauss values on faces 1,2,3 and 4
gpf(3).gp = [1 1 1   ;a 0  -a]';  gpf(4).gp = [a 0 -a  ;-1 -1 -1]';         % 
nm = [-1 0;0 1;1 0;0 -1];kloc = 0;tfdof = nfcs*neDoF;                       % normal for the faces and the total number of face DOF
k.v1 = zeros(tfdof,1); k.v2=k.v1; k.v3 = k.v1; k.v4=k.v1; k.r1 = k.v1; k.c1=k.v1; % Zeroed vectors for row and column steering matrix:
k.r2 = k.v1; k.c2=k.v1; k.r3 = k.v1; k.c3=k.v1; k.r4 = k.v1; k.c4=k.v1;     %
for fn = 1:nfcs
    co_pve = coord(etpl(ftlp(fn,1),:),:);                            	    % Positive element coordinates
    co_nve = coord(etpl(ftlp(fn,2),:),:);                            	    % Negative element coordinates
    ke_1 = zeros(nov);ke_2 = ke_1;ke_3 = ke_1;ke_4 = ke_1;                  % Assigning space for temporary local surface matrices: matrix 1
    n = nm(ftlp(fn,3),:);                                                   % Normal
    gp_p = gpf(ftlp(fn,3)).gp; gp_n = gpf(ftlp(fn,4)).gp;
    for gp = 1:ngp                                                          % Gauss point loop
        Jp = fdNr(gp_p(gp,:))*co_pve; Bp=Jp\fdNr(gp_p(gp,:));               % Jacobian for positive element and Global shape function derivatives positive
        Jn = fdNr(gp_n(1+ngp-gp,:))*co_nve; Bn=Jn\fdNr(gp_n(1+ngp-gp,:));   % Jacobian for negative element and Global shape function derivatives negative
        Np = fNR(gp_p(gp,:)); Nn = fNR(gp_n(1+ngp-gp,:));                   % shape function
        W = wp(gp)*dx/2; pen = 1^2/dx;                                      % Integral weight
        ke_1 = ke_1 - ((Bp'*n'*Np)/2)*W - ((Np'*n*Bp)/2)*W + 10*pen*(Np'*Np)*W; % Computing the surface stiffness matrix 1
        ke_2 = ke_2 + ((Bp'*n'*Nn)/2)*W - ((Np'*n*Bn)/2)*W - 10*pen*(Np'*Nn)*W; % Computing the surface stiffness matrix 2
        ke_3 = ke_3 - ((Bn'*n'*Np)/2)*W + ((Nn'*n*Bp)/2)*W - 10*pen*(Nn'*Np)*W; % Computing the surface stiffness matrix 3
        ke_4 = ke_4 + ((Bn'*n'*Nn)/2)*W + ((Nn'*n*Bn)/2)*W + 10*pen*(Nn'*Nn)*W; % Computing the surface stiffness matrix 4
    end
    kloc = max(kloc) + (1:neDoF);                                           % Steering matrix 1 index
    k.v1(kloc) = reshape(ke_1,neDoF,1);                                     % Storage of local stiffness matrix 1
    k.v2(kloc) = reshape(ke_2,neDoF,1);                                     % Storage of local stiffness matrix 2
    k.v3(kloc) = reshape(ke_3,neDoF,1);                                     % Storage of local stiffness matrix 3
    k.v4(kloc) = reshape(ke_4,neDoF,1);                                     % Storage of local stiffness matrix 4
    ed_p = ed(ftlp(fn,1),:);ed_n = ed(ftlp(fn,2),:);                        % Steering vector rows for matrix 1
    k.r1(kloc) = reshape(ed_p.'*ones(1,nov),neDoF,1);                       %
    k.c1(kloc) = reshape(ones(nov,1)*ed_p  ,neDoF,1);                       %
    k.r2(kloc) = reshape(ed_p.'*ones(1,nov),neDoF,1);                       %
    k.c2(kloc) = reshape(ones(nov,1)*ed_n  ,neDoF,1);                       %
    k.r3(kloc) = reshape(ed_n.'*ones(1,nov),neDoF,1);                       %
    k.c3(kloc) = reshape(ones(nov,1)*ed_p  ,neDoF,1);                       %
    k.r4(kloc) = reshape(ed_n.'*ones(1,nov),neDoF,1);                       %
    k.c4(kloc) = reshape(ones(nov,1)*ed_n  ,neDoF,1);                       %
end

% Compute solution --------------------------------------------------------
K = sparse([krow;k.r1;k.r2;k.r3;k.r4],[kcol;k.c1;k.c2;k.c3;k.c4],[kval; k.v1; k.v2;k.v3; k.v4]);
fd = ed;
for i = 1:nnls
    for j = 1:4
        x = coord(etpl(i,j),1);
        y = coord(etpl(i,j),2);
        if abs(x)==1 || abs(y)==1
            fd(i,j)=-1;
        end
    end
end
fd=unique(fd);
fd(fd==-1)=[];
u = zeros(max(ed(:)),1); u(fd)=K(fd,fd)\F(fd);                              % solve for diffusion
plot_u(etpl,coord,u,dx)
end

function [etpl,ftlp,nnls,nls,nfcs,coord] = makemesh(dx)
% making the coordinates  --------------------------------------------------
[X,Y] = meshgrid(-1:dx:1);nnx = length(-1:dx:1);                            % create 1D mesh stencil
x=reshape(X,[],1);y=reshape(Y,[],1);coord = [x,y];                          % create coordinates
nnds = nnx^2; nnls = (nnx-1)^2;nls = 1:nnls;                                % number of nodes and elements
n1 = 1:nnds;n1((end-nnx):end)=[];n1(nnx:nnx:end)=[];                        % create the list for node 1

% make the element topology -----------------------------------------------
etpl = n1'.*ones(1,4); etpl(:,2)=etpl(:,1)+1; etpl(:,3)=etpl(:,1)+nnx+1; etpl(:,4)=etpl(:,3)-1; % creating the element topology matrix

% make the face topology --------------------------------------------------
ftlp_in=ones(size(etpl));                                                   % Creating index to keep track of faces already visited
ftlp=zeros(size(etpl,1)*3,4);                                               % Allocating memory for etpl_face
face_count=0; fst=1:4;                                                      % Row counter for ftlp and local face list
etpl_face_nodes=[etpl(:,[1 2]),etpl(:,[2 3]),etpl(:,[3 4]),etpl(:,[4 1])];  % Matrix of nodes for faces 1,2 and 3 anticlockwise
for el = 1:size(ftlp_in,1)                                                  % Looping through all elements
    for loc_face = 1:size(ftlp_in,2)                                        % Looping through element faces
        if ftlp_in(el,loc_face)==1                                          % If the face has not been used before
            pos_nodes=etpl_face_nodes(el,(loc_face*2-1:loc_face*2));        % Positive element nodes
            face_count=face_count+1;                                        % Face topology row count
            ftlp(face_count,[1,3])= [el,loc_face];                          % Adding positive element characteristics to the matrix
            el_neg   = nls(sum(etpl_face_nodes(:,2:2:end)==pos_nodes(1),2)==1); % Finding possible negative element
            for i = 1:length(el_neg)                                        % Looping through possible negative elements
                neg_face=fst(etpl_face_nodes(el_neg(i),2:2:end)==pos_nodes(1)); % Negative element face
                neg_global_node=etpl_face_nodes(el_neg(i),(neg_face*2)-1);  % Corressponding node to negative element face
                if neg_global_node==pos_nodes(2)                            % Confirming that the positive and element elements have an adjacent face
                    ftlp(face_count,[2,4])=[el_neg(i),neg_face];            % Adding negative element contributions to etpl_face
                    ftlp_in(el_neg(i),neg_face)=0;                          % Face is marked as visited
                end
            end
        end
    end
end
ftlp(ftlp(:,2)==0,:) = [];nfcs = size(ftlp,1);
end






