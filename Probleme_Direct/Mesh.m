classdef Mesh
    properties
        domain
        nodes
        boundnodes
        intnodes
        triangles
        tricenters
        edges
        Nnodes
        Ntriangles
        Nedges
        int1
        intx
        intxx
        basis
        stats
    end
    methods
        
        % Constructors
        function obj = Mesh(domain,dx,s)
            % Class MESH
            %
            % Constructor :
            %
            % obj = Mesh(domain,dx) Class constuctor.
            % Create a mesh from a domain object (see Domain) a maximum
            % triangle obj = Mesh(domain,dx,'save') allows matlab to save
            % the mesh object and load it for the next calls which save time.
            %
            % METHOD LIST :
            %
            %
            % obj.disp
            % obj.plot
            %
            % u = obj.P0(input,order)
            % u = obj.P1(input,order)
            % u = obj.P0bound(input,order)
            % z = obj.P1togrid(u,x,y)
            % u = obj.gridtoP1(u_grid,x,y)
            % uP0 = obj.P1toP0(uP1)
            % p = obj.P0prod(u,v)
            % n = obj.L2norm(u)
            % v = obj.norm(u,p)
            %
            % surf(u)
            % contour(u)
            % quiver(u)
            %
            % gu = obj.grad(u)
            % gsu = obj.grads(u)
            % divu = obj.div(u)
            %
            % u = obj.solvesystem(A,F)
            %
            % u = obj.solve(a,c,f,bclist)
            % u = obj.solvepde(a,c,f,bclist)
            % K = obj.stiffness(a)
            % M = obj.mass(c)
            % F = obj.source(f)
            % Fb = obj.boundsource(E,f)
            % [Mb,Fb,Id,Vd] = obj.bcmatrices(bclist)
            % [A,F,Id,Vd] = obj.pde(a,c,f,bclist)
            % [A,F] = obj.system(a,c,f,bclist)
            %
            % u = obj.solve_vect(a,c,f,bclist)
            % u = obj.solvepde_vect(a,c,f,bclist)
            % K = obj.stiffness_vect(a)
            % A = obj.mass_vect(a)
            % F = obj.source_vect(f)
            % F = obj.boundsource_vect(E,f)
            % [Mb,Fb,Id,Vd] = obj.bcmatrices_vect(bclist)
            % [A,F,Id,Vd] = obj.pde_vect(a,c,f,bclist)
            % [A,F] = obj.system_vect(a,c,f,bclist)
            %
            % E = obj.buildbasis
            % [int1,intx,intxx] = obj.integrals
            % [A,F] = obj.assemble(A,F,Id,Vd)
            % u = obj.eval(input,x,y,order)
            % m = obj.displace(u)
            
            
            
            
            % If already exists then load it else save it !
            if nargin == 3 && strcmp(s,'save')
                a = dir('mesh*.mat');
                N = length(a);
                for i = 1 : N
                    load(a(i).name,'saved_dx');
                    load(a(i).name,'saved_domain');
                    if saved_dx == dx && saved_domain == domain
                        disp('loading existing mesh')
                        load(a(i).name,'saved_mesh');
                        obj = saved_mesh;
                        return
                    end
                end
                obj = Mesh(domain,dx);
                saved_mesh = obj;
                saved_dx = dx;
                saved_domain = domain;
                name = ['mesh' num2str(N+1) '.mat' ];
                save(name,'saved_mesh','saved_dx','saved_domain')
                return

            end
            disp('buiding mesh')
            
            time=cputime;
            [gm,Ibound]=domain.matrix(dx);
            
            [p,e,t]=initmesh(gm,'Hmax',dx,'Hgrad',1.5,'Jiggle','mean','JiggleIter',10);
            p=p';
            e=e';
            t=t';
            
            obj.domain=domain;
            obj.nodes=p;
            obj.edges=[e(:,[1 2]) Ibound(e(:,5))];
            obj.triangles=t(:,1:3);
            
            obj.tricenters=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            
            obj.boundnodes = union(obj.edges(:,1),obj.edges(:,2));
            obj.intnodes = setdiff(1:size(p,1), obj.boundnodes)';
            
            
            q=mean(pdetriq(p',t'));
            a=2*sqrt(pdetrg(p',t'))/(5^(1/4));
            r=mean(a);
            rmin=min(a);
            rmax=max(a);
            [i1,ix,ixx] = obj.integrals;
            
            obj.Nnodes = size(p,1);
            obj.Ntriangles = size(t,1);
            obj.Nedges = size(e,1);
            
            obj.int1 = i1;
            obj.intx = ix;
            obj.intxx = ixx;
            obj.basis = obj.buildbasis;
            w=whos('obj');
            mem=w.bytes/1000;
            
            time=cputime-time;
            obj.stats=[dx rmin rmax r rmax/rmin q mem time];
            
           
        end
        function m = displace(obj,u)
            m=obj;
            m.nodes = m.nodes + u;
        end
        
        
        % Refiners
        function [m,P] = naiverefine(obj)
            time=cputime;
            m = obj;
            Nn = obj.Nnodes;
            Nt = obj.Ntriangles;
            m.nodes = [obj.nodes ; obj.tricenters];
            m.intnodes = [obj.intnodes ; Nn + (1:Nt)'];
            
            m.triangles = [];
            for i = 1 : obj.Ntriangles
                t = obj.triangles(i,:);
                m.triangles = [m.triangles ; [t(1) t(2) Nn+i ; t(2) t(3) Nn+i ; t(3) t(1) Nn+i]];
            end
            p = m.nodes;
            t = m.triangles;
            m.tricenters=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            m.Nnodes = Nn+Nt;
            m.Ntriangles = 3*Nt;
            
            q=mean(pdetriq(p',t'));
            a=2*sqrt(pdetrg(p',t'))/(5^(1/4));
            r=mean(a);
            rmin=min(a);
            rmax=max(a);
            [i1,ix,ixx] = m.integrals;
            

            m.int1 = i1;
            m.intx = ix;
            m.intxx = ixx;
            m.basis = m.buildbasis;
            w=whos('m');
            mem=w.bytes/1000;
            
            I = (1:3*Nt)';
            J = ones(3,1)*(1:Nt);
            J  = J(:);
            V = ones(3*Nt,1);
            P = sparse(I,J,V,3*Nt,Nt);
            
            
            
            
            time=cputime-time;
            m.stats=[obj.stats(1) rmin rmax r rmax/rmin q mem time];
            
            
            

            
        end
        function m = refine(obj,domain)
            disp('refining mesh')
            m = obj;
            
            dx = obj.stats(1);
            time=cputime;
            [gm,Ibound]=domain.matrix(dx);
            
            p = obj.nodes;
            e = obj.edges;
            t = obj.triangles;
            
            ne = size(e,1);
            nt = size(t,1);
            e2 = [e(:,1:2) zeros(ne,1) ones(ne,1) e(:,1) ones(ne,1) zeros(ne,1)];
            t2 = [t ones(nt,1)];
            
            [p,e,t]=refinemesh(gm,p',e2',t2');
            p=jigglemesh(p,e,t,'Iter',2);
            p=p';
            e=e';
            t=t';
            m.nodes=p;
            m.edges=[e(:,[1 2]) Ibound(e(:,5))];
            m.triangles=t(:,1:3);
            
            m.tricenters=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            
            m.boundnodes = union(m.edges(:,1),m.edges(:,2));
            m.intnodes = setdiff(1:size(p,1), m.boundnodes)';
            
            
            q=mean(pdetriq(p',t'));
            a=2*sqrt(pdetrg(p',t'))/(5^(1/4));
            r=mean(a);
            rmin=min(a);
            rmax=max(a);
            [i1,ix,ixx] = m.integrals;
            
            m.Nnodes = size(p,1);
            m.Ntriangles = size(t,1);
            m.Nedges = size(e,1);
            
            m.int1 = i1;
            m.intx = ix;
            m.intxx = ixx;
            m.basis = m.buildbasis;
            w=whos('obj');
            mem=w.bytes/1000;
            
            time=cputime-time;
            m.stats=[dx rmin rmax r rmax/rmin q mem time];
            
           
        end

                   
        % Display
        function disp(obj)
            % MESH.DISP
            %       obj.disp displays data of the Mesh object.
            disp('Mesh object')
            disp(['Nodes :            ' num2str(obj.Nnodes)])
            disp(['Triangles :        ' num2str(obj.Ntriangles)])
            disp(['dx :               ' num2str(obj.stats(1))])
            disp(['Min resolution :   ' num2str(obj.stats(2))])
            disp(['Max resolution :   ' num2str(obj.stats(3))])
            disp(['Mean resolution :  ' num2str(obj.stats(4))])
            disp(['Uniformity :       ' num2str(obj.stats(5))])
            disp(['Mean quality :     ' num2str(obj.stats(6))])
            disp(['Memory used(ko) :  ' num2str(obj.stats(7))])
            disp(['Time(s) :          ' num2str(obj.stats(8))])
            disp(' ')
        end
        function plot(obj)
            % MESH.PLOT
            %       obj.plot plots the Mesh object.
            t = obj.triangles;
            t = [t ones(length(t),1)];
            pdemesh(obj.nodes',obj.edges',t');
            hold on
            obj.domain.plot;
            axis image
        end
        
        % Functionnal spaces
        
        function u = P0(obj,input,order)
            % MESH.P0:
            %       u = mesh.P0(input) creates a P0 function from input. u
            %       = mesh.P0(input,order) creates a P0 tensor of a given
            %       order.
            switch nargin
                case 1
                    u = zeros(obj.Ntriangles,1);
                case 2
                    if ischar(input)
                        order = 0;
                    else
                        order = log2(size(input,2));
                    end
                    x = obj.tricenters(:,1);
                    y = obj.tricenters(:,2);
                    u = obj.eval(input,x,y,order);
                case 3
                    x = obj.tricenters(:,1);
                    y = obj.tricenters(:,2);
                    u = obj.eval(input,x,y,order);
                otherwise
                    error('Wrong inputs number')
            end
        end
        function u = P1(obj,input,order)
            % MESH.P1:
            %       u = mesh.P1(input) creates a P1 function from input. u
            %       = mesh.P1(input,order) creates a P1 tensor of a given
            %       order.
            
            if nargin == 2
                order = 0;
            end
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            u = obj.eval(input,x,y,order);
            
            
        end 
        function u = P0bound(obj,input,order)
                    if nargin == 2
                if ischar(input)
                    order = 0;
                else
                    order = log2(length(input));
                end
                    end
            e1 = obj.edges(:,1);
            e2 = obj.edges(:,2);
            x = (obj.nodes(e1,1) + obj.nodes(e2,1))/2;
            y = (obj.nodes(e1,2) + obj.nodes(e2,2))/2;
            u = obj.eval(input,x,y,order);
        end   
        
        % Interpolation
        function z = P1togrid(obj,u,x,y)
            z = tri2grid(obj.nodes',obj.triangles',u,x,y);
        end
        function u = gridtoP1(obj,x,y,u_grid)
            [X,Y] = meshgrid(x,y);
            u = interp2(X,Y,u_grid,obj.nodes(:,1),obj.nodes(:,2));
        end
        function u = gridtoP0(obj,X,Y,u_grid)
            u = interp2(X,Y,u_grid,obj.tricenters(:,1),obj.tricenters(:,2));
        end
        function uP0 = P1toP0(obj,uP1)
            N=size(uP1,2);
            nT=obj.Ntriangles;
            uP0=zeros(nT,N);
            for i = 1 : N
                uP0(:,i)=pdeintrp(obj.nodes',obj.triangles',uP1(:,i));
            end
        end
        function uP1 = P0toP1(obj,uP0)
            N=size(uP0,2);
            nN=obj.Nnodes;
            uP1=zeros(nN,N);
            for i = 1 : N
                uP1(:,i)=pdeprtni(obj.nodes',obj.triangles',uP0(:,i)');
            end
        end
        
        function p = P0prod(obj,u,v)
            a = obj.int1;
            p = sum(u.*v.*a);
        end 
        function n = L2norm(obj,u)
            % MESH.L2NORM
            %       n = obj.L2norm(u) compute the L2 norm of P1 function u.
           Nt = obj.Ntriangles;
                M = obj.mass(ones(Nt,1));
                n = sqrt(u'*(M*u));
                
               
        end
        function v = norm(obj,u,p)
            v = (sum(abs(u).^p,2)).^(1/p); 
        end
        
        % Drawers
        function surf(obj,u)
            % MESH.SURF
            %       obj.surf(u) plots a 3D colored graph of u.
            N = size(u,1);
            t = obj.triangles;
            t = [t ones(length(t),1)];
            
            if obj.Ntriangles == N
                if isreal(u)
                    pdesurf(obj.nodes',t',u');
                    view([0 90])
                    colormap jet
                    colorbar
                    %axis image
                else
                    subplot(1,2,1);
                    pdesurf(obj.nodes',t',real(u'));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                    
                    subplot(1,2,2);
                    pdesurf(obj.nodes',t',imag(u'));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                end
            elseif  obj.Nnodes == N
                if isreal(u)
                    pdesurf(obj.nodes',t',u);
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                else
                    subplot(1,2,1);
                    pdesurf(obj.nodes',t',real(u));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                    
                    subplot(1,2,2);
                    pdesurf(obj.nodes',t',imag(u));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                end
            end
        end
        function contour(obj,u)
            % MESH.SURF
            %       obj.surf(u) plots a 3D colored graph of u.
            N = size(u,1);
            
            if obj.Ntriangles == N
                if isreal(u)
                    pdecont(obj.nodes',obj.triangles',u');
                    colorbar
                else
                    subplot(1,2,1);
                    pdecont(obj.nodes',obj.triangles',real(u'));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                    
                    subplot(1,2,2);
                    pdecont(obj.nodes',obj.triangles',imag(u'));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                end
            elseif  obj.Nnodes == N
                T = obj.triangles;
                T = [T ones(size(T,1),1)]';
                if isreal(u)
                    pdecont(obj.nodes',T,u);
                    colorbar
                else
                    subplot(1,2,1);
                    pdesurf(obj.nodes',T,real(u));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                    
                    subplot(1,2,2);
                    pdecont(obj.nodes',T,imag(u));
                    view([0 90])
                    colormap jet
                    colorbar
                    axis image
                end
            end
            
        end
        function quiver(obj,u)
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            quiver(x,y,u(:,1),u(:,2));
            axis equal
        end
        
        % Differential operators
        function gu = grad(obj,u)
            Nn = obj.Nnodes;
            Nt = obj.Ntriangles;
            if size(u,1) == Nn
                Nc = size(u,2);
                gu = zeros(Nt,2*Nc);
                for i = 1 : Nc
                    [gu1,gu2]=pdegrad(obj.nodes',obj.triangles',u(:,i));
                    gu(:,i) = gu1';
                    gu(:,Nc+i) = gu2';
                end

            else
                disp('Mesh.grad: Input must be P1.')
                gu=[];
            end 
        end
        function gsu = grads(obj,u)
            gu = obj.grad(u);
            gsu = gu;
            gsu(:,2:3) = ((gu(:,2) + gu(:,3))/2)*[1 1];
        end
        function divu = div(obj,u)
            gu1=obj.grad(u(:,1));
            gu2=obj.grad(u(:,2));
            divu=gu1(:,1)+gu2(:,2);
        end
        
        % Solvers
        function u = solvesystem(obj,A,F)
            Nn = obj.Nnodes;
           U = A\F;
           Nu = length(U);
           switch Nu
               case Nn
                   u = U;
               case 2*Nn
                   u = reshape(U,Nn,2);
               case Nn+1
                   u = U(1:Nn);
               case 2*Nn+3
                   u = reshape(U(1:2*Nn,:),Nn,2);
           end
               
        end
        function u = solvepde_vect(obj,a,c,f,bclist)
            a = obj.P0(a,4);
            c = obj.P0(c,2);
            f = obj.P0(f,1);
            A = obj.intgugv_vect(a);
            C = obj.intuv_vect(c);
            F = obj.intv_vect(f);
            [Ab,Fb,Id,Vd] = obj.bcmatrices_vect(bclist);
            [Ad,Fd] = obj.dirichlet(A+C+Ab,F+Fb,Id,Vd);
            u = obj.solvesystem(Ad,Fd);
        end
        
        % Scalar PDE matrices
        function K = stiffness(obj,a)
            Nt = obj.Ntriangles;
            if size(a,2) ~= 4 || size(a,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 order 2 tensor !');
                K = [];
                return
            end
            
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);
            
            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);
            
            a = a(rep9,:);
            
            i1 = obj.int1(rep9);
            
            agej = [a(:,1).*alphaj(:,1) + a(:,2).*alphaj(:,2) a(:,3).*alphaj(:,1) + a(:,4).*alphaj(:,2)];
            geiagej = i1.*(alphai(:,1).*agej(:,1) + alphai(:,2).*agej(:,2));
            K = sparse(I,J,geiagej);
        end
        function A = intugv(obj,b)
            Nt = obj.Ntriangles;
            if size(b,2) ~= 2 || size(b,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 vector function !');
                A = [];
                return
            end
            
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);

            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);
            
            betamat = obj.basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';
            betaj = beta(rep33,:);
            
            b = b(rep9,:);
            
            i1 = obj.int1(rep9);
            ix = obj.intx(rep9,:);
            
            
            V1 = (alphai(:,1).*b(:,1) + alphai(:,2).*b(:,2)).*(alphaj(:,1).*ix(:,1) + alphaj(:,2).*ix(:,2));
            V2 = (alphai(:,1).*b(:,1) + alphai(:,2).*b(:,2)).*betaj.*i1;
           
            
        
            aeiej = V1+V2;
           
            A = sparse(I,J,aeiej);
        end
        function M = mass(obj,c)
            Nt = obj.Ntriangles;
            if nargin == 1
                M = obj.mass(ones(Nt,1));
                return
            end
            if size(c,2) ~= 1 || size(c,1)~= Nt
                disp('MESH.mass error : input must be an P0 scalar function !');
                M = [];
                return
            end
            
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);

            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);

            betamat = obj.basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';
            betai = beta(rep3,:);
            betaj = beta(rep33,:);
            
            c = c(rep9,:);
            
            i1 = obj.int1(rep9);
            ix = obj.intx(rep9,:);
            ixx = obj.intxx(rep9,:);
            
            V1 = alphai(:,1).*(ixx(:,1).*alphaj(:,1) + ixx(:,2).*alphaj(:,2)) +... 
                 alphai(:,2).*(ixx(:,3).*alphaj(:,1) + ixx(:,4).*alphaj(:,2));
            V2 = (alphai(:,1).*ix(:,1) + alphai(:,2).*ix(:,2)).*betaj +...
                 (alphaj(:,1).*ix(:,1) + alphaj(:,2).*ix(:,2)).*betai;
            V3 = betai.*betaj.*i1;
            
        
            aeiej = c.*(V1 + V2 + V3);
           
            M = sparse(I,J,aeiej);
        end 
        function Mb = boundmass(obj,E,g)
             Ne = obj.Nedges;
             Nn = obj.Nnodes;
            if size(g,2) ~= 1 || size(g,1)~= Ne
                error('input 2 must be a P0bound scalar function !');
            end
            
            ebound = obj.edges;

            im = ismember(ebound(:,3),E);
            g = g(im,:);
            e = ebound(im,1:2);
            
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            L = sqrt((x(e(:,2)) - x(e(:,1))).^2 + (y(e(:,2)) - y(e(:,1))).^2);
            e = e';
            I = e([1 1 2 2],:);
            J = e([1 2 1 2],:);
            
            V = 1/6*[2;1;1;2]*(L.*g)';
            Mb = sparse(I,J,V,Nn,Nn);
        end
        
        function F = source(obj,f)
            Nt = obj.Ntriangles;
            if size(f,2) ~= 1 || size(f,1)~= Nt
                disp('MESH.source error : input must be a P0 scalar function !');
                F = [];
                return
            end
            tri = obj.triangles;
            tp=tri';
            rep3 = floor((3:3*Nt+2)/3);
            
            I = tp(:);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            
            betamat = obj.basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';

            f = f(rep3,:);
            
            i1 = obj.int1(rep3);
            ix = obj.intx(rep3,:);
            
            fei = f.*(alpha(:,1).*ix(:,1) + alpha(:,2).*ix(:,2) + beta.*i1);
            F = full(sparse(I,ones(3*Nt,1),fei));

        end  
        function Fb = boundsource(obj,E,f)
            Ne = obj.Nedges;
            Nn = obj.Nnodes;
            if size(f,2) ~= 1 || size(f,1)~= Ne
                disp('MESH.boundsource error : input must be a P0bound scalar function !');
                Fb = [];
                return
            end
            ebound = obj.edges;

            im = ismember(ebound(:,3),E);
            f = f(im,:);
            e = ebound(im,1:2);
            Ne = length(e);
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            L = sqrt((x(e(:,2)) - x(e(:,1))).^2 + (y(e(:,2)) - y(e(:,1))).^2);
            ep = e';
            
            I = ep(:);
            J = ones(2*Ne,1);

            rep2 = floor((2:2*Ne+1)/2);

            f = f(rep2,:);
            
            fei = f.*L(rep2)/2;
            Fb = full(sparse(I,J,fei,Nn,1));
        end
        
        % Vectorial PDE matrices
        function K = stiffness_vect(obj,a)
            Nt = obj.Ntriangles;
            if size(a,2) ~= 16 || size(a,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 order 4 tensor !');
                K = [];
                return
            end
            
            A1 = obj.stiffness(a(:,1:4));
            A2 = obj.stiffness(a(:,5:8));
            A3 = obj.stiffness(a(:,9:12));
            A4 = obj.stiffness(a(:,13:16));
            K = [A1 A2 ; A3 A4];
            
        end
        function A = mass_vect(obj,a)
            Nt = obj.Ntriangles;
            if size(a,2) ~= 4 || size(a,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 order 2 tensor !');
                A = [];
                return
            end
            
            A1 = obj.mass(a(:,1));
            A2 = obj.mass(a(:,2));
            A3 = obj.mass(a(:,3));
            A4 = obj.mass(a(:,4));
            A = [A1 A2 ; A3 A4];
            
        end       
        function F = source_vect(obj,f)
            Nt = obj.Ntriangles;
            if size(f,2) ~= 2 || size(f,1)~= Nt
                disp('MESH.intv_vect error : input must be a P0 vector function !');
                F = [];
                return
            end
            F1 = obj.source(f(:,1));
            F2 = obj.source(f(:,2));
            F = [F1;F2];

        end 
        function F = boundsource_vect(obj,E,f)
            Ne = obj.Nedges;
            if size(f,2) ~= 2 || size(f,1)~= Ne
                disp('MESH.intboundv_vect error : input must be a P0 vector function !');
                F = [];
                return
            end
            F1 = obj.boundsource(E,f(:,1));
            F2 = obj.boundsource(E,f(:,2));
            F = [F1;F2];
         end
        
        % Boundary Conditions applyers
        function [Mb,Fb,Id,Vd] = bcmatrices(obj,bclist)
            Nbclist = length(bclist);
            Nn = obj.Nnodes;
            e = obj.edges(:,3);
            if Nbclist == 0
                Id = -1;
                Vd = [];
                Mb = sparse(Nn,Nn);
                Fb = zeros(Nn,1);
                return
            end
            onlyneumann = 1;
            Id = [];
            Vd = zeros(Nn,1);
            Mb = sparse(Nn,Nn);
            Fb = zeros(Nn,1);
            for k = 1 : Nbclist
                bc = bclist{k};
                if bc{2} == 0
                    onlyneumann = 0;
                    ed=bc{1};
                    id=[];
                    for i = 1 : length(ed)
                        id=[id ; obj.edges(e==ed(i),1) ; obj.edges(e==ed(i),2)];
                    end
                    id=unique(id);
                    Id=union(Id,id);
                    Vdk=obj.P1(bc{4},0);
                    Vd(id)=Vdk(id);
                end
                if bc{2} == 1
                    g = obj.P0bound(bc{3});
                    if sum(abs(g))~=0
                        onlyneumann = 0;
                        Mb = Mb + obj.boundmass(bc{1},g);
                    end
                    h = obj.P0bound(bc{4},0);
                    Fb = Fb + obj.boundsource(bc{1},h);
                    
                        
                end    
                
            end
            if onlyneumann
                    Id = -1;
            end
        end   
         function [Mb,Fb,Id,Vd] = bcmatrices_vect(obj,bclist)
            Nbclist = length(bclist);
            Nn = obj.Nnodes;
            e = obj.edges(:,3);
            if Nbclist == 0
                Id = -1;
                Vd = [];
                Mb = sparse(2*Nn,2*Nn);
                Fb = zeros(2*Nn,1);
                return
            end
            onlyneumann = 1;
            Id = [];
            Vd = zeros(2*Nn,1);
            Mb = sparse(2*Nn,2*Nn);
            Fb = zeros(2*Nn,1);
            for k = 1 : Nbclist
                bc = bclist{k};
                if bc{2} == 0
                    onlyneumann = 0;
                    ed=bc{1};
                    id=[];
                    for i = 1 : length(ed)
                        id=[id ; obj.edges(e==ed(i),1) ; obj.edges(e==ed(i),2)];
                    end
                    id=unique(id);
                    Id=union(Id,id);
                    h = obj.P1(bc{4},1);
                    Vd(id)=h(id,1);
                    Vd(Nn+id)=h(id,2);
                end
                
                if bc{2} == 1
                    g = obj.P0bound(bc{3},2);
                    if sum(abs(g))~=0
                        onlyneumann = 0;
                        Mb = Mb + obj.boundmass_vect(bc{1},g);
                    end
                    h = obj.P0bound(bc{4},1);
                    Fb = Fb + obj.boundsource_vect(bc{1},h);
                    
                        
                end    
                
            end
            Id = [Id ; Nn+Id];
            
            if onlyneumann
                    Id = -1;
            end
        end   
        
        % Scalar pde builders and solvers
        function [A,F,Id,Vd] = pde(obj,a,c,f,bclist)
            Nn = obj.Nnodes;
            K = sparse(Nn,Nn);
            M = sparse(Nn,Nn);
            S = zeros(Nn,1);
            
            
            if nnz(a)
                K = obj.stiffness(obj.P0(a,2));
            end
            if nnz(c)
                M = obj.mass(obj.P0(c));
            end
             if nnz(f)
                S = obj.source(obj.P0(f));
             end
             
             [Mb,Sb,Id,Vd] = obj.bcmatrices(bclist);
             A = K+M+Mb;
             F = S + Sb;
        end     
        function [A,F] = system(obj,a,c,f,bclist)
            [A,F,Id,Vd] = obj.pde(a,c,f,bclist);
            [A,F] = obj.assemble(A,F,Id,Vd);
        end
        function u = solve(obj,a,c,f,bclist)
            if nargin == 3
                A = a;
                F = c;
            else
                [A,F] = obj.system(a,c,f,bclist);
            end
            u = A\F;
        end
        
        % Vectorial pde builders and solvers
        function [A,F,Id,Vd] = pde_vect(obj,a,c,f,bclist)
            Nn = obj.Nnodes;
            K = sparse(2*Nn,2*Nn);
            M = sparse(2*Nn,2*Nn);
            S = zeros(2*Nn,1);
            
            if nnz(a)
                K = obj.stiffness_vect(obj.P0(a,4));
            end
            if nnz(c)
                M = obj.mass_vect(obj.P0(c,2));
            end
             if nnz(f)
                S = obj.source_vect(obj.P0(f,1));
             end
             
             [Mb,Sb,Id,Vd] = obj.bcmatrices_vect(bclist);
             A = K+M+Mb;
             F = S + Sb;
        end
        function [A,F] = system_vect(obj,a,c,f,bclist)
            [A,F,Id,Vd] = obj.pde_vect(a,c,f,bclist);
            [A,F] = obj.assemble(A,F,Id,Vd);
        end
        function u = solve_vect(obj,a,c,f,bclist)
            if nargin == 3
                A = a;
                F = c;
            else
                [A,F] = obj.system_vect(a,c,f,bclist);
            end
            Nn = obj.Nnodes;
            U = A\F;
            u = [U(1:Nn) U(Nn+1:end)];
        end

        
        % Tools
        function E = buildbasis(obj)
            Node1=obj.nodes(obj.triangles(:,1),:);
            Node2=obj.nodes(obj.triangles(:,2),:);
            Node3=obj.nodes(obj.triangles(:,3),:);
            
            E = zeros(obj.Ntriangles,9);
            for j = 1 : obj.Ntriangles
                M=[[Node1(j,:) ; Node2(j,:) ; Node3(j,:)] [1;1;1]];
                N=inv(M);
                E(j,:)=N(:)';
            end
        end
        function [int1,intx,intxx] = integrals(obj)
            p=obj.nodes;
            t=obj.triangles;

            p1=p(t(:,1),:)';
            p2=p(t(:,2),:)';
            p3=p(t(:,3),:)';
            M=[p2-p1;p3-p1];
            det=abs(M(1,:).*M(4,:)-M(2,:).*M(3,:));
            
            int1 = det'/2;
            
            intx=[det.*(1/2*p1(1,:)+1/6*(M(1,:)+M(3,:))) ; det.*(1/2*p1(2,:)+1/6*(M(2,:)+M(4,:)))]';
            
            I1=1/2*[p1(1,:).^2;p1(1,:).*p1(2,:);p1(1,:).*p1(2,:);p1(2,:).^2];
            I2=1/6*[p1(1,:).*(M(1,:)+M(3,:));p1(1,:).*(M(2,:)+M(4,:));p1(2,:).*(M(1,:)+M(3,:));p1(2,:).*(M(2,:)+M(4,:))];
            I3=[I2(1,:);I2(3,:);I2(2,:);I2(4,:)];
            
            L1=1/24*(2*M(1,:).^2+2*M(1,:).*M(3,:)+2*M(3,:).^2);
            L2=1/24*(2*M(2,:).*M(1,:)+M(4,:).*M(1,:)+M(3,:).*M(2,:)+2*M(4,:).*M(3,:));
            L4=1/24*(2*M(2,:).^2+2*M(2,:).*M(4,:)+2*M(4,:).^2);
            
            I4=[L1;L2;L2;L4];
            I=I1+I2+I3+I4;
            intxx=[det.*I(1,:);det.*I(2,:);det.*I(3,:);det.*I(4,:)]';
            
        end
        function [A,F] = assemble(obj,A,F,Id,Vd)
            Nn = obj.Nnodes;
            if Id == -1
                if length(F) == Nn
                    if abs(sum(F))> 1e-12
                        disp('Mesh.assemble warning : sum(F) must be zero for full Neumann boundary conditions !')
                        F = F - sum(F)/Nn;
                    end
                    A(1,:)=0;
                    A(1,1)=1;
                    F(1)=0;
                elseif length(F) == 2*Nn
                    F1 = F(1:Nn);
                    F2 = F(Nn+1:2*Nn);
                    if abs(sum(F))> 2e-12
                        F1 = F1 - sum(F1)/Nn;
                        F2 = F2 - sum(F2)/Nn;
                    end
                    A([1 2 Nn+1],:)=0;
                    A(1,1) = 1;
                    A(2,2) = 1;
                    A(Nn+1,Nn+1) = 1;
                    F1([1 2 Nn+1]) = 0;
                    
                end
            elseif ~isempty(Id)
                NId = length(Id);
                A(Id,:) = 0;
                A(Id,Id) = speye(NId,NId);
                F(Id) = Vd(Id);
            end
        end
        function u = eval(obj,input,x,y,order)
            N = length(x);

            w = whos('input');
            switch w.class
                case 'char'
                    mat = eval(input);
                case 'cell'
                    l = length(input);
                    mat = zeros(N,l);
                    for i = 1 : l
                        if ischar(input{i})
                            mat(:,i) = eval(input{i});
                        elseif size(input{i},1) == 1 && size(input{i},2) == 1
                            mat(:,i) = input{i}*ones(N,1);
                        elseif size(input{i},1) == N && size(input{i},2) == 1
                            mat(:,i) = input{i};
                        else
                            disp('Mesh.eval : wrong format in input cell !')
                            return
                        end
                    end
                case 'double'
                    Nrow = size(input,1);
                    switch Nrow
                        case 1
                            mat = ones(N,1)*input;
                        case N
                            mat = input;
                        otherwise
                            disp('Mesh.eval : wrong format in input matrix !')
                            mat = [];
                            return
                    end
               case 'function_handle'
                   mat = input(x,y);
            end
            
            
            Ncol = size(mat,2);
            if nargin == 2
                u = mat;
            else
                switch order
                    case 0
                        if Ncol == 1
                            u = mat;
                        else
                            disp('Mesh.eval : mismatch dimension for order 0!')
                            u = [];
                            return
                        end
                    case 1
                        switch Ncol
                            case 1
                                u = mat*[1 1];
                            case 2
                                u = mat;
                            otherwise
                                disp('Mesh.eval : mismatch dimension for order 1!')
                                u = [];
                                return
                        end
                        
                    case 2
                        switch Ncol
                            case 1
                                u = mat*[1 0 0 1];
                            case 2
                                u = [mat(:,1) 0 0 mat(:,2)];
                            case 3
                                u = [mat(:,1) mat(:,2) mat(:,2) mat(:,3)];
                            case 4
                                u = mat;
                            otherwise
                                disp('Mesh.eval : mismatch dimension for order 2!')
                                u = [];
                                return
                        end
                    case 4
                        switch Ncol
                            case 1
                                I = [1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1];
                                u = mat(:,1)*I;

                            case 2
                                I = [1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1];
                                T = [1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1];
                                D = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
                                u = (2*mat(:,1))*(I+T)/2 + mat(:,2)*D;
                                
                            case 4
                                s11 = mat(:,1);
                                s12 = mat(:,2);
                                s21 = mat(:,3);
                                s22 = mat(:,4);
                                u = [s11.*s11 s11.*s21  s21.*s11 s21.*s21 ...  
                                     s11.*s12 s11.*s22  s21.*s12 s21.*s22 ...
                                     s12.*s11 s12.*s21  s22.*s11 s22.*s21 ...  
                                     s12.*s12 s12.*s22  s22.*s12 s22.*s22 ...
                                      ];
                           case 8
                                s11 = mat(:,1);
                                s12 = mat(:,2);
                                s21 = mat(:,3);
                                s22 = mat(:,4);
                                t11 = mat(:,5);
                                t12 = mat(:,6);
                                t21 = mat(:,7);
                                t22 = mat(:,8);

                                u = [s11.*t11 s11.*t21  s21.*t11 s21.*t21 ...  
                                     s11.*t12 s11.*t22  s21.*t12 s21.*t22 ...
                                     s12.*t11 s12.*t21  s22.*t11 s22.*t21 ...  
                                     s12.*t12 s12.*t22  s22.*t12 s22.*t22 ...
                                      ];
                                
                                
                            case 16
                                u = mat;

                                
                            otherwise
                                disp('Mesh.eval error : mismatch dimension for order 4!')
                                u = [];
                                return
                        end
                    otherwise
                        disp('Mesh.eval error : order is not treated !')
                        u = [];
                        return
                end
            end
        end
        
        function I = boundary(obj,blist)
           if nargin == 1
               I = obj.boundnodes;
           end
           if nargin == 2
               N = length(blist);
               I = [];
               for i = 1 : N
                   ind = obj.edges(:,3)==blist(i);
                    Ii = union(obj.edges(ind,1),obj.edges(ind,2));
                    I = union(I,Ii);
               end
           end
        end
        
    end
end


function m=refine(obj,g,f,level)
disp('Refine is obsolete')
% MESH.REFINE
%       m=obj.refine(g) refines the Mesh object. g is the
%       initial Geometry object.
%
%       m=obj.refine(g,f) refines only the area specified by
%       the boolean function f. Matrix f is 1xnT, 1xnP or expr.
%
%       m=obj.refine(g,f,level) refines around the highest
%       value of the non negative map f. Matrix f is 1xnT, 1xnP
%       or expr. level is a double in [0,1].
time=cputime;
nt=obj.stats(3);
m=obj;
GM=g.matrix;
p=m.nodes;
e=m.edges;
t=m.triangles;

switch nargin
    case 2
        F=ones(1,nt);
        s=1/2;
    case 3
        F=obj.tensor(f,0,1);
        s=1/2;
    case 4
        F=obj.tensor(f,0,1);
        s=level;
end
maxf=max(F);
it=find(F>((1-s)*maxf));
[p,e,t]=refinemesh(GM,p',e',t');
p=jigglemesh(p,e,t,'Iter',2);

m.nodes=p;
m.edges=e;
m.triangles=t;
w=whos('obj');
mem=w.bytes/1000;
q=mean(pdetriq(p,t));
a=pdetrg(p,t);
r=mean(a);
rmin=min(a);
rmax=max(a);
time=cputime-time;
obj.Nnodes = size(p,2);
obj.Ntriangles = size(t,2);
m.stats=[dx rmin rmax r rmax/rmin q mem time+obj.stats(10)];
end




