classdef special_quad
    % special quadrature class. Contains definitions and helper functions
    % for the Cauchy and hypersinglar integrals
    
    properties
        L;
        Nup;
        x;
        w;
        xup;
        wup;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function sq = special_quad(Nup)
            sq.Nup = Nup;
            sq.L = legendre.matrix(16);
            [sq.x, sq.w] = legendre.gauss(16);    
            [sq.xup, sq.wup] = legendre.gauss(Nup);            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Ic = cauchy_integral(o, fsrc, ztar, zsrc, p0)
            
            cf = o.L*fsrc;
            cz = o.L*zsrc;

            P = legendre.vec(15, o.xup);
            fsrc_up = P*cf;
            zsrc_up = P*cz;
            
%             poly_coeff = o.vandernewtonT(zsrc_up,fsrc_up,32);
%             rel = min(abs(poly_coeff))/max(abs(poly_coeff));
%             if max(abs(poly_coeff)) > 1 && rel > 1e-10
%                 a = 0;
%             end
%             assert(max(abs(poly_coeff)) > 1 && rel > 1e-10, 'Cauchy integral: Polynomial coefficient too large (%e)', abs(coeff(end)));
            
            p = zeros(o.Nup,1);
            p(1) = p0;
            
            if abs(ztar) > 1.2
                p(end) = sum(o.xup.^(o.Nup - 1)./(o.xup - ztar).*o.wup);
                for k = o.Nup-1:-1:2
                    p(k) = (p(k+1) - (1 - (-1)^k)/k)/ztar;
                end
            else
                for k=2:o.Nup
                    p(k) = ztar*p(k-1) + (1 - (-1)^(k-1))/(k-1);
                end
            end
            
            w = o.vandernewton(zsrc_up,p,o.Nup);
            Ic = -(sum(fsrc_up.*w));
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Ih = hypersingular_integral(o, fsrc, ztar, zsrc, p0, b, a)
            
            p = zeros(o.Nup,1);
            r = zeros(o.Nup,1);
            
            % upsample density and source points
            cf = o.L*fsrc;
            cz = o.L*zsrc;

            P = legendre.vec(15, o.xup);
            fsrc_up = P*cf;
            zsrc_up = P*cz;
            
%             poly_coeff = o.vandernewtonT(zsrc_up,fsrc_up,32);
%             rel = min(abs(poly_coeff))/max(abs(poly_coeff));
%             if max(abs(poly_coeff)) > 1 && rel > 1e-10
%                 a = 0;
%             end
%             assert(max(abs(poly_coeff)) > 1 && rel > 1e-10, 'Hyper integral: Polynomial coefficient too large (%e)', abs(coeff(end)));
            
            p(1) = p0;
            r(1) = -1/(1 + ztar) - 1/(1 - ztar);
            
            if abs(ztar) > 1.2
                % compute p_32 and r_32 numerically, run recursion backwards
                p(end) = sum(o.xup.^(o.Nup - 1)./(o.xup - ztar).*o.wup);
                r(end) = sum(o.xup.^(o.Nup - 1)./(o.xup - ztar).^2.*o.wup);
                
                for k = o.Nup - 1:-1:2
                    p(k) = (p(k+1) - (1 - (-1)^k)/k)/ztar;
                    r(k) = (r(k+1) - p(k))/ztar;
                end

            else
                %run forward recursion
                for k=2:o.Nup
                    p(k) = ztar*p(k-1) + (1 - (-1)^(k-1))/(k-1);
                    r(k) = ztar*r(k-1) + p(k-1);
                end
            end
            
            w = o.vandernewton(zsrc_up, r, o.Nup);
            Ih = sum(2*fsrc_up/(b-a).*w);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Is = supersingular_integral(o, fsrc, ztar, zsrc, p0, b, a)
            
            p = zeros(o.Nup,1);
            r = zeros(o.Nup,1);
            s = zeros(o.Nup,1);
            
            p(1) = p0;
            r(1) = -1/(1 + ztar) - 1/(1 - ztar);
            s(1) = 1/(2*(1+ztar)^2) - 1/(2*(1-ztar)^2);
            
            % upsample density and source points
            cf = o.L*fsrc;
            cz = o.L*zsrc;

            P = legendre.vec(15, o.xup);
            fsrc_up = P*cf;
            zsrc_up = P*cz;
            
%             poly_coeff = o.vandernewtonT(zsrc_up,fsrc_up,32);
%             rel = min(abs(poly_coeff))/max(abs(poly_coeff));
%             if max(abs(poly_coeff)) > 1 && rel > 1e-10
%                 a = 0;
%             end
%             assert(max(abs(poly_coeff)) > 1 && rel > 1e-10, 'Super integral: Polynomial coefficient too large (%e)', abs(coeff(end)));
            
            if abs(ztar) > 1.2
                % compute p_32, r_32, and s_32 numerically, run recursion backwards
                p(end) = sum(o.xup.^(o.Nup - 1)./(o.xup - ztar).*o.wup);
                r(end) = sum(o.xup.^(o.Nup - 1)./(o.xup - ztar).^2.*o.wup);
                s(end) = sum(o.xup.^(o.Nup - 1)./(o.xup - ztar).^3.*o.wup);
                for k = o.Nup - 1:-1:2
                    p(k) = (p(k+1) - (1 - (-1)^k)/k)/ztar;
                    r(k) = (r(k+1) - p(k))/ztar;
                    s(k) = (s(k+1) - r(k))/ztar;
                end
            else
                for k=2:o.Nup
                    p(k) = ztar*p(k-1) + (1 - (-1)^(k-1))/(k-1);
                    r(k) = ztar*r(k-1) + p(k-1);
                    s(k) = ztar*s(k-1) + r(k-1);
                end
            end
            
            w = o.vandernewton(zsrc_up,s,o.Nup);
            Is = -sum(4*fsrc_up/(b-a)^2.*w);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function p0 = compute_exact_log(o, nz, nzsrc)
            
            
            lg1 = log(1-nz);
            lg2 = log(-1-nz);
            
            %Check if the point nz is between the panel and the real axis
            if real(nz) > -1 && real(nz) < 1
                if imag(nz) > 0 %above real axis, check if enclosed by axis and panel
                    furthercheck = 0;
                    
                    if sum(imag(nzsrc) > imag(nz)) ~= 0
                        furthercheck = 1;
                    end
                    
                    if furthercheck
                        %interpol. nzpan to poly and check value for
                        %real(nz)
                        tmpT = real(nzsrc);
                        tmpb = imag(nzsrc);
                        
                        p = o.vandernewtonT(tmpT,tmpb,16);
                        
                        kk = (0:15)';
                        test = sum(p.*real(nz).^kk);
                        if abs(test - imag(nz)) > 1e-14
                            if test > imag(nz) %Correct value of integral
                                lg1 = lg1 - pi*1i;
                                lg2 = lg2 + pi*1i;
                            end
                        end
                    end
                else
                    if imag(nz) < 0 %below the real axis, check enclosed
                        furthercheck = 0;
                        
                        if sum(imag(nzsrc) < imag(nz)) ~= 0
                            furthercheck = 1;
                        end
                        
                        if furthercheck
                            tmpT = real(nzsrc);
                            tmpb = imag(nzsrc);
                            
                            p = o.vandernewtonT(tmpT,tmpb,16);
                            
                            kk = (0:15)';
                            test = sum(p.*real(nz).^kk);
                            
                            if abs(test - imag(nz)) > 1e-14
                                if test < imag(nz) %Correct value of integral
                                    lg1 = lg1 + pi*1i;
                                    lg2 = lg2 - pi*1i;
                                end
                            end
                        end
                    end
                end
            end
            
            p0 = lg1-lg2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function b = vandernewton(~,T,b,n)
            
            for k=1:n-1
                for i=n:-1:k+1
                    b(i) = b(i) - T(k)*b(i-1);
                end
            end
            
            for k=n-1:-1:1
                for i=k+1:n
                    b(i) = b(i)/(T(i)-T(i-k));
                end
                for i=k:n-1
                    b(i) = b(i) - b(i+1);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a] = vandernewtonT(~,T,b,n)
            x = T;
            c = b;
            for k=1:n-1
                for i=n:-1:k+1
                    c(i) = (c(i)-c(i-1))/(x(i)-x(i-k));
                end
            end
            a = c;
            for k=n-1:-1:1
                for i=k:n-1
                    a(i) = a(i)-x(k)*a(i+1);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [accurate, testsum, err] = sq_necessary(~, exact_int, ztar, zsrc, zpsrc, wsrc)
            testsum = sum(zpsrc.*wsrc./(zsrc-ztar));
            err = abs(exact_int-testsum);
            accurate = err < 1e-14;
        end
    end
    
end

