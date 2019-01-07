function [ vA1, vT1, alpha1deg, vA2, vT2, alpha2deg ] = velocity_evolution_ConstantAngle( triangle_mid, r_mid, is_stator, first, second, r )
    
    
    %% side 1
    if first
        if strcmp(is_stator, 'stator')
            alpha1 = triangle_mid.alpha1;
        else alpha1 = triangle_mid.beta1;
        end
        
        sinalpha1sq = sin(alpha1)^2;
        tanalpha1 = tan(alpha1);
        vA1mid = triangle_mid.v1A; % v1A = w1A
        K1 = vA1mid * r_mid ^ sinalpha1sq;
        
        u = @(r) triangle_mid.u / r_mid * r;
        vA1 = @(r) K1 ./ (r.^sinalpha1sq);
        vT1 = @(r) vA1(r) * tanalpha1;
        alpha1deg = @(r) ones(size(r)) * rad2deg(alpha1);
        
        if(nargin == 6)
            %if strcmp(is_stator, 'stator')
                vA1 = vA1(r);
                vT1 = vT1(r);
                alpha1deg = alpha1deg(r);
            %else
            %    vA1 = vA1(r); % WA1
            %    vT1 = vT1(r) - u(r); % WT1
            %    alpha1deg = atand(vT1./vA1); % beta1
            %end
        end
        
        if ~second
            vA2 = [];
            vT2 = [];
            alpha2deg = [];
        end
    end
    
    %% side 2
    if second
        if strcmp(is_stator, 'stator')
            alpha2 = triangle_mid.alpha2;
        else alpha2 = triangle_mid.beta2;
        end
            
        sinalpha2sq = sin(alpha2)^2;
        tanalpha2 = tan(alpha2);
        vA2mid = triangle_mid.v2A;
        K2 = vA2mid * r_mid ^ sinalpha2sq;
        
        u = @(r) triangle_mid.u / r_mid * r;
        vA2 = @(r) K2 ./ (r.^sinalpha2sq);
        vT2 = @(r) vA2(r) * tanalpha2;
        alpha2deg = @(r)  ones(size(r)) * rad2deg(alpha2);
        
        
        if(nargin == 6)
            %if strcmp(is_stator, 'stator')
                vA2 = vA2(r); % or WA2
                vT2 = vT2(r); % or WT2
                alpha2deg = alpha2deg(r);
            %else
            %    vA2 = vA2(r); % WA2
            %    vT2 = vT2(r) - u(r); % WT2
            %    alpha2deg = atand(vT2./vA2); % beta2
            %end
        end
        
        if ~first
            vA1 = [];
            vT1 = [];
            alpha1deg = [];
        end
    end
    
end

