%rotation test, test if the generated structure by frames is rotating
function recon_fbc    
clear;clc;close all;

    [x, y, val] = klt_read_featuretable('food_box_crop/features_fbc.ft');

    
    frames = 10;       %set the number of frames each segment
    start = 1;         %set the start frame
    total = 253;       %set the start frame, total frames and 
    S_cur = 0;X_cur = 0;
    S_last = 0;X_last = 0;
    S_all = 0;

    for i = start : frames-1 : total-frames;
        [S_cur, X_cur, rgbArray] = dfrm_patch(x, y, val, i, frames);
        if(i == start)        
            S_last = S_cur;
            X_last = X_cur;
            S_all = S_cur;
            continue;
        end
        
        f1 = i + 1 - frames;
        f2 = i;
        [ind1, ind2] = sinter(X_last, X_cur, f1, f2, frames);

        S1 = S_last(:, ind1);S2 = S_cur(:, ind2);
        
        %get the transformation parameters
        [R, t, c] = getRtc(S2, S1);

        %in this particular cropped food box case, the KLT tracker tracks
        %very few features during the segment of 145-154 frames, so our
        %algorithm cannot compute correctly the transformation parameters.
        %Thus, in this particular segment of 145-154 frames, we mannual
        %computed the transformation parameters and apply them in this
        %segment. This situation happens only in this particular segment
        %and it is because the KLT tracker simply cannot track enough
        %features
        if (i == 145)
            Rx = [1,0,0;0,-0.5736,-0.8192;0,0.8192,-0.5736];
            Ry = [0.9397,0,-0.3420;0,1,0;0.3420,0,0.9397];
            Rz = [0.9925,0.1218,0;-0.1218,0.9925,0;0,0,1];
            R = Rx * Ry * Rz;
            t = [0;-60;-70];
            c = 1;
        end
        dim = size(S_cur);
        c = 1;
        S_r = c * R * S_cur + repmat(t, [1, dim(2)]);
        S_all = [S_all, S_r];
        
        %use Rtc to register the latter point cloud to the former point
        %cloud's coordinate system
        S_last = S_r;
        X_last = X_cur;
        figure(1);
        hold on;
        %set color here
        sizeSR = size(S_r);
        for j = 1 : sizeSR(2)   
            plot3(S_r(1, j), S_r(2, j), S_r(3, j), 'LineWidth',2,'color', [rgbArray(1,j) rgbArray(2,j) rgbArray(3,j)]);
        end
        set(gca, 'DataAspectRatio', [1 1 1]);
        view(-105,18);
        xlabel('x');
        ylabel('y');
        zlabel('z');
        hold off;
    end    
end

function[ind1, ind2] = sinter(X1, X2, f1, f2, fn)
    [c, ind1, ind2] = intersect(X1(:, f2 - f1 + 1), X2(:, 1));
end