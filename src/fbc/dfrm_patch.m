%result structures, start frame1, start frame2, frame number
%raw x, y

function[S1, x, rgbArray] = dfrm_patch(X, Y, val, f1, fn)
    I = imread(['food_box_crop/food_box_8_1_',num2str(f1),'_crop.png']);
    X1 = X(:, f1:f1+fn - 1);
    Y1 = Y(:, f1:f1+fn - 1);
    val1 = val(:, f1:f1+fn - 1);
   
    dim = size(X1);
    
    %filter: only keep those features in fn frames
    [m1, n1] = find(val1(:, 2:end) ~= 0);
    m1 = unique(m1);
    
    m1 = setdiff([1: dim(1)], m1);
    t1 = size(m1);
    
    X1 = X1(m1, :);
    Y1 = Y1(m1, :);
    
    dist = (X1(:, 1) - X1(:, fn)) .^ 2 + ...,
    (Y1(:, 1) - Y1(:, fn)) .^ 2;
    moved1 = find(dist > 3);
    
    X1 = X1(moved1, :);
    Y1 = Y1(moved1, :);
    
    [S1, x, ind] = screen_hao(X1, Y1);
    y = Y1(ind, :);
    
    %just get 1 pixel
    rgbArray = [];
    sizeInd = size(ind);
    for ci = 1:sizeInd(2)
        pat = getPatch(I, round(x(ci)), round(y(ci)), 1);
        pat = reshape(pat, 3, 1);
        pat = double(pat/255);
        rgbArray = [rgbArray,pat]; 
    end
    
end

function[s, x] = subdfrm(X, Y)
    dim = size(X);

    x_h = mean(X);
    y_h = mean(Y);
    
    X_r = X - repmat(x_h, dim(1), 1);
    Y_r = Y - repmat(y_h, dim(1), 1);
    
    W = zeros(dim(2) * 2, dim(1));
    W(1:2:dim(2) * 2, :) = X_r';
    W(2:2:dim(2) * 2, :) = Y_r';
    a = size(W);
    
    [U, S, V] = svd(W);
    %V = V';
    U_ = U(:, 1:3);
    V_ = V(:, 1:3);
    S_ = S(1:3, 1:3);
    Msvd = U_ * sqrt(S_);
    Ssvd = sqrt(S_) * V_';
    
    W_ = Msvd * Ssvd;
    
    W_diff = W_ - W;
    diff_x = W_diff(1:2:dim(2) * 2, :);
    diff_y = W_diff(2:2:dim(2) * 2, :);
    
    dist = sqrt(diff_x .^ 2 + diff_y .^ 2);
    dist_max = max(dist);
    
    ind = find(dist_max < sqrt(1.38));
    
    x = X(ind, :);
    y = Y(ind, :);
    
    W_n = W(:, ind);

    %calculate svd a second time
    W_c = [W_n(1: 2: dim(2) * 2, :);W_n(2: 2: dim(2) * 2, :)];
    %2nd calc
    [U2, S2, V2] = svd(W_c);
    %V2 = V2';
    U2_ = U2(:, 1:3);
    V2_ = V2(:, 1:3);
    S2_ = S2(1:3, 1:3);
    Msvd2 = U2_ * sqrt(S2_);
    Ssvd2 = sqrt(S2_) * V2_';
    
    [m, s] = orthonormalize(Msvd2, Ssvd2);
    
    v = getv(m);
    
end

function[ind1, ind2] = sinter(X1, X2, f1, f2, fn)
    [c, ind1, ind2] = intersect(X1(:, f2 - f1 + 1), X2(:, 1));
end

function[v] = getv(R)
    dimR = size(R);
    Rr = -R(1 : dimR(1)/2, :);%result of cross(n, u)
    Ru = R((dimR(1)/2 + 1) : dimR(1), :);
    Rn = cross(Ru, Rr, 2);
    r1 = Rn(1, :);r2 = Rn(2,:);r3 = Rn(3,:);
    r12 = r1 - r2;
    r13 = r2 - r3;
    v = cross(r12, r13);
end
