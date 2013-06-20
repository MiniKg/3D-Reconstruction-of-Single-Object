function [patch] = getPatch(image,x,y,k)

s = size(image);
width = (k - 1)/2;

if x - width <= 0 || x + width >= s(2) || y - width <= 0 || y + width >= s(1)
    patch = image(y,x,:);
else
    patch = image(y-width:y+width,x-width:x+width,:);
end

end