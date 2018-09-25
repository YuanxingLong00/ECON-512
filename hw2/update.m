function [xNew,yNew] =update(x, y)
xNew=1/( 1- exp(2-x)/( 1+ exp(2-x)+exp(2-y) ));
yNew=1/( 1- exp(2-y)/( 1+ exp(2-x)+exp(2-y) ));
end
