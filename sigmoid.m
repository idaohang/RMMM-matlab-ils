function y=sigmoid(x,height,center,width,bottom)

y=height./(1+exp(x/width-center))+bottom;
end