function count = mandelbrot_count(maxIterations, xGrid, yGrid) %#codegen

z0 = xGrid + 1i*yGrid;
count = ones(size(z0));

% Map computation to GPU.
coder.gpu.kernelfun;

z = z0;
for n = 0:maxIterations
    z = z.*z + z0;
    inside = abs(z)<=2;
    count = count + inside;
end
count = log(count);