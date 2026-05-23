clear;clc;close all;

ip=imageprocessor();

img=ip.load_image('');
stars=ip.detect_stars(img);
mu = mean(img(:));
sig = std(img(:));
fprintf('mu=%.4f, sigma=%.4f\n',mu,sig);
fprintf('threshold T= %.4f\n',mu+5*sig);

fprintf('Detected Stars: %d\n', size(stars,1));
fprintf('\n%-6s %-10s %-10s %-12s\n','Star','u[px]', 'v[px]','brightness');
fprintf('%s\n',repmat('-', 1, 42)) %ripulisce la formattazione della tabella in console

for i=1:size(stars,1)
    fprintf('%-6d %-10.2f %-10.2f %-12.4f\n',i,stars(i,1),stars(i,2),stars(i,3));
end
ip.visualize(img,stars);
