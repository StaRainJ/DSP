% VL_DEMO_SIFT_MATCH  Demo: SIFT: basic matching
close all;
pfx = fullfile(vl_root,'figures','demo') ;
randn('state',0) ;
rand('state',0) ;
figure(1) ; clf ;

% --------------------------------------------------------------------
%                                                    Create image pair
% --------------------------------------------------------------------
% 
% Ia = imread(fullfile(vl_root,'data','bk2.jpg')) ;
% Ib = imread(fullfile(vl_root,'data','bk3.jpg')) ;
% Ia = imread(fullfile(vl_root,'data','fox3.jpg')) ;
% Ib = imread(fullfile(vl_root,'data','fox2.jpg')) ;
% Ia = imread(fullfile(vl_root,'data','1.jpg')) ;
% Ib = imread(fullfile(vl_root,'data','1s.jpg')) ;
% aa=size(Ia);bb=size(Ib);
% Ia=imresize(Ia,[max(aa(1),bb(1)),max(aa(2),bb(2))]);
% Ib=imresize(Ib,[max(aa(1),bb(1)),max(aa(2),bb(2))]);
Ia = imread(fullfile(vl_root,'data/image','25_l.JPG')) ;
Ib = imread(fullfile(vl_root,'data/image','25_r.JPG')) ;
thresh=1.5;
th = pi/4 ;
sc = 4 ;
c = sc*cos(th) ;
s = sc*sin(th) ;
A = [c -s; s c] ;
T = [- size(Ia,2) ; - size(Ia,1)]  / 2 ;

tform = maketform('affine', [A, A * T - T ; 0 0  1]') ;
% Ib = imtransform(Ia,tform,'size',size(Ia), ...
%                  'xdata', [1 size(Ia,2)], ...
%                  'ydata', [1 size(Ia,1)], ...
%                  'fill', 255);
  


% --------------------------------------------------------------------
%                                           Extract features and match
% --------------------------------------------------------------------

[fa,da] = vl_sift(im2single(rgb2gray(Ia))) ;
[fb,db] = vl_sift(im2single(rgb2gray(Ib))) ;

[matches, scores] = vl_ubcmatch(da,db,thresh) ;

[drop, perm] = sort(scores, 'descend') ;
matches = matches(:, perm) ;
scores  = scores(perm) ;
xa = fa(1,matches(1,:)) ;
ya = fa(2,matches(1,:)) ;
xb = fb(1,matches(2,:)) + size(Ia,2) ;
yb = fb(2,matches(2,:)) ;

X=[xa;ya]';
Y=[xb-size(Ia,2);yb]';

figure(1) ; clf ;
imagesc(cat(2, Ia, Ib)) ;
hold on ;
%%
h = line([xa ; xb], [ya ; yb]) ;
set(h,'linewidth', 1, 'color', 'b') ;

