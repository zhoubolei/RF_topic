function im = genOptTopicIm_color(A, bgim)

CV = [1 0 0; 1 0 1; 0 1 1; 0 0 1;];

C1 = repmat(reshape(CV(1,:), [1 1 3]), [48 72 1]);
C2 = repmat(reshape(CV(2,:), [1 1 3]), [48 72 1]);
C3 = repmat(reshape(CV(3,:), [1 1 3]), [48 72 1]);
C4 = repmat(reshape(CV(4,:), [1 1 3]), [48 72 1]);

% C1 = repmat(reshape(CV(1,:), [1 1 3]), [360 295 1]);
% C2 = repmat(reshape(CV(2,:), [1 1 3]), [360 295 1]);
% C3 = repmat(reshape(CV(3,:), [1 1 3]), [360 295 1]);
% C4 = repmat(reshape(CV(4,:), [1 1 3]), [360 295 1]);

% C1 = repmat(reshape(CV(1,:), [1 1 3]), [96 110 1]);
% C2 = repmat(reshape(CV(2,:), [1 1 3]), [96 110 1]);
% C3 = repmat(reshape(CV(3,:), [1 1 3]), [96 110 1]);
% C4 = repmat(reshape(CV(4,:), [1 1 3]), [96 110 1]);

% C1 = repmat(reshape(CV(1,:), [1 1 3]), [144 118 1]);
% C2 = repmat(reshape(CV(2,:), [1 1 3]), [144 118 1]);
% C3 = repmat(reshape(CV(3,:), [1 1 3]), [144 118 1]);
% C4 = repmat(reshape(CV(4,:), [1 1 3]), [144 118 1]);

% C1 = repmat(reshape(CV(1,:), [1 1 3]), [75 61 1]);
% C2 = repmat(reshape(CV(2,:), [1 1 3]), [75 61 1]);
% C3 = repmat(reshape(CV(3,:), [1 1 3]), [75 61 1]);
% C4 = repmat(reshape(CV(4,:), [1 1 3]), [75 61 1]);

% C1 = repmat(reshape(CV(1,:), [1 1 3]), [48 72 1]);
% C2 = repmat(reshape(CV(2,:), [1 1 3]), [48 72 1]);
% C3 = repmat(reshape(CV(3,:), [1 1 3]), [48 72 1]);
% C4 = repmat(reshape(CV(4,:), [1 1 3]), [48 72 1]);

% C1 = repmat(reshape(CV(1,:), [1 1 3]), [24 35 1]);
% C2 = repmat(reshape(CV(2,:), [1 1 3]), [24 35 1]);
% C3 = repmat(reshape(CV(3,:), [1 1 3]), [24 35 1]);
% C4 = repmat(reshape(CV(4,:), [1 1 3]), [24 35 1]);


%bgim = rgb2gray(bgim);
bgim = im2double(bgim);
%bgim = repmat(bgim, [1 1 3]);
%B = reshape(A, [36*2 48*2 4]);
%B = reshape(A, [144 118 4]);
%B = reshape(A, [360 295 4]);
%B = reshape(A, [96 110 4]);
%B = reshape(A, [75 61 4]);
sumI=zeros(480,720,3);
for i=1:size(A,1)
    curA=A(i,:);
     B = reshape(curA, [48 72 4]);
%B = reshape(A, [24 35 4]);
S = sum(B, 3);

b = max(S(:));

B = B/b;

I = C1.* repmat(B(:,:,1), [1 1 3]) + C2.* repmat(B(:,:,2), [1 1 3]) + C3.* repmat(B(:,:,3), [1 1 3]) + C4.* repmat(B(:,:,4), [1 1 3]);
% I = imresize(I, [360 480], 'bilinear');
%I = imresize(I, [720 591], 'bilinear');
%I = imresize(I, [477 547], 'bilinear');
%I = imresize(I, [1500 1228]);
I = imresize(I, [480 720]);
%I = imresize(I, [240 350]);
sumI=sumI+I;
end
im = 1*bgim + 1.2*sumI ;