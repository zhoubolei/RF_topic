clear
%% read tracklets from files
trks=[];
for i=1:16
    curfile=['..\vc_RTM\release\filttracks_grand_' num2str(i) '.txt'];
    trks=[trks readTraks_ss(curfile)];
    i
end
% read the background image
bg=imread('bgcrowdlarge.png');
bg=im2double(bg);

%% visualize all the tracklets.
nTopic=30;
nTrks=length(trks);
curImg=displayTrkImgColor(trks,bg);
figure,imshow(curImg);
