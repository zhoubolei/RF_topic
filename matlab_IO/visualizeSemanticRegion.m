clear
clc

sizeImg=[480 720];
curFile='\';
topics=readTopic([curFile 'classqq.txt'],sizeImg); 
% read the topic distribution learned from Random Field Topic Model.

bg=imread('bgcrowdlarge.png');
bg=imresize(bg,[480 720]);

bg=im2double(bg);
for i=1:size(topics,1)
    curImg=genOptTopicIm_color(topics(i,:),bg);
    imwrite(curImg,[num2str(i) '.jpg'], 'jpg');% output the semantic regions
end

