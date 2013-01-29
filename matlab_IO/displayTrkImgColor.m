function curImg=displayTrkImgColor(trks,bg)
%DISPLAYTRKIMG Summary of this function goes here
%   Detailed explanation goes here
curImg=bg;
%curImg=repmat(curImg,[1 1 3]);curImg=curImg*0.5;
hold off
%imshow(curImg);
for i=1:length(trks)
   color=rand(3,1);
   for j=1:length(trks(i).y)
       
        curImg(trks(i).y(j),trks(i).x(j),:)=color(:,1);
        curImg(trks(i).y(j)+1,trks(i).x(j)+1,:)=color(:,1);
        curImg(trks(i).y(j)+1,trks(i).x(j),:)=color(:,1);

   end

end

end
