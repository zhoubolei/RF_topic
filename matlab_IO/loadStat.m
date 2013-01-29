function trkStat= loadStat(fileName,nTopics,nTrks)
fid = fopen(fileName, 'r');

trkStat=zeros(nTrks,nTopics);
for i=1:nTrks   
    index= fscanf(fid,'%d');
    A = fscanf(fid,'(%d)');
    
        
    trkStat(index,:)=A;
   
    

end
fclose(fid);