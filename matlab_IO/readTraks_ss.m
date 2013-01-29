function trk= readTraks_ss(fileName)

fid = fopen(fileName, 'r');

i = 0;
len = fscanf(fid,'%d');
while (len) > 0    
    i = i + 1;
    trk(i).ss=fscanf(fid,'[%d,%d]');
    A = fscanf(fid,'(%d,%d,%d)');
     
    trk(i).x = A(1:3:end);
    trk(i).y = A(2:3:end);
    trk(i).t = A(3:3:end);
   
    len = fscanf(fid,'%d');
end
fclose(fid);