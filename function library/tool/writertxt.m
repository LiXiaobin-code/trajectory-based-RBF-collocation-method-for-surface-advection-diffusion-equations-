function   writertxt(A,mypath)
[m,n]=size(A);
fid=fopen(mypath,'wt');
for i=1:m
        for j=1:n  
            if j==n
                fprintf(fid,'%5.5g \\\\ \n',A(i,j));
            elseif j==1
                fprintf(fid,'%5.6g\t & ',A(i,j));
            else
                fprintf(fid,'%5.4e\t & ',A(i,j));
            end
        end 
end      
fclose(fid);
end

