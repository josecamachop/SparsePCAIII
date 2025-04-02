function CS = cossim(A,B) 

    pcs = 1:size(A,2); 
    CS = 0; 
    for i = pcs 
        CS = CS + (A(:,i)'*B(:,i))/(norm(A(:,i),2)*norm(B(:,i),2)); 
    end 
    CS = 1/length(pcs) * CS;
    
end