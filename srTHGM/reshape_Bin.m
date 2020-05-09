function A = reshape_Bin(Org,d,n)

    if(d==1)
        Org = Org';
    end
    
    B = reshape(Org,n,[]);
    Bc = mean(B);
    
    A = reshape(Bc,[],size(Org,2));
    
    if(d==1)
        A = A';
    end
end