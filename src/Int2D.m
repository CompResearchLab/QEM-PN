function I = Int2D(f,dA)
    f(:,1)=f(:,1)/2;
    f(:,end)=f(:,end)/2;
    f(1,:)=f(1,:)/2;
    f(end,:)=f(end,:)/2;
    f(1,1)=f(1,1)/2;
    f(end,end)=f(end,end)/2;
    f(1,end)=f(1,end)/2;
    f(end,1)=f(end,1)/2;
    I=sum(f*dA,"all");
end

