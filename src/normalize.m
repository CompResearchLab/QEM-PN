function A=normalize(Z,dA);
      SqrZ=Z.^2;
      A=1/sqrt(Int2D(SqrZ,dA));
end

