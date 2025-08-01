
function [L, det, error] = cholesky_regressaoLinear(matrix)
   
    n = leght(matrix);
    det=1;
    L=zeros(n,n);

    for j=1:n
        soma=0;
        for k = 1:j-1
            soma = soma + L(j,k)^2;
        end
        t=matrix(j,j)-soma;
        det=det*t;
        error = t<=0;
        if error==1
            return;
            %prompt = "A matriz não é definida positiva";
            %error(prompt);
        else
            L(j,j)=sqrt(t);
            r=1/L(j,j);
        end
        for i=j+1:n
            soma = 0;
            for k=1:j-1
                soma = soma+L(i,k)*L(j,k);
            end
            L(i,j)=(matrix(i,j)-soma)*r;
        end
    end
end
