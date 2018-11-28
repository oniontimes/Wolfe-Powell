% 用Wolfe-Powell 求解步长
% 返回函数 f 在 x 处，方向 d 时的步长 alpha
% 因为在其他函数，如SteepestDescent或ConjugateGrandient调用WolfePowell前就在那些函数内部已设定过
% rho∈（0,1/2） sigma∈（rho，1）

function alpha = WolfePowell(alpha_,rho,sigma,f,x,d)    

phi0 = f(x); %函数在x处值
t = sym('t',size(x)); %设置符号变量t 
grd = gradient(f(t)); %相应的梯度grd
phi0_ = dot(double(subs(grd,t,x)),d); %计算df(x)*d

a1 = 0;
a2 = alpha_;
phi1 = phi0;
phi1_ = phi0_;
alpha = (a1+a2)/2;

n = 0;
while(n < 1000)    %限制迭代上限n<1000,避免时间太长
    phi = f(x+alpha*d);
    if phi <= phi0 + rho*alpha*phi0_
        phi_ = dot(double(subs(grd,t,x+alpha*d)),d);
        if phi_ >= sigma*phi0_
            break;
        else
            alpha_new = alpha - (a1-alpha)*phi_/(phi1_-phi_);
            a1 = alpha;
            alpha = alpha_new;
            phi1 = phi;
            phi1_ = phi_;
        end
    else
        alpha_new = a1 + 0.5*(a1-alpha)*(a1-alpha)*phi1_/((phi1-phi)-(a1-alpha)*phi1_);
        a2 = alpha;
        alpha = alpha_new;
    end
    n = n + 1;
end

end







        
