% ��Wolfe-Powell ��ⲽ��
% ���غ��� f �� x �������� d ʱ�Ĳ��� alpha
% ��Ϊ��������������SteepestDescent��ConjugateGrandient����WolfePowellǰ������Щ�����ڲ����趨��
% rho�ʣ�0,1/2�� sigma�ʣ�rho��1��

function alpha = WolfePowell(alpha_,rho,sigma,f,x,d)    

phi0 = f(x); %������x��ֵ
t = sym('t',size(x)); %���÷��ű���t 
grd = gradient(f(t)); %��Ӧ���ݶ�grd
phi0_ = dot(double(subs(grd,t,x)),d); %����df(x)*d

a1 = 0;
a2 = alpha_;
phi1 = phi0;
phi1_ = phi0_;
alpha = (a1+a2)/2;

n = 0;
while(n < 1000)    %���Ƶ�������n<1000,����ʱ��̫��
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







        
