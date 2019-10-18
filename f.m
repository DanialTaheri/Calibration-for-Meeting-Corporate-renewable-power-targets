function [x] = f(sigmaSq, beta, b, Pt, Pt_1)
N=3;
tildea=2*b*beta/sigmaSq;
betatilde= 2*b*(1-beta)/sigmaSq;
P=zeros(N,length(Pt));
P1=zeros(N,length(Pt_1));
h=zeros(N,length(Pt));
h1=zeros(N,length(Pt_1));
   lambda=zeros(N,1);
   ply_func_scl=zeros(N,1);
   S=zeros(N,length(Pt));
   x=zeros(length(Pt),1);
for n=1:N
    for m=0:n
        h(n,:)=h(n,:)+(-1)^m * (factorial(n)/(factorial(m)*factorial(n-m)))* ...
            (gamma(tildea+betatilde+n+m-1)/gamma(tildea+m)) * Pt'.^m;
       h1(n,:)=h1(n,:)+(-1)^m * (factorial(n)/(factorial(m)*factorial(n-m)))* ...
            (gamma(tildea+betatilde+n+m-1)/gamma(tildea+m)) * Pt_1'.^m;
          %display(h);
    end
   % display(h);
   lambda(n)=-b*n-0.5* sigmaSq*n*(n-1); 
    ply_func_scl(n)=((gamma(tildea+n)*(2*n+tildea+betatilde-1)*gamma(tildea)*gamma(betatilde))/...
       ( factorial(n)*gamma(tildea+betatilde+n-1)*gamma(tildea+betatilde)*gamma(betatilde+n)))^0.5;
  % ply_func_scl_n_1=((gamma(tildea+n-1)*(2*(n-1)+tildea+betatilde-1)*gamma(tildea)*gamma(betatilde))/...
     %  ( factorial(n-1)*gamma(tildea+betatilde+n-2)*gamma(tildea+betatilde)*gamma(betatilde+(n-1))))^0.5;
    P(n,:)= ply_func_scl(n).* h(n,:);
    P1(n,:)=ply_func_scl(n).* h1(n,:); 
    S(n,: )=exp(lambda(n))*P(n,:).*P1(n,:);
end
    

for j=1:length(Pt)
    for i=1:N
x(j,1)=x(j,1)+S(i,j);
    end
end
 %   display(x(j,1));
  % if x(j,1)>10000000
   %    display('error');
   %end
    %  if x(j,1)<-10000000
     %  display('error');
   %end

   
 end

