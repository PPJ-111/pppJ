# pppJ
verbosity=0.;
macro UgradV(u1,u2,v1,v2) [[u1,u2]'*[dx(v1),dy(v1)],[u1,u2]'*[dx(v2),dy(v2)]]  //

//macro Grad(u1,u2) [dx(u1),dy(u1),dx(u2),dy(u2)]  //

macro f1(t) ( 
 
cos(sin(t))*(100*x^2*y^4*cos(t)^2*(2*x - 1)^2*(x - 1)^2*(y - 1)^4 + 100*x^4*y^2*cos(t)^2*(2*y - 1)^2*(x - 1)^4*(y - 1)^2)*(x*cos(sin(t)) + y*sin(sin(t)) + 2) - (10*x^2*y*cos(t)*(2*y - 1)*(x*sin(sin(t))*cos(t) - y*cos(sin(t))*cos(t))*(x - 1)^2*(y - 1) + 10*x^2*y*sin(t)*(2*y - 1)*(x - 1)^2*(y - 1)*(x*cos(sin(t)) + y*sin(sin(t)) + 2))*(x*cos(sin(t)) + y*sin(sin(t)) + 2) - 40*x^2*cos(t)*(x - 1)^2*(y - 1) - 20*x^2*cos(t)*(2*y - 1)*(x - 1)^2 - 40*x^2*y*cos(t)*(x - 1)^2 - 20*x^2*y*cos(t)*(2*y - 1)*(y - 1) - 20*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) - 40*x*y*cos(t)*(2*x - 2)*(2*y - 1)*(y - 1) - 10*x*y^2*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2*(10*x^2*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) + 20*x^2*y*cos(t)*(x - 1)^2*(y - 1) + 10*x^2*y*cos(t)*(2*y - 1)*(x - 1)^2)*(x*cos(sin(t)) + y*sin(sin(t)) + 2)^2 + 10*x^2*y*cos(t)*(2*y - 1)*(20*x*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) + 10*x^2*y*cos(t)*(2*x - 2)*(2*y - 1)*(y - 1))*(x - 1)^2*(y - 1)*(x*cos(sin(t)) + y*sin(sin(t)) + 2)^2 - 5*x^2*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1)*(x*cos(sin(t)) + y*sin(sin(t)) + 2)^2*(20*x*y*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 - 20*x*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) + 10*x*y^2*cos(t)*(2*x - 1)*(2*y - 2)*(x - 1) - 10*x^2*y*cos(t)*(2*x - 2)*(2*y - 1)*(y - 1))
 
 ) //

macro  f2(t)  ( 
  
 (10*x*y^2*cos(t)*(2*x - 1)*(x*sin(sin(t))*cos(t) - y*cos(sin(t))*cos(t))*(x - 1)*(y - 1)^2 + 10*x*y^2*sin(t)*(2*x - 1)*(x - 1)*(y - 1)^2*(x*cos(sin(t)) + y*sin(sin(t)) + 2))*(x*cos(sin(t)) + y*sin(sin(t)) + 2) + sin(sin(t))*(100*x^2*y^4*cos(t)^2*(2*x - 1)^2*(x - 1)^2*(y - 1)^4 + 100*x^4*y^2*cos(t)^2*(2*y - 1)^2*(x - 1)^4*(y - 1)^2)*(x*cos(sin(t)) + y*sin(sin(t)) + 2) + 40*y^2*cos(t)*(x - 1)*(y - 1)^2 + 20*y^2*cos(t)*(2*x - 1)*(y - 1)^2 + 40*x*y^2*cos(t)*(y - 1)^2 + 20*x*y^2*cos(t)*(2*x - 1)*(x - 1) + 20*x*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 + 40*x*y*cos(t)*(2*x - 1)*(2*y - 2)*(x - 1) - 10*x^2*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1)*(10*y^2*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 + 20*x*y^2*cos(t)*(x - 1)*(y - 1)^2 + 10*x*y^2*cos(t)*(2*x - 1)*(y - 1)^2)*(x*cos(sin(t)) + y*sin(sin(t)) + 2)^2 + 10*x*y^2*cos(t)*(2*x - 1)*(20*x*y*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 + 10*x*y^2*cos(t)*(2*x - 1)*(2*y - 2)*(x - 1))*(x - 1)*(y - 1)^2*(x*cos(sin(t)) + y*sin(sin(t)) + 2)^2 + 5*x*y^2*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2*(x*cos(sin(t)) + y*sin(sin(t)) + 2)^2*(20*x*y*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 - 20*x*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) + 10*x*y^2*cos(t)*(2*x - 1)*(2*y - 2)*(x - 1) - 10*x^2*y*cos(t)*(2*x - 2)*(2*y - 1)*(y - 1))
 
 )//

macro u1e(t)  (  10*x^2*(x-1)^2*y*(y-1)*(2*y-1)*cos(t)  )  //

macro u2e(t)  (  -10*x*(x-1)*(2*x-1)*y^2*(y-1)^2*cos(t)  )  //
macro sige(t) (   2+x*cos(sin(t))+y*sin(sin(t))         )   //

macro g(t) (  

(20*x*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) + 10*x^2*y*cos(t)*(2*x - 2)*(2*y - 1)*(y - 1))*((x*cos(sin(t)))/2 + (y*sin(sin(t)))/2 + 1) - (20*x*y*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 + 10*x*y^2*cos(t)*(2*x - 1)*(2*y - 2)*(x - 1))*((x*cos(sin(t)))/2 + (y*sin(sin(t)))/2 + 1) - x*sin(sin(t))*cos(t) + y*cos(sin(t))*cos(t) + 10*x^2*y*cos(sin(t))*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) - 10*x*y^2*sin(sin(t))*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2
 
)//

macro s(t)(  20*x*y*cos(t)*(2*y - 1)*(x - 1)^2*(y - 1) - 20*x*y*cos(t)*(2*x - 1)*(x - 1)*(y - 1)^2 - 10*x*y^2*cos(t)*(2*x - 1)*(2*y - 2)*(x - 1) + 10*x^2*y*cos(t)*(2*x - 2)*(2*y - 1)*(y - 1) 

  )  //


real nn,nnn,t,dt,T=0.2,nu=1.;
real[int] L2error(4), L2RRerror(4),h(4);

//func re=2+x*cos(sin(t))+y*sin(sin(t));
//func pe=sin(x)*sin(y)*sin(t);

for(int n=1;n<4;n++)
{
nn=2^(n+3);
dt=1/(nn^1.5);//
//t=0;
nnn=T/dt;

mesh Th=square(nn,nn);
fespace Wh(Th,P2);   Wh sig0,sig1,sig2,ah;   //
fespace Vh(Th,P2);  Vh v1,v2,uold1,uold2,uh11,uh12,u1,u2;//
fespace Ph(Th,P1);   Ph p1,p,q;
// v1 v2 q ah 是test function
uold1=  10*x^2*(x-1)^2*y*(y-1)*(2*y-1);
uold2= -10*x*(x-1)*(2*x-1)*y^2*(y-1)^2 ;//真解不好？？？？？？？？？
sig0=2.0+x; //开根号？
//rh0=2+x*cos(sin(t))+y*sin(sin(t))
t=dt;//////

//算sig1
solve SIG(sig1,ah,solver=sparsesolver)=
                  int2d(Th)(sig1*ah/dt)
                - int2d(Th)(sig0*ah/dt)
                + int2d(Th)(dx(sig1)*uold1*ah+dy(sig1)*uold2*ah)
                - int2d(Th)(g(t)*ah)      //g(t)要重新算
                + int2d(Th)(0.5*(sig1*dx(uold1)*ah+sig1*dy(uold2)*ah) );
               // + on(1,2,3,4,rh1=re);

              
solve NS([uh11,uh12,p1],[v1,v2,q],solver=UMFPACK) =
      int2d(Th)(  sig1*sig1*uh11*v1/dt+sig1*sig1*uh12*v2/dt)//TT

    + int2d(Th)(  nu *(dx(uh11)*dx(v1) + dy(uh11)*dy(v1)+ dx(uh12)*dx(v2) + dy(uh12)*dy(v2))) //TT
    + int2d(Th)(  sig1*sig1*  UgradV(uold1,uold2,uh11,uh12)  '*[v1,v2] ) //TT

    //+ int2d(Th)( rh1*uold1*dx(uh11)*v1 +rh1*uold2*dy(uh11)*v1 +rh1*uold1*dx(uh12)*v2+ rh1*uold2*dy(uh12)*v2 )
    //+ int2d(Th)(  0.5* ( rh1*uh11*v1/dt+rh1*uh12*v2/dt-rh0*uh11*v1/dt-rh0*uh12*v2/dt ) ) //T
   
    + int2d(Th)(  0.5*2.0*sig1*dx(sig1)*uold1*uh11*v1 + 0.5*2.0*sig1*dy(sig1)*uold2*uh11*v1
                + 0.5*dx(uold1)*sig1*sig1*uh11*v1 + 0.5*dy(uold2)*sig1*sig1*uh11*v1
                + 0.5*2.0*sig1*dx(sig1)*uold1*uh12*v2 + 0.5*2.0*sig1*dy(sig1)*uold2*uh12*v2
                + 0.5*dx(uold1)*sig1*sig1*uh12*v2 + 0.5*dy(uold2)*sig1*sig1*uh12*v2   )//TT
    - int2d(Th)((dx(v1)+dy(v2))*p1)//散度 TT
    + int2d(Th)((dx(uh11)+dy(uh12))*q)//TT
    + int2d(Th)((1e-10)*p1*q)// TT
    //- int2d(Th)(0.5*g(t)*(uh11*v1+uh12*v2))
    - int2d(Th)(f1(t)*v1+f2(t)*v2)//f1(t) f2(t) 要重新算
    - int2d(Th)(s(t)*q)// s(t) 要重新算
    - int2d(Th)(sig1*sig0*uold1*v1/dt+sig1*sig0*uold2*v2/dt)//TT
    + on(1,2,3,4,uh11=0,uh12=0);    

//解出来 sig1 uh11 uh12 
 

problem SIG1(sig2,ah,solver=sparsesolver)=
                  int2d(Th)(1.5*sig2*ah/dt)
                - int2d(Th)(2.0*sig1*ah/dt)
                + int2d(Th)(0.5*sig0*ah/dt)
                + int2d(Th)(2.0*uh11*dx(sig2)*ah+2.0*uh12*dy(sig2)*ah)
                - int2d(Th)(uold1*dx(sig2)*ah+uold2*dy(sig2)*ah)
                - int2d(Th)(g(t)*ah) //g(t)重新算
                + int2d(Th)(sig2*dx(uh11)*ah+sig2*dy(uh12)*ah)
                - int2d(Th)(0.5*sig2*dx(uold1)*ah+0.5*sig2*dy(uold2)*ah);


problem NS1([u1,u2,p],[v1,v2,q],solver=UMFPACK) =
      int2d(Th)( 1.5*sig2*sig2*u1*v1/dt+1.5*sig2*sig2*u2*v2/dt )
    - int2d(Th)( 2.0*sig2*sig1*uh11*v1/dt+2.0*sig2*sig1*uh12*v2/dt)
    + int2d(Th)( 0.5*sig2*sig0*uold1*v1/dt+0.5*sig2*sig0*uold2*v2/dt)//TT

    - int2d(Th)((dx(v1)+dy(v2))*p)//散度 TT
    + int2d(Th)((dx(u1)+dy(u2))*q)//TT
    + int2d(Th)(  nu *(dx(u1)*dx(v1) + dy(u1)*dy(v1)+ dx(u2)*dx(v2) + dy(u2)*dy(v2))) //TT
   
    + int2d(Th)(  2.0*sig2*sig2*  UgradV(uh11,uh12,u1,u2)  '*[v1,v2] ) //TT
    - int2d(Th)(  sig2*sig2*  UgradV(uold1,uold2,u1,u2)  '*[v1,v2] ) //TT

   // + int2d(Th)(  0.75*rh2*u1*v1/dt+0.75*rh2*u2*v2/dt)
   // - int2d(Th)(  rh1*u1*v1/dt+rh1*u2*v2/dt)
   // + int2d(Th)(  0.25*rh0*u1*v1/dt+0.25*rh0*u2*v2/dt)

    + int2d(Th)(  2.0*sig2*dx(sig2)*uh11*u1*v1 +2.0*sig2*dx(sig2)*uh11*u2*v2
                + 2.0*sig2*dy(sig2)*uh12*u1*v1 +2.0*sig2*dy(sig2)*uh12*u2*v2
                - 0.5*2.0*sig2*dx(sig2)*uold1*u1*v1 - 0.5*2.0*sig2*dx(sig2)*uold1*u2*v2
                - 0.5*2.0*sig2*dy(sig2)*uold2*u1*v1 - 0.5*2.0*sig2*dy(sig2)*uold2*u2*v2    )//TT

    + int2d(Th)(  sig2*sig2*dx(uh11)*u1*v1 + sig2*sig2*dx(uh11)*u2*v2
                + sig2*sig2*dy(uh12)*u1*v1 + sig2*sig2*dy(uh12)*u2*v2
                - 0.5*sig2*sig2*dx(uold1)*u1*v1 - 0.5*sig2*sig2*dx(uold1)*u2*v2
                - 0.5*sig2*sig2*dy(uold2)*u1*v1 - 0.5*sig2*sig2*dy(uold2)*u2*v2    )//TT
    + int2d(Th)((1e-10)*p*q)
   // - int2d(Th)( 0.5*g(t)*(u1*v1+u2*v2))
    - int2d(Th)(f1(t)*v1+f2(t)*v2)
    - int2d(Th)(s(t)*q)
    + on(1,2,3,4,u1=0,u2=0);


for (int m=1;m<=nnn-1;m++)
  {

       t=t+dt;      
       SIG1;
       NS1;
       sig0=sig1;
       sig1=sig2; 
       uold1=uh11;
       uold2=uh12;     
       uh11=u1;
       uh12=u2;
       //t=t+dt;
      // rh0=rh1;
     //cout<<"t="<<t<<endl;  
    }
  L2error[n-1] =sqrt(int2d(Th)((u1e(t)-uh11)^2+(u2e(t)-uh12)^2));
  L2RRerror[n-1] =sqrt(int2d(Th)(sige(t)-sig1)^2);

  h[n-1]=nn;
    cout << " n = "<< n <<endl;

 }
 for(int n=1;n<4;n++)
  {
   cout << " L2error " << n << " = "<< L2error[n-1] <<endl;
   cout << " L2RRerror " << n << " = "<< L2RRerror[n-1] <<endl;

 
    }
 for(int n=1;n<3;n++)
   {
   cout <<" L2 convergence rate "<< n << " = "<< log(L2error[n-1]/L2error[n])/log(h[n]/h[n-1]) <<endl;
   cout <<" L2RR convergence rate "<< n << " = "<< log(L2RRerror[n-1]/L2RRerror[n])/log(h[n]/h[n-1]) <<endl;

    }
