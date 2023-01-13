#include "AMReX_Reduce.H"
#include <ReactorBDF.H>
#define NEQNS (NUM_SPECIES+1)
#define NEQNS2 NEQNS*NEQNS
#define KSPSIZE NEQNS
#define NOPC 0
#define GJPC 1
#define SGSPC 2
#define NEARZERO 1e-15

//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void noprecond(amrex::Real MinvX[NEQNS],
                                              amrex::Real A[NEQNS][NEQNS],
                                              amrex::Real X[NEQNS])
{
    for(int i=0;i<NEQNS;i++)
    {
        MinvX[i]=X[i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void GJprecond(amrex::Real MinvX[NEQNS],
                                               amrex::Real A[NEQNS][NEQNS],
                                               amrex::Real X[NEQNS])
{
    //solve Dy = X
    for(int i=0;i<NEQNS;i++)
    {
        MinvX[i]=X[i]/A[i][i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void SGSprecond(amrex::Real MinvX[NEQNS],
                                               amrex::Real A[NEQNS][NEQNS],
                                               amrex::Real X[NEQNS])
{
    amrex::Real sum;
    amrex::Real y_d[NEQNS],y_dd[NEQNS];
    amrex::Real diag_inv[NEQNS];
    amrex::Real diag[NEQNS];

    //solve (D+L)D^-1(D+U)y = X
    //MinvX = y

    for(int i=0;i<NEQNS;i++)
    {
        diag[i]=A[i][i];
        diag_inv[i]=1.0/diag[i];
    }


    //solve (D+L) y' = X
    for(int i=0;i<NEQNS;i++)
    {
        sum  = 0.0;
        for(int j=0;j<i;j++)
        {
            sum += A[i][j]*y_d[j];
        }
        y_d[i] = (X[i]-sum)*diag_inv[i];
    } 

    //solve D^(-1)y'' = y' = (D+L)^(-1) X
    for(int i=0;i<NEQNS;i++)
    {
        y_dd[i]  = diag[i] * y_d[i];
    }

    //solve (D+U) y = y'' = D (D+L)^(-1) X
    for(int i=NEQNS-1;i>=0;i--)
    {
        sum=0.0;
        for(int j=NEQNS-1;j>i;j--)
        {
            sum += A[i][j]*MinvX[j];
        }

        MinvX[i] = (y_dd[i]-sum)*diag_inv[i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void precond(amrex::Real MinvX[NEQNS],
                                            amrex::Real A[NEQNS][NEQNS],
                                            amrex::Real X[NEQNS],int type)
{
    if(type==NOPC)
    {
        noprecond(MinvX,A,X);   
    }
    else if(type==GJPC)
    {
        GJprecond(MinvX,A,X);
    }
    else if(type==SGSPC)
    {
        SGSprecond(MinvX,A,X);
    }
    else
    {
        amrex::Abort("Unknown preconditioner type in BDF");
    }

}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void findAX(amrex::Real AX[NEQNS],
                                           amrex::Real A[NEQNS][NEQNS],
                                           amrex::Real X[NEQNS])
{
    for(int i=0;i<NEQNS;i++)
    {
        AX[i]=0.0;
        for(int j=0;j<NEQNS;j++)
        {
            AX[i] += A[i][j]*X[j];
        }
    }
}
//========================================================================
#ifndef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
printvec(std::string str,amrex::Real v[NEQNS])
{
    amrex::Print()<<"\n"<<str<<"\t";
    for(int i=0;i<NEQNS;i++)
    {
        amrex::Print()<<v[i]<<"\t";
    }
    amrex::Print()<<"\n";

}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
printmat(std::string str,amrex::Real mat[NEQNS][NEQNS])
{
    amrex::Print()<<"\n"<<str<<"\n";

    for(int i=0;i<NEQNS;i++)
    {
        for(int j=0;j<NEQNS;j++)
        {
            amrex::Print()<<mat[i][j]<<"  ";
        }
        amrex::Print()<<"\n";
    }
    amrex::Print()<<"\n";

}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
print_hessmat(std::string str,amrex::Real hmat[KSPSIZE+1][KSPSIZE])
{
    amrex::Print()<<"\n"<<str<<"\n";

    for(int i=0;i<KSPSIZE+1;i++)
    {
        for(int j=0;j<KSPSIZE;j++)
        {
            amrex::Print()<<hmat[i][j]<<"  ";
        }
        amrex::Print()<<"\n";
    }
    amrex::Print()<<"\n";

}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
print_arnoldivecs(std::string str,amrex::Real arnoldivecs[KSPSIZE+1][NEQNS])
{
    amrex::Print()<<"\n"<<str<<"\n";

    for(int i=0;i<(KSPSIZE+1);i++)
    {
        for(int j=0;j<NEQNS;j++)
        {
            amrex::Print()<<arnoldivecs[i][j]<<"  ";
        }
        amrex::Print()<<"\n";
    }
    amrex::Print()<<"\n";

}
#endif
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void scalevector(amrex::Real v1[NEQNS],
                                                     amrex::Real factor)
{
    for(int i=0;i<NEQNS;i++)
    {
        v1[i] = v1[i]*factor;
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
multiplyvectors(amrex::Real v1[NEQNS],
                amrex::Real v2[NEQNS],
                amrex::Real v12[NEQNS])
{
    for(int i=0;i<NEQNS;i++)
    {
        v12[i]=v1[i]*v2[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real findnorm(amrex::Real v1[NEQNS])
{
    amrex::Real norm=0.0;
    for(int i=0;i<NEQNS;i++)
    {
        norm=norm+v1[i]*v1[i];
    }

    return(sqrt(norm));
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real 
innerproduct(amrex::Real v1[NEQNS],amrex::Real v2[NEQNS])
{
    amrex::Real innerprod=0.0;
    for(int i=0;i<NEQNS;i++)
    {
        innerprod=innerprod+v1[i]*v2[i];
    }

    return(innerprod);
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void addvectors(amrex::Real v1[NEQNS],
                                               amrex::Real v2[NEQNS],
                                               amrex::Real v12[NEQNS],
                                               amrex::Real a,amrex::Real b)
{
    for(int i=0;i<NEQNS;i++)
    {
        v12[i]=a*v1[i]+b*v2[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void copyvector(amrex::Real v1[NEQNS],
                                               amrex::Real v2[NEQNS]) //(dest,source,size)
{
    for(int i=0;i<NEQNS;i++)
    {
        v1[i]=v2[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void getkspvector(amrex::Real v1[NEQNS],
                        amrex::Real kspvectors[KSPSIZE+1][NEQNS],int vecnum)
{
    for(int i=0;i<NEQNS;i++)
    {
        v1[i]=kspvectors[vecnum][i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void setkspvector(amrex::Real vec[NEQNS],
                          amrex::Real kspvectors[KSPSIZE+1][NEQNS],int vecnum)
{
    for(int i=0;i<NEQNS;i++)
    {
        kspvectors[vecnum][i]=vec[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool isit_NaN(amrex::Real val)
{
    amrex::Real valbyval;
    amrex::Real valplus3;

    valbyval=val/val;
    valplus3=val+3.0;

    if(val != val)
    {
        return(true);
    }
    else if((valbyval != valbyval) && (valplus3 == val))
    {
        return(true);
    }
    else
    {
        return(false);
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool arnoldi(int j,
                                                 amrex::Real mat[NEQNS][NEQNS],
                                                 amrex::Real kspvectors[KSPSIZE+1][NEQNS],
                                                 amrex::Real Hessmat[KSPSIZE+1][KSPSIZE],
                                                 int precond_type,
                                                 amrex::Real tol,
                                                 bool &nanflag)
{
    int i,n;
    amrex::Real Avj[NEQNS],vj[NEQNS],vi[NEQNS];
    amrex::Real wj[NEQNS],tempvec[NEQNS];
    amrex::Real MinvAvj[NEQNS];
    

    bool lucky; //when norm becomes 0,
    //KSP vectors is no longer linearly independent.
    //we would have gotten the best solution.

    lucky   = false;
    nanflag = false;

    getkspvector(vj,kspvectors,j);
    findAX(Avj,mat,vj);
    precond(MinvAvj,mat,Avj,precond_type);
    copyvector(Avj,MinvAvj);
    //Avj is now M^-1 A vj
    //remember we are solving M^-1 A X = M^-1 b

    for(i=0;i<=j;i++)
    {
        getkspvector(vi,kspvectors,i);
        Hessmat[i][j]=innerproduct(Avj,vi);
    }

    copyvector(wj,Avj);

    for(i=0;i<=j;i++)
    {
        getkspvector(vi,kspvectors,i);
        addvectors(wj,vi,tempvec,1.0,-Hessmat[i][j]);
        copyvector(wj,tempvec);
    }

    Hessmat[j+1][j] = findnorm(wj);

    if(Hessmat[j+1][j] > tol)
    {
        for(int i=0;i<NEQNS;i++)
        {
            wj[i]=wj[i]/Hessmat[j+1][j];
        }
        setkspvector(wj,kspvectors,j+1);
    }
    else
    {
        if( !isit_NaN(Hessmat[j+1][j]) )
        {
#ifndef AMREX_USE_GPU
            amrex::Print()<<"lucky\n";
#endif
            lucky=true;
        }
        else
        {
            amrex::Abort("NaN or Inf detected in Hessenberg matrix\n");
            nanflag=true;
        }
    }		

    return(lucky);
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
givens_rotate(amrex::Real r,amrex::Real s,amrex::Real &cos,
              amrex::Real &sin)
{
    if(r==0)
    {
        cos=0.0;
        sin=1.0;
    }
    else
    {
        cos=r/std::sqrt(r*r+s*s);
        sin=s/std::sqrt(r*r+s*s);
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real 
triangularize(amrex::Real Hessmat[KSPSIZE+1][KSPSIZE],
              amrex::Real cos_arr[KSPSIZE],
              amrex::Real sin_arr[KSPSIZE],
              amrex::Real beta_e1[KSPSIZE+1],
                int k)
{

  amrex::Real temp,cos,sin;

  //apply for ith column
  for(int i=0;i<=(k-1); i++)
  {
    temp   =  cos_arr[i] * Hessmat[i][k] + sin_arr[i] * Hessmat[i+1][k];
    Hessmat[i+1][k] = -sin_arr[i] * Hessmat[i][k] + cos_arr[i] * Hessmat[i + 1][k];
    Hessmat[i][k]   = temp;
  }

  //update the next sin cos values for rotation
  givens_rotate(Hessmat[k][k], Hessmat[k+1][k], cos, sin);

  //eliminate H(i + 1, i)
  Hessmat[k][k] = cos * Hessmat[k][k] + sin * Hessmat[k+1][k];
  Hessmat[k+1][k] = 0.0;

  cos_arr[k]=cos;
  sin_arr[k]=sin;

  beta_e1[k+1] = -sin*beta_e1[k];
  beta_e1[k]   =  cos*beta_e1[k];

  amrex::Real error=amrex::Math::abs(beta_e1[k+1]);

  return(error);
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
uppertrisolve(amrex::Real Hessmat[KSPSIZE+1][KSPSIZE],
         amrex::Real y[KSPSIZE],amrex::Real beta_e1[KSPSIZE+1],
         int k)
{
    for(int i=k-1;i>=0;i--)	
    {
        y[i] = beta_e1[i];
        for(int j=i+1;j<=k;j++)
        {
            y[i]=y[i] - Hessmat[i][j] * y[j];
        }
        y[i] = y[i]/Hessmat[i][i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
bool performgmres(amrex::Real Mat[NEQNS][NEQNS],
                  amrex::Real b[NEQNS],
                  amrex::Real  x0[NEQNS], //initial soln for GMRES
                  amrex::Real x[NEQNS],   //current solution for GMRES
                  int precond_type,
                  int restart_it,
                  amrex::Real tol,int printflag)
{
    amrex::Real beta;
    bool arnoldistopped;
    bool success;
    bool nanflag;

    amrex::Real v1[NEQNS]={0.0};
    amrex::Real tempvec[NEQNS]={0.0};
    amrex::Real v[NEQNS]={0.0};
    amrex::Real y[KSPSIZE]={0.0};
    amrex::Real r[NEQNS]={0.0};
    amrex::Real r0[NEQNS]={0.0};
    amrex::Real Ax0[NEQNS]={0.0};
    amrex::Real Ax[NEQNS]={0.0};
    amrex::Real Minvr[NEQNS]={0.0};
    amrex::Real kspvectors[KSPSIZE+1][NEQNS]={0.0};
    amrex::Real Hessmat[KSPSIZE+1][KSPSIZE]={0.0};
    amrex::Real cos_arr[KSPSIZE]={0.0};
    amrex::Real sin_arr[KSPSIZE]={0.0};
    amrex::Real beta_e1[KSPSIZE+1]={0.0};
    amrex::Real residnorm,bnorm;
    
    success = true;

    //finding r0
    findAX(Ax0,Mat,x0);
    addvectors(b,Ax0,r0,1.0,-1.0);  //b-Ax0
    precond(Minvr,Mat,r0,precond_type); //Minv * r0
    copyvector(r0,Minvr);

    //initial residual is r0=M^-1(b-Ax0)
    //we are solving M^-1 A X = M^-1 b
    copyvector(r,r0);
    copyvector(x,x0); 

#ifndef AMREX_USE_GPU
    if(printflag) amrex::Print()<<"initial norm of residual:"
        <<findnorm(r0)<<"\n";
#endif

    for(int it=0;it<restart_it;it++)
    {
#ifndef AMREX_USE_GPU
        if(printflag) amrex::Print()<<"restart iteration:"
            <<it<<"\n";
#endif

        copyvector(x0,x);
        beta = findnorm(r);

        for(int i=0;i<NEQNS;i++)
        {
            v1[i]=r[i]/beta;
        }
        beta_e1[0]=beta;

        setkspvector(v1,kspvectors,0);
        
        int optkspsize=KSPSIZE;
        for(int kspdim=0;kspdim<KSPSIZE;kspdim++)
        {
            copyvector(x,x0);
            // finds the ksp vector at kspdim+1
            arnoldistopped = arnoldi(kspdim,Mat,kspvectors,Hessmat,
                                     precond_type,tol*0.01,nanflag);

            residnorm=triangularize(Hessmat,cos_arr,sin_arr,beta_e1,
                                    kspdim);
            
            if(nanflag)
            {
                success=!nanflag;
                amrex::Abort("nan detected in arnoldi iteration\n");
                break;
            }


#ifndef AMREX_USE_GPU
            if(printflag) amrex::Print()<<"norm of residual:"
                <<residnorm<<"\t"<<tol<<"\t"<<kspdim<<"\n";
#endif

            if(arnoldistopped) 
            {
                optkspsize=kspdim+1;
#ifndef AMREX_USE_GPU
                if(printflag) amrex::Print()<<"lucky condition:"
                    <<kspdim<<"\n";
#endif
                break;
            }

            if(residnorm <= tol)
            {
                optkspsize=kspdim+1;
                break;
            }
        }

        uppertrisolve(Hessmat,y,beta_e1,optkspsize);
        for(int i=0;i<optkspsize;i++)
        {
            getkspvector(v,kspvectors,i);
            addvectors(x,v,tempvec,1.0,y[i]);
            copyvector(x,tempvec);
        }

        if( (residnorm <= tol) || (!success) || (arnoldistopped) )
        {
            break;
        }
    }


    return(success);
}
//========================================================================
