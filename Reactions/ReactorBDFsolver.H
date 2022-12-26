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
                                              amrex::Real X[NEQNS],
                                              amrex::Real dt)
{
    for(int i=0;i<NEQNS;i++)
    {
        MinvX[i]=X[i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void GJprecond(amrex::Real MinvX[NEQNS],
                                               amrex::Real A[NEQNS][NEQNS],
                                               amrex::Real X[NEQNS],
                                               amrex::Real dt)
{
    amrex::Real sum,diag;
    amrex::Real y_d[NEQNS],y_dd[NEQNS];
    amrex::Real dt_inv=1.0/dt;

    //solve Dy = X
    //MinvX = y

    for(int i=0;i<NEQNS;i++)
    {
        diag=dt_inv-A[i][i];
        MinvX[i]=X[i]/diag;
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void SGSprecond(amrex::Real MinvX[NEQNS],
                                               amrex::Real A[NEQNS][NEQNS],
                                               amrex::Real X[NEQNS],
                                               amrex::Real dt)
{
    amrex::Real sum;
    amrex::Real y_d[NEQNS],y_dd[NEQNS];
    amrex::Real diag_inv[NEQNS];
    amrex::Real diag[NEQNS];
    amrex::Real dt_inv=1.0/dt;

    //solve (D+L)D^-1(D+U)y = X
    //MinvX = y

    for(int i=0;i<NEQNS;i++)
    {
        diag[i]=dt_inv-A[i][i];
        diag_inv[i]=1.0/diag[i];
    }


    //solve (D+L) y' = X
    for(int i=0;i<NEQNS;i++)
    {
        sum  = 0.0;
        for(int j=0;j<i;j++)
        {
            sum += -A[i][j]*y_d[j];
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
            sum += -A[i][j]*MinvX[j];
        }

        MinvX[i] = (y_dd[i]-sum)*diag_inv[i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void precond(amrex::Real MinvX[NEQNS],
                                            amrex::Real A[NEQNS][NEQNS],
                                            amrex::Real X[NEQNS],
                                            amrex::Real dt,int type)
{
    if(type==NOPC)
    {
        noprecond(MinvX,A,X,dt);   
    }
    else if(type==GJPC)
    {
        GJprecond(MinvX,A,X,dt);
    }
    else if(type==SGSPC)
    {
        SGSprecond(MinvX,A,X,dt);
    }
    else
    {
        amrex::Abort("Unknown preconditioner type in BDF");
    }

}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void findAX(amrex::Real AX[NEQNS],
                                           amrex::Real A[NEQNS][NEQNS],
                                           amrex::Real X[NEQNS],
                                           amrex::Real dt)
{
    //The linear system for first order backward Euler 
    //is of the form
    //[I/dt - df/du] du = (un-u)/dt + f(u)

    //The linear system for second order Crank Nicholson
    //is of the form
    //[I/dt - 0.5*df/du] du = (un-u)/dt + 0.5*f(u) + 0.5*f(un)

    //time stepping contribution
    for(int i=0;i<NEQNS;i++)
    {
        AX[i]=X[i]/dt;
    }

    //Jacobian contribution
    for(int i=0;i<NEQNS;i++)
    {
        for(int j=0;j<NEQNS;j++)
        {
            AX[i] -= A[i][j]*X[j];
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
                                                 amrex::Real dt,
                                                 int precond_type,
                                                 int luckykspdim,
                                                 bool &nanflag)
{
    int i,n;
    amrex::Real Avj[NEQNS],vj[NEQNS],vi[NEQNS];
    amrex::Real wj[NEQNS],tempvec[NEQNS];
    amrex::Real MinvAvj[NEQNS];
    luckykspdim = KSPSIZE;
    

    bool lucky; //when norm becomes 0,
    //KSP vectors is no longer linearly independent.
    //we would have gotten the best solution.

    lucky   = false;
    nanflag = false;

    getkspvector(vj,kspvectors,j);
    //printvec("vj:",vj);
    findAX(Avj,mat,vj,dt);
    //printvec("Avj:",Avj);
    precond(MinvAvj,mat,Avj,dt,precond_type);
    //printvec("MinvAvj:",MinvAvj);
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

    if(Hessmat[j+1][j] > NEARZERO)
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
            lucky=true;
            luckykspdim=j+1;
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
leastsqminimize(amrex::Real Hessmat[KSPSIZE+1][KSPSIZE],
                amrex::Real y[KSPSIZE],amrex::Real beta,
                int currentkspsize,bool lucky)
{
    amrex::Real c,s,h_up,h_down,dtr;
    amrex::Real val1,val2;

    //local copy
    amrex::Real Hmat[KSPSIZE+1][KSPSIZE];
    for(int i=0;i<(KSPSIZE+1);i++)
    {
        for(int j=0;j<KSPSIZE;j++)
        {
            Hmat[i][j]=Hessmat[i][j];
        }
    }    

    amrex::Real beta_e1[KSPSIZE+1]={0.0};
    beta_e1[0] = beta;

    //convert H into QR
    for(int i=0;i<currentkspsize;i++)
    {
        h_up   = Hmat[i][i] ;
        h_down = Hmat[i+1][i];

        dtr = sqrt(h_up*h_up + h_down*h_down);

        c=h_up/dtr; s=h_down/dtr;

        for(int j=0;j<currentkspsize;j++)
        {
            h_up   = Hmat[i][j];
            h_down = Hmat[i+1][j];

            //perform rotations
            //ith row
            Hmat[i][j] =  c*h_up + s*h_down;
            //(i+1)th row
            Hmat[i+1][j] = -s*h_up + c*h_down;

        }

        val1 =  c*beta_e1[i] + s*beta_e1[i+1];
        val2 = -s*beta_e1[i] + c*beta_e1[i+1];

        beta_e1[i]=val1; 
        beta_e1[i+1]=val2;

    }

    // ||Hm y - beta e1|| = || QR y - Q Q^T beta e1||
    // || Q ( Ry - Q^T beta e1) || = || Ry - Q^T beta e1||

    //solve least squares problem
    y[currentkspsize-1] = beta_e1[currentkspsize-1]/
    Hmat[currentkspsize-1][currentkspsize-1];

    for(int i=currentkspsize-2;i>=0;i--)	
    {
        y[i] = beta_e1[i];

        for(int j=i+1;j<currentkspsize;j++)
        {
            y[i]=y[i] - Hmat[i][j] * y[j];
        }

        y[i] = y[i]/Hmat[i][i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
bool performgmres(amrex::Real Mat1d[NEQNS2],
                  amrex::Real ydot[NEQNS],
                  amrex::Real  x0[NEQNS], //initial soln for GMRES
                  amrex::Real x[NEQNS],   //current solution
                  amrex::Real dsoln_n[NEQNS],   //for explicit term
                  amrex::Real dt_dim,
                  int precond_type,
                  int restart_it,
                  amrex::Real tol,int printflag)
{
    amrex::Real beta;
    bool arnoldistopped;
    bool success;
    bool nanflag;
    int luckykspdim;

    amrex::Print()<<"kspsize:"<<KSPSIZE<<"\n";

    amrex::Real v1[NEQNS]={0.0};
    amrex::Real tempvec[NEQNS]={0.0};
    amrex::Real v[NEQNS]={0.0};
    amrex::Real y[KSPSIZE]={0.0};
    amrex::Real r[NEQNS]={0.0};
    amrex::Real r0[NEQNS]={0.0};
    amrex::Real Ax0[NEQNS]={0.0};
    amrex::Real Ax[NEQNS]={0.0};
    amrex::Real Minvr[NEQNS]={0.0};
    amrex::Real b[NEQNS]={0.0};
    amrex::Real Mat[NEQNS][NEQNS]={0.0};
    amrex::Real kspvectors[KSPSIZE+1][NEQNS]={0.0};
    amrex::Real Hessmat[KSPSIZE+1][KSPSIZE]={0.0};
    amrex::Real residnorm,bnorm;

    success = true;


    //set RHS
    for(int i=0;i<NEQNS;i++)
    {
        b[i]=-(dsoln_n[i]/dt_dim-ydot[i]);
    }
    bnorm=findnorm(b);
    //bnorm=1.0;
    amrex::Real dt=dt_dim*bnorm;
    
    //set matrix
    for(int i=0;i<NEQNS;i++)
    {
        for(int j=0;j<NEQNS;j++)
        {
              Mat[i][j]=Mat1d[i*NEQNS+j]/bnorm;
        }
        b[i]=b[i]/bnorm;
    }

    printmat("Jacobian",Mat);
    //finding r0
    findAX(Ax0,Mat,x0,dt);
    addvectors(b,Ax0,r0,1.0,-1.0);  //b-Ax0
    precond(Minvr,Mat,r0,dt,precond_type); //Minv * r0
    copyvector(r0,Minvr);

    //initial residual is r0=M^-1(b-Ax0)
    //we are solving M^-1 A X = M^-1 b
    copyvector(r,r0);
    printvec("r0",r0);
    printvec("r",r);
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

        setkspvector(v1,kspvectors,0);

        for(int kspdim=0;kspdim<KSPSIZE;kspdim++)
        {
            copyvector(x,x0);
            // finds the ksp vector at kspdim+1
            arnoldistopped = arnoldi(kspdim,Mat,kspvectors,Hessmat,
                                     dt,precond_type,luckykspdim,nanflag);

            if(nanflag)
            {
                success=!nanflag;
                amrex::Abort("nan detected in arnoldi iteration\n");
                break;
            }

            leastsqminimize(Hessmat,y,beta,kspdim+1,arnoldistopped);

            for(int i=0;i<(kspdim+1);i++)
            {
                getkspvector(v,kspvectors,i);
                addvectors(x,v,tempvec,1.0,y[i]);
                copyvector(x,tempvec);
            }

            //finding new residual
            findAX(Ax,Mat,x,dt);
            addvectors(b,Ax,r,1.0,-1.0);
            precond(Minvr,Mat,r,dt,precond_type); //Minv * r0
            copyvector(r,Minvr);

            residnorm = findnorm(r);
#ifndef AMREX_USE_GPU
            if(printflag) amrex::Print()<<"norm of residual:"
                <<residnorm<<"\t"<<tol<<"\t"<<kspdim<<"\n";
#endif

            if(arnoldistopped) 
            {
#ifndef AMREX_USE_GPU
                if(printflag) amrex::Print()<<"lucky condition:"
                    <<luckykspdim<<"\n";
#endif
                break;
            }

            if(residnorm <= tol)
            {
                break;
            }

        }

        if( (residnorm <= tol) || (!success) || (arnoldistopped) )
        {
            break;
        }


        //amrex::Print()<<"******************\n";
    }

    return(success);
}
//========================================================================
