double rhoY2rho(amrex::Array4<amrex::Real> rhoY) 
{
   double rho = 0;
   for (int i=0; i<NUM_SPECIES; i++)
   {
       rho += rhoY.dataPtr()[i];
   }
   return rho;
}
void rhoY2Y(amrex::Array4<amrex::Real> rhoY, double* Y) 
{
   double rho = rhoY2rho(rhoY);
   for (int i=0; i<NUM_SPECIES; i++)
   {
       Y[i] = rhoY.dataPtr()[i]/rho;
   }
}


