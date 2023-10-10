field calculviscosite(Float CArrh,Float TrefArrh,Float taumu,Float nmu,Float amu,field &TKh,field &uh) {
  
  /* Récupération de la géométrie */
  const geo &omega=uh.get_geo();

  /* Recuperation de l'espace d'approximation de uh */
  const space &Xh=uh.get_space();
  
   /* Espace d'approximation pour le calcul des taux de cisaillement et de la viscosité */
  string Pkd="P" +itos(Xh.degree()-1)+"d";
  space Sh(omega,Pkd);
  
  /* Détermination du tenseur des taux de déformation par formulation faible */
  field Dh=rateofstraintensor(uh);

  /* Interpolation de la température sur l'espace Sh */
  field TSh=interpolate(Sh,TKh);
  
  /* Renvoie de la viscosite */
  return interpolate(Sh,compose(eta(CArrh,TrefArrh,taumu,nmu,amu),2.*ddot(Dh,Dh),TSh));
}
