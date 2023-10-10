field rateofstraintensor(field &uh) {

  /* Récupération de la géométrie */
  const geo &omega=uh.get_geo();

  /* Recuperation de l'espace d'approximation de uh */
  const space &Xh=uh.get_space();
  
   /* Espace d'approximation pour le calcul des taux de cisaillement et de la viscosité */
  string Pkd="P" +itos(Xh.degree()-1)+"d";
  space Th(omega,Pkd,"tensor");
  
  /* Declaration des fonctions test et d'essai */
  trial gamma(Th);
  test tau(Th);

  /* Forme bilinéaire correspondant à la matrice masse inverse */
  integrate_option iopt;
  iopt.invert = true;
  form invm=integrate(ddot(gamma,tau),iopt);

  /* Determination du terme source */
  field lh=integrate(ddot(D(uh),tau));

  /* Calcul de la solution correspondant au  tenseur des taux de deformation */
  return invm*lh;
}
