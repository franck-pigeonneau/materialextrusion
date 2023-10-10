/* 
  Fonction  résolvant le problème de Stokes.
 */

void Stokessolver(Float CArrh,Float TrefArrh,Float taumu,Float nmu,Float amu,field &TKh,field& uh, field& ph,quadrature_option &qopt) {

  /* Récupération des espaces de la vitesse et de la pression */
  const space &Xh=uh.get_space();
  const space &Qh=ph.get_space();
  const geo  &omega=uh.get_geo();
  
  /* Déclaration des fonctions test et d'essai*/
  trial u(Xh),p(Qh);
  test  v(Xh),q(Qh);

  /* Calcul de la viscosite */
  field etah=calculviscosite(CArrh,TrefArrh,taumu,nmu,amu,TKh,uh);
  
  /* Ecriture des formes varioationnelles */
  form a  = integrate(2.*etah*ddot(D(u),D(v)),qopt);
  form b  = -integrate (div(u)*q);
  form mp = integrate (p*q);

  /* Résolution du problème de Stokes */
  problem_mixed stokes(a,b);
  stokes.solve(field(Xh,0.),field(Qh,0.),uh,ph);
}
