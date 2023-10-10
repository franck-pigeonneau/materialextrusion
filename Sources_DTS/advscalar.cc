field advscalar(Float adt,string Pkd,Float c0,field &uh,field &fh,quadrature_option &qopt) {
  
  /* Récupération de la géométrie */
  const geo &omega=uh.get_geo();
  
  /* Définition de coefficient de decentrement */
  Float alphadv=1.;

  /* Definition de l'espace d'approximation sur omega de la concentration */
  space Xh(omega,Pkd);

  /* Recuperation de la dimension du probleme et de l'ordre polynomial
     d'approximation */
  size_t d=omega.dimension();
  size_t k=Xh.degree();

  /* Determination du coefficient beta pour la penalisation des sauts */
  Float beta=(k+1)*(k+d)/d;

  /* Declaration des fonctions d'essai et de test */
  trial c(Xh);
  test ct(Xh);

  /* Determination de la forme integrale du terme convectif */
  form a=integrate(dot(uh,grad_h(c))*ct+adt*c*ct)+
              integrate("internal_sides",-dot(uh,normal())*jump(c)*average(ct)+
			0.5*alphadv*abs(dot(uh,normal()))*jump(c)*jump(ct))+
              integrate("left",max(0,-dot(uh,normal()))*c*ct);

  /* Determination des formes lineaires correspondant au second membre */
  field lh=integrate(fh*ct)+integrate("left",max(0,-dot(uh,normal()))*c0*ct);

  /****************************************/
  /* Resolution du probleme a(c,ct)=lh(ct) */
  /****************************************/
  
  solver sa (a.uu());
  field ch(Xh);
  ch.set_u()=sa.solve(lh.u());

  /* Retour à la fonction appelante */
  return ch;

}
