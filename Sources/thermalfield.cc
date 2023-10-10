field thermalfield(Float Pe,Float Gn,Float Nucyl,Float Nunoz,string PkT,Float adt,Float CArrh,
		   Float TrefArrh,Float taumu,Float nmu,Float amu,field &TKh,field &uh,field &fh,
		   quadrature_option &qopt) {
  
  /* Récupération de la géométrie */
  const geo &omega=uh.get_geo();
  
  /* Définition de coefficient de decentrement */
  Float alphadv=1.;

  /* Definition de l'espace d'approximation sur omega de la temperature */
  space Th(omega,PkT);

  /* Recuperation de la dimension du probleme et de l'ordre polynomial
     d'approximation */
  size_t d=omega.dimension();
  size_t k=Th.degree();

  /* Determination du coefficient beta pour la penalisation des conditions
     aux limites */
  Float beta=(k+1)*(k+d)/d;

  /* Declaration des fonctions d'essai et de test */
  trial theta(Th);
  test thetat(Th);

  /* Determination de la forme intégrale correspondant au terme instationnaire */
  form ainsta=adt*integrate(theta*thetat);

  /* Determination de la forme integrale du terme convectif */
  form aconv=integrate(dot(uh,grad_h(theta))*thetat)+
    integrate("internal_sides",-dot(uh,normal())*jump(theta)*average(thetat)+
	      0.5*alphadv*abs(dot(uh,normal()))*jump(theta)*jump(thetat));

  /* Determination de la forme integrale de -div[grad(theta)] */
  form adif=integrate(dot(grad_h(theta),grad_h(thetat)))+ 
    integrate("internal_sides",beta*penalty()*jump(theta)*jump(thetat)-
	      jump(theta)*average(dot(grad_h(thetat),normal()))-
	      jump(thetat)*average(dot(grad_h(theta),normal())))+
    integrate("left",beta*penalty()*jump(theta)*jump(thetat)-
	      jump(theta)*average(dot(grad_h(thetat),normal()))-
	      jump(thetat)*average(dot(grad_h(theta),normal())))+
    integrate("cylinder",Nucyl*theta*thetat)+integrate("nozzle",Nunoz*theta*thetat);

  /* Constitution de la forme globale */
  form a=Pe*ainsta+Pe*aconv+adif;

  /* Determination des formes lineaires correspondant aux second membre */

  /* 1. Imposition faible de la temperature en entree */
  //field lh=integrate("left",beta*penalty()*thetat-dot(grad_h(thetat),normal()));

  /* 2. ajout des termes correspondant à la condition de Fourier */
  field lh=integrate("cylinder",Nucyl*thetat)+integrate("nozzle",Nunoz*thetat);

  /* 3. Ajout du terme source correspondant à la partie explicite de la dérivée temporelle */
  lh+=Pe*integrate(fh*thetat);

  /* 4. Ajout de la dissipation visqueuse */
  field etah=calculviscosite(CArrh,TrefArrh,taumu,nmu,amu,TKh,uh);
  field Dh=rateofstraintensor(uh);
  lh+=Gn*integrate(2.*etah*ddot(Dh,Dh)*thetat);

  /****************************************/
  /* Resolution du probleme a(T,Tt)=lh(Tt) */
  /****************************************/
  
  problem pa (a);
  field thetah(Th);
  pa.solve(lh,thetah);

  /* Retour à la fonction appelante */
  return thetah;

}
