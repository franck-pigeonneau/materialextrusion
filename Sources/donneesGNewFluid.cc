int donneesGNewFluid(string &nomfichin,string &fichgeo,string &Pku,string &PkT,Float &tol,size_t &maxit,
		     size_t &bdfk,Float &dt,Float &CArrh,Float &TrefArrh,Float &Tinlet,Float &Textr,Float &Pe,Float &Gn,Float &Nucyl,
		     Float &Nunoz,Float &taumu,Float &nmu,Float &amu,Float &u0,string &nomfieldu,string &nomfieldT,string &nomfichdat,
		     size_t &fbranch,size_t &ffield) {

  /* Ouverture du fichier en lecture */
  ifstream entree;
  entree.open (nomfichin,ios::in);
  if (!entree) {
    cout << "Ouverture du fichier " << nomfichin << " impossible" << endl;
    return -1;
  }

  /************************************/
  /* Lecture des parametres de calcul */
  /************************************/

  string nomvar;
  entree >> nomvar >> fichgeo
	 >> nomvar >> Pku
         >> nomvar >> PkT
	 >> nomvar >> tol
	 >> nomvar >> maxit
    	 >> nomvar >> bdfk
    	 >> nomvar >> dt
	 >> nomvar >> CArrh
	 >> nomvar >> TrefArrh
	 >> nomvar >> Tinlet
	 >> nomvar >> Textr
	 >> nomvar >> Pe
	 >> nomvar >> Gn
	 >> nomvar >> Nucyl
	 >> nomvar >> Nunoz
    	 >> nomvar >> taumu
	 >> nomvar >> nmu
	 >> nomvar >> amu
	 >> nomvar >> u0
	 >> nomvar >> nomfieldu
	 >> nomvar >> nomfieldT
	 >> nomvar >> nomfichdat
	 >> nomvar >> fbranch
	 >> nomvar >> ffield;

  /************************/
  /* Fermeture du fichier */
  /************************/

  entree.close();

  /**************************************/
  /* Affichage des parametres de calcul */
  /**************************************/

  dout << "fichgeo=" << fichgeo << endl
       << "Pku=" << Pku << endl
       << "PkT=" << PkT << endl
       << "tol=" << tol << endl
       << "maxit=" << maxit << endl
       << "bdfk=" << bdfk << endl
       << "dt=" << dt << endl
       << "CArrh=" << CArrh << endl
       << "TrefArrh=" << TrefArrh << endl
       << "Tinlet=" << Tinlet << endl
       << "Textr=" << Textr << endl
       << "Pe=" << Pe << endl
       << "Gn=" << Gn << endl
       << "Nucyl=" << Nucyl << endl
       << "Nunoz=" << Nunoz << endl
       << "taumu=" << taumu << endl
       << "nmu=" << nmu << endl
       << "amu=" << amu << endl
       << "u0=" << u0 << endl
       << "nomfieldu=" << nomfieldu << endl
       << "nomfieldT=" << nomfieldT << endl
       << "nomfichdat=" << nomfichdat << endl
       << "fbranch=" << fbranch << endl
       << "ffield=" << ffield << endl;

  /*****************************/
  /* Fin normal de la fonction */
  /*****************************/
  return 0;
}
