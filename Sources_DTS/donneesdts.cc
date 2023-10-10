int donneesdts(string &nomfichin,string &Pkd,size_t &bdfk,Float &dt,Float &injtime,size_t &itmax,Float &tol,string &nomflow,
	               string &nomfield,string &nomdata,size_t &ffield,size_t &fbranch) {

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
  entree >> nomvar >> Pkd
	 >> nomvar >> bdfk
	 >> nomvar >> dt
	 >> nomvar >> injtime
	 >> nomvar >> itmax
	 >> nomvar >> tol
	 >> nomvar >> nomflow
	 >> nomvar >> nomfield
	 >> nomvar >> nomdata
    	 >> nomvar >> ffield
    	 >> nomvar >> fbranch;

  /************************/
  /* Fermeture du fichier */
  /************************/

  entree.close();
  
  /**************************************/
  /* Affichage des parametres de calcul */
  /**************************************/

  dout << "Pkd=" << Pkd << endl
       << "bdfk=" << bdfk << endl
       << "dt=" << dt << endl
       << "injtime=" << injtime << endl
       << "itmax=" << itmax << endl
       << "tol=" << tol << endl
       << "nomflow=" << nomflow << endl
       << "nomfield=" << nomfield << endl
       << "nomdata=" << nomdata << endl
       << "ffield=" << ffield << endl
       << "fbranch=" << fbranch << endl;

  /*****************************/
  /* Fin normal de la fonction */
  /*****************************/
  return 0;
}
