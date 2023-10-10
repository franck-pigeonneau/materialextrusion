/************************************************************************
  Ce programme C++-RHEOLEF réalise la résolution temporelle des équations
  de Stokes et de la chaleur par une discretisation par éléments finis.

  La viscosité est une fonction du taux de cisaillement et de la température
  selon une loi de puissance, soit une loi de Carreau.

  Auteur : Franck Pigeonneau.

  Date de création : 22/03/2019.

  Version Rheolef: 7.0
************************************************************************/

/* Introduction du fichier d'entête de la librairie rheolef */
#include "rheolef.h"

/* Utilisation des classes cpp */
#include <iostream> // flux entrée/sortie
#include <fstream> // flux d'entrée/sortie dans des fichiers
#include <string> // Manipulation des chaines de caractères

/* Introduction des identificateurs des bibliotheques standard et rheolef */
using namespace rheolef;
using namespace std;

/* Introduction des fonctions */
#include "bdfk.h"
#include "donneesGNewFluid.cc"
#include "viscopowerlaw.cc"
#include "rateofstraintensor.cc"
#include "calculviscosite.cc"
#include "stokessolver.cc"
#include "thermalfield.cc"

/*********************************/
/* Début du programme principal */
/*********************************/

int main(int argc, char**argv) {
  environment rheolef (argc, argv);

  /* Récupération de la précision machine */
  Float eps = numeric_limits<Float>::epsilon();

  /*******************************************/
  /* Rappel de la forme d'appel du programme */
  /*******************************************/

  if (argc < 2) {
    cout << "Appel de la fonction incorrecte" << endl;
    exit(0);
  }

  /*******************************************/
  /* Recuperation du nom du fichier d'entrée */
  /*******************************************/

  string nomfichin=argv[1];

  /*********************************************/
  /* Lecture des données nécessaires au calcul */
  /*********************************************/

  int info;
  string fichgeo;
  string Pku;
  string PkT;
  Float tol;
  size_t maxit;
  size_t bdfk;
  Float dt;
  Float CArrh;
  Float TrefArrh;
  Float Tinlet;
  Float Textr;
  Float Pe;
  Float Gn;
  Float Nucyl;
  Float Nunoz;
  Float taumu;
  Float nmu;
  Float amu;
  Float u0;
  string nomfieldu;
  string nomfieldT;
  string nomfichdat;
  size_t fbranch;
  size_t ffield;

  info=donneesGNewFluid(nomfichin,fichgeo,Pku,PkT,tol,maxit,bdfk,dt,CArrh,
			TrefArrh,Tinlet,Textr,Pe,Gn,Nucyl,Nunoz,taumu,nmu,
			amu,u0,nomfieldu,nomfieldT,nomfichdat,fbranch,ffield);

  if (info!=0) {
    dout << "Probleme de lecture des donnees." << endl;
    exit(0);
  }
  
  /* Lecture du domaine de calcul */
  geo omega(fichgeo);

  /* Declaration de l'espace d'approximation de la vitesse */
  space Xh(omega,Pku,"vector");

  /* Declaration de l'espace d'approximation de la temperature */
  space Yh(omega,PkT);

  /* Blockage des degres de liberté sur l'entree et sur les parois du cylindre et de la buse */
  Xh.block("cylinder");
  Xh.block("nozzle");
  Xh.block("left");
  
  /* Declaration de l'espace d'approximation de la pression */
  string Pkp="P" +itos(Xh.degree()-1);
  space Qh(omega,Pkp);

  /* Espace d'approximation pour le calcul des taux de cisaillement et de la viscosité */
  string Pkd="P" +itos(Xh.degree()-1)+"d";
  space Th(omega,Pkd,"tensor");
  space Sh(omega,Pkd);

  /* Option de la quadrature */
  quadrature_option_type qopt;
  qopt.set_family(quadrature_option_type::gauss);
  qopt.set_order(2*Xh.degree()+1);
  quadrature_option_type qopt2;
  qopt2.set_family(quadrature_option_type::gauss);
  qopt2.set_order(2*Yh.degree()+1);

  /* Formes bilinéaires correspondant aux matrices masses*/
  trial u(Xh),p(Qh),theta(Yh);
  test v(Xh),q(Qh),alpha(Yh);
  
  form mu=integrate(dot(u,v),qopt);
  form mp=integrate(p*q,qopt);
  form mt=integrate(theta*alpha,qopt2);
  
  /****************************************************/
  /* Initialisation des champs de vitesse et de pression */
  /***************************************************/

  dout << "Initialisation de la solution." << endl;

  vector<field> uh(bdfk+1);
  uh[0]=field(Xh,0.);
  field ph(Qh,0.);
  
  /* Initialisation de la solution en Stokes */
  uh[0][0]["left"]=u0;
  uh[0]["cylinder"]=0.;
  uh[0]["nozzle"]=0.;
  Float thetainit=0.;
  Float Tinit=Tinlet+thetainit*(Textr-Tinlet);
  field TKh(Yh,Tinit);
  Stokessolver(CArrh,TrefArrh,taumu,nmu,amu,TKh,uh[0],ph,qopt);
  uh[1]=uh[0];
  
  /*****************************************/
  /* Initialisation du champ de temperature */
  /*****************************************/
  
  vector<field> thetah(bdfk+1);
  thetah[0]=field(Yh,thetainit);
  thetah[1]=thetah[0];

  /* Branch definition */
  branch eventu("t","u","p");
  branch eventT("t","T");

  /* Sauvagarde des donnees initiales */
  odiststream outbranchu(nomfieldu,"branch");
  odiststream outbranchT(nomfieldT,"branch");
  outbranchu << eventu(0.,uh[0],ph);
  outbranchT << eventT(0.,thetah[0]);

  /* Ouverture du fichier de sauvegarde de la dérivée temporelle */
  ofstream outdat(nomfichdat,ios::out);
  outdat << "#t      ||dudt||      ||dTdt||" << endl;

  /*********************************************************************************/
  /* Calcul de la solution du problème de Stokes pour un fluide de Newton généralisé */
  /*********************************************************************************/
  
  dout << "Début du calcul ..." << endl;

  Float minL2=1.e1;
  
  for (size_t n=1;(n<maxit)&&(minL2>tol);n++) {

    /* Determination of the time coefficient of the temporal derivatives */
    size_t ktime=min(bdfk,n);
    Float adt= (dt>0.) ? ABDF[ktime][0]/dt : 0.;

    /* Determination of the explicit part of the temporal derivatives */
    field STh(Yh,0.);
    if (dt!=0.) {
      for (size_t i=1;i<=ktime;i++) {
	STh+=ABDF[ktime][i]*thetah[i];
      }
      STh/=dt;
    }

    /* Extrapolation of theta */
    thetah[0]=BBDF[ktime][1]*thetah[1];
    for (size_t i=2;i<=ktime;i++) {
      thetah[0]+=BBDF[ktime][i]*thetah[i];
    }

    /* Calcul de la nouvelle solution de uh et ph*/
    TKh=Tinlet+(Textr-Tinlet)*thetah[0];
    Stokessolver(CArrh,TrefArrh,taumu,nmu,amu,TKh,uh[0],ph,qopt);

    /* Calcul de la nouvelle solution de thetah */
    thetah[0]=thermalfield(Pe,Gn,Nucyl,Nunoz,PkT,adt,CArrh,TrefArrh,taumu,nmu,amu,TKh,uh[0],STh,qopt2);

    /* Temporal derivations of the velocity and temperature */
    field duhdt=ABDF[ktime][0]*uh[0];
    field dthetahdt=ABDF[ktime][0]*thetah[0];
    for (size_t i=1;i<=ktime;i++) {
      duhdt-=ABDF[ktime][i]*uh[i];
      dthetahdt-=ABDF[ktime][i]*thetah[i];
    }
    if (dt!=0.) {
      duhdt/=dt;
      dthetahdt/=dt;
    }
    
    /* Calcul de la norme L2 de la dérivée temporelle de u et de T */
    Float normL2dudt=sqrt(mu(duhdt,duhdt));
    Float normL2dTdt=sqrt(mt(dthetahdt,dthetahdt));
    
    /* Calcul d'eccart quadratique des termes instationnaires */
    minL2=sqrt(pow(normL2dudt,2)+pow(normL2dTdt,2));
    
    /* Sauvegarde de thetah et uh de aux ktime pas de temps précédents */
    for (size_t i=min(bdfk,ktime+1);i>=1;i--) {
      thetah[i]=thetah[i-1];
      uh[i]=uh[i-1];
    }

    /* Sauvegarde des champs dans les fichiers branch */
    float t=n*dt;
    if (n % fbranch == 0) {
      outbranchT << eventT(t,thetah[0]);
      outbranchu << eventu(t,uh[0],ph);
    }
    
    if (n % ffield == 0) {
      /* Sauvegarde de la temperature dans le fichier nomfieldtemp */
      TKh=Tinlet+(Textr-Tinlet)*thetah[0];
      field etah=calculviscosite(CArrh,TrefArrh,taumu,nmu,amu,TKh,uh[0]);
      odiststream Tfield (nomfieldT,"field");
      Tfield << catchmark("T") << thetah[0]
	     << catchmark("eta") << etah;
      Tfield.close();

      /* Sauvegarde de la vitesse et de la pression dans le fichier nomfieldu */
      odiststream motionfield (nomfieldu,"field");
      field Dh=rateofstraintensor(uh[0]);
      field gammadot=interpolate(Sh,sqrt(2.*ddot(Dh,Dh)));
      field phimuh=Gn*interpolate(Sh,2.*ddot(Dh,Dh)*etah);
      motionfield << catchmark("u") << uh[0]
		  << catchmark("p") << ph
		  << catchmark("gammadot") << gammadot
		  << catchmark("phimu") << phimuh;
      motionfield.close();
    }

    outdat << t << "  " << normL2dudt << " " << normL2dTdt << endl;
    dout << t << "  " << normL2dudt << " " << normL2dTdt << endl;
  }

  /* Sauvegarde de la temperature dans le fichier nomfieldtemp */
  odiststream Tfield (nomfieldT,"field");
  TKh=Tinlet+(Textr-Tinlet)*thetah[0];
  field etah=calculviscosite(CArrh,TrefArrh,taumu,nmu,amu,TKh,uh[0]);
  Tfield << catchmark("T") << thetah[0]
	 << catchmark("eta") << etah;
  Tfield.close();

  /* Sauvegarde de la vitesse et de la pression dans le fichier nomfieldu */
  odiststream motionfield (nomfieldu,"field");
  field Dh=rateofstraintensor(uh[0]);
  field gammadot=interpolate(Sh,sqrt(2.*ddot(Dh,Dh)));
  field phimuh=Gn*interpolate(Sh,2.*ddot(Dh,Dh)*etah);
  motionfield << catchmark("u") << uh[0]
              << catchmark("p") << ph
	      << catchmark("gammadot") << gammadot
              << catchmark("phimu") << phimuh;
  motionfield.close();

  /* Fermeture des fichiers branch et dat */
  outbranchu.close();
  outbranchT.close();
  outdat.close();

  /****************************/
  /* Fin normal du programme */
  /****************************/
  return 0;
}
