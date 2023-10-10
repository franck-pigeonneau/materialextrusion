/*
  Cette fonction calcul la viscosité suivant une loi de puissance. Le variable x
  correspond au taux de cissaillement au carré, soit 2DijDij.
 */

struct eta {
  Float operator() (const Float &x,const Float &T) const {
    /* Calcul du taux de cisaillement sachant que x en entrée est égal 2D:D */
    Float gammadot=sqrt(x);
    
    /* Calcul du facteur de décallage */
    Float aT=exp(CArrh*(1./T-1./TrefArrh));

    /* Renvoie de la viscosite */
    return aT*pow(1.+pow(tau*aT*gammadot,a),(n-1.)/a);
  }
  eta (const Float &CArrh1, const Float &TrefArrh1,const Float &tau1, const Float &n1, const Float &a1) : CArrh(CArrh1), TrefArrh(TrefArrh1), tau(tau1), n(n1), a(a1) {}
Float CArrh,TrefArrh,tau,n,a;
};
