//
// helper functions
//

//template <class T> const T& min (const T& a, const T& b) {
//  return (a<b)?a:b;   
//}

//template <class T> const T& max (const T& a, const T& b) {
//  return (a>b)?a:b;   
//}

class UnitConversion {
  double m,l,t,th,v,d,s,scl;

  public:

  UnitConversion(){
    m=0.0;
    l=0.0;
    t=0.0;
    th=0.0; 
    d=0.0;
    v=0.0;
    s=0.0;
    scl=0.;
  } // constructor
  
  void set_mass(double tmp)  {m =tmp;}
  void set_length(double tmp){l =tmp;}
  void set_time(double tmp)  {t =tmp;}
  void set_tk(double tmp)    {th=tmp;}

  void set_scale(double scl_in, int MKS2CGS){
   scl=scl_in;

   if(MKS2CGS==1){// assumes conversion from MKS -- FARGO OUTPUT!!!
     l*=1e2;
     d*=1e-3;
     s*=1e-1;
     m*=1e3;
     v*=1e2;
   }
   d /=pow(scl,3.);
   s /=pow(scl,2.);
   l /=scl;
   th/=scl;
   v /=scl;
  }

  void make(){
    d=m/pow(l,3.);
    s=d*l;
    v=l/t;
  }

  double mass()  {return m; }
  double length(){return l; }
  double time()  {return t; }
  double tk()    {return th;}
  double vel()   {return v; }
  double den()   {return d; }
  double sig()   {return s; }
  double scale() {return scl;}
};

