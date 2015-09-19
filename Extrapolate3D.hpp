//
// nora bolig
//
class Extrapolate3D{
  long nx,ny,nz;
  double mstar,adindx,dx,dy,dz,mu;
  std::vector< std::vector< std::vector<double> > > den3D;
  std::vector< std::vector<double > > h;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
//
//
//
  public:

  Extrapolate3D(){
  
  }//constructor

  void set_nx(long n){nx=n;}
  void set_ny(long n){ny=n;}
  void set_nz(long n){nz=n;}

  void set_dx(double tmp){dx=tmp;}
  void set_dy(double tmp){dy=tmp;}
  void set_dz(double tmp){dz=tmp;}
  void set_mstar(double tmp){mstar=tmp;}
  void set_adindx(double tmp){adindx=tmp;}
  void set_mu(double tmp){mu=tmp;}

  void set_x(std::vector<double> tmp){
    long n = tmp.size();
    x.resize(n);
    x=tmp;
  }
  void set_y(std::vector<double> tmp){
    long n = tmp.size();
    y.resize(n);
    y=tmp;
  }

  void initialize(long nx_in, long ny_in, long nz_in, double dx_in, double dy_in, 
                  double dz_in, double mstar_in, double adindx_in, double mu_in,
                  std::vector<double> x_in, std::vector<double>  y_in){
    set_nx(nx_in);
    set_ny(ny_in);
    set_nz(nz_in);
    set_mstar(mstar_in);
    set_adindx(adindx_in);
    set_dx(dx_in);
    set_dy(dy_in);
    set_dz(dz_in);
    set_mu(mu_in);
    set_x(x_in);
    set_y(y_in);

    h.resize(nx);
    for(long ix = 0; ix<nx; ++ix){h[ix].resize(ny);}

    z.resize(nz);
    for(long iz = 0; iz<nz; ++iz){z[iz] = (0.5-double(nz)*0.5 + double(iz))*dz;}

  }

  void set_scaleheight(const std::vector< std::vector<double> > &tkxy){

    for(long ix=0;ix<nx;++ix){
      for(long iy=0;iy<ny;++iy){
        double r = sqrt(x[ix]*x[ix]+y[iy]*y[iy]);
        double o = sqrt(mstar/pow(r,3.));
        double c = sqrt(adindx*tkxy[ix][iy]/mu);
        h[ix][iy] = c/o;
//        maxh=std::max<double>(maxh,h[ix][iy]); 
//        std::cout << x[ix]<<" "<<y[iy]<<" "<<h[ix][iy]<<" "<<h[ix][iy]/r<<" "<<h[ix][iy]/dz<<std::endl;
      }
    }
//    std::cout << "# CHECK MAX SCALE HEIGHT "<<maxh<<" and H/dz "<<maxh/dz<<std::endl;

   
  }

  void CreateIso3D(const std::vector< std::vector<double> > &sden){
    if(den3D.size()<nx){
      den3D.resize(nx);
      for(long ix=0;ix<nx;++ix){
        den3D[ix].resize(ny);
        for(long iy=0;iy<ny;++iy){
          den3D[ix][iy].resize(nz);
        }
      }
    }

    double m=0;
    double vol=dx*dy*dz;
    for(long ix=0;ix<nx;++ix){
      for(long iy=0;iy<ny;++iy){
        double rho0 = sden[ix][iy]/(2.5*h[ix][iy]);
        double hloc = h[ix][iy];
        for(long iz=0;iz<nz;++iz){
           if(hloc>0){den3D[ix][iy][iz] = rho0*exp(-0.5*pow(z[iz]/hloc,2.));}
           else{den3D[ix][iy][iz]=0.0;}
           m+=den3D[ix][iy][iz]*vol;
        }
      }
    }
    std::cout << "#3D MASS FROM EXTRAPOLATION = "<<m<<std::endl;
    
  }

  void print_scaleheight_grid(){ 
    for(long ix=0;ix<nx;++ix){
        for(long iy=0;iy<ny;++iy){
          double r = sqrt(x[ix]*x[ix]+y[iy]*y[iy]);
          std::cout << x[ix]<<" "<<y[iy]<<" "<<h[ix][iy]<<" "<<h[ix][iy]/r<<" "<<h[ix][iy]/dz<<std::endl;
        }
    }
  }

  void writeArrayRadMC3D(std::string fname, double conversion){

    std::ofstream myfile;
    
    myfile.open(fname.c_str(),std::ios::out);
    myfile << 1 << std::endl;
    long nall = nx*ny*nz;

    myfile << nall << std::endl;
    myfile << 1 << std::endl;

    myfile<<std::scientific;
    for(long iz=0;iz<nz;++iz){
      for(long iy=0;iy<ny;++iy){
        for(long ix=0;ix<nx;++ix){
          myfile << den3D[ix][iy][iz]*conversion <<std::endl;
        }
      }
    }
    myfile.close();
  }

  void writeArray2DRadMC3D(std::string fname, double conversion,
        const std::vector< std::vector<double> > &array){

    std::ofstream myfile;
    
    myfile.open(fname.c_str(),std::ios::out);
    myfile << 1 << std::endl;
    long nall = nx*ny*nz;

    myfile << nall << std::endl;
    myfile << 1 << std::endl;

    myfile<<std::scientific;
    for(long iz=0;iz<nz;++iz){
      for(long iy=0;iy<ny;++iy){
        for(long ix=0;ix<nx;++ix){
          myfile << array[ix][iy]*conversion <<std::endl;
        }
      }
    }
    myfile.close();
  }

  void writeFixedGridRadMC3D(std::string fname, double conversion){
    std::ofstream myfile;
    myfile.open(fname.c_str(),std::ios::out);
    myfile << 1 << std::endl;
    myfile << 0 << std::endl;
    myfile << 1 << std::endl;
    myfile << 0 << std::endl;
    myfile << "1 1 1"<<std::endl;
    myfile << nx << " " << ny << " " << nz <<std::endl;
    myfile << std::scientific;
    for(long ix=0;ix<nx;++ix){ myfile << x[ix]*conversion << " ";}
    myfile << std::endl;
    for(long iy=0;iy<ny;++iy){ myfile << y[iy]*conversion << " ";}
    myfile << std::endl;
    for(long iz=0;iz<nz;++iz){ myfile << z[iz]*conversion << " ";}
    myfile << std::endl;
    myfile.close();
  }

  void writeSlice3D(std::string fname,char dir, long slice){
    std::ofstream handle;

    std::cout << "# Writing slice for density information. In code units.\n";

    handle.open(fname.c_str(),std::ios::out);
    if(dir=='y'){
       for(long iy=0;iy<ny;++iy){
         for(long iz=0;iz<nz;++iz){
           handle << y[iy] << " " << z[iz] <<" "<<den3D[slice][iy][iz]<<std::endl;
         }
       }
    }else{
       for(long ix=0;ix<nx;++ix){
         for(long iz=0;iz<nz;++iz){
           handle << x[ix] << " " << z[iz] <<" "<<den3D[ix][slice][iz]<<std::endl;
         }
       }

    } 
    handle.close();

  }

  void dealloc(){
    x.clear();
    y.clear();
    z.clear();
    
    for(long ix = 0; ix<nx;++ix){
      h[ix].clear();
      for(long iy = 0; iy<ny;++iy) den3D[ix][iy].clear();
      den3D[ix].clear();
    }
    h.clear();
    den3D.clear();
  }

};



