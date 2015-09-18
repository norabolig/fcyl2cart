////////////////////////////////////////////////////////////////////////
//
// Written by nora bolig. 17 Sept. 2015
//
// Takes 2D Fargo output and XY Cartesian images and/or produces
// input files for Rad3DMC.
//
// See Readme.md for details.
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include "UnitConversion.hpp"
#include "Extrapolate3D.hpp"
#include "CommandLine.hpp"
using namespace std;


#define FNAME "fname_default.par" // This is defaulted.  Use this name for a list of default file names.

int main(int argc, char* argv[]) {

  int sdble = sizeof(double); // define commonly-used names. 
  double pi = acos(-1.0); 
  ifstream myfile;

  CommandLine cl; // set command line class

  myfile.open(FNAME,ios::in); // read default file names
  if(!myfile){ 
    cout << "\nDefault file name file 'fname_default.par' not found\n\n";
    return 0;
  }
  while(!myfile.eof()){
    string key;
    string value;
    myfile >> key >> value;
    cl.putValue(key,value);
  }
  myfile.close();

  bool ok=0; // now we can get the command line options.
  ok=cl.readCL(argc,argv);
  if(!ok){
    cl.printOptions();
    return 0;
  }
  //cl.printValues(); // uncomment if you want to see what file names will be used.

  {
    string options[]={"xy","rad3dmc","both"};
    string * first = find(options,options+3,cl.getValue("-product"));
    if (first==options+3){
      cout <<"\n"<< cl.getValue("-product") << " is not a valid option for -product (xy | rad3dmc | both)\n";
      cout << " Excluding the -product flag defaults to 'both'.\n\n";
      return 0;
    } 
  }

//
// Basic control.  Read in the param.par file
// and use that to determine grid size, grid spacing, etc.
// param.par should be produced by FARGO.
//
  long int nrad,nsec;
  double rmin,rmax,adindx;
  char gridSpacing;
  myfile.open(cl.getValue("-param").c_str(),ios::in);
  
  if(!myfile){
    cout << "\nParameter file "<<cl.getValue("-param")<<" not found\n\n";
    return 0;
  }

  while(!myfile.eof()){
     string pname, pvalue;
     myfile >> pname >> pvalue;
 
     if(      pname == "Nrad"){
         nrad = atol(pvalue.c_str());
     }else if(pname == "Nsec"){
         nsec = atol(pvalue.c_str());
     }else if(pname == "Rmin"){
         rmin = atof(pvalue.c_str());
     }else if(pname == "Rmax"){
         rmax = atof(pvalue.c_str());
     }else if(pname == "RadialSpacing"){
         gridSpacing=char(pvalue[0]);
     }else if(pname == "AdiabaticIndex"){
         adindx = atof(pvalue.c_str());
     }   
     
  }
  myfile.close();


//
// read in units file.
// This should be the units file produced by FARGO.
//
  UnitConversion units; //initialize units class.  Used to get and modify units.
  {
    double mass,length,time,tk;
    myfile.open(cl.getValue("-units").c_str(),ios::in);

    if(!myfile){
      cout << "\nUnits file "<<cl.getValue("-units")<<" not found\n\n";
      return 0;
    }

    myfile  >> mass >> length >> time >> tk;
    myfile.close();
    units.set_mass(mass);
    units.set_length(length);
    units.set_time(time);
    units.set_tk(tk);
    units.make();
  }

// 
// Set additional scaling
// The actual sim can be scaled in a way that is different from the units.par file. 
// This file is user-defined and will scale the units accordingly.  
// The star's mass is also defined here.
//
  double mstar=0.0;
  double mu=0.0;
  double dustToGas=0.0;
  {
    double rscale=0.0;
    int    MKS2CGS=0;
    myfile.open(cl.getValue("-scale").c_str(),ios::in);

    if(!myfile){
      cout << "\nScaling file "<<cl.getValue("-scale")<<" not found\n\n";
      return 0;
    }

    myfile >> mstar >> rscale >> mu >> dustToGas >> MKS2CGS;
    myfile.close();
    units.set_scale(rscale,MKS2CGS);
  }

//
// Calculate radial units.
// The used_radius.dat is not necessary.
//
  vector<double>radB(nrad);
  if (gridSpacing=='L'){
    double dlogr = (log10(rmax)-log10(rmin))/double(nrad);
    double logrmin = log10(rmin);
    for(long int r=0;r<nrad+1;++r){
      double base=10.0;
      radB[r]=pow(base,logrmin + dlogr * double(r));
    }
  }else{
    double dr = (rmax-rmin)/double(nrad);
    for(long int r=0;r<nrad+1;++r){ radB[r]=rmin + double(r)*dr;}  
    cout << "#Grid spacing not set to 'L', assuming linear."<<endl;
  }


//
// Now read 2D FARGO data from FARGO output.
//

  vector< vector<double> > density(nrad,vector<double>(nsec)); //density is really surface density here.

  long int ntot=nrad*nsec;
  double *tmp;
  tmp = new double [ntot]; // will be deleted below

  myfile.open(cl.getValue("-gas_in").c_str(),ios::in|ios::binary|ios::ate);

  if(!myfile){
    cout << "\n2D density file "<<cl.getValue("-gas_in")<<" not found\n\n";
    return 0;
  }

  streampos size=myfile.tellg();
  myfile.seekg(0,ios::beg);
  myfile.read((char *)tmp,sdble*ntot);
  myfile.close();
  cout << "#File content memory size is " << size << " and "<<sdble*ntot<<endl;  // basic file checks.
//
// binary read does not like double pointers, so convert 1D array to 2D. 
//
  { //keep scope local
    long iter=0;
    for(long irad=0;irad<nrad;++irad){
      for(long isec=0;isec<nsec;++isec){
        density[irad][isec]=tmp[iter++];
      }
    }
  }
//
// rinse, repeat for temperature file
//
  vector< vector<double> > temperature(nrad,vector<double>(nsec));

  myfile.open(cl.getValue("-tk_in").c_str(),ios::in|ios::binary|ios::ate);

  if(!myfile){
    cout << "\n2D temperature  file "<<cl.getValue("-tk_in")<<" not found\n\n";
    return 0;
  }

  size=myfile.tellg();
  myfile.seekg(0,ios::beg);
  myfile.read((char *)tmp,sdble*ntot);
  myfile.close();
  cout << "#File content memory size is " << size << " and "<<sdble*ntot<<endl;  
//
// make 2D array
//
  {
    long iter=0;
    for(long irad=0;irad<nrad;++irad){
      for(long isec=0;isec<nsec;++isec){
        temperature[irad][isec]=tmp[iter++];
      }
    }
  }
  delete[] tmp; // clean up
//
// Get dimensions for desired interpolation
// Now we need 2D and 3D Cartesian grid parameters.
// Get from user-supplied function. Even if 3D is not used, the code is expecting it.
// May change that later.
//
  long int nx,ny,nz,nradcell,nseccell;
  double dx,dy,dz;
  myfile.open(cl.getValue("-grid").c_str(),ios::in);
  myfile >> nx >> ny >> nz >> dx >> dy >> dz >> nradcell >> nseccell;
  myfile.close();
  cout << "#Got the following values for nx, ny, dx, and dy, respectively\n";
  cout << "#"<<nx << " "<<ny<<" "<< nz << " "<< dx << " " << dy<<" "<<dz<<endl;

//
// make density and temperature grids
//
  vector<double> x(nx);
  vector<double> y(ny);
  for(long int ix=0;ix<nx;++ix){x[ix] = (0.5-double(nx)*0.5 + double(ix))*dx;}
  for(long int iy=0;iy<ny;++iy){y[iy] = (0.5-double(ny)*0.5 + double(iy))*dy;}

  vector< vector<double> > sdenxy(nx, vector<double>(ny));
  vector< vector<double> > tkxy(nx, vector<double>(ny));

//
// fill grid values with zeros
//
  for(long ix=0;ix<nx;++ix){
    for(long iy=0;iy<ny;++iy){
         sdenxy[ix][iy]  =0.0;
         tkxy[ix][iy] =0.0;
     }
  }


//
// Cycle through cylindrical grid, dividing cells by mass, and integrating Cartesian grid.
// This effetively integrates the fraction of mass that corresponds to the area intersected
// by each cell.  This way, the mass is conservative.  Temperature is not a conserved
// quantity, so we use a mass-weighted temperature as a proxy for the energy density, which
// is conserved.  This should not be a problem as long as the adiabatic index is fixed.
//
  { // keep scope local
    double dsec = 2.0*pi/double(nsec);
    double mtot_sum=0.0;
    double m_sum=0.0;
    for(long r=0;r<nrad;++r){
      for(long s=0;s<nsec;++s){
        double dr = radB[r+1]-radB[r]; // define dr regardless of spacing type
        double drcell = dr/double(nradcell); // subdivide cell for integration
        double sect= dsec*double(s); // get lower bound of azimuthal sector
        mtot_sum += (density[r][s])*(pow(radB[r+1],2.0)-pow(radB[r],2.0))*dsec*0.5;
        for(long rcel=0;rcel<nradcell;++rcel){
          //
          // get mass of a "particle" used for integration. Masses are calculated
          // based on area of subcell.
          //
          double m = density[r][s]*(pow(radB[r]+(double(rcel)+1.0)*drcell,2.0)
                                   -pow(radB[r]+double(rcel)*drcell,2.0))*dsec/double(nseccell)*0.5;

          double mt = m*temperature[r][s];
          double radius = radB[r]+drcell*(double(rcel)+0.5);
          for(long scel=0;scel<nseccell;++scel){
             //
             // Now assign masses to Cartesian grid cells based on bucket 
             //
             double sectcel = sect + (double(scel)+0.5)*dsec/double(nseccell); 
             double xcel = radius*cos(sectcel);
             double ycel = radius*sin(sectcel);
             long ixcel = long(xcel/dx +0.5*(double(nx)));
             long iycel = long(ycel/dy +0.5*(double(ny)));
             sdenxy[ixcel][iycel]+=m;
             tkxy[ixcel][iycel]+=mt;
             m_sum+=m;
          }
        }
      } 
    } // loop over radial points
    cout << "# COMPARE MASSES (before after) : "<<mtot_sum<<" "<<m_sum<<endl;
  } // keep scope local

  for(long ix=0;ix<nx;++ix){
    for(long iy=0;iy<ny;++iy){  
      // 
      // remove mass from energy-like term to get temperature.
      // turn mass into surface density.
      // keep in code units.
      //
      if(sdenxy[ix][iy]>0.0)tkxy[ix][iy]/=sdenxy[ix][iy]; 
      sdenxy[ix][iy]/=(dx*dy);
    }
  }
 
//
// spit out mages if wanted.
// units for surface density and temperature are converted tp physical
// positions are outputted as scaled code units (as usually 1 code unit has some special meaning). 
//
  if(cl.getValue("-product") == "xy" || cl.getValue("-product") == "both"){
    cout << "# Writing out 2D cyl2cart files\n";
    ofstream outfile;
    outfile.open(cl.getValue("-gas_out").c_str(),ios::out);
    for(long ix=0;ix<nx;++ix){
      for(long iy=0;iy<ny;++iy){  
        outfile << x[ix]*units.scale() <<" " << y[iy]*units.scale() << " " << sdenxy[ix][iy]*units.sig() << endl;
      }
    }
    outfile.close();

    outfile.open(cl.getValue("-tk_out").c_str(),ios::out);
    for(long ix=0;ix<nx;++ix){
      for(long iy=0;iy<ny;++iy){  
        outfile << x[ix]*units.scale() << " " << y[iy]*units.scale() << " " << tkxy[ix][iy]*units.tk() << endl;
      }
    }
    outfile.close();
  }
//
// Check scale heights.
//
  if(cl.getValue("-product") == "rad3dmc" || cl.getValue("-product") == "both"){
    cout << "# Calculating 3D structure\n";

    Extrapolate3D xtr3d;

    xtr3d.initialize(nx,ny,nz,dx,dy,dz,mstar,adindx,mu,x,y);

    xtr3d.set_scaleheight(tkxy);

    xtr3d.CreateIso3D(sdenxy);

  //xtr3d.writeSlice3D("test.txt",'x',(nx/2)); // example for producing a slice SLICE

    {
      double lconv=units.length();
      xtr3d.writeFixedGridRadMC3D(cl.getValue("-amr_grid"),lconv);
    }

    {
      double dconv=units.den()*dustToGas;
      cout << std::scientific<<"#CHECK DENSITY DUST CONVERSION "<<dconv<<std::endl;
      xtr3d.writeArrayRadMC3D(cl.getValue("-dust_density"),dconv);
    }

    {
      double tconv=units.tk();
      xtr3d.writeArray2DRadMC3D(cl.getValue("-dust_tk"),tconv,tkxy);
    }

    xtr3d.dealloc();
  }

}
