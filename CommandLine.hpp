//
//
//
#include <map>

class CommandLine{

  std::map<std::string,std::string> entries;

  public:

  CommandLine(){
    entries.insert (std::pair<std::string,std::string>("-gas_in","gas.dat"))  ;
    entries.insert (std::pair<std::string,std::string>("-gas_out","gout.dat"))  ;
    entries.insert (std::pair<std::string,std::string>("-tk_in","tk.dat"))  ;
    entries.insert (std::pair<std::string,std::string>("-tk_out","tout.dat"))  ;
    entries.insert (std::pair<std::string,std::string>("-units","units.dat"))  ;
    entries.insert (std::pair<std::string,std::string>("-scale","scale.par"))  ;
    entries.insert (std::pair<std::string,std::string>("-param","param.par"))  ;
    entries.insert (std::pair<std::string,std::string>("-grid","cart.par"))  ;
    entries.insert (std::pair<std::string,std::string>("-amr_grid","amr_grid.inp"))  ;
    entries.insert (std::pair<std::string,std::string>("-dust_density","dust_density,inp"))  ;
    entries.insert (std::pair<std::string,std::string>("-dust_tk","dust_temperature.inp"))  ;
    entries.insert (std::pair<std::string,std::string>("-product","both"));
    entries.insert (std::pair<std::string,std::string>("-h"," "));
  };


  void printOptions(){
    std::cout<<"Option List:\n";
    std::cout<<"  -gas_in       File name for surface density input\n";
    std::cout<<"  -gas_out      File name for gridded surface density output\n";
    std::cout<<"  -tk_in        File name for temperature input (2D)\n";
    std::cout<<"  -tk_out       File name for gridded temperature output\n";
    std::cout<<"  -units        File name for units file from FARGO\n";
    std::cout<<"  -scale        File name for scaling file used to rescale units\n";
    std::cout<<"  -param        File name for parameters file from FARGO\n";
    std::cout<<"  -grid         File name for grid setup for making rad3dmc files\n";
    std::cout<<"  -amr_grid     File for amr_grid.inp in rad3dmc (don't change unless you are prepared to accept the consequences\n";
    std::cout<<"  -dust_density File for dust_density,inp in rad3dmc (don't change unless you are prepared to accept the consequences\n";
    std::cout<<"  -dust_tk      File for dust_temparture.inp in rad3dmc (don't change unless you are prepared to accept the consequences\n";
    std::cout<<"  -product      Define products of program. Options include 'xy', 'rad3dmc', or 'both' (default).\n";
    std::cout<<"  -h            Displays this list\n";

  }

  void printValues(){
    std::cout <<" Printing all keys and values:\n";
    std::map<std::string,std::string>::iterator ikey;
    for(ikey = entries.begin(); ikey != entries.end();++ikey){
      std::cout << ikey->first<<" "<<ikey->second<<std::endl;
    }
   
  }

  bool readCL(int argc,char* argv[]){

    bool OK=1;
    std::map<std::string,std::string>::iterator ikey;

    for(int a=1;a<argc;a+=2){
       ikey = entries.find(std::string(argv[a]));
       if(ikey != entries.end()){
           if(ikey->first == "-h"){return OK=false;}
           entries.at(std::string(argv[a])) = std::string(argv[a+1]);
           //std::cout << ikey->second<<std::endl;
       }else{
           std::cout << "Invalid option: '"<<argv[a]<<"'"<<std::endl;
           return OK=false;
       }

    }
    return OK;

  }

  void putValue(std::string key, std::string value){
     std::map<std::string,std::string>::iterator ikey;
     ikey = entries.find(key);
     if(ikey !=entries.end()) ikey->second=value;
  }

  std::string getValue(std::string key){
     std::map<std::string,std::string>::iterator ikey;
     ikey = entries.find(key);
     if(ikey !=entries.end()){ return ikey->second;}
     else{
       std::cout<<" WARNING WARNING WARNING: Invalid Key in Map. Returning foo.bar.\n";
       return "foo.bar";
     }
  }
  

};
