#include <cmath>
#include <memory>
#include <time.h>
#include <string>
#include <iomanip>

#include "System.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "include/mpi.hpp"
#include "VerletListAdress.hpp"
#include "FixedTupleListAdress.hpp"

#include "analysis/Temperature.hpp"
#include "analysis/Pressure.hpp"
#include "analysis/PressureTensor.hpp"

#include "integrator/VelocityVerlet.hpp"
#include "integrator/Adress.hpp"
#include "integrator/LangevinThermostat.hpp"

#include "storage/DomainDecompositionAdress.hpp"
#include "storage/CellGrid.hpp"
#include "storage/NodeGrid.hpp"

#include "esutil/RNG.hpp"

#include "interaction/VerletListAdressInteractionTemplate.hpp"
#include "interaction/FixedPairListInteractionTemplate.hpp"
#include "interaction/LennardJones.hpp"
#include "interaction/LennardJonesCapped.hpp"
#include "interaction/FENE.hpp"
#include "interaction/Morse.hpp"
#include "interaction/Tabulated.hpp"

#include "include/esconfig.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "main/espressopp_common.hpp"
#include "boost/make_shared.hpp"
#include "io/DumpXYZ.hpp"


static void trim (std::string& s) {
  //trim left
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
    return !std::isspace(ch);
  }));
  
  //trim right
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

static void erase_char(std::string& s, char c) {
  for(std::string::iterator it=s.begin();it!=s.end();++it) {
    if(*it == c) { 
      it = s.erase(it);
      it--;
    }
  }
}

static std::vector<std::string> split_string(std::string s) {
  std::vector<std::string> result;
  size_t posS = 0;
  size_t posE = 0;
  std::string token;
  do {
    posE = s.find(" ", posS);
    token = s.substr(posS, posE-posS);
    if(token.find_first_not_of(' ') != std::string::npos)
    {
      result.push_back(token);
    }
    posS = posE+1;
  } while(posE != std::string::npos);

  return result;
}
static void read_espresso_data(std::string input_file,
                        espressopp::real& Lx, espressopp::real& Ly, espressopp::real& Lz,
                        std::vector<espressopp::real>& x,
                        std::vector<espressopp::real>& y,
                        std::vector<espressopp::real>& z,
                        std::vector<int>& type,
                        std::vector<int>& q,
                        std::vector<espressopp::real>& vx,
                        std::vector<espressopp::real>& vy,
                        std::vector<espressopp::real>& vz,
                        std::vector<espressopp::real>& fx,
                        std::vector<espressopp::real>& fy,
                        std::vector<espressopp::real>& fz,
                        std::vector<std::pair<int,int>>& bondpairs) 
{
  std::ifstream input (input_file);
  std::string line;
  
  std::vector<std::string> props;
  
  bool variable = false;
  bool particles = false;
  bool bonds = false;
  
  if (input.is_open()) {
    while (getline(input,line)) {
      //Trim line
      trim(line);
      
      //skip empty lines
      if (line == "") continue;
      
      if (line.find("variable") != std::string::npos) {
        variable = true;
        continue;
      }
      
      if (variable) {
        if(line.compare("}") == 0) {
          variable = false;
          continue;
        }
        erase_char(line, '{');
        erase_char(line, '}');

        std::vector<std::string> tmp = split_string(line);
        if (tmp[0].compare("box_l") == 0) {
          Lx = std::atof(tmp[1].c_str());
          Ly = std::atof(tmp[2].c_str());
          Lz = std::atof(tmp[3].c_str());
        }
        continue;
      }
      
      if (line.find("particles") != std::string::npos) {
        particles = true;
        erase_char(line, '}');
        line = line.substr(12);
        std::vector<std::string> tmp = split_string(line);
        for(auto it=tmp.begin();it!=tmp.end();++it) {
          if (it->compare("id") == 0) props.push_back("id");
          if (it->compare("pos") == 0) props.push_back("pos");
          if (it->compare("type") == 0) props.push_back("type");
          if (it->compare("q") == 0) props.push_back("q");
          if (it->compare("v") == 0) props.push_back("v");
          if (it->compare("f") == 0) props.push_back("f");
        }
        continue;
      }
      
      if(particles) {
        trim(line);
        if (line.compare("}") == 0) {
          particles = false;
          continue;
        }
        
        erase_char(line, '{');
        erase_char(line, '}');
        std::vector<std::string> tmp = split_string(line);
        
        int index = -1;
        for(auto it=props.begin();it!=props.end();++it) {
          index++;
          
          if (it->compare("id") == 0) continue;
          if (it->compare("pos") == 0) {
            x.push_back(std::stod(tmp[index]));
            y.push_back(std::stod(tmp[index+1]));
            z.push_back(std::stod(tmp[index+2]));
            index += 2;
            continue;
          }
          if (it->compare("type") == 0) {
            type.push_back(std::stoi(tmp[index]));
            continue;
          }
          if (it->compare("q") == 0) {
            q.push_back(std::stoi(tmp[index]));
            continue;
          }
          if (it->compare("v") == 0) {
            vx.push_back(std::stod(tmp[index]));
            vy.push_back(std::stod(tmp[index+1]));
            vz.push_back(std::stod(tmp[index+2]));
            index += 2;                        
            continue;
          }
          if (it->compare("f") == 0) {
            fx.push_back(std::stod(tmp[index]));
            fy.push_back(std::stod(tmp[index+1]));
            fz.push_back(std::stod(tmp[index+2]));
            index += 2;                        
            continue;
          }          
        }
      }
      
      if (line.find("bond") != std::string::npos) {
        bonds = true;
        continue;
      }
      
      if (bonds) {
        trim(line);
        if (line.compare("}") == 0) {
          bonds = false;
          continue;
        }
        
        erase_char(line, '{');
        erase_char(line, '}');
        std::vector<std::string> tmp = split_string(line);
        
        if (tmp.size() > 2) {
          int first = std::stoi(tmp[0]);
          for(int i=2; i<tmp.size();i+=2) {
            bondpairs.push_back(std::make_pair(first, std::stoi(tmp[i])));
          }
        }
      }
    
    }
  }
  input.close();  
}


static void write_tab_file(boost::shared_ptr<espressopp::interaction::Morse> potMorse,
                                  std::string tabMorse, int N, espressopp::real low, 
                                  espressopp::real high, espressopp::real body) {
  std::ofstream output(tabMorse);
  espressopp::real delta = (high-low) / (N-1);

  if(output.is_open()) {
    for(int i=0;i<N;++i) {
      espressopp::real r = low + i * delta;
      espressopp::real energy = potMorse->computeEnergy(r);
      espressopp::real force;
      espressopp::Real3D dist(r,0,0);
      if (body == 2) { //This is for 2-body potentials
        force = potMorse->computeForce(dist)[0];
      }
      else { // This is for 3- and 4-body potentials
        force = potMorse->computeForce(r);
      }
      output << std::fixed << std::setw(15) << std::setprecision(8) << r << " ";
      output << std::fixed << std::setw(15) << std::setprecision(8) << energy << " ";
      output << std::fixed << std::setw(15) << std::setprecision(8) << force << "\n";

    }
    output.close();
  }
  else {
    std::cout << "ERROR: write_tab_file(): could not open file!\n";
  }
}












